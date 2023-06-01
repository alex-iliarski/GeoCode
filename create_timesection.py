from segpy.reader import create_reader
from segpy.writer import write_segy, TraceHeaderRev1
from obspy.io.segy.segy import _read_segy, SEGYFile

import numpy as np
import math

"""
Docs:
    Segpy https://segpy.readthedocs.io/en/latest/
    
    NumPy https://docs.scipy.org/doc/numpy-1.13.0/reference/index.html
    
    Obspy https://docs.obspy.org/index.html
"""


if __name__ == '__main__':

    #read in data from file into SEGYFile object
    #file is assumed to be in the data/ folder
    filename = input('Enter the .sgy file: ')

    #keep taking inputs until valid input or user quits
    valid_data = False
    while not valid_data:
        try:
            infile = _read_segy('data/' + filename)
            valid_data = True
        except:
            response = input('\nFile could not be read.\nEnter a new file or enter \'quit\' to exit the program: ')
            if response == 'quit':
                quit()
            else:
                filename = response

    print(infile.__str__())

    #CONSTANTS AND VARIABLES
    trace_lst = infile.traces #list of all input traces
    num_traces = len(trace_lst)  # number of total input traces
    trace_length = trace_lst[0].data.size  # number of points in a trace

    V = 2 #average velocity of propogation of reflected wave (meters/second)
    t = 0.2  # angular coefficient of lines T_0i (seconds/meter) (VOTR)
    t = 0.3 # (VOTR_2) 
    V_stk = (2*V) / ((4-((t**2)*(V**2)))**0.5) #stacking velocity

    delta_x = trace_lst[1].header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
    delta_x /= 2  #half of the physical distance between each successive data collection

    sample_interval = (infile.binary_file_header.sample_interval_in_microseconds/1000) #sample interval in milliseconds

    traces_per_seismogram = 120 #number of traces in a seismogram

    num_seismograms = math.ceil(num_traces/traces_per_seismogram) #number of seismograms in the input file, rounded up.

    # CREATE TIMESECTION
    timesection_result_traces = np.zeros((num_seismograms, trace_length)) # result traces default to zero, then get added to
    # num_fails = 0

    offset_array = np.zeros(traces_per_seismogram) #stores the offset for a particular trace in a seismogram
    for j in range(traces_per_seismogram):
        offset_array[j] = (4 * ((j*delta_x)**2)) / (V_stk**2) #formula for offset

    T0_squared_array = np.zeros(trace_length) #stores the time value for each point in a trace, squared.
        #(4ms, 8ms, 12ms, ... ) -> (16, 64, 144, ... )
    for k in range(trace_length):
        T0_squared_array[k] = ((k+1) * sample_interval) ** 2

    for trace in range(num_traces): #loop through each trace in the input .sgy file
        trace_num_within_seismogram = trace % traces_per_seismogram
        offset = offset_array[trace_num_within_seismogram] #offset for the given trace

        for point in range(trace_length): #loop through each point of the trace
            T = math.sqrt(T0_squared_array[point] + offset)
            index = math.floor(T / sample_interval) - 1  # index we grab number from in the input trace and index we add to in resultant trace

            if index < trace_length:
                timesection_result_traces[trace // traces_per_seismogram][index] += trace_lst[trace].data[index]

            # if index >= trace_length:
            #     num_fails += 1
            # else:
            #     timesection_result_traces[trace // traces_per_seismogram][index] += trace_lst[trace].data[index]

    # print('total number of points: ' + str(trace_length*traces_per_seismogram*num_seismograms))
    # print('number of failed points: ' + str(num_fails))

    #repopulate infile with transformed data
    for i in range(num_seismograms):
        trace_lst[i].data = np.asarray(timesection_result_traces[i], dtype = np.int16)
    infile.traces = trace_lst[:num_seismograms] #remove excess traces

    #write transformed data to new .sgy file
    outfile_name = 'data/' + filename[:-4] + '_TIMESECTION.sgy'
    outfile = open(outfile_name, 'wb')
    # infile.textual_file_header.replace(b'C40  END EBCDIC       ', b'C40 END EBCDIC        ') #warning handle
    infile.write(outfile_name)
    outfile.flush()
    outfile.close()



    # for i in range(num_seismograms): #loop through each seismogram
    #     for j in range(traces_per_seismogram): #loop through each trace in a seismogram
    #         # x = j * delta_x #distance source point
    #         # offset = (4 * (x**2)) / (V_stk**2)
    #         offset = offset_array[j] #for a given trace, the offset will be the same for all points in that trace
    #
    #         input_trace_number = i * traces_per_seismogram + j
    #
    #         for k in range(trace_length): #loop through each point in a trace
    #             # T = math.sqrt(((k+1) * sample_interval)**2 + offset)
    #             T = math.sqrt(T0_squared_array[k] + offset)
    #             index = math.floor(T/sample_interval) - 1 #index we grab number from in the input trace and index we add to in resultant trace
    #
    #             if index > 1248 or input_trace_number > 11998:
    #                 num_fails += 1
    #             else:
    #                 # print(i)
    #                 # print(index)
    #                 # print(input_trace_number)
    #                 result_traces[i][index] += trace_lst[input_trace_number].data[index]


    # T0 = 500 #0.5s or 500ms ? # double travel time of the reflected wave along the normal to the first boundary at the coordinate system origin
    # delta_T0 = 1000 #1s or 1000ms? # step between values T_0i

    # #create empty numpy matrix used to hold and transform data
    # data = np.empty([num_traces, trace_length]) #make into integer? not float?
    #
    # #populate numpy data matrix
    # for i in range(num_traces):
    #     current_trace = trace_lst[i].data
    #     data[i] = np.array(current_trace)
    # data *= 3

    # infile = open('data/CUTE.sgy', 'rb')
    # outfile = open('data/CUTE_COPY.sgy', 'wb')
    # in_data = create_reader(infile)
    # # how to transform data?
    # write_segy(outfile, in_data)

    # outfile = open('data/CUTE_COPY.sgy', 'wb')
    # infile.write('data/CUTE_COPY.sgy') #copies contents of infile to new file
    # outfile.flush()
    # outfile.close()

    # with open('data/CUTE.sgy', 'rb') as segy_in_file:
    #     # The seg_y_dataset is a lazy-reader, so keep the file open throughout.
    #     seg_y_dataset = create_reader(segy_in_file)  # Non-standard Rev 1 little-endian
    #     print(seg_y_dataset.num_traces())
    #     # # Write the seg_y_dataset out to another file, in big-endian format
    #     # with open('seismic_big.sgy', 'wb') as segy_out_file:
    #     #     write_segy(segy_out_file, s  bbb  eg_y_dataset, endian='>')  # Standard Rev 1 big-endian
    #
    #     data = []
    #     for i in range(0, seg_y_dataset.num_traces()):
    #         trace = seg_y_dataset.trace_samples(i)
    #         data.append(trace)
    #
    #     data_matrix = np.array(data, np.int32)
    #     ## each row is a trace
    #     # print(data_matrix)
    #     segy_in_file.close()

    # with open('data/CUTE.sgy', 'rb') as infile:
    #     outfile = open('data/CUTE_COPY.sgy', 'wb')
    #     in_data = create_reader(infile)
    #     # how to transform data?
    #     write_segy(outfile, in_data)

    # # Open the input data file
    # infile = open('data/CUTE.sgy', 'rb')
    # in_data = create_reader(infile)
    # print(in_data.num_traces())
    #
    # data = []
    # for i in range(0, in_data.num_traces()):
    #     trace = in_data.trace_samples(i)
    #     data.append(trace)
    #
    # data_matrix = np.array(data, np.int32)
    #
    # outfile = open('data/CUTE_COPY.sgy', 'wb')
    # h = in_data.trace_header(0)
    #
    # write_segy(outfile, in_data)

    # output_segy_spec: segyio.spec = segyio.spec()
    # output_segy_spec.tracecount = in_data.num_traces()
    # # output_segy_spec.samples = in_data.trace_samples(0)
    # output_segy_spec.format = SegySampleFormat.SIGNED_INTEGER_4_BYTE
    # output_segy_spec.samples = data
    #
    # result_file = segyio.create("data/result.sgy", output_segy_spec)

    # write_segy(result_file, in_data)

    # with open('data/CUTE.sgy', 'rb') as infile:
    #     with open('data/CUTE_COPY.sgy', 'wb') as outfile:
    #         in_data = create_reader(infile)
    #         # how to transform data?
    #         write_segy(outfile, in_data)

