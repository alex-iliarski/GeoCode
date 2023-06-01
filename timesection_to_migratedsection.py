from segpy.reader import create_reader
from segpy.writer import write_segy, TraceHeaderRev1
from obspy.io.segy.segy import _read_segy, SEGYFile

import numpy as np
import math
from numba import jit

"""
Docs:
    Segpy https://segpy.readthedocs.io/en/latest/

    NumPy https://docs.scipy.org/doc/numpy-1.13.0/reference/index.html

    Obspy https://docs.obspy.org/index.html
    
    Numba https://numba.pydata.org/
"""


if __name__ == '__main__':

    # read in timesection data from file into SEGYFile object
    # file is assumed to be in the data/ folder
    filename = input('Enter the .sgy file: ')

    # keep taking inputs until valid input or user quits
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

    # CONSTANTS AND VARIABLES
    input_traces = infile.traces #list of all input traces
    num_traces = len(input_traces)  # number of total input traces
    trace_length = input_traces[0].data.size  # number of points in a trace

    # create numpy matrix to store timesection input data in, faster
    timesection_traces = np.zeros((num_traces, trace_length), dtype = np.int16)
    for i in range(num_traces):
        timesection_traces[i] = np.asarray(input_traces[i].data, dtype = np.int16)


    velocity = 2 #average velocity of propogation of reflected wave (meters/second)
    velocity_squared = velocity ** 2 #used in creation of migrated time section formula

    # delta_x = timesection_traces[1].header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
    # delta_x /= 2  # half of the physical distance between each successive data collection (meters)
    delta_x = 15
    total_distance = num_traces * delta_x #total physical distance encompassed by this data

    sample_interval = (infile.binary_file_header.sample_interval_in_microseconds / 1000)  # sample interval in milliseconds

    #CREATE MIGRATED SECTION
    migrated_result_traces = np.zeros((num_traces, trace_length)) # result traces default to zero, then get added to

    distances_array = np.zeros(num_traces) #array of the physical distance of each measurement at a trace number (ex: 0m, 15m, 30m, etc)
    for i in range(num_traces): #populate distances_array
        distances_array[i] = i * delta_x

    time_squared_array = np.zeros(trace_length, dtype = np.int32) # array of time squared values for each point of a trace (milliseconds)
    for point in range(trace_length): #populate time_squared_array
        time_squared_array[point] = int ((sample_interval * point) ** 2) #time squared at that point

    num_x_before_and_after = 100 #number of x's looped through before and after each y in the following loop
        # must be less than half the number of traces
    num_x_before_and_after_plus_one = num_x_before_and_after + 1
    #TODO: THIS NUMBER CAN BE MODIFIED

    # "trace number" refers to the index of a particular trace

    # the i-th value of this list corresponds to the value of the expression ((4 * ((i * delta_x) ** 2)) / velocity_squared)
        # where i*delta_x is also the difference in trace numbers of trace x and trace y
        # the maximum difference between x and y is (num_x_before_and_after + 1)
    expression_values = np.zeros(num_x_before_and_after_plus_one)
    for i in range(num_x_before_and_after_plus_one): #populate expression_values
        expression_values[i] = ((4 * ((i * delta_x) ** 2)) / velocity_squared)

    # table storing the index values for a particular (x - y) value and time value according to the formula:
        # result_T = math.sqrt(((sample_interval * point_num) ** 2) + ((4 * ((x - y) ** 2)) / velocity_squared))
        # index = math.floor(result_T / sample_interval)
    index_matrix = np.zeros((num_x_before_and_after_plus_one, trace_length), dtype=np.int32)
    for x_trace_minus_y_trace in range(num_x_before_and_after_plus_one): #populate index_matrix
        exp_val = expression_values[x_trace_minus_y_trace] #value of exprssion for particular difference in x and y trace nums
        for point_num in range(trace_length):
            t_squared = time_squared_array[point_num] #time squared value for a particular point number in a trace
            index_matrix[x_trace_minus_y_trace][point_num] = math.floor((math.sqrt(t_squared + exp_val)) / sample_interval)

    sample_num_at_100_ms = math.floor(100 / sample_interval) # the number of the point within a trace that occurs at 100 milliseconds
        # to avoid certain unnecessary calculations, we start at 100 ms instead of 0 ms

    # use njit to make the process faster
    # @jit(nopython=True)
    def populate_migrated_section(timesection_traces, num_traces, sample_num_at_100_ms, trace_length):
        for y_trace_num in range(num_traces): #loop through each distance, y
            #assign variables to loop through a particular amount of traces before and after trace y
            y_minus = y_trace_num - num_x_before_and_after
            y_plus = y_trace_num + num_x_before_and_after + 1

            # to keep inside boundaries of list
            if y_minus < 0:
                y_minus = 0
            if y_plus > num_traces:
                y_plus = num_traces

            # print(y_trace_num)

            # loop through num_x_before_and_after trace numbers before and after y, with y fixed.
            for x_trace_num in range(num_traces)[y_minus:y_plus]:

                x_y_difference = abs(x_trace_num - y_trace_num) #difference in trace numbers of trace x and trace y (absolute value)

                for point_num in range(sample_num_at_100_ms, trace_length): #loop through each point of trace y, starting at 100ms point

                    index = index_matrix[x_y_difference][point_num] #index where we take number from?? TODO:: CORRECT?
                    if index < trace_length:
                        migrated_result_traces[y_trace_num][point_num] += timesection_traces[y_trace_num][index]

    populate_migrated_section(timesection_traces, num_traces, sample_num_at_100_ms, trace_length)

    #repopulate infile with migrated section data
    for i in range(num_traces):
        input_traces[i].data = np.asarray(migrated_result_traces[i], dtype = np.int16)

    #write transformed data to new .sgy file
    outfile_name = 'data/' + filename[:-4] + '_MIGRATEDSECTION.sgy'
    outfile = open(outfile_name, 'wb')
    # infile.textual_file_header.replace(b'C40  END EBCDIC       ', b'C40 END EBCDIC        ') #warning handle
    infile.write(outfile_name)
    outfile.flush()
    outfile.close()

# NUMBA, WELD