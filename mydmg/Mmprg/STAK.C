#include        <math.h>
#include	<stdio.h>
#include	<time.h>
#include	<process.h>
#include 	<stdlib.h>
#include 	<stdlib.h>
#include        <disp.h>
#include        <malloc.h>


void main()

{

    short mv[3000],mlv[1000],mr[3000],im[3620],om[3120];
    int mkt[3000];

 static long ln[2500],r5,rnd,nr,k1=1L,ic;

static char  nam[12];

 int K1 = 76;

 long mm[20],*rez;

 short cod,ier[2];

 long samp,si,krb,rog,dog,nw,l,l1,l2,ktn,ni,ng,itn,itk,nwb,pr,ti,bt,ll,

tn1,tk1,nx,mxi,j,tn,tk,sam,mx,my,i,n,kn,ngt,rx,ii,kwb,dv,clt,cnc,cpr

,mt[8],amp,r1,r2,r3;

static float xt,xtn,xtk,xto,xsi,zv[3000],xrx,xv,xvv,xkn,xr2,xr1;

   rez=(long*)calloc(6000000, sizeof(long));

         pr=0;

       clt=1;


      inifila(0,"LINE",ier);


    if(*ier != 0L) goto S999;


    rdfila(0,ln,10000);

    wtfila(ier);


      if(*ier != 0L) goto S999;


    *mt=1;


    bits(ln,mm,mt);



    gdopfl(ln,mkt);



      mvic((ln+150),1,nam,1,10);


    wrfini(1,nam,ier);


    if(*ier != 0) goto S999;



    *mt=2;

    *(mt+1)=0;
    

    bits(ln,mm,mt);


    mvii((ln+300),1,mlv,1,2000);
    samp = *(ln+2);
    si = *(ln+3);
    xsi=si;
    sam = samp*si-si;

    bt=3000;


    tn = *(ln+176);

    if(tn >= sam || tn < 100L) tn = 100L;

    xtn=(float)tn;


    tk = *(ln+177);
    if(tk > sam) tk = sam;

    xtk=(float)tk;


    tn1=tn/si+1;

    tk1=tk/si;

    kwb = samp**(ln+71)+240L;
    krb = samp**(ln+70)+240L;

    rog = *(ln+112);
    dog = *(ln+113);

    rnd=(dog-rog)+1;

    rdfrd(0,im,3600);
    wrfwr(0,im,3600);

    nw = 0L;



    *mt=5;


    bits(ln,mm,mt);



        for(ll=1; ll<=rnd; ll++) {

            for(i=1L; i<=samp; i++) *(rez+(short)i-1+((short)ll-1)*bt) = 0L;

        }





    *mt=6;

    *(mt+1)=0;


    *(mt+2)=(dog-rog)*clt/100;


    bits(ln,mm,mt);

    cnc=kbhit();

/*

                printf(" cnc %ld ic %ld  rog %ld dog %ld nr %ld\n",
                cnc,ic,rog,dog,nr);
                getch();
*/

    if(cnc ==0)goto S1003;

    r1=getch();



    if(r1==27)goto S1001;

S1003:



    for(l=rog; l<=dog; l++) {



            *(ln+49) = (*(mkt+(short)l-1)-1L)*krb+3600L;

                  kn = *(mkt+(short)(l+1L)-1)-*(mkt+(short)l-1);

                if(kn <= 0L) goto S490;

            if(*(mt+1)==1) goto S120;


            dvelin(l,mv,mr,mlv,ln);

        for(i=1L; i<=samp; i++) *(zv+(short)i-1) = (float)(
          (double)*(mv+(short)i-1)/1000.0);


S120:

        for(ni=1L; ni<=kn; ni++) {


                readts(im,om,ln,ier);

                if(*ier != 0L) goto S999;

                *(ln+49) += krb;

                ngt = *(ln+270);

                if(ngt != l) goto S999;

                rx = *(ln+271);


                  xrx=(float)rx;



       if(*(ln+229) ==0)goto S999;

/*
                printf(" ni %ld kn %ld  rog %ld dog %ld rx %ld\n",
                ni,kn,rog,dog,rx);
                getch();
*/

        for(ii=tn1; ii<=tk1; ii++) {


        xto=(float)(ii-1)*si;


        if(*(mt+1)!=1)goto S110;

        xt=(float)(double)xto;

        goto S130;

S110:
            xv=(float)*(zv+(short)ii-1);

            xvv=(float)(double)xv*(double)xv;


            xt=(float)sqrt((double)xto*(double)xto+(double)xrx*(double)xrx/

            (double)xvv);

//    printf("xt %f xto %f xv %f xrx %f xtn %f xtk %f\n",
//       xt,xto,xv,xrx,xtn,xtk);
//       r5=getch();

S130:
            if(xt < xtn ||xt > xtk)goto S10;

            ti=(double)xt/(double)xsi+1;

            r1=*(om+ti-1);


           *(rez+(short)ii-1+((short)l-1)*bt)+=r1;

S10:
         continue;

            }

/*
                printf(" 2 %ld ni %ld kn %ld  rog %ld dog %ld rx %ld\n",
                2,ni,kn,rog,dog,rx);
                getch();
*/



        }




/*
                printf(" 3 %ld ni %ld kn %ld  rog %ld dog %ld rx %ld l %ld\n",
                3,ni,kn,rog,dog,rx,l);
                getch();
*/

        xkn=(float)kn;

        xr1=(float)sqrt((double)xkn);

        xr2=(float)sqrt((double)xr1);

        xkn=(float)(double)xr1*(double)xr2;



        for(ii=tn1; ii<=tk1; ii++) {


            xr1=*(rez+(short)ii-1+((short)l-1)*bt);

            xr2=(float)(double)xr1/(double)xkn;


            r1=xr2;


            *(rez+(short)ii-1+((short)l-1)*bt)=r1;


                  }
/*

                printf(" kn %ld  rog %ld dog %ld l %ld\n",
                kn,rog,dog,l);
                getch();


                printf(" l %ld kn %ld  rog %ld dog %ld l %ld\n",
                l,kn,rog,dog,l);

*/

S490:
		continue;

                    }





        amp=0;

        for(l=rog; l<=dog; l++) {

//        printf("  slav-results-lN.O.CDP= %ld\n",l);

                n = l-rog+1L;

            for(i=1L; i<=samp; i++) {         //13


                 r2=(float)(long)(*(rez+(short)i-1+((short)n-1)*bt));


                if(r2 > amp)amp=r2;

                }

                }

                if(amp == 0)amp=1.;



        for(l=rog; l<=dog; l++) {              //12

            for(i=1L; i<=120L; i++) *(im+(short)i-1) = 0;
            for(i=1L; i<=samp; i++) *(im+(short)(i+120L)-1) = 0;

                n = l-rog+1L;

            for(i=1L; i<=samp; i++) {         //13


                 xr1=(float)(*(rez+(short)i-1+((short)n-1)*bt));



                 r3=xr1*16000./amp;



                 if(r3 < -16000)r3=-16000;

                 if(r3 > 16000)r3=16000;


              *(im+(short)(i+120L)-1) =r3;

S530:
		continue;

                    }



            *(ln+74) = 0;
            *(ln+75) = 0;

	    ic = l-rog+1L;




            dpastty(ic,ic,k1,(im),
            (im+14),(im+44),ln);
            nr=(ic-1L)*kwb+3600;


            wrfwr(nr,im,kwb);
            wrfwt(ier);
            if(*ier != 0L) goto S999;


    if(l != rog)goto S1002;


    *mt=9;

    *(mt+3)=rog;


    bits(ln,mm,mt);


     
S1002:


    *mt=10;

    *(mt+3)=l;


    bits(ln,mm,mt);




S7:
        continue;

    }



S998:

    pr=1;


    *mt=7;


    bits(ln,mm,mt);



    goto S1001;


S999:


             if(pr!=0)goto S1001;


    *mt=8;


    bits(ln,mm,mt);


S1001:

             wrclose();

             rdfclose();



   disp_open();


   disp_move(0,0);

   disp_eeop();

   disp_close();



}
