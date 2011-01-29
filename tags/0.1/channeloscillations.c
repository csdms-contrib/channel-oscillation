#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1
#define PI 3.141592653589793

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float gasdev(idum)
int *idum;
{
        static int iset=0;
        static float gset;
        float fac,r,v1,v2;
        float ran3();

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

main()
{    FILE *fp1;
     int idum,*iup,*idown,count,printed,outputinterval,length_y,lattice_size_x,i;
     float courantwave,courantdiff,maxdelw,del,*h2,*w2,h0,*flux,*deltaw,*deltah;
     float w0,simulationlength,deltax,elapsedtime,timestep,channelslope;
     float *height,*heightold,*width,*widthold,maxheightchange; 

     fp1=fopen("channeloscillations","w");
     lattice_size_x=100;
     deltax=100;          /* m */
     idum=-76;
     channelslope=0.01;   /* m/m */
     w0=100;              /* m */
     h0=1;                /* m */
     iup=ivector(1,lattice_size_x);
     idown=ivector(1,lattice_size_x);
     height=vector(1,lattice_size_x);
     heightold=vector(1,lattice_size_x);
     deltah=vector(1,lattice_size_x);
     width=vector(1,lattice_size_x);
     widthold=vector(1,lattice_size_x);
     deltaw=vector(1,lattice_size_x); 
     flux=vector(1,lattice_size_x);
     h2=vector(1,lattice_size_x);
     w2=vector(1,lattice_size_x);
     simulationlength=5000000.0; /* yr */ 
     outputinterval=100000.0;
     for (i=1;i<=lattice_size_x;i++)
      {width[i]=w0;
       widthold[i]=width[i];
       height[i]=-0.05*sin((i-1)*2*PI*3/(lattice_size_x-1))+0.03*gasdev(&idum);
       heightold[i]=height[i];
       flux[i]=0;h2[i]=0;w2[i]=0;
       deltaw[i]=0;deltah[i]=0;
       iup[i]=i+1;idown[i]=i-1;}
     iup[lattice_size_x]=1;idown[1]=lattice_size_x;
     timestep=1;
     elapsedtime=0.0;
     printed=0;
     while (elapsedtime<simulationlength)
      {printed=elapsedtime/outputinterval;
       maxdelw=0;
       for (i=1;i<=lattice_size_x;i++)
        flux[i]=((heightold[idown[i]]-heightold[i])/deltax+channelslope)/widthold[i];
       for (i=1;i<=lattice_size_x;i++)
        {del=0.5*timestep/deltax*(flux[idown[i]]-flux[i]);
         h2[i]=0.5*(heightold[idown[i]]+heightold[i])+del;
         w2[i]=widthold[i]-0.5*timestep*0.5*
          (flux[idown[i]]+flux[i]-2*channelslope/w0)/(h0-h2[i]);}
       for (i=1;i<=lattice_size_x;i++)
        flux[i]=((h2[i]-h2[iup[i]])/deltax+channelslope)/w2[i];
       for (i=1;i<=lattice_size_x;i++) 
        {del=timestep/deltax*(flux[i]-flux[iup[i]]);
         height[i]+=del;
         width[i]-=0.5*timestep*(flux[iup[i]]+flux[i]-
          2*channelslope/w0)/(h0-heightold[i]);
         if (fabs(widthold[idown[i]]-widthold[iup[i]])/widthold[i]>maxdelw)
          maxdelw=fabs(widthold[idown[i]]-widthold[iup[i]])/widthold[i];
         if (width[i]<10) width[i]=10;
         if (width[i]>300) width[i]=300;}
       elapsedtime+=timestep;
       if ((maxdelw>0.5*deltax/timestep)||(timestep>deltax*deltax))
         {elapsedtime-=timestep;
          timestep/=2.0;
          for (i=1;i<=lattice_size_x;i++)
           {height[i]=heightold[i];
            width[i]=widthold[i];}}
        else 
          if ((maxdelw<0.5*deltax*deltax/timestep/10)&&(timestep<deltax*deltax/10))
           timestep*=1.2;
       for (i=1;i<=lattice_size_x;i++)
         {heightold[i]=height[i];
          widthold[i]=width[i];}
       if ((int)(elapsedtime/outputinterval)>printed)
        {for (i=1;i<=lattice_size_x;i++)
          fprintf(fp1,"%d %f %f\n",i,height[i]/h0,width[i]/w0);
         printed=(int)(elapsedtime/outputinterval);}}
     fclose(fp1);
}  
