/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "sion.h"
#include "mpilinktest_mapping.h"
#include "mpilinktest_datastructures.h"
#include "ppmwrite.h"
#include "postscript_driver.h"
#include "pingponganalysis_tools.h"

#define FEHLERINDATENno
#define FILENAME_LENGTH 1024
#define MAXENTRIES 10000
#define MB *1024*1024

// int length_of_message = (128*1024);

double opt_minbw=1;
double opt_maxbw=10000;
double opt_mintime=-1e20;
double opt_maxtime=1e20;
int do_range=0;


int main(int argc, char **argv)
{
  FILE *fp, *outfp;
  char infilename[FILENAME_LENGTH];

  int i,j,rank,f_rank,blknum,t;
  int nullcount;
  sion_int64 chunksize=0;
  char *localbuffer;
  sion_int64 *chunksizes=NULL;
  int  *globalranks=NULL;

  /* options */
  int  verbose=0;
  sion_int64 opt_localsize=-1;
  sion_int32 opt_fsblocksize=-1;
  int opt_steps=500;
  int opt_doaccesspatternppm=0;
  int opt_doaccesspatternascii=0;
  int opt_dobwpatternppm=0;
  int opt_dognuplotfile=0;
  int opt_dopsreport=0;
  int opt_dobadlinkreport=0;

  /* for file infomation */
  int ntasks,numtasks,nfiles;
  sion_int32 fsblksize;
  int sid,outsid,size,blocks;
  sion_int64  globalskip;
  sion_int64  start_of_varheader;
  sion_int64 *sion_chunksizes;
  sion_int64 *sion_globalranks;
  sion_int64 *sion_blockcount;
  sion_int64 *sion_blocksizes;
  sion_int64 *siondump_sum_bytes_per_task;
  sion_int64  siondump_sum_bytes;
  sion_int64  siondump_filesize;
  char filename_ppm[256];
  char distfilename_ppm[256];
  char bwfilename_ppm[256];
  FILE *outfile_ppm,*outdistfile,*accessfile,*stepfile1,*stepfile2, *outbwfile_ppm, *gnuplotfile, *psreportfile;
  double opt_factor=1.0;
  char helpstr[256];
  
  int opt_dodistpattern=0;
  int minentriesX[MAXENTRIES];
  int minentriesY[MAXENTRIES];
  double minentriesVAL[MAXENTRIES];
  int minentriesfound=0;
  int maxentriesX[MAXENTRIES];
  int maxentriesY[MAXENTRIES];
  double maxentriesVAL[MAXENTRIES];
  int maxentriesfound=0;
  int maxtorus[MAPPING_DIM],maxdist=0;
  int distance[MAPPING_DIM];
  int dist_right[MAPPING_DIM];
  int dist_left[MAPPING_DIM];
  int collectpnumR,pos;
  int number_of_messages;
  int do_alltoall,do_serial,do_mix,opt_donodelist;
  int opt_mask=0;
  int opt_torus=0;
  int opt_oldformat=0;
  char opt_maskfile[FILENAME_LENGTH];
  double opt_mintimecol=-1.0,opt_maxtimecol=-1.0;
  int mapping_dim;

  /* analysis data */
  analysis_fields_t fields;

  char       location[MAX_LOCATION_LENGTH];
  sion_int32 mappingout[MAPPING_DIM];
  double     statdata_dp[10];
  sion_int32 statdata_int[10];
  sion_int32 version_int[3];
  char       linkteststr[5];

  /* derived analysis data */
  double  ptimings_min
	, ptimings_max
	, ptimings_avg
	, stepwidth
	, atimings_min
	, atimings_max
	, atimings_avg
	/** Variables inserted by kondil */
	, scaleTime_min
	, scaleTime_max;
  int     help, doswap;
  opt_donodelist=0;
  

  /* parse command line */
  i=1;
  if (argc < 2)
    usage(argv[0]);

  while( i < argc ) {
    if( argv[i][0] == '-' ) {
      switch( argv[i][1] ) {
      case 'a':
        opt_doaccesspatternppm=1;
        break;
      case 'A':
        opt_doaccesspatternascii=1;
        break;
      case 'b':
        opt_dobwpatternppm=1;
        break;
      case 'B':
        opt_dobadlinkreport=1;
        break;
      case 'd':
        opt_dodistpattern=1;
        break;
      case 'D':
        opt_donodelist=1;
        break;
      case 'f':
        opt_factor = atof(argv[++i]);
        break;
      case 'g':
        opt_dognuplotfile=1;
        break;
      case 'l':
        opt_minbw = atof(argv[++i]);
	do_range=1;
        break;
      case 'L':
        opt_maxbw = atof(argv[++i]);
	do_range=1;
        break;
      case 'm':
        opt_mask = 1;
	strcpy(opt_maskfile,argv[++i]);
        break;
      case 'O':
        opt_oldformat=1;
        break;
      case 'p':
        opt_dopsreport=1;
        break;
      case 'q':
        opt_fsblocksize = atoi(argv[++i]);
        break;
      case 'Q':
        opt_fsblocksize = atoi(argv[++i]) MB;
        break;
      case 's':
        opt_steps = atoi(argv[++i]);
        break;
      case 't':
        opt_mintimecol = atof(argv[++i]);
        break;
      case 'T':
        opt_maxtimecol = atof(argv[++i]);
        break;
      case 'v':
        verbose++;
        break;
      case 'V':
	fprintf(stderr,"FZJ Linktest version (pingponganalysis): %d.%dp%d\n",VERSION,VERSIONSUB,VERSIONPATCHLEVEL);
	exit(1);
      default:
        usage(argv[0]);
      }
    }
    i++;
  }


  /* STEP: open binary SION file */
  strcpy(infilename,argv[argc-1]);
  printf("pingponganalysis: infilename:                %-30s\n",infilename);
  chunksizes=NULL;globalranks=NULL; /* will be allocated by sion_open */
  sid=sion_open(infilename,"br,ansi",&ntasks,&nfiles,&chunksizes,&fsblksize,&globalranks,&fp);

  /* STEP: print header information */
  printf("pingponganalysis: sid:                       %d\n",sid);
  printf("pingponganalysis: filename:                  %-30s\n",infilename);
  printf("pingponganalysis: fsblksize:                 %lu bytes (%6.2f MB)\n",(unsigned long)fsblksize,fsblksize/1024.0/1024.0);
  printf("pingponganalysis: current endianness:        %s\n",(sion_get_endianness())?"big":"little");
  printf("pingponganalysis: file endianness:           %s\n", (sion_get_file_endianness(sid)) ? "big" : "little");
  
  doswap=sion_get_file_endianness(sid)!=sion_get_endianness();
  printf("pingponganalysis: big/little  swapping:      --> %d\n", doswap);
  
  /* STEP: get more information from SION file and print it */
  sion_get_locations(sid,&size,&blocks,&globalskip,&start_of_varheader,&sion_chunksizes,
		     &sion_globalranks,&sion_blockcount,&sion_blocksizes);
  printf("pingponganalysis: start_of_varheader:        %lld bytes  (%6.2f MB)\n",start_of_varheader,start_of_varheader/1024.0/1024.0);
  printf("pingponganalysis: max number of blocks:      %d\n",blocks);
  printf("pingponganalysis: ntasks:                    %d\n",ntasks);
  if(blocks>1) {
    fprintf(stderr, "wrong sion file max #blocks > 1 (%d)\n", blocks);
    return(1);
  } 
  
  /* STEP: read general data of measurement */
  printf("pingponganalysis: id:                      %s %d sid=%d fp=%x\n",LINKTEST_ID,strlen(LINKTEST_ID),sid,fp);
  sion_seek(sid,0,0,0);

  if(!opt_oldformat) {
    freaddata(linkteststr,(long) strlen(LINKTEST_ID),fp);
    freaddata(version_int,(long) 3 * sizeof(sion_int32),fp);
    sion_swap(version_int,version_int,sizeof(sion_int32),3,doswap);
    printf("pingponganalysis: data version:            %d.%dp%d\n",version_int[0],version_int[1],version_int[2]);
    
    if(version_int[1]>=4) {
      mapping_dim=6;
    } else {
      mapping_dim=4;
    }
  } else {
    mapping_dim=4;
    version_int[0]=1;
    version_int[1]=2;
    version_int[2]=0;
  }

  freaddata(location,(long) MAX_LOCATION_LENGTH,fp);

  for(i=0;i<MAPPING_DIM;i++) mappingout[i]=0;
  freaddata(&mappingout[0],(long) mapping_dim * sizeof(sion_int32),fp);
  sion_swap(mappingout,mappingout,sizeof(sion_int32),mapping_dim,doswap);

  freaddata(&statdata_int[0],(long) 10 * sizeof(sion_int32),fp);
  sion_swap(statdata_int,statdata_int,sizeof(sion_int32),10,doswap);

  freaddata(&statdata_dp[0],(long) 10 * sizeof(double),fp);
  sion_swap(&statdata_dp[0],&statdata_dp[0],sizeof(double),10,doswap);

  
  /* STEP: set internal variables */
  number_of_messages=statdata_int[0];  length_of_message=statdata_int[1];
  do_alltoall=statdata_int[2];         collectpnumR=statdata_int[3];
  do_serial=statdata_int[4];           do_mix=statdata_int[5];
  if(collectpnumR>1) {
    printf("WARNING: collectpnum>1 (%d) not supported by this version, existing ...\n",collectpnumR);
    exit(0);
  }
  ptimings_min=statdata_dp[0];
  ptimings_max=statdata_dp[1];  
  ptimings_avg=statdata_dp[2];

  // kondil
  if(opt_mintimecol>0)
    scaleTime_min = opt_mintimecol;
  else
    scaleTime_min = ptimings_min;
  
  if(opt_maxtimecol>0) 
    scaleTime_max = opt_maxtimecol;
  else
    scaleTime_max = ptimings_max;
  
//   if(opt_mintimecol>0) ptimings_min=opt_mintimecol;
//   if(opt_maxtimecol>0) ptimings_max=opt_maxtimecol;

  printf("pingponganalysis: opt_min max timecol:       %11.2f us %11.2f us\n",opt_mintimecol*1000*1000.0, opt_maxtimecol*1000*1000.0);

    // kondil
  stepwidth=(scaleTime_max - scaleTime_min)/opt_steps;
//   stepwidth=(ptimings_max-ptimings_min)/opt_steps;

  printf("pingponganalysis: histogram:                 %11.2f us  to %11.2f us, stepwidth=%11.2f us, #steps=%d\n",
	 ptimings_min*1000*1000.0, ptimings_max*1000*1000.0,
	 stepwidth*1000*1000.0,opt_steps);
  
  if(do_alltoall) {
    atimings_min=statdata_dp[3];
    atimings_max=statdata_dp[4];
    atimings_avg=statdata_dp[5];
  }

  /* STEP: calculate min max limits for scan */
  opt_maxtime=mbs2time(opt_minbw);
  opt_mintime=mbs2time(opt_maxbw);
  printf("pingponganalysis: reporting link outside range:       %11.2f us   < t < %11.2f us\n",opt_mintime*1000*1000.0, opt_maxtime*1000*1000.0);
  printf("pingponganalysis: reporting link outside range:       %11.2f MB/s < t < %11.2f MB/s\n",opt_minbw, opt_maxbw);

  printf("pingponganalysis: ->stat_data_int:              %10d %10d %10d %10d %10d\n",
         statdata_int[0],statdata_int[1],statdata_int[2],statdata_int[3],statdata_int[4]);
  printf("pingponganalysis: ->statdata_dp:               %e %e %e %e %e\n",
         statdata_dp[0],statdata_dp[1],statdata_dp[2],statdata_dp[3],statdata_dp[4]);
  printf("pingponganalysis: length_of_message:              %d\n",length_of_message);
  printf("pingponganalysis: collectpnumR:                   %d\n",collectpnumR);
  printf("pingponganalysis: alltoall:                       %d\n",do_alltoall);
  numtasks=ntasks;
  printf("pingponganalysis: matsize:                        %d\n",numtasks);
  printf("pingponganalysis: opt_mask:                       %d\n",opt_mask);
  if(opt_mask) {
    printf("pingponganalysis: opt_maskfile:                   %s\n",opt_maskfile);
  }

  /* STEP: allocate fields for reading data of one task */
  allocate_fields(&fields, numtasks, opt_steps);

  if(opt_mask) {
    read_maskfile(&fields,numtasks,opt_maskfile);
  }

  /* STEP: init output files */
  ppminitsmooth(1); 
  sprintf(filename_ppm,"accesspattern.ppm");
  sprintf(distfilename_ppm,"distpattern.ppm");
  sprintf(bwfilename_ppm,"bwpattern.ppm");
  if(opt_doaccesspatternppm) {
    outfile_ppm=ppmopen(numtasks,numtasks,0.0,1.0*numtasks,filename_ppm);
  }
  if(opt_doaccesspatternascii) {
    accessfile=fopen("accesspattern.dat","w");
  }
  if(opt_dodistpattern) {
    outdistfile=fopen("distancepattern.dat","w");
  }

  if(opt_dobwpatternppm) {
    outbwfile_ppm=ppmopen(numtasks,numtasks,0.0,1.0*numtasks,bwfilename_ppm);
  }
  if(opt_dognuplotfile) {
    gnuplotfile=fopen("gnuplotinput.dat","w");
  }
  if(opt_dopsreport) {
    int offy;
#define FIRSTROW 5
#define SECONDROW 300
#define FONTSIZE 9

    psreportfile=fopen("report_fzjlinktest.ps","w");
    PSsetmsglen(length_of_message);
    PSheader(psreportfile);
    textr(psreportfile,infilename,PSTITELX,PSTITELY,14);
    rect(psreportfile,PSDESCX,PSDESCY,PSDESCW,PSDESCH,0,0.1);
    offy=PSDESCY+PSDESCH;
    sprintf(helpstr,"length_of_message:   %d bytes (%8.2f KBytes)",length_of_message,length_of_message/1024.0);
    text(psreportfile,helpstr,PSDESCX+FIRSTROW,offy-11,FONTSIZE);
    sprintf(helpstr,"number_of_tasks:     %d ",numtasks);
    text(psreportfile,helpstr,PSDESCX+SECONDROW,offy-11,FONTSIZE);
    offy-=12;
    sprintf(helpstr,"number_of_messages:  %d ",number_of_messages);
    text(psreportfile,helpstr,PSDESCX+FIRSTROW,offy-11,FONTSIZE);
    sprintf(helpstr,"Execution order:     %s ",(do_serial)?"Serial":"Parallel");
    text(psreportfile,helpstr,PSDESCX+SECONDROW,offy-11,FONTSIZE);
    offy-=12;
    sprintf(helpstr,"Alltoall:            %d ",do_alltoall);
    text(psreportfile,helpstr,PSDESCX+FIRSTROW,offy-11,FONTSIZE);
    sprintf(helpstr,"Mixing PE rank:     %s ",(do_mix)?"Yes":"No");
    text(psreportfile,helpstr,PSDESCX+SECONDROW,offy-11,FONTSIZE);
    offy-=12;
    sprintf(helpstr,"Min Value:           %6.1fus (%7.2f MB/s)",ptimings_min*1000*1000,time2mbs(ptimings_min));
    text(psreportfile,helpstr,PSDESCX+FIRSTROW,offy-11,FONTSIZE);
    if(do_alltoall) {
      sprintf(helpstr,"Alltoall Min Value:  %8.1fus (1 Byte)",atimings_min*1000*1000);
      text(psreportfile,helpstr,PSDESCX+SECONDROW,offy-11,FONTSIZE);
    }
    offy-=12;
    sprintf(helpstr,"Max Value:           %6.1fus (%7.2f MB/s)",ptimings_max*1000*1000,time2mbs(ptimings_max));
    text(psreportfile,helpstr,PSDESCX+FIRSTROW,offy-11,FONTSIZE);
    if(do_alltoall) {
      sprintf(helpstr,"Alltoall Max Value:  %8.1fus (1 Byte)",atimings_max*1000*1000);
      text(psreportfile,helpstr,PSDESCX+SECONDROW,offy-11,FONTSIZE);
    }
    offy-=12;
    sprintf(helpstr,"Avg Value:           %6.1fus (%7.2f MB/s)",ptimings_avg*1000*1000,time2mbs(ptimings_avg));
    text(psreportfile,helpstr,PSDESCX+FIRSTROW,offy-11,FONTSIZE);
    if(do_alltoall) {
      sprintf(helpstr,"Alltoall Avg Value:  %8.1fus (1 Byte)",atimings_avg*1000*1000);
      text(psreportfile,helpstr,PSDESCX+SECONDROW,offy-11,FONTSIZE);
    }
    
  }

  /* STEP: scan all data for mappingdata -> calculate size of torus */
  for(i=0;i<MAPPING_DIM;i++) maxtorus[i]=0;
  for(rank=0;rank<size;rank++) {
    sion_seek(sid,rank,0,0);
    if(!opt_oldformat) {
      freaddata(linkteststr,(long) strlen(LINKTEST_ID),fp);
      freaddata(version_int,(long) 3 * sizeof(sion_int32),fp);
      sion_swap(version_int,version_int,sizeof(sion_int32),3,doswap);
    }
    freaddata(location,(long) MAX_LOCATION_LENGTH,fp);
    strcpy(fields.alllocation+rank*MAX_LOCATION_LENGTH,location);
    freaddata(&mappingout[0],(long) mapping_dim * sizeof(sion_int32),fp);
    sion_swap(mappingout,mappingout,sizeof(sion_int32),mapping_dim,doswap);

    if((rank)%32==0) {
      printf("pingponganalysis: ->rank: %3d, -> location=>%s< %dx%dx%dx%d pos=%d\n",
	     rank,fields.alllocation+rank*MAX_LOCATION_LENGTH,
	     mappingout[0],mappingout[1],mappingout[2],mappingout[3],pos);
    }
      
    for(i=0;i<MAPPING_DIM;i++) {
      fields.mapping[rank*MAPPING_DIM+i]=mappingout[i];
      if (mappingout[i]>maxtorus[i]) maxtorus[i]=mappingout[i];
    }
  }
  maxdist=0;
  printf("pingponganalysis: -> size:");
  for(i=0;i<MAPPING_DIM;i++) {
    maxtorus[i]++;
    maxdist+=(int) (maxtorus[i]/2.0);
    printf(" %d",maxtorus[i]);
  }
  printf(" maxdist=%d\n",maxdist);

  
  if(opt_dopsreport) {
    PSplotscala (psreportfile, PSSCALENELEM, scaleTime_min, scaleTime_max);
//     PSplotscala (psreportfile, PSSCALENELEM, ptimings_min,ptimings_max);

    /* PSplotrowimage_start(psreportfile,fields.ptimings,0,numtasks,ptimings_min,ptimings_max,fields.alllocation+rank*MAX_LOCATION_LENGTH,  */
    /*  			 do_range,opt_mintime,opt_maxtime);  */
    
  }

  /* STEP: read data for each rank and process data */
  for(rank=0;rank<size;rank++) {
   
    sion_seek(sid,rank,0,0);
    
    /* read data from block */
    freaddata(linkteststr,(long) strlen(LINKTEST_ID),fp);

    freaddata(version_int,(long) 3 * sizeof(sion_int32),fp);
    sion_swap(version_int,version_int,sizeof(sion_int32),3,doswap);
    
    freaddata(location,(long) MAX_LOCATION_LENGTH,fp);
    strcpy(fields.alllocation+rank*MAX_LOCATION_LENGTH,location);
    freaddata(&mappingout[0],(long) mapping_dim * sizeof(sion_int32),fp);
    sion_swap(mappingout,mappingout,sizeof(sion_int32),mapping_dim,doswap);
    if(rank%32==0) {
      printf("pingponganalysis: ->rank:                    %d \n",rank);
    }
    
    /* read parameter and statistic data */
    if(rank==0) {
      freaddata(&statdata_int[0],(long) 10 * sizeof(sion_int32),fp);
      sion_swap(statdata_int,statdata_int,sizeof(sion_int32),10,doswap);
      freaddata(&statdata_dp[0],(long) 10 * sizeof(double),fp);
      sion_swap(&statdata_dp[0],&statdata_dp[0],sizeof(double),10,doswap);
    }
    
    /* timings, first run  */
    freaddata(&fields.ptimings[0],(long) numtasks * sizeof(double),fp);
    sion_swap(&fields.ptimings[0],&fields.ptimings[0],sizeof(double),numtasks,doswap);

    /*
      for(j=0;j<size;j++) {
      printf("WF: %d   %2d->%10.8f\n",rank,j,fields.ptimings[j]);
      }*/
    
    /* analysis of ptimings */
    for(t=0;t<numtasks;t++) {
      
      if(opt_mask) {
	if(!fields.mask[rank*numtasks+t]) fields.ptimings[t]=-1;
      }
      
      if(opt_dognuplotfile) {
	if(fields.ptimings[t]>=0) {
	  fprintf(gnuplotfile,"%10d %10d %18.8f %18.8f\n",rank,t,fields.ptimings[t],time2mbs(fields.ptimings[t]));
	} else {
	  fprintf(gnuplotfile,"%10d %10d %18.8f %18.8f\n",rank,t,0.0,0.0);
	}
      }
      
      if(rank != t ) {
	if(fields.ptimings[t]<fields.timings_min) fields.timings_min=fields.ptimings[t];
	if(fields.ptimings[t]>fields.timings_max) fields.timings_max=fields.ptimings[t];
        
	fields.nodeslist_fromavg[rank]+=fields.ptimings[t];
	fields.nodeslist_toavg[t]+=fields.ptimings[t];
	if(fields.ptimings[t]>fields.nodeslist_frommax[rank]) fields.nodeslist_frommax[rank]=fields.ptimings[t];
	if(fields.ptimings[t]>fields.nodeslist_tomax[t])      fields.nodeslist_tomax[t]=     fields.ptimings[t];
	  
	if(time2mbs(fields.ptimings[t])<1000)          fields.nodeslist_fromcntmax[rank]++;
	if(time2mbs(fields.ptimings[t])<1000)          fields.nodeslist_tocntmax[t]++;

	if((do_range) && (fields.ptimings[t]<opt_mintime) && (minentriesfound<MAXENTRIES)) {
	  minentriesX[minentriesfound]=rank;
	  minentriesY[minentriesfound]=t;
	  minentriesVAL[minentriesfound]=fields.ptimings[t];
	  /*           printf("[%d: %dx%d %e]\n",minentriesfound,(rank),t,ptimings[t]); */
	  minentriesfound++;
	}

	if((do_range) && (fields.ptimings[t]>opt_maxtime) && (maxentriesfound<MAXENTRIES)) {
	  maxentriesX[maxentriesfound]=rank;
	  maxentriesY[maxentriesfound]=t;
	  maxentriesVAL[maxentriesfound]=fields.ptimings[t];
	  maxentriesfound++;
	}
	  
	if(fields.ptimings[t]>=0) {
	  help= (int) ((fields.ptimings[t]-ptimings_min)/stepwidth);
	  if(help>=opt_steps) {
	    /* 	    printf("r=%d t=%d %d>%d val=%e stepwidth=%e\n",(rank),t,help,opt_steps,ptimings[t],stepwidth);  */
	    help=opt_steps-1;
	  }
	  if(help<0) {
	    printf("r=%d t=%d %d<0 val=%e stepwidth=%e\n",rank,t,help,fields.ptimings[t],stepwidth);
	    help=0;
	  }
	  /* 	  printf("r=%d t=%4d help=%4d val=%10.4f(%7.1fMB/s) stepwidth=%10.4f(%7.1fMB/s)\n", */
	  /* 		 rank,t,help,ptimings[t]*1000*1000,time2mbs(ptimings[t]),stepwidth*1000*1000,time2mbs(stepwidth)); */
	  fields.ptimings_stepcnt[help]++;
	}
      }
    }

    if(opt_dopsreport) {
      /* PSplotrow(psreportfile,fields.ptimings,rank,numtasks,ptimings_min,ptimings_max,fields.alllocation+rank*MAX_LOCATION_LENGTH, */
      /* 		do_range,opt_mintime,opt_maxtime); */
      
      /* PSplotrow(psreportfile,fields.ptimings,rank,numtasks,scaleTime_min,scaleTime_max,fields.alllocation+rank*MAX_LOCATION_LENGTH, */
      /* 		do_range,opt_mintime,opt_maxtime); */
      
      PSplotrowimage(psreportfile,fields.ptimings,rank,numtasks,ptimings_min,ptimings_max,fields.alllocation+rank*MAX_LOCATION_LENGTH,  
		     do_range,opt_mintime,opt_maxtime);  
      
		     
		     
		     
		     
      /* PSplotrowimage_row(psreportfile,fields.ptimings,rank,numtasks,ptimings_min,ptimings_max,fields.alllocation+rank*MAX_LOCATION_LENGTH,  */
      /* 			 do_range,opt_mintime,opt_maxtime);  */
    }
    
    
    /* scan accesspattern */
    freaddata(&fields.accesspattern[0],(long) numtasks * sizeof(sion_int32),fp);
    sion_swap(&fields.accesspattern[0],&fields.accesspattern[0],sizeof(sion_int32),numtasks,doswap);
    if(opt_doaccesspatternascii) {
      for(t=0;t<numtasks;t++) {
	fprintf(accessfile,"%8d %8d %8d\n",rank,t,fields.accesspattern[t]);  
      }
    }
    if(opt_doaccesspatternppm) {
      for(t=0;t<numtasks;t++) {
	if(fields.accesspattern[t]<0) {
	  printf("pingponganalysis: accesspattern[%d]=%d < 0, rank=%d\n",t,fields.accesspattern[t],rank);
	  fields.accesspattern[t]=0;
	}
	fields.accesspatterndp[t]=fields.accesspattern[t];
      }
      ppmwriteblock(outfile_ppm,fields.accesspatterndp,numtasks,1,0.0,1.0*numtasks);
    }
    
    /* analysis of ptimings in respect to distance */
    if(opt_dodistpattern) {
      int maxdim=6;
      for(i=MAPPING_DIM-1;i>=0;i--) {
	if(maxtorus[i]==1) maxdim--;
      }
      /* BG/P */
      if ((maxdim==4) && (maxtorus[3]==4)) {
	maxdim=3; 
      }
      
      for(t=0;t<numtasks;t++) {
	fprintf(outdistfile, "%4d %4d ",rank,t);
	  
	fields.distance[t]=0;
	for(i=0;i<maxdim;i++) {
	  dist_right[i]=abs(fields.mapping[rank*MAPPING_DIM+i]-fields.mapping[t*MAPPING_DIM+i]);
	  dist_left[i] =maxtorus[i]-abs(fields.mapping[rank*MAPPING_DIM+i]-fields.mapping[t*MAPPING_DIM+i]);
	  if(maxtorus[i]>1) distance[i]=(dist_right[i]<dist_left[i])?dist_right[i]:dist_left[i];
	  else              distance[i]=0;
	  fields.distance[t]+=distance[i];
	}
	
	for(i=0;i<MAPPING_DIM;i++) fprintf(outdistfile, "%2d ",fields.mapping[rank*MAPPING_DIM+i]);
	fprintf(outdistfile, "  ");
	for(i=0;i<MAPPING_DIM;i++) fprintf(outdistfile, "%2d ",fields.mapping[t*MAPPING_DIM+i]);
	fprintf(outdistfile, "    ");
	for(i=0;i<MAPPING_DIM;i++) fprintf(outdistfile, "%2d ",distance[i]);
	fprintf(outdistfile, "  ");
	fprintf(outdistfile, "%2d %18.8e s\n",fields.distance[t],fields.ptimings[t]);
	
      }
    }
    
    /* Blue Gene/P torus report  */
    if(opt_torus) {
      for(t=0;t<numtasks;t++) {
      }
    }
    
    if(opt_dobwpatternppm) {
      for(t=0;t<numtasks;t++) {
	if(fields.ptimings[t]<0) fields.ptimings[t]=0.0;
	/*           if(distance[t]>0) { */
	/*             ptimings[t]*=(double) opt_factor * 1.0/distance[t]; */
	/*           } */
      }
      ppmwriteblock(outbwfile_ppm,fields.ptimings,numtasks,1,0.0,opt_factor*fields.timings_max);
    }
    
    /* all to all timings  */
    if(do_alltoall==1) {
      freaddata(&fields.atimings[0],(long) numtasks * sizeof(double),fp);
      sion_swap(&fields.atimings[0],&fields.atimings[0],sizeof(double),numtasks,doswap);
    } 
    
    /* single link tests, only store on task 0  */
    if(rank==0) {
      freaddata(&fields.stimings1[0],(long) numtasks * sizeof(double),fp);
      sion_swap(&fields.stimings1[0],&fields.stimings1[0],sizeof(double),numtasks,doswap);
      
      freaddata(&fields.stimings2[0],(long) numtasks * sizeof(double),fp);
      sion_swap(&fields.stimings2[0],&fields.stimings2[0],sizeof(double),numtasks,doswap);
	
      freaddata(&fields.fromlist[0],(long) numtasks * sizeof(sion_int32),fp);
      sion_swap(&fields.fromlist[0],&fields.fromlist[0],sizeof(sion_int32),numtasks,doswap);
	
      freaddata(&fields.tolist[0],(long) numtasks * sizeof(sion_int32),fp);
      sion_swap(&fields.tolist[0],&fields.tolist[0],sizeof(sion_int32),numtasks,doswap);

#ifdef FEHLERINDATEN
      freaddata(&fromlist[0],(long) numtasks * sizeof(sion_int32),fp);
      freaddata(&tolist[0],(long) numtasks * sizeof(sion_int32),fp);
#endif
    } 

  } /* rank */

  
    if(opt_dopsreport) {
      
      /* PSplotrowimage(psreportfile,fields.ptimings,rank,numtasks,scaleTime_min,scaleTime_max,fields.alllocation+rank*MAX_LOCATION_LENGTH,  */
      /*  		     do_range,opt_mintime,opt_maxtime);  */
      
//  PSplotrowimage(psreportfile,fields.ptimings,rank,numtasks,ptimings_min,ptimings_max,fields.alllocation+rank*MAX_LOCATION_LENGTH, 
//        		     do_range,opt_mintime,opt_maxtime); 
    }
  


  /* STEP: report timings outside valid range (min,max) */
  printf("pingponganalysis: timings_min:          %e s %10.6f MB/s\n",fields.timings_min,time2mbs(fields.timings_min));
  printf("pingponganalysis: timings_max:          %e s %10.6f MB/s\n",fields.timings_max,time2mbs(fields.timings_max));
  for(t=0;t<minentriesfound;t++) {
    printf("pingponganalysis:  suspected link #%05d: %6d <-> %6d: %e s %16.10f MB/s (%s <-> %s)\n",t,
           minentriesX[t],minentriesY[t],minentriesVAL[t],time2mbs(minentriesVAL[t]),
           fields.alllocation+minentriesX[t]*MAX_LOCATION_LENGTH,fields.alllocation+minentriesY[t]*MAX_LOCATION_LENGTH);
    
  }

  for(t=0;t<maxentriesfound;t++) {
    printf("pingponganalysis:  FAST      link #%05d: %6d <-> %6d: %e s %16.10f MB/s (%s <-> %s)\n",t,
           maxentriesX[t],maxentriesY[t],maxentriesVAL[t],time2mbs(maxentriesVAL[t]),
           fields.alllocation+maxentriesX[t]*MAX_LOCATION_LENGTH,fields.alllocation+maxentriesY[t]*MAX_LOCATION_LENGTH);
    
  }

  /* STEP: write some histograms */
  stepfile1=fopen("distribution.dat","w");
  for(t=0;t<opt_steps;t++) fprintf(stepfile1,"%5d %10e %10ld\n",t,t*stepwidth+ptimings_min,
                                   fields.ptimings_stepcnt[t] ) ;
  fclose(stepfile1);
  printf(" wrote %s\n","distribution.dat");

  stepfile1=fopen("mapping.dat","w");
  for(t=0;t<numtasks;t++) {
    fprintf(stepfile1,"%6d",t);
    for(i=0;i<MAPPING_DIM;i++) {
      fprintf(stepfile1," %4d",fields.mapping[t*MAPPING_DIM+i]);
    }
    fprintf(stepfile1," %s\n",fields.alllocation+t*MAX_LOCATION_LENGTH);

  }
  close(stepfile1);
  printf(" wrote %s\n","mapping.dat");

  stepfile2=fopen("distribution_bw.dat","w");
  for(t=0;t<opt_steps;t++) fprintf(stepfile2,"%5d %10e %10ld\n",t,time2mbs(t*stepwidth+ptimings_min),
                                   fields.ptimings_stepcnt[t] ) ;
  fclose(stepfile2);
  printf(" wrote %s\n","distribution_bw.dat");

  if(opt_dopsreport) {
    PSplothist(psreportfile,fields.ptimings_stepcnt,opt_steps,scaleTime_min,scaleTime_max,stepwidth);
    PSplothist(psreportfile,fields.ptimings_stepcnt,opt_steps,ptimings_min,ptimings_max,stepwidth);
//     PSplothist(psreportfile,fields.ptimings_stepcnt,opt_steps,ptimings_min,ptimings_max,stepwidth);
  }

  if(opt_dobadlinkreport) {
    stepfile1=fopen("bad_links.dat","w");
    fprintf(stepfile1,"# %5s %5s %5s %16s %16s %16s %16s\n","Nr.","From","To","1st (s)","1st (MB/s)","2nd (s)","2nd (MB/s)");
    for(t=0;t<numtasks;t++) {
      fprintf(stepfile1," %5d %5d %5d %16.8e %16.10f %16.8e %16.10f %s -> %s\n",t,
	      fields.fromlist[t],fields.tolist[t],fields.stimings2[t],time2mbs(fields.stimings2[t]),fields.stimings1[t],time2mbs(fields.stimings1[t]),
	      fields.alllocation+fields.fromlist[t]*MAX_LOCATION_LENGTH,fields.alllocation+fields.tolist[t]*MAX_LOCATION_LENGTH);
    }
    fclose(stepfile1);
    printf(" wrote %s\n","bad_links.dat");
  }

  if(opt_donodelist) {
    stepfile1=fopen("nodelist.dat","w");
    fprintf(stepfile1,"# %5s %15s %25s %25s %25s %25s %10s %10s\n","Nr.","Host",
	    "from_avg_bw (MB/s)","to_avg_bw (MB/s)",
	    "from_min_bw (MB/s)","to_min_bw (MB/s)",
	    "fromcntmax<2000","tocntmax<2000");
    for(t=0;t<numtasks;t++) {
      fprintf(stepfile1," %5d %15s %25.10f %25.10f %25.10f %25.10f %10d %10d\n",t,
	      fields.alllocation+t*MAX_LOCATION_LENGTH,
	      time2mbs(fields.nodeslist_fromavg[t]/(numtasks-1.0)),
	      time2mbs(fields.nodeslist_toavg[t]/(numtasks-1.0)),
	      time2mbs(fields.nodeslist_frommax[t]),
	      time2mbs(fields.nodeslist_tomax[t]),
	      fields.nodeslist_fromcntmax[t],
	      fields.nodeslist_tocntmax[t]
	      );
    }
    fclose(stepfile1);
    printf(" wrote %s\n","nodelist.dat");
  }

  /* STEP: close all files */
  sion_close(sid);
  printf(" after sion_close\n");
  if(opt_doaccesspatternppm) {
    ppmclose(outfile_ppm);
    printf(" after close outfile_ppm\n");
  }
  if(opt_dodistpattern) {
    fclose(outdistfile);
    printf(" after close outdistfile\n");
  }
  if(opt_dobwpatternppm) {
    ppmclose(outbwfile_ppm);
    printf(" after close outdistfile_ppm\n");
  }
  if(opt_doaccesspatternascii) {
    fclose(accessfile);
    printf(" after close accessfile\n");
  }
  if(opt_dognuplotfile) {
    fclose(gnuplotfile);
    printf(" after close gnuplotfile\n");
  }
  if(opt_dopsreport) {
    PSfooter(psreportfile);
    fclose(psreportfile);
    printf(" after close psreportfile\n");
  }
  
  return(0);
}


