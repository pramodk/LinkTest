/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <string.h>

#ifdef LINKTEST_BGP
/* #include <spi/kernel_interface.h> */
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#endif

#ifdef LINKTEST_BGQ
#include <firmware/include/personality.h>
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#include <hwi/include/common/uci.h>
#endif

#include <assert.h>
#include <stdlib.h>

#include "mpilinktest_mapping.h"
#include "mpilinktest_datastructures.h"


int write_mapping (char *filename) {
  int num_procs, my_rank;
  int i,j;
  int *pos;
  int mypos[MAPPING_DIM];
  double a;

#ifdef LINKTEST_BGP
  char location[BGPPERSONALITY_MAX_LOCATION];
  FILE *outfile;
  _BGP_Personality_t personality;

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if(my_rank==10) {
    pos=(int *) malloc(num_procs*MAPPING_DIM*sizeof(int));
    for(i=0;i<MAPPING_DIM;i++) pos[i]=0;
    if(!pos) {
      fprintf(stderr,"cannot allocate memory for mapping data\n");
      return(0);
    }
  } else {
    pos=NULL;
  }
  
  /*
   * BlueGene runtime: get and print BG personality
   */

  Kernel_GetPersonality(&personality, sizeof(personality));
  BGP_Personality_getLocationString(&personality, location);
  
  for(i=0;i<MAPPING_DIM;i++) mypos[i]=0;
  mypos[0] = personality.Network_Config.Xcoord;
  mypos[1] = personality.Network_Config.Ycoord;
  mypos[2] = personality.Network_Config.Zcoord;
  mypos[3] = Kernel_PhysicalProcessorID();

  MPI_Gather(mypos,MAPPING_DIM,MPI_INT,pos,MAPPING_DIM,MPI_INT,10,MPI_COMM_WORLD);

  if(my_rank==10) {
      outfile=fopen(filename,"w");
      if(outfile) {
	  fprintf(outfile,"BGL-Mapping: type=MPI_COMM_WORLD\n");
	  fprintf(outfile,"MPI-SIZE:  %d\n",num_procs);
	  fprintf(outfile,"VN-MODE:   %d\n",personality.Kernel_Config.ProcessConfig);
	  fprintf(outfile,"PART-SIZE: %d %d %d %d\n",	 personality.Network_Config.Xnodes,
		  personality.Network_Config.Ynodes,
		  personality.Network_Config.Znodes,
		  personality.Kernel_Config.ProcessConfig+1);
	  for(i=0;i<num_procs;i++) {
 	    fprintf(outfile,"%04d\n",i);
	    for(j=0;j<MAPPING_DIM;j++) {
	      fprintf(outfile," %4d",pos[i*MAPPING_DIM+j]);
	    }
 	    fprintf(outfile,"\n");
	  }
	  fclose(outfile);
      } else {
	  fprintf(stderr,"cannot open output file for mapping data: %s\n",filename);
	  return(0);
      }
      
  }

#endif
  return(1);
}

int  get_mapping (char **locationptr, int **mappingptr) {
  int num_procs, my_rank;
  int i,j;
  int *pos;
  int mypos[MAPPING_DIM];
  double a;
#ifdef LINKTEST_BGP
  char location[BGPPERSONALITY_MAX_LOCATION];
  _BGP_Personality_t personality;
/* #NODIST# */
#elif defined(LINKTEST_BGQ)
  char location[MAX_LOCATION_LENGTH];
  Personality_t personality;
/* #NODIST# */
#else
  char location[MAX_LOCATION_LENGTH];
#endif
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  *mappingptr=(int *) malloc(MAPPING_DIM*sizeof(int));
  if(!*mappingptr) {
    fprintf(stderr,"cannot allocate memory for mapping data\n");
    return(0);
  }

  *locationptr=(char *) malloc(MAX_LOCATION_LENGTH);
  if(!*locationptr) {
    fprintf(stderr,"cannot allocate memory for mapping data\n");
    return(0);
  }
  
  /*
   * BlueGene runtime: get and print BG personality
   */
 
#ifdef LINKTEST_BGP
  Kernel_GetPersonality(&personality, sizeof(personality));
  BGP_Personality_getLocationString(&personality, location);
  
  for(i=0;i<MAPPING_DIM;i++) (*mappingptr)[i]=0;
  (*mappingptr)[0] = personality.Network_Config.Xcoord;
  (*mappingptr)[1] = personality.Network_Config.Ycoord;
  (*mappingptr)[2] = personality.Network_Config.Zcoord;
  (*mappingptr)[3] = Kernel_PhysicalProcessorID();
#elif defined(LINKTEST_BGQ)
   Kernel_GetPersonality(&personality, sizeof(personality));
  (*mappingptr)[0] = personality.Network_Config.Acoord;
  (*mappingptr)[1] = personality.Network_Config.Bcoord;
  (*mappingptr)[2] = personality.Network_Config.Ccoord;
  (*mappingptr)[3] = personality.Network_Config.Dcoord;
  (*mappingptr)[4] = personality.Network_Config.Ecoord;
  (*mappingptr)[5] = Kernel_ProcessorID();
  BG_UniversalComponentIdentifier uci = personality.Kernel_Config.UCI;
  unsigned int row, col, mp, nb, cc;
  bg_decodeComputeCardOnNodeBoardUCI(uci, &row, &col, &mp, &nb, &cc);
  sprintf(location, "R%x%x-M%d-N%02x-J%02x <%d,%d,%d,%d,%d,%d>", row, col, mp, nb, cc,
          (*mappingptr)[0], (*mappingptr)[1], (*mappingptr)[2],
          (*mappingptr)[3], (*mappingptr)[4], (*mappingptr)[5]);
#else
  for(i=0;i<MAPPING_DIM;i++) (*mappingptr)[i]=0;
  (*mappingptr)[0] = my_rank;
  gethostname(location,256);
#endif
  strncpy(*locationptr,location,MAX_LOCATION_LENGTH-1); 
  (*locationptr)[MAX_LOCATION_LENGTH-1]='\0'; 
/*   printf("WF: in get_mapping rank=%d >%s< max_location_length=%d bgppersonality_max_location=%d\n",my_rank,*locationptr,MAX_LOCATION_LENGTH,BGPPERSONALITY_MAX_LOCATION); */


  return(1);
}

    
