#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>

#include "ppmwrite.h"

#define SHADES 256
#define MAXCOL2 (SHADES*SHADES*SHADES)
#define MAXCOL (7*SHADES)

#define BUFLEN 2048

static int color_r[MAXCOL];
static int color_g[MAXCOL];
static int color_b[MAXCOL];

static double gmaxval=0;
static double gminval=0;



int ppminit( int d ) {
  int ix,iy,iz;
  int count = d==1 ? 0 : MAXCOL-1;
  for (ix=0; ix<SHADES; ++ix) {
    for (iy=0; iy<SHADES; ++iy) {
      for (iz=0; iz<SHADES; ++iz) {
        color_r[count] = ix;
        color_g[count] = iy;
        color_b[count] = iz;
        count += d;
      }
    }
  }
  return(1);
}

/*
 * SmoothColorTable
 */

int ppminitsmooth( int d) {
  int S = SHADES - 1;
  int i;

  int count = d==1 ? 0 : MAXCOL-1;

    color_r[count] = 0;
    color_g[count] = 0;
    color_b[count] = 0;
    count += d;

  for (i=1; i<=S; ++i) {
    color_r[count] = 0;
    color_g[count] = 0;
    color_b[count] = i;
    count += d;
  }

  for (i=1; i<=S; ++i) {
    color_r[count] = 0;
    color_g[count] = i;
    color_b[count] = S;
    count += d;
  }

  for (i=1; i<=S; ++i) {
    color_r[count] = 0;
    color_g[count] = S;
    color_b[count] = S-i;
    count += d;
  }

  for (i=1; i<=S; ++i) {
    color_r[count] = i;
    color_g[count] = S;
    color_b[count] = 0;
    count += d;
  }

  for (i=1; i<=S; ++i) {
    color_r[count] = S;
    color_g[count] = S-i;
    color_b[count] = 0;
    count += d;
  }

  for (i=1; i<=S; ++i) {
    color_r[count] = S;
    color_g[count] = 0;
    color_b[count] = i;
    count += d;
  }

  for (i=1; i<=S; ++i) {
    color_r[count] = S;
    color_g[count] = i;
    color_b[count] = S;
    count += d;
  }

  /*
  for (i=0; i<MAXCOL; ++i) {
    printf("COL[%5d]: %3d,%3d,%3d\n",i,color_r[i],color_g[i],color_b[i]);
  }
  */
}

/*
 * ppmwrite
 */
void ppmwrite( int *a, int nx, int ny, int minval, int maxval, char *filename) {
  FILE *outfile;
  double distance=maxval-minval;
  double factor = (double) distance / MAXCOL;
  char buf[BUFLEN];
  int b,ix,iy,iz;

  gminval=minval;
  gmaxval=maxval;

  if((outfile=fopen(filename,"wt"))==NULL)
    printf("\tUnable to open output file %s.\n",filename);

  b = sprintf(buf, "P6 %5d %5d %d\n", nx, ny, SHADES-1);

  for (iy=0; iy<ny; ++iy) {
    for (ix=0; ix<nx; ++ix) {
/*        int idx = a[ix+iy*nx]; */
      int idx = a[ix*ny+iy];
      if ( idx == maxval ) {
        for (iz=0; iz<3; ++iz) {
          buf[b++] = 0;
          if ( b >= BUFLEN ) {
            fwrite(buf, (size_t) BUFLEN, 1, outfile);
            b = 0;
          }
        }
      } else {
        idx = (int) (idx / factor);
        for (iz=0; iz<3; ++iz) {
	  if(iz==0) buf[b++] = color_r[idx];
	  if(iz==1) buf[b++] = color_g[idx];
	  if(iz==2) buf[b++] = color_b[idx];
/*  	  printf("WF: %3d %3d %1d -> %4d %d\n",ix,iy,iz,idx,buf[b-1]); */
          if ( b >= BUFLEN ) {
            fwrite(buf, (size_t) BUFLEN, 1, outfile);
            b = 0;
          }
        }
      }
    }
  }
  fwrite(buf, (size_t) b, 1, outfile);
  fclose(outfile);
}


FILE* ppmopen( int nx, int ny, int minval, int maxval, char *filename) {
  FILE *outfile;
  double distance=maxval-minval;
  double factor = (double) distance / MAXCOL;
  char buf[BUFLEN];
  int b,ix,iy,iz;

  gminval=minval;
  gmaxval=maxval;

  if((outfile=fopen(filename,"wt"))==NULL) {
    printf("\tUnable to open output file %s.\n",filename);
    return(NULL);
  }
  b = sprintf(buf, "P6 %5d %5d %d\n", nx, ny, SHADES-1);
  fwrite(buf, (size_t) b, 1, outfile);
  return(outfile);
}

/*
 * ppmwrite
 */
void ppmwriteblock(FILE* outfile, double *a, int nx, int ny, double minval, double maxval) {
/*    double distance=maxval-minval; */
  double distance=maxval;
  double factor = (double) distance / MAXCOL;
  char buf[BUFLEN];
  int b,ix,iy,iz;
  gminval=minval;
  gmaxval=maxval;
  b=0;
  for (iy=0; iy<ny; ++iy) {
    for (ix=0; ix<nx; ++ix) {
/*        int idx = a[ix+iy*nx]; */
      int idx = (int) a[ix+iy*nx];
/*         printf("WF: %3d %3d -> %4d %e dist=%e fact=%e\n",ix,iy,idx,a[ix+iy*nx],distance,factor);   */
      if ( a[ix+iy*nx] >= maxval ) {
        for (iz=0; iz<3; ++iz) {
          buf[b++] = 0;
          if ( b >= BUFLEN ) {
            fwrite(buf, (size_t) BUFLEN, 1, outfile);
            b = 0;
          }
        }
      } else {
        idx = (int) (a[ix+iy*nx] / factor);
        for (iz=0; iz<3; ++iz) {
	  if(iz==0) buf[b++] = color_r[idx];
	  if(iz==1) buf[b++] = color_g[idx];
	  if(iz==2) buf[b++] = color_b[idx];
/*     	  printf("WF: %3d %3d %1d -> %4d %d\n",ix,iy,iz,idx,buf[b-1]);   */
          if ( b >= BUFLEN ) {
            fwrite(buf, (size_t) BUFLEN, 1, outfile);
            b = 0;
          }
        }
      }
    }
  }
  fwrite(buf, (size_t) b, 1, outfile);
}

void ppmclose( FILE *outfile) {
  FILE *psoutfile;
  int i,x,y;
  if((psoutfile=fopen("./colortab.ps","wt"))==NULL) {
    printf("\tUnable to open file ./colortab.ps.\n");
    return;
  }
 
  fprintf(psoutfile, "%!PS-Adobe-2.0\n");

  fprintf(psoutfile, "/rect {\n");
  fprintf(psoutfile, "  0 8 rlineto 8 0 rlineto 0 -8 rlineto -8 0 rlineto fill\n");
  fprintf(psoutfile, "} bind def\n");
  fprintf(psoutfile, "/Helvetica findfont 6 scalefont setfont\n");
	  
  for (i=0; i<MAXCOL; ++i) {
    x= (int) (i/100.0) * 40;
    y=i%100 * 8;
    fprintf(psoutfile, "%d %d moveto %3d %3d %3d setrgbcolor rect %d %d moveto (%4d:%14.8f) show\n",
	    x,y,color_r[i],color_g[i],color_b[i],x+20,y, i, 1.0*i/MAXCOL*gmaxval*1000*1000);
  }
 

  fclose(psoutfile);

  fclose(outfile);
}

