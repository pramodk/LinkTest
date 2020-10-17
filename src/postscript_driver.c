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

#include "postscript_driver.h"

static int length_of_message_ps = (128*1024);

int PSsetmsglen(int len) {
  length_of_message_ps=len;
}

double time2mbs_ps( double time) {
  if(time>0) {
    return((length_of_message_ps * sizeof(char))/1024.0/1024.0/time );
  } else {
    return(0.0);
  }
}

double mbs2time_ps( double mbs) {
  if(mbs>0) {
    return( (length_of_message_ps*sizeof(char))/(mbs*1024*1024));
  } else {
    return(0.0);
  }
}


int PSplotrowimage (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
		    int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;
  int r,b,g;

  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;

  /* frectgrey(OUT,PSMATX,posy,PSMATW,bheight,linewidth,0.9);  */


  fprintf(OUT, "gsave\n");
  fprintf(OUT, "%f %f translate\n",(double) PSMATX, posy);
  fprintf(OUT, "%f %f scale\n",(double) PSMATW,bheight);
  fprintf(OUT, "%d %d 8 [%d 0 0 %d 0 %d]\n",ntasks,1,ntasks,-1,1);
  fprintf(OUT, "{<");


  for(t=0;t<ntasks;t++) {
    
    if(do_range) {
      printit=((ptimings[t]<opt_mintime) || (ptimings[t]>opt_maxtime));
    }


    posx=PSMATX+t*bwidth;
    val=(ptimings[t]-minval)/(maxval-minval);

    /* val=  1.0*t/ntasks; */

    nval2rgb(val,&r,&g,&b);

    if(printit) {
      if(ptimings[t]==0) {
	r=b=g=0;
      }
    } else {
      r=b=g=0;
    }

    
    fprintf( OUT, "%02X%02X%02X", (unsigned char) r , (unsigned char) g , (unsigned char) b );

    if(t%10==0) fprintf(OUT, "\n");

  }

  fprintf(OUT, ">}\n");

  fprintf(OUT, "false 3 colorimage\n");
  fprintf(OUT, "grestore\n");

  righttext(OUT,location,PSMATX+PSTEXTXOFF,posy+2,textheight);
  righttextrot(OUT,location,PSMATX+rank*bwidth+bwidth,PSMATY+PSTEXTYOFF,textheight);
  sprintf(helpstr,"%d",rank);
  text(OUT,helpstr,PSMATX+PSMATW+2,posy+2,textheight);
  textrot(OUT,helpstr,PSMATX+rank*bwidth+bwidth,PSMATY+PSMATH+2,textheight);
  return(1);
}

int PSplotrowimageClip (FILE *OUT, double *ptimings, int rank, int rankF,int rankT, int ntasks, double minval, double maxval,
			char *locationF, char *locationT, int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;
  int r,b,g;

  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;


  fprintf(OUT, "gsave\n");
  fprintf(OUT, "%f %f translate\n",(double) PSMATX, posy);
  fprintf(OUT, "%f %f scale\n",(double) PSMATW,bheight);
  fprintf(OUT, "%d %d 8 [%d 0 0 %d 0 %d]\n",ntasks,1,ntasks,-1,1);
  fprintf(OUT, "{<");


  for(t=0;t<ntasks;t++) {
    
    if(do_range) {
      printit=((ptimings[t]<opt_mintime) || (ptimings[t]>opt_maxtime));
    }


    posx=PSMATX+t*bwidth;
    val=(ptimings[t]-minval)/(maxval-minval);

    nval2rgb(val,&r,&g,&b);

    if(printit) {
      if(ptimings[t]==0) {
	r=b=g=0;
      }
    } else {
      r=b=g=0;
    }

    
    fprintf( OUT, "%02X%02X%02X", (unsigned char) r , (unsigned char) g , (unsigned char) b );

    if(t%10==0) fprintf(OUT, "\n");

  }

  fprintf(OUT, ">}\n");

  fprintf(OUT, "false 3 colorimage\n");
  fprintf(OUT, "grestore\n");

  righttext(OUT,locationF,PSMATX+PSTEXTXOFF,posy+2,textheight);
  righttextrot(OUT,locationT,PSMATX+rank*bwidth+bwidth,PSMATY+PSTEXTYOFF,textheight);

  sprintf(helpstr,"%d",rankF);
  text(OUT,helpstr,PSMATX+PSMATW+2,posy+2,textheight);

  sprintf(helpstr,"%d",rankT);
  textrot(OUT,helpstr,PSMATX+rank*bwidth+bwidth,PSMATY+PSMATH+2,textheight);
  return(1);
}


int PSplotrow (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
	       int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;

  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;

  /* frectgrey(OUT,PSMATX,posy,PSMATW,bheight,linewidth,0.9);  */


  fprintf(OUT, "[");


  for(t=0;t<ntasks;t++) {
    
    if(do_range) {
      printit=((ptimings[t]<opt_mintime) || (ptimings[t]>opt_maxtime));
    }


    posx=PSMATX+t*bwidth;
    val=(ptimings[t]-minval)/(maxval-minval);
    col=nval2col(val);
    if(printit) {
      if(ptimings[t]>0) {
	fprintf(OUT, "%4.3f ",col);
      } else {
	fprintf(OUT, "%d ",-1);
      }
    } else {
      fprintf(OUT, "%d ",-1);
    }
    if(t%10==0) fprintf(OUT, "\n");

  }

  fprintf(OUT, "]\n");

  fprintf(OUT, "%d   %8.6f %8.6f   %8.6f %8.6f %f  frow\n",ntasks-1,(double) PSMATX,posy,(double) (PSMATX+PSMATW),posy+bheight,linewidth);

  righttext(OUT,location,PSMATX+PSTEXTXOFF,posy+2,textheight);
  righttextrot(OUT,location,PSMATX+rank*bwidth+bwidth,PSMATY+PSTEXTYOFF,textheight);
  sprintf(helpstr,"%d",rank);
  text(OUT,helpstr,PSMATX+PSMATW+2,posy+2,textheight);
  textrot(OUT,helpstr,PSMATX+rank*bwidth+bwidth,PSMATY+PSMATH+2,textheight);
  return(1);
}

int PSplotrow_ori (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
	       int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;

  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;

  frectgrey(OUT,PSMATX,posy,PSMATW,bheight,linewidth,0.9); 

  for(t=0;t<ntasks;t++) {
    
    if(do_range) {
      printit=((ptimings[t]<opt_mintime) || (ptimings[t]>opt_maxtime));
    }


    posx=PSMATX+t*bwidth;
    val=(ptimings[t]-minval)/(maxval-minval);
    col=nval2col(val);
    if(printit) {
      if(ptimings[t]>0) {
	frect(OUT,posx,posy,bwidth,bheight,linewidth,col);
      } else {
/* 	frectgrey(OUT,posx,posy,bwidth,bheight,linewidth,0.8); */
      }
      if(linewidth>0) {
	rect(OUT,posx,posy,bwidth,bheight,linewidth,0.1);
      }
    }
  }
  righttext(OUT,location,PSMATX+PSTEXTXOFF,posy+2,textheight);
  righttextrot(OUT,location,PSMATX+rank*bwidth+bwidth,PSMATY+PSTEXTYOFF,textheight);
  sprintf(helpstr,"%d",rank);
  text(OUT,helpstr,PSMATX+PSMATW+2,posy+2,textheight);
  textrot(OUT,helpstr,PSMATX+rank*bwidth+bwidth,PSMATY+PSMATH+2,textheight);
  return(1);
}

int PSplotrowClip (FILE *OUT, double *ptimings, int rank, int rankF,int rankT, int ntasks, double minval, double maxval, char *locationF, char *locationT, 
	           int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;

  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;

  frectgrey(OUT,PSMATX,posy,PSMATW,bheight,linewidth,0.9); 

  for(t=0;t<ntasks;t++) {
    
    if(do_range) {
      printit=((ptimings[t]<opt_mintime) || (ptimings[t]>opt_maxtime));
    }


    posx=PSMATX+t*bwidth;
    val=(ptimings[t]-minval)/(maxval-minval);
    col=nval2col(val);
    if(printit) {
      if(ptimings[t]>0) {
	frect(OUT,posx,posy,bwidth,bheight,linewidth,col);
      } else {
/* 	frectgrey(OUT,posx,posy,bwidth,bheight,linewidth,0.8); */
      }
      if(linewidth>0) {
	rect(OUT,posx,posy,bwidth,bheight,linewidth,0.1);
      }
    }
  }
  righttext(OUT,locationF,PSMATX+PSTEXTXOFF,posy+2,textheight);

  righttextrot(OUT,locationT,PSMATX+rank*bwidth+bwidth,PSMATY+PSTEXTYOFF,textheight);

  sprintf(helpstr,"%d",rankF);
  text(OUT,helpstr,PSMATX+PSMATW+2,posy+2,textheight);

  sprintf(helpstr,"%d",rankT);
  textrot(OUT,helpstr,PSMATX+rank*bwidth+bwidth,PSMATY+PSMATH+2,textheight);
  return(1);
}


int PSheader(FILE *OUT) {
  fprintf( OUT, "%%!PS-Adobe-3.0\n");
  fprintf( OUT, "%%%%Title: Communication analysis FZJ Linktest \n");
  fprintf( OUT, "%%%%Creator: fzjlinktest Version 1.0, Forschungszentrum Juelich, JSC\n");
  fprintf( OUT, "%%%%DocumentFonts: Times-Roman, Times-Bold Courier-Bold\n");
  fprintf( OUT, "%%%%BoundingBox: 0 0 610 790\n");
  fprintf( OUT, "%%%%Orientation: Portrait\n");
  fprintf( OUT, "%%%%EndComments\n");
/*   fprintf( OUT, "40 60 translate\n"); */
/*   fprintf( OUT, "0.45 0.56 scale\n"); */
  fprintf( OUT, "%%%%BeginProlog\n");
  fprintf( OUT, "/fr {\n");
  fprintf( OUT, "  /col exch def /lw exch def /y2 exch def /x2 exch def /y1 exch def /x1 exch def \n");
  fprintf( OUT, "  gsave lw setlinewidth col 0.8 0.8 sethsbcolor x1 y1 moveto x2 y1 lineto\n");
  fprintf( OUT, "        x2 y2 lineto x1 y2 lineto closepath fill\n");
  fprintf( OUT, "  grestore\n");
  fprintf( OUT, "} bind def\n");

  fprintf( OUT, "/frectgrey {\n");
  fprintf( OUT, "  /col exch def /lw exch def /y2 exch def /x2 exch def /y1 exch def /x1 exch def \n");
  fprintf( OUT, "  gsave lw setlinewidth col col col setrgbcolor x1 y1 moveto x2 y1 lineto\n");
  fprintf( OUT, "        x2 y2 lineto x1 y2 lineto closepath fill\n");
  fprintf( OUT, "  grestore\n");
  fprintf( OUT, "} bind def\n");

  fprintf( OUT, "/rect {\n");
  fprintf( OUT, "  /col exch def /lw exch def /y2 exch def /x2 exch def /y1 exch def /x1 exch def \n");
  fprintf( OUT, "  gsave lw setlinewidth col 0.8 0.8 sethsbcolor x1 y1 moveto x2 y1 lineto\n");
  fprintf( OUT, "        x2 y2 lineto x1 y2 lineto closepath stroke\n");
  fprintf( OUT, "  grestore\n");
  fprintf( OUT, "} bind def\n");

  fprintf( OUT, "/frow {\n");
  fprintf( OUT, "  /lw exch def /y2 exch def /x2 exch def /y1 exch def /x1 exch def /num exch def /data exch def\n");
  fprintf( OUT, "  gsave\n");
  fprintf( OUT, "  x2 x1 sub num div /distw exch def\n");
  fprintf( OUT, "  0 1 num {\n");
  fprintf( OUT, "	  /ind exch def\n");
  /* fprintf( OUT, "	  x1 =\n"); */
  fprintf( OUT, "	  /x1 x1 distw add def\n");
  fprintf( OUT, "	  x1 y1\n");
  fprintf( OUT, "	  x1 y2 lw data ind get fr\n");
  fprintf( OUT, "      } for \n");
  fprintf( OUT, "  grestore\n");
  fprintf( OUT, "} bind def\n");

  fprintf( OUT, "/rightText {				%% (str) rightText -\n");
  fprintf( OUT, "  	dup stringwidth	        	%% (str) wx wy\n");
  fprintf( OUT, "  	pop     			%% (str) wx\n");
  fprintf( OUT, "  	neg 0 rmoveto			%% (str) \n");
  fprintf( OUT, "  	show\n");
  fprintf( OUT, "} def\n");
  fprintf( OUT, "%%%%EndProlog\n");
  fprintf( OUT, "%%%%Page: 1 1\n");
  fprintf( OUT, "%%%%%%BeginPageSetup\n");
  fprintf( OUT, "%%%%EndPageSetup\n");

  return(1);
}

int PSfooter(FILE *OUT) {
  fprintf( OUT, "%% Footer \n");
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "/Times-Bold findfont 6 scalefont setfont 0 0 0 setrgbcolor\n");
  fprintf(OUT, "340 20 moveto\n");
  fprintf(OUT, "(Report generated by FZJ Linktest Result Analyzer, Forschungszentrum Juelich GmbH) show\n");
  fprintf(OUT, "grestore\n");
  fprintf( OUT, "showpage\n");
  fprintf( OUT, "%%%%PageTrailer\n");
  fprintf( OUT, "%%%%Trailer\n");
  fprintf( OUT, "%%%%Pages: 1\n");
  fprintf( OUT, "%%%%EOF\n");
   
  return(1);
}

double nval2col (double val) {
  double col=1.0-(0.2+0.8*val); 
/*   double col=0.2+0.8*val;  */
/*     double col=1.0-0.9*val; */
    return(col);
}

double nval2rgb (double val, int *rr, int *gg, int *bb) {
  double col=1.0-(0.2+0.8*val); 
  double h,s,v;
  double res1, val1, val2, val3;
  int sel; 
  int r,g,b;
  h=(int) (col*360.0);
  s=0.8*100;
  v=0.8*100;
  
  if (s == 0) {
    r = b = g = (int) v;
  } else {
    h /= 60;
    sel = (int) h;
    res1 = h - sel;
    val1 = v/100.0 * (100.0 - s);
    val2 = v/100.0 * (100.0 - s * res1);
    val3 = v/100.0 * (100.0 - s * (1 - res1));
    if (sel == 0) {
      r = (int) v;
      g = (int) val3;
      b = (int) val1;
    } else if (sel == 1) {
      r = (int) val2;
      g = (int) v;
      b = (int) val1;
    } else if (sel == 2) {
      r = (int) val1;
      g = (int) v;
      b = (int) val3;
    } else if (sel == 3) {
      r = (int) val1;
      g = (int) val2;
      b = (int) v;
    } else if (sel == 4) {
      r = (int) val3;
      g = (int) val1;
      b = (int) v;
    } else {
      r = (int) v;
      g = (int) val1;
      b = (int) val2;
    }
  }

  *rr=(int) (255.0/100.0*r);
  *gg=(int) (255.0/100.0*g);
  *bb=(int) (255.0/100.0*b);

  /* printf( "nval2rgb %8.4f -> %3d %3d %3d -> %3d %3d %3d\n",val,h,s,v,*rr,*gg,*bb);  */
  return(1);
}


int PSplotscala (FILE *OUT, int nelem, double minval, double maxval) {
  int t;
  double col,rval,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  bheight=(double) PSSCALEH/nelem;
  bwidth =(double) PSSCALEW;
  posx=PSSCALEX;
  for(t=0;t<nelem;t++) {
    posy=PSSCALEY+PSSCALEH-(t+1)*bheight;
    rval=minval+(maxval-minval)*(double) t/ (double) nelem;
    val=(rval-minval)/(maxval-minval);
    col=nval2col(val);
    frect(OUT,posx,posy,bwidth,bheight,0,col);
    rect(OUT,posx,posy,bwidth,bheight,0.1,0.1);
    sprintf(helpstr,"% 7.1fus (%7.2f MB/s)",rval*1000*1000,time2mbs_ps(rval));
    text(OUT,helpstr,posx+bwidth+1,posy+2,bheight);
  }
  return(1);
}

int PSplothist (FILE *OUT, long *count, int numsteps, double minval, double maxval, double stepwidth) {
  int s;
  long c,maxcount;
  double maxcountr,cr;
  double col,val,rval;
  double posx,posy,bheight,bwidth;
  char helpstr[256];

  bwidth =(double) PSHISTW/ (double)numsteps;
  maxcount=0;
  for(s=0;s<numsteps;s++) {
    if(count[s]>maxcount) maxcount=count[s];
  }
  if(maxcount==0) return(0);

  rect(OUT,PSHISTX,PSHISTY,PSHISTW,PSHISTH,0.1,0.1);

  maxcountr=log(maxcount);
  posy=PSHISTY;

  for(c=1;c<maxcount;c*=2) {
      bheight=log((double) c)/maxcountr * PSHISTH;
      sprintf(helpstr,"% 10d",(int) c);
      text(OUT,helpstr,PSHISTX-40,posy+bheight,8);
      line(OUT,PSHISTX,posy+bheight,PSHISTX+PSHISTW,posy+bheight,0.1,0.1);
  }
  for(s=0;s<=10;s++) {
    posx=PSHISTX+PSHISTW-(s+0.0)*PSHISTW/10.0;
    rval=minval+(maxval-minval)*(double) s/ (double) 10.0;
    sprintf(helpstr,"%7.2fus",rval*1000*1000);
    text(OUT,helpstr,posx-PSHISTW/10*0.4,posy-10,6);
    sprintf(helpstr,"%7.2fMB/s",time2mbs_ps(rval));
    text(OUT,helpstr,posx-PSHISTW/10*0.4,posy-15,6);
    line(OUT,posx,PSHISTY,posx,PSHISTY-5,0.1,0.1);
  }

  for(s=0;s<numsteps;s++) {
    posx=PSHISTX+PSHISTW-((double)s+1.0)*bwidth;
    rval=minval+(maxval-minval)*(double) s/ (double) numsteps; 
/*     rval=minval+s*stepwidth; */
    val=(rval-minval)/(maxval-minval);
    col=nval2col(val);
    if(count[s]>0) {
      bheight=log((double) count[s])/maxcountr * PSHISTH;
      frect(OUT,posx,posy,bwidth,bheight,0.1,col);
    }
  }

  return(1);
}

int frect (FILE *OUT, double x1, double y1, double w, double h, double lw, double color ) {
  double x2=x1+w;
  double y2=y1+h;
  fprintf(OUT, "%6.2f %6.2f %6.2f %6.2f %6.2f %5.3f fr\n",x1,y1,x2,y2,lw,color);
  return(1);
}

int frectgrey (FILE *OUT, double x1, double y1, double w, double h, double lw, double color ) {
  double x2=x1+w;
  double y2=y1+h;
  fprintf(OUT, "%f %f %f %f %f %f frectgrey\n",x1,y1,x2,y2,lw,color);
  return(1);
}

int rect (FILE *OUT, double x1, double y1, double w, double h, double lw, double color ) {
  double x2=x1+w;
  double y2=y1+h;
  fprintf(OUT, "%f %f %f %f %f %f rect\n",x1,y1,x2,y2,lw,color);
  return(1);
}

int line(FILE *OUT, double x1, double y1, double x2, double y2, double lw, double color ) {
  fprintf(OUT, "gsave %f setlinewidth %f 0.8 0.8 sethsbcolor %f %f moveto %f %f lineto stroke grestore\n",lw,color,x1,y1,x2,y2);
  return(1);
}

int text (FILE *OUT, char *strtext, double x, double y, double size) {
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "/Courier-Bold findfont %f scalefont setfont 0 0 0 setrgbcolor\n",size);
  fprintf(OUT, "%f %f moveto\n",x,y);
  fprintf(OUT, "(%s) show\n",strtext);
  fprintf(OUT, "grestore\n");
  return(1);
}

int righttext (FILE *OUT, char *strtext, double x, double y, double size) {
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "/Courier-Bold findfont %f scalefont setfont 0 0 0 setrgbcolor\n",size);
  fprintf(OUT, "%f %f moveto\n",x,y);
  fprintf(OUT, "(%s) rightText\n",strtext);
  fprintf(OUT, "grestore\n");
  return(1);
}

int textr (FILE *OUT, char *strtext, double x, double y, double size) {
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "/Times-Bold findfont %f scalefont setfont 0 0 0 setrgbcolor\n",size);
  fprintf(OUT, "%f %f moveto\n",x,y);
  fprintf(OUT, "(%s) show\n",strtext);
  fprintf(OUT, "grestore\n");
  return(1);
}

int textrot (FILE *OUT, char *strtext, double x, double y, double size) {
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "/Courier-Bold findfont %f scalefont setfont 0 0 0 setrgbcolor\n",size);
  fprintf(OUT, "%f %f moveto\n",x,y);
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "90 rotate\n");
  fprintf(OUT, "(%s) show\n",strtext);
  fprintf(OUT, "grestore\n");
  fprintf(OUT, "grestore\n");
  return(1);
}

int righttextrot (FILE *OUT, char *strtext, double x, double y, double size) {
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "/Courier-Bold findfont %f scalefont setfont 0 0 0 setrgbcolor\n",size);
  fprintf(OUT, "%f %f moveto\n",x,y);
  fprintf(OUT, "gsave\n");
  fprintf(OUT, "90 rotate\n");
  fprintf(OUT, "(%s) rightText\n",strtext);
  fprintf(OUT, "grestore\n");
  fprintf(OUT, "grestore\n");
  return(1);
}


int PSplotrowimage_start (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			  int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;
  int r,b,g;

  rank=0;
  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;

  frectgrey(OUT,PSMATX,posy,PSMATW,bheight,linewidth,0.9); 


  fprintf(OUT, "gsave\n");
  fprintf(OUT, "%f %f translate\n",(double) PSMATX, posy);
  fprintf(OUT, "%f %f scale\n",(double) PSMATW,bheight);
  fprintf(OUT, "%d %d 8 [%d 0 0 %d 0 %d]\n",ntasks,ntasks,ntasks,-1,1);
  fprintf(OUT, "{<");
  return(1);
}


int PSplotrowimage_row (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;
  int r,b,g;

  for(t=0;t<ntasks;t++) {
    
    if(do_range) {
      printit=((ptimings[t]<opt_mintime) || (ptimings[t]>opt_maxtime));
    }

    val=(ptimings[t]-minval)/(maxval-minval);

    nval2rgb(val,&r,&g,&b);

    if(printit) {
      if(ptimings[t]==0) {
	r=b=g=0;
      }
    } else {
      r=b=g=0;
    }

    
    fprintf( OUT, "%02X%02X%02X", (unsigned char) r , (unsigned char) g , (unsigned char) b );

    if(t%10==0) fprintf(OUT, "\n");

  }

  return(1);
}

int PSplotrowimage_end (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			int do_range, double opt_mintime, double opt_maxtime) {

  fprintf(OUT, ">}\n");

  fprintf(OUT, "false 3 colorimage\n");
  fprintf(OUT, "grestore\n");

  return(1);
}

int PSplotrowimage_text (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			int do_range, double opt_mintime, double opt_maxtime) {
  int t;
  double col,val;
  double posx,posy,bheight,bwidth;
  char helpstr[256];
  double linewidth;
  double textheight;
  int printit=1;
  int r,b,g;

  bheight=(double) PSMATH/ntasks;
  bwidth =(double) PSMATW/ntasks;
  textheight=0.6*bheight;
  if(textheight>10) textheight=10;
  posy=PSMATY+rank*bheight;
  linewidth=0.5;
  if(ntasks>32) linewidth=0.0;
  if(ntasks>1024) linewidth=0.0;

  righttext(OUT,location,PSMATX+PSTEXTXOFF,posy+2,textheight);
  righttextrot(OUT,location,PSMATX+rank*bwidth+bwidth,PSMATY+PSTEXTYOFF,textheight);
  sprintf(helpstr,"%d",rank);
  text(OUT,helpstr,PSMATX+PSMATW+2,posy+2,textheight);
  textrot(OUT,helpstr,PSMATX+rank*bwidth+bwidth,PSMATY+PSMATH+2,textheight);

  return(1);
}
