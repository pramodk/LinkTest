/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef POSTSCRIPT_DRIVER
#define POSTSCRIPT_DRIVER

#define PSMATX   80
#define PSMATY  320
#define PSMATW  420
#define PSMATH  420
#define PSTEXTXOFF  -2
#define PSTEXTYOFF -2
#define PSSCALEX  520
#define PSSCALEY  320
#define PSSCALEH  420
#define PSSCALEW   20
#define PSSCALENELEM   96
#define PSTITELX   80
#define PSTITELY  760
#define PSDESCX    40
#define PSDESCY    40
#define PSDESCH    75
#define PSDESCW   540
#define PSHISTX    60
#define PSHISTY   160
#define PSHISTH   100
#define PSHISTW   520

int PSsetmsglen(int len);
int PSheader(FILE *OUT);
int PSfooter(FILE *OUT);
double nval2col (double val);
double nval2rgb (double val, int *rr, int *gg, int *bb);

int PSplotrowimage_start (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			  int do_range, double opt_mintime, double opt_maxtime);
int PSplotrowimage_row (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			int do_range, double opt_mintime, double opt_maxtime);
int PSplotrowimage_end (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			int do_range, double opt_mintime, double opt_maxtime);
int PSplotrowimage_text (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
			int do_range, double opt_mintime, double opt_maxtime);
int PSplotrowimage (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location, 
		    int do_range, double opt_mintime, double opt_maxtime);


int PSplotrow (FILE *OUT, double *ptimings, int rank, int ntasks, double minval, double maxval, char *location,
	       int do_range, double opt_mintime, double opt_maxtime);
int PSplotrowClip (FILE *OUT, double *ptimings, int rank, int rankF,int rankT, int ntasks, double minval, double maxval, char *locationF, char *locationT, 
	           int do_range, double opt_mintime, double opt_maxtime);
int PSplotrowimageClip (FILE *OUT, double *ptimings, int rank, int rankF,int rankT, int ntasks, double minval, double maxval,
			char *locationF, char *locationT, int do_range, double opt_mintime, double opt_maxtime);


int PSplotscala (FILE *OUT, int nelem, double minval, double maxval);
int PSplothist (FILE *OUT, long *count, int numsteps, double minval, double maxval, double stepwidth);
int frect (FILE *OUT, double x1, double y1, double w, double h, double lw, double color );
int frectgrey (FILE *OUT, double x1, double y1, double w, double h, double lw, double color );
int rect (FILE *OUT, double x1, double y1, double w, double h, double lw, double color );
int line(FILE *OUT, double x1, double y1, double x2, double y2, double lw, double color );
int text (FILE *OUT, char *strtext, double x, double y, double size);
int righttext (FILE *OUT, char *strtext, double x, double y, double size);
int textr (FILE *OUT, char *strtext, double x, double y, double size);
int textrot (FILE *OUT, char *strtext, double x, double y, double size);
int righttextrot (FILE *OUT, char *strtext, double x, double y, double size);
#endif
