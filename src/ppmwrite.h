int ppminit(int direction);
int ppminitsmooth(int direction);
void ppmwrite( int *a, int nx, int ny, int minval, int maxval, char *filename);
FILE* ppmopen( int nx, int ny, int minval, int maxval, char *filename);
void ppmwriteblock(FILE* outfile, double *a, int nx, int ny, double minval, double maxval);
void ppmclose( FILE *outfile);



