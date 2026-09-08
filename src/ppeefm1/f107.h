/*
 * f107.h
 * Patrick Alken
 */

#ifndef INCLUDED_f107_h
#define INCLUDED_f107_h

#include <time.h>

typedef struct
{
  double *f107_data;
  time_t t0;         /* timestamp of first date */
  long t0_d;         /* timestamp in days */
  long offset_d;     /* days offset between epoch and 1-Jan-2000 */
  size_t n;          /* number of f107 data points */
} f107_workspace;

#define F107_MAX_DATA      10000
#define F107_MAX_BUFFER    512

#define FILE_F107          "F107.txt"

#define F107_MIN(a,b)      ((a) < (b) ? (a) : (b))

/*
 * Prototypes
 */

f107_workspace *f107_alloc(const char *filename);
void f107_free(f107_workspace *w);
int f107_get(double fday, double *result, f107_workspace *w);
int f107_get2(time_t ts, double *result, f107_workspace *w);
int f107a_get(double fday, double *result, f107_workspace *w);

extern int putenv(char *string);

#endif /* INCLUDED_f107_h */
