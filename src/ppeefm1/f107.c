/*
 * f107.c
 * Patrick Alken
 *
 * This module reads data from a SPIDR data file containing
 * F10.7 index values and then returns the correct value for
 * a given day.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>

#include "common.h"
#include "f107.h"

/*
read_data()
  Read F107 data from a SPIDR data file
*/

static size_t
read_data(const char *filename, f107_workspace *w)
{
  FILE *fp;
  size_t n;
  int ret;
  char buf[F107_MAX_BUFFER];
  char date[11]; /* date in yyyy-mm-dd format */
  char hour[6];  /* hour in hh:mm format */
  double f107;   /* f107 value */

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  n = 0;
  while (fgets(buf, F107_MAX_BUFFER, fp) != NULL)
    {
      /* ignore comments */
      if (*buf == '#')
        continue;

      /*
       * file format is:
       *
       * yyyy-mm-dd hh:mm value qualifier description
       *
       * with one data value per day
       */
      ret = sscanf(buf,
                   "%10s %5s %lf%*[^\n]",
                   date,
                   hour,
                   &f107);
      if (ret <= 0)
        continue;

      if (w->t0 == 0)
        {
          struct tm tm0;
          int year, month, day;

          ret = sscanf(date, "%d-%d-%d", &year, &month, &day);
          if (ret <= 0)
            {
              fprintf(stderr, "read_data: F107 data file is corrupted\n");
              return 0;
            }

          /*
           * this is the first data line - record the first
           * date timestamp
           */
          tm0.tm_sec = 0;
          tm0.tm_min = 0;
          tm0.tm_hour = 0;
          tm0.tm_mday = day;
          tm0.tm_mon = month - 1;
          tm0.tm_year = year - 1900;
          tm0.tm_isdst = 0;

          w->t0 = mktime(&tm0);

          /* convert to days */
          w->t0_d = (long) (w->t0 / 86400);
        }

      w->f107_data[n] = f107;

      if (++n >= F107_MAX_DATA)
        {
          fprintf(stderr, "read_data: F107_MAX_DATA not large enough\n");
          return n;
        }
    }

  fclose(fp);

  return n;
} /* read_data() */

/*
f107_alloc()
  Allocate a f107 workspace and read in data from filename
*/

f107_workspace *
f107_alloc(const char *filename)
{
  f107_workspace *w;
  size_t i;

  w = (f107_workspace *) calloc(1, sizeof(f107_workspace));
  if (!w)
    return 0;

  w->f107_data = (double *) malloc(sizeof(double) * F107_MAX_DATA);
  for (i = 0; i < F107_MAX_DATA; ++i)
    w->f107_data[i] = 0.0;

  w->t0 = 0;

  /* calculate offset in days between 1-Jan-1970 and 1-Jan-2000 */
  {
    struct tm tm0;
    time_t offset_s;
    char tz_str[] = "TZ=GMT";

    /* use GMT time */
    putenv(tz_str);

    tm0.tm_sec = 0;
    tm0.tm_min = 0;
    tm0.tm_hour = 0;
    tm0.tm_mday = 1;
    tm0.tm_mon = 0;
    tm0.tm_year = 100;
    tm0.tm_isdst = 0;

    offset_s = mktime(&tm0);
    w->offset_d = offset_s / 86400;
  }

  w->n = read_data(filename, w);

  return (w);
} /* f107_alloc() */

void
f107_free(f107_workspace *w)
{
  if (!w)
    return;

  if (w->f107_data)
    free(w->f107_data);

  free(w);
} /* f107_free() */

/*
f107_get()
  Obtain the F10.7 value for a given fday

Inputs: fday   - floating point days since 01-Jan-2000 00:00:00 UTC
        result - where to store F10.7 value
        w      - f107 workspace

Return: 0 on success, 1 on failure
*/

int
f107_get(double fday, double *result, f107_workspace *w)
{
  long fdayi;
  size_t idx;

  /* days since 01-Jan-2000 */
  fdayi = (long) fday;

  /* days since 01-Jan-1970 */
  fdayi += w->offset_d;

  idx = (size_t) (fdayi - w->t0_d);

  if ((fdayi < w->t0_d) || (idx >= w->n))
    {
      fprintf(stderr, "f107_get: data not available for given fday\n");
      return 1;
    }

  *result = w->f107_data[idx];

  return 0;
} /* f107_get() */

/*
f107_get2()
  Obtain the F10.7 value for a given fday

Inputs: ts     - UT timestamp
        result - where to store F10.7 value
        w      - f107 workspace

Return: 0 on success, 1 on failure
*/

int
f107_get2(time_t ts, double *result, f107_workspace *w)
{
  double fday;
  long fdayi;
  size_t idx;

  fday = time2fday(ts);

  /* days since 01-Jan-2000 */
  fdayi = (long) fday;

  /* days since 01-Jan-1970 */
  fdayi += w->offset_d;

  idx = (size_t) (fdayi - w->t0_d);

  if ((fdayi < w->t0_d) || (idx >= w->n))
    {
      fprintf(stderr, "f107_get: data not available for given t\n");
      return 1;
    }

  *result = w->f107_data[idx];

  return 0;
} /* f107_get() */

/*
f107a_get()
  Obtain the F10.7A value for a given day (the 81 day average of
F10.7)

Inputs: fday   - floating point days since 01-Jan-2000 00:00:00 UTC
        result - where to store F10.7A value
        w      - f107 workspace

Return: 0 on success, 1 on failure
*/

int
f107a_get(double fday, double *result, f107_workspace *w)
{
  long fdayi;
  size_t idx;
  size_t mini, maxi;
  double sum;

  /* days since 01-Jan-2000 */
  fdayi = (long) fday;

  /* days since 01-Jan-1970 */
  fdayi += w->offset_d;

  idx = (size_t) (fdayi - w->t0_d);

  if ((fdayi < w->t0_d) || (idx >= w->n))
    {
      fprintf(stderr, "f107_get: data not available for given fday\n");
      return 1;
    }

  maxi = F107_MIN(idx + 40 , w->n - 1);
  if ((int) maxi - 80 < 0)
    mini = 0;
  else
    mini = maxi - 80;

  sum = 0.0;
  for (idx = mini; idx <= maxi; ++idx)
    sum += w->f107_data[idx];

  *result = sum / 81.0;

  return 0;
} /* f107a_get() */
