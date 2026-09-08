/*
 * rteef.c
 *
 * Usage:
 *
 * rteef [start UT] [end UT] [longitude(degrees)]
 *
 * where start and end UT are strings of the form:
 *
 * Month DD YYYY HH:MM
 *
 * and Month is the complete month or 3 letter abbreviation
 */

#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#include "common.h"
#include "f107.h"
#include "eefm_calc.h"
#include "lunar.h"
#include "rteef.h"

static int read_rt_coeffs(rteef_workspace *w);
static int read_ief_data(const char *filename, time_t t_start, time_t t_end,
                         rteef_workspace *w);

/* output file */
#define FILE_OUTPUT         "eef_data"

/* time window in hours to read in ACE data prior to desired EEF time */
#define ACE_WINDOW_PRIOR    (1.0)

/*
 * maximum number of days to compute EEF at one time - the user specifies
 * start and end times but the end time cannot exceed this number of
 * days from the start
 */
#define MAX_TIME_INTERVAL   (5.0)

#define MAX_BUFFER          4096

/*
rteef_alloc()
  Allocate a rteef workspace
*/

rteef_workspace *
rteef_alloc(void)
{
  rteef_workspace *w;
  int s;
  size_t i;

  w = calloc(1, sizeof(rteef_workspace));
  if (!w)
    {
      fprintf(stderr, "rteef_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  s = read_rt_coeffs(w);

  /*
   * To find out how much ACE data to read in, take the total window
   * size and multiply by 12 since data occurs at 5 minute intervals
   */
  w->ndata_ief_max = (ACE_WINDOW_PRIOR + MAX_TIME_INTERVAL*24.0) * 12.0;

  w->data_ief = calloc(w->ndata_ief_max, sizeof(double));
  w->eef_rt = calloc(w->ndata_ief_max, sizeof(double));
  w->eef_mean = calloc(w->ndata_ief_max, sizeof(double));
  w->dbuffer = calloc(w->nb, sizeof(double));

  if (!w->data_ief || !w->eef_rt || !w->eef_mean || !w->dbuffer)
    {
      fprintf(stderr, "rteef_alloc: calloc failed: %s\n", strerror(errno));
      rteef_free(w);
      return 0;
    }

  for (i = 0; i < w->ndata_ief_max; ++i)
    {
      w->eef_rt[i] = 0.0;
      w->eef_mean[i] = 0.0;
    }

  w->eefm_calc_p = eefm_calc_alloc(EEF_CALC_CHAMP);
  w->f107_p = f107_alloc(FILE_F107);

  return w;
} /* rteef_alloc() */

void
rteef_free(rteef_workspace *w)
{
  if (w->data_ief)
    free(w->data_ief);

  if (w->eef_rt)
    free(w->eef_rt);

  if (w->eef_mean)
    free(w->eef_mean);

  if (w->dbuffer)
    free(w->dbuffer);

  if (w->eefm_calc_p)
    eefm_calc_free(w->eefm_calc_p);

  if (w->f107_p)
    f107_free(w->f107_p);

  free(w);
} /* rteef_free() */

/*
rteef_compute()
  Compute the real time EEF estimate for a given time

Inputs: t_start - timestamp start (UT)
        t_end   - timestamp end (UT)
        longitude - longitude of EEF in degrees
        w       - workspace

Return: success or failure

Notes:
*/

int
rteef_compute(time_t t_start, time_t t_end, double longitude,
              rteef_workspace *w, const double *ief_override, size_t nief)
{
  size_t j, k;
  struct tm *tm_p;
  char ief_filename[100];
  size_t na = w->na;
  time_t t, t_ace_start;
  int s;
  double ut, lt, season, euvac, f107, f107a, tau;
  double phi = longitude * M_PI / 180.0;

  /* read in ACE data 1 hour prior to desired start time */
  t_ace_start = t_start - 3600;

  if (ief_override)
    {
      /* synthetic IEF Ey: fill the buffer directly, skip the ACE file */
      if (nief > w->ndata_ief_max)
        nief = w->ndata_ief_max;
      for (j = 0; j < nief; ++j)
        w->data_ief[j] = ief_override[j];
      w->ndata_ief = nief;
    }
  else
    {
      tm_p = gmtime(&t_ace_start);
      sprintf(ief_filename, "ace%d.dat", tm_p->tm_year + 1900);
      s = read_ief_data(ief_filename, t_ace_start, t_end, w);
      if (s)
        return s;
    }

  na += 1;
  t = t_ace_start;

  for (j = 0; j < w->ndata_ief; j++)
    {
      for (k = 0; k < (w->nb - 1); k++)
        w->dbuffer[k] = w->dbuffer[k + 1];

      w->dbuffer[w->nb-1] = 0.0;

      for (k = 0; k < w->nb; k++)
        w->dbuffer[k] += w->data_ief[j] * w->b[k];

      for (k = 2; k < na; k++)
        w->dbuffer[k - 1] -= w->dbuffer[0] * w->a[k - 1];

      w->eef_rt[j] = w->dbuffer[0];

      /* add climatological part (skipped for a synthetic transfer-function run) */
      if (ief_override)
        {
          t += 300;
          continue;
        }

      tm_p = gmtime(&t);

      ut = (double) tm_p->tm_hour +
           (double) tm_p->tm_min / 60.0 +
           (double) tm_p->tm_sec / 3600.0;

      lt = ut + longitude / 15.0;
      if (lt < 0.0)
        lt += 24.0;
      if (lt > 24.0)
        lt -= 24.0;

      if (lt >= EEF_BASIS_LT_MIN && lt <= EEF_BASIS_LT_MAX)
        {
          season = (double) tm_p->tm_yday;

          f107_get2(t, &f107, w->f107_p);
          f107a_get(time2fday(t), &f107a, w->f107_p);
          euvac = 0.5 * (f107 + f107a);

          tau = lunartime(t, phi);

          w->eef_mean[j] = eefm_calc_mean(phi, lt, season, euvac,
                                          tau, w->eefm_calc_p);
        }

      /* advance t by 5 minutes */
      t += 300;
    }

  return 0;
} /* rteef_compute() */

/*******************************************************
 *             INTERNAL ROUTINES                       *
 *******************************************************/

/*
read_rt_coeffs()
  Read real-time coeffs from TF.COF file and store in w->a and w->b
*/

static int
read_rt_coeffs(rteef_workspace *w)
{
  int i;
  int na, nb;
  FILE *fp;

  fp = fopen("TF.COF", "r");
  if (!fp)
    {
      fprintf(stderr, "TF.COF not found\n");
      return 1;
    }

  fscanf(fp, "%d", &na);
  for (i = 0; i < na; ++i)
    fscanf(fp, "%lf", &(w->a[i]));

  fscanf(fp, "%d", &nb);
  for (i = 0; i < nb; ++i)
    fscanf(fp, "%lf", &(w->b[i]));

  w->na = na;
  w->nb = nb;

  fclose(fp);

  return 0;
} /* read_rt_coeffs() */

/*
read_ief_data()
  Read IEF (ACE) data file

Inputs: filename - ACE data file
        w        - workspace
*/

static int
read_ief_data(const char *filename, time_t t_start, time_t t_end,
              rteef_workspace *w)
{
  FILE *fp;
  char buffer[MAX_BUFFER];
  int year, doy, hour, min, month, day;
  double Ey;
  time_t t;
  struct tm tms;
  size_t cnt;
  double lastEy = 0.0;
  char tz_str[] = "TZ=GMT";

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "%s not found\n", filename);
      return 1;
    }

  putenv(tz_str);

  /* allow for prior window */
  t_start -= ACE_WINDOW_PRIOR * 3600.0;

  cnt = 0;

  while (fgets(buffer, MAX_BUFFER, fp) != NULL)
    {
      int n;

      n = sscanf(buffer,
                 "%d %d %d %d %lf\n",
                 &year,
                 &doy,
                 &hour,
                 &min,
                 &Ey);
      if (n < 5)
        continue;

      doy2md(year, doy, &month, &day);

      tms.tm_sec = 0;
      tms.tm_min = min;
      tms.tm_hour = hour;
      tms.tm_mon = month - 1;
      tms.tm_mday = day;
      tms.tm_year = year - 1900;
      tms.tm_isdst = 0;

      t = mktime(&tms);

      if (t < t_start)
        continue;

      if (t > t_end)
        break; /* done */

      if (Ey > 500.0)
        w->data_ief[cnt] = lastEy;
      else
        w->data_ief[cnt] = Ey;

      if (++cnt == w->ndata_ief_max)
        {
          fprintf(stderr, "read_ief_data: too much data?\n");
          break;
        }

      lastEy = Ey;
    }

  fclose(fp);

  w->ndata_ief = cnt;

  return 0;
} /* read_ief_data() */

int
main(int argc, char *argv[])
{
  size_t i;
  FILE *fp;
  time_t t_start, t_end, t;
  rteef_workspace *rteef_p;
  double longitude;
  char *s;
  struct tm tm_start, tm_end;
  char buffer[MAX_BUFFER];
  char tz_str[] = "TZ=GMT";

  putenv(tz_str);

  if (argc > 1)
    strncpy(buffer, argv[1], MAX_BUFFER);
  else
    {
      fprintf(stderr, "Please enter starting universal time (Ex: Aug 01 2004 12:30)\n");
      fgets(buffer, MAX_BUFFER, stdin);
    }

  s = strstr(buffer, "\n");
  if (s)
    *s = '\0';

//  strncat(buffer, " UTC", MAX_BUFFER);
    
  s = strptime(buffer, "%b %d %Y %H:%M", &tm_start);
    
  if (!s)
    {
      fprintf(stderr, "Error in starting time: %s\n", buffer);
      exit(1);
    }

  tm_start.tm_sec = 0;
  t_start = mktime(&tm_start);

  if (argc > 2)
    strncpy(buffer, argv[2], MAX_BUFFER);
  else
    {
      fprintf(stderr, "Please enter ending universal time (Ex: Aug 02 2004 12:30)\n");
      fgets(buffer, MAX_BUFFER, stdin);
    }

  s = strstr(buffer, "\n");
  if (s)
    *s = '\0';

//  strncat(buffer, " UTC", MAX_BUFFER);
    
  s = strptime(buffer, "%b %d %Y %H:%M", &tm_end);
  if (!s)
    {
      fprintf(stderr, "Error in ending time: %s\n", buffer);
      exit(1);
    }

  tm_end.tm_sec = 0;
  t_end = mktime(&tm_end);

  if (argc > 3)
    longitude = atof(argv[3]);
  else
    {
      fprintf(stderr, "Please enter longitude in degrees: ");
      fscanf(stdin, "%lf", &longitude);
    }

  fprintf(stderr, "Starting time: %s", asctime(&tm_start));
  fprintf(stderr, "Ending time: %s", asctime(&tm_end));
  fprintf(stderr, "Longitude: %f degrees\n", longitude);

  rteef_p = rteef_alloc();

  fprintf(stderr, "Computing real time EEF...");
  rteef_compute(t_start, t_end, longitude, rteef_p, NULL, 0);
  fprintf(stderr, "done\n");

  /* write output file */

  fprintf(stderr, "Writing output file %s...", FILE_OUTPUT);

  fp = fopen(FILE_OUTPUT, "w");
  if (!fp)
    {
      perror("fopen");
      exit(1);
    }

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT)\n", i++);
  fprintf(fp, "# Field %zu: Mean electric field (mV/m)\n", i++);
  fprintf(fp, "# Field %zu: Real-time electric field (mV/m)\n", i++);
  fprintf(fp, "# Field %zu: Total electric field (mV/m)\n", i++);

  /*
   * Since we read in 1 hour of ACE data prior to desired start time,
   * begin indexing at 12 (there are 12 5 minute intervals in 1 hour)
   */
  i = 12;
  for (t = t_start; t <= t_end; t += 300)
    {
      fprintf(fp, "%ld %.12e %.12e %.12e\n",
              t,
              rteef_p->eef_mean[i],
              rteef_p->eef_rt[i],
              rteef_p->eef_mean[i] + rteef_p->eef_rt[i]);
      ++i;
    }

  fclose(fp);

  fprintf(stderr, "done\n");

  rteef_free(rteef_p);

  return 0;
} /* main() */
