/*
 * lunar.c
 * Patrick Alken
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "lunar.h"

static double normalize(double x);

/*
lunartime()
  Compute the local lunar time for a given UT and longitude

Return: local lunar time in hours (0 to 24.8333)
*/

double
lunartime(time_t t, double longitude)
{
  double tau = lunartime_r(t, longitude);

  /* convert to hours */
  tau = tau / 2.0 / M_PI * 24.83333333333;

  return (tau);
} /* lunartime() */

/*
lunartime_r()
  Compute the local lunar time for a given UT and longitude

Return: local lunar time in radians (0 to 2pi)
*/

double
lunartime_r(time_t t, double longitude)
{
  double nu = lunarphase(t);
  double lt = get_localtime2(t, longitude);
  double tau;

  /* convert LT to an angle between 0 and 2pi */
  lt = lt / 23.9444444 * 2.0 * M_PI;
  if (lt > 2.0 * M_PI)
    lt -= 2.0 * M_PI;

  /* compute local lunar time in radians */
  tau = lt - nu;

  if (tau < 0.0)
    tau += 2.0 * M_PI;

  return (tau);
} /* lunartime_r() */

/*
lunarphase()
  Compute the lunar phase for a given timestamp

Inputs: t - timestamp UTC

Return: lunar phase in radians between 0 and 2pi. A phase of 0 corresponds
to a new moon and a phase of pi corresponds to a full moon
*/

double
lunarphase(time_t t)
{
  lunardata data;
  struct tm *tm_p;
  int year, day, month;

  tm_p = gmtime(&t);
  year = tm_p->tm_year + 1900;
  month = tm_p->tm_mon + 1;
  day = tm_p->tm_mday;

  lunarcalc(year, month, day, &data);

  return data.phase;
} /* lunarphase() */

int
lunarcalc(int year, int month, int day, lunardata *data)
{
  int yy, mm;
  int k1, k2, k3;
  int jd;
  double ip, ag;
  double dp, di, np, rp;
  double lat, lon;

  yy = year - ((12 - month) / 10);

  mm = month + 9;
  if (mm >= 12)
    mm -= 12;

  k1 = (int) (365.25 * (yy + 4712));
  k2 = (int) (30.6 * mm + 0.5);
  k3 = (int) ((int) ((yy / 100.0) + 49.0) * 0.75) - 38;

  jd = k1 + k2 + day + 59;
  if (jd > 2299160)
    jd -= k3;

  /* calculate illumination (synodic) phase */
  ip = normalize((jd - 2451550.1) / 29.530588853);
  ag = ip * 29.53;  /* moon's age in days */
  ip *= 2.0 * M_PI; /* convert phase to radians */

  /* calculate distance from anomalistic phase */
  dp = normalize((jd - 2451562.2) / 27.55454988);
  dp *= 2.0 * M_PI;
  di = 60.4 - 3.3 * cos(dp) - 0.6 * cos(2*ip-dp) - 0.5*cos(2*ip);

  /* calculate latitude from nodal phase */
  np = 2.0 * M_PI * normalize((jd - 2451565.2) / 27.212220817);
  lat = 5.1 * sin(np);

  /* calculate longitude from sidereal motion */
  rp = normalize((jd - 2451555.8) / 27.321582241);
  lon = 360.0 * rp + 6.3 * sin(dp) + 1.3 * sin(2*ip-dp) + 0.7 * sin(2*ip);

  data->phase = ip;
  data->age = ag;
  data->eclat = lat * M_PI / 180.0;
  data->eclon = lon * M_PI / 180.0;
  data->dist = di;

  return 0;
} /* lunarcalc() */

static double
normalize(double x)
{
  double y = x - (int)x;

  if (y < 0.0)
    y += 1.0;

  return y;
} /* normalize() */

/********************************************************
 * Stefan's Lunar code
 ********************************************************/

#define JD0 2437077
#define JD2000 2451544.5

#define M2_PHASE0 75.8
#define M2_PERIOD 12.42060122

#define N2_PHASE0 325.14
#define N2_PERIOD 12.65834824

double lunarphase_m2(time_t t)
{
  double fday = time2fday(t);
   double phase,cumphase,jd_diff;
   jd_diff = fday + JD2000 - JD0;
   cumphase = M2_PHASE0 + ((jd_diff*24.0/M2_PERIOD) - ((int) 
(jd_diff*24.0/M2_PERIOD)))*360.0;
   phase = cumphase - ((int) (cumphase/360.0))*360.0;
   phase = phase * M_PI / 180;
   return phase;
}

double lunarphase_n2(time_t t)
{
  double fday = time2fday(t);
double phase,cumphase,jd_diff;
   jd_diff = fday + JD2000 - JD0;
   cumphase = N2_PHASE0 + ((jd_diff*24.0/N2_PERIOD) - ((int) 
(jd_diff*24.0/N2_PERIOD)))*360.0;
   phase = cumphase - ((int) (cumphase/360.0))*360.0;
   phase = phase * M_PI / 180;
   return phase;
}

/*
lunartime_m2_r()
  Compute the local lunar time for a given UT and longitude

Return: local lunar time in radians (0 to 2pi)
*/

double
lunartime_m2_r(time_t t, double longitude)
{
  double m2 = lunarphase_m2(t);
  double lt = get_localtime2(t, longitude);
  double tau;

  /*tau = M_PI - (longitude - m2 - M_PI);*/
  tau = longitude - m2 - M_PI;

  while (tau < 0.0)
    tau += 2.0 * M_PI;
  while (tau > 2.0 * M_PI)
    tau -= 2.0 * M_PI;

  return (tau);
} /* lunartime_m2_r() */
