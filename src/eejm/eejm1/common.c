/* common.c
 * 
 * Copyright (C) 2006, 2007 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"
#include "coord.h"
#include "Gdef.h"

#include <gsl/gsl_math.h>

/*
 * time stamp of 01-Jan-2000 00:00:00 UTC
 *
 * command used: date -d "1/1/2000 00:00:00 UTC" +"%s"
 */
static time_t t0_fday = 946684800;

static char tz_str[] = "TZ=GMT";

time_t
fday2timet(double fday)
{
  time_t t;

  /* find time in seconds since 1-Jan-2000 00:00:00 UTC */
  t = (time_t) (fday * 86400.0);

  /* add time to the epoch 1-Jan-1970 00:00:00 UTC */
  t += t0_fday;

  return (t);
} /* fday2timet() */

double
time2fday(time_t t)
{
  double fday;

  fday = (double) (t - t0_fday);
  fday /= 86400.0;

  return (fday);
} /* time2fday() */

/*
get_season()
  Determine the day of year for a given fday

Inputs: t - timestamp in seconds since 1-Jan-1970 00:00:00 UTC

Return: floating point days since Jan 1, in the range 0 to 365
*/

double
get_season(time_t t)
{
  struct tm *tm_p;
  double days;

  tm_p = gmtime(&t);

  days = (double) tm_p->tm_yday;
  days += tm_p->tm_hour / 24.0 +
          tm_p->tm_min / 1440.0 +
          tm_p->tm_sec / 86400.0;

  return (days);
} /* get_season() */

/*
get_year()
  Determine the year for a given fday

Inputs: fday - floating point days since 1-Jan-2000 00:00:00 UTC

Return: floating point year
*/

double
get_year(double fday)
{
  time_t t;             /* fday in seconds */
  struct tm *tm_p;
  double year;
  double days;

  /* find time in seconds since 1-Jan-2000 00:00:00 UTC */
  t = (time_t) (fday * 86400.0);

  /* add time to the epoch 1-Jan-1970 00:00:00 UTC */
  t += t0_fday;

  /*
   * 't' now holds the representation in seconds of 'fday' since
   * the epoch - find the day of the year from the corresponding
   * tm struct
   */
  tm_p = gmtime(&t);

  days = (double) tm_p->tm_yday;
  days += tm_p->tm_hour / 24.0 +
          tm_p->tm_min / 1440.0 +
          tm_p->tm_sec / 86400.0;

  year = (double) tm_p->tm_year + 1900.0 + days / 365.0;

  return (year);
} /* get_year() */

/*
get_localtime()
  Calculate the localtime for a given fday

Inputs: t         - timestamp
        longitude - longitude in radians of where to compute the
                    local time

Return: local time in floating point hours
*/

double
get_localtime2(time_t t, double longitude)
{
  struct tm *tmp;
  double lt;

  /* convert utc to local time by adding longitude offset */
  lt = (double) t + longitude / (15.0 * M_PI / 180.0 / 3600.0);

  t = (time_t) lt;

  tmp = gmtime(&t);

  lt = (double) tmp->tm_hour +
       (double) tmp->tm_min / 60.0 +
       (double) tmp->tm_sec / 3600.0;

  return (lt);
} /* get_localtime() */

double
get_localtime(double fday, double longitude)
{
  time_t t;
  struct tm *tmp;
  double lt;

  /* find time in seconds since the epoch 1-Jan-1970 00:00:00 UTC */
  t = (time_t) (fday * 86400.0) + t0_fday;

  /* convert utc to local time by adding longitude offset */
  lt = (double) t + longitude / (15.0 * M_PI / 180.0 / 3600.0);

  t = (time_t) lt;

  tmp = gmtime(&t);

  lt = (double) tmp->tm_hour +
       (double) tmp->tm_min / 60.0 +
       (double) tmp->tm_sec / 3600.0;

  return (lt);
} /* get_localtime() */

/*
get_zenith()
  Calculate the solar zenith angle for given fday, latitude
and longitude

Inputs: fday      - floating point days since 1-Jan-2000 00:00:00 UTC
        latitude  - latitude (radians)
        longitude - longitude (radians)

Return: zenith angle in radians
*/

double
get_zenith(double fday, double latitude, double longitude)
{
  double gse_lon, gse_lat;
  double sun_lon = 0.0;
  double sun_lat = 0.0;
  double delta;
  double az;

  trans(GEO2GSE,
        fday,
        longitude,
        latitude,
        &gse_lon,
        &gse_lat);
  gse_lon *= 180.0 / M_PI;
  gse_lat *= 180.0 / M_PI;
  my_delaz(&gse_lat, &gse_lon, &sun_lat, &sun_lon, &delta, &az);

  return (delta * M_PI / 180.0);
} /* get_zenith() */

/*
get_fday()
  Convert parameterized time to fday

Inputs: year       - year
        longitude  - longitude in radians
        local_time - local time in hours
        season     - day of year (0 - 365)
*/

double
get_fday(int year, double longitude, double local_time, int season)
{
  struct tm tminfo;
  time_t offset;
  double fday;
  double utc;
  double lt = local_time * 3600.0;
  double s = season * 86400.0;

  /* use GMT time */
  putenv(tz_str);

  /*
   * get time offset since epoch of 00:00:00 UTC of the specified year
   */
  tminfo.tm_sec = 0;
  tminfo.tm_min = 0;
  tminfo.tm_hour = 0;
  tminfo.tm_mday = 1;
  tminfo.tm_mon = 0;
  tminfo.tm_year = year - 1900;
  tminfo.tm_isdst = 0;

  offset = mktime(&tminfo);

  /* get offset since 01-Jan-2000 00:00:00 UTC and convert to days */
  fday = (double) (offset - t0_fday);

  /* convert local time to utc by adding longitude offset */
  utc = lt - longitude / (15.0 * M_PI / 180.0 / 3600.0);

  /* add day of year and time */
  fday += s;
  fday += utc;

  fday /= 86400.0;
  
  return (fday);
} /* get_fday() */
