/* common.h
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

#ifndef INCLUDED_common_h
#define INCLUDED_common_h

#include <time.h>

/*
 * Prototypes
 */

time_t fday2timet(double fday);
double time2fday(time_t t);
double time2lunar(time_t t);
double get_season(time_t t);
double get_year(double fday);
double get_localtime(double fday, double longitude);
double get_localtime2(time_t t, double longitude);
double get_zenith(double fday, double latitude, double longitude);
double get_fday(int year, double longitude, double local_time,
                int season);
int doy2md(int year, int doy, int *month, int *day);

#ifdef USE_CYGWIN
extern int putenv(const char *string);
#else
extern int putenv(char *string);
#endif

#endif /* INCLUDED_common_h */
