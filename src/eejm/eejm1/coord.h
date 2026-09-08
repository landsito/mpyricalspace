/*
 * coord.h
 * Patrick Alken
 */

#ifndef INCLUDED_coord_h
#define INCLUDED_coord_h

/*
 * Prototypes
 */

void trans(int action, double fday, double lon_rad, double lat_rad,
           double *long_rad_s, double *lat_rad_s);
void my_delaz(double *lat1, double *lon1, double *lat2, double *lon2,
              double *delta, double *az);

#endif /* INCLUDED_coord_h */
