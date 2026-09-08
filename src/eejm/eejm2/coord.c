/* The GEO,GSE and GSM routines have been tested to be correct
   to about 0.1 degrees */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "Gdef.h"

#define r_erde          6371.2

#define WGS84_A 6378.1370       /* in km */
#define WGS84_B 6356.752314     /* in km */
#define WGS84_E sqrt(WGS84_A*WGS84_A - WGS84_B*WGS84_B)
#define WGS84_E2 (WGS84_E  * WGS84_E)
#define WGS84_E4 (WGS84_E2 * WGS84_E2)

double vec_norm(double x, double y, double z){ return sqrt(x*x + y*y + z*z); }
double sqr(double x){ return x*x; }
int imin(int x,int y){ if (x < y) return x; else return y; }
int imax(int x,int y){ if (x > y) return x; else return y; }
double dmin(double x,double y){ if (x < y) return x; else return y; }
double dmax(double x,double y){ if (x > y) return x; else return y; }
int iround(double x){ return (int) (x+0.5); }

double gps2fday(double gps){return (gps-630720013.0) / 86400.0;}
double fday2gps(double fday){return fday*86400.0 + 630720013.0;}

int indx(int n, int m)
{
  int i;

  if (m == 0) {
    i = (n - 1) * (n + 1) + 1;
    return i-1;
  }
  if (m >= 0)
    i = (n - 1) * (n + 1) + m * 2;
  else
    i = (n - 1) * (n + 1) - m * 2 + 1;
  return i-1;
}  /*indx*/


int fday2season(double fday)    /* 1=spring, 2=summer, 3=fall, 4=winter */
{
  double mday, test, rest;
  int season;

  mday = fday + 36525;          /* add 100 years to avoid negative days*/
  test = (mday - 35.0) / 365.25;
  rest = test - floor(test);
  season = (int) (rest*4.0) + 1;
  return season;
}


void rotate_clockwise(double alfa_rad,double bx,double by,double *bxs,double *bys)
{
  double nx,ny;
  nx = cos(alfa_rad)*bx - sin(alfa_rad)*by;
  ny = sin(alfa_rad)*bx + cos(alfa_rad)*by;
  *bxs = nx;
  *bys = ny;
}


void star_camera_correction(double az_rad,double  scx_deg,double  scy_deg,double  scz_deg,
                            double  x,double  y,double  z,
                            double  *sx,double  *sy,double  *sz)
{
  double cx, cy, cz, nx, ny, nz;

  rotate_clockwise(az_rad, x, y, &cx, &cy);
  cz = z;
  rotate_clockwise(scx_deg * M_PI / 180.0, cy, cz, &ny, &nz);
  nx = cx;
  rotate_clockwise(scy_deg * M_PI / 180.0, nz, nx, &nz, &nx);
  rotate_clockwise(scz_deg * M_PI / 180.0, nx, ny, &nx, &ny);
  rotate_clockwise(-az_rad, nx, ny, &nx, &ny);
  *sx = nx - x;
  *sy = ny - y;
  *sz = nz - z;
}


void my_delaz(double *lat1,double *lon1,double *lat2,double *lon2,double *delta,double *az)
{

  /* output: 0 <= az < 360.0 */
  /* tested to be correct at high latitudes */

  double theta1,theta2,dlon,dlat,arg;

  theta1 = (90.0-*lat1) * M_PI/180.0;
  theta2 = (90.0-*lat2) * M_PI/180.0;
  dlon = (*lon2-*lon1) * M_PI/180.0;
  if (dlon >  M_PI) dlon -= 2*M_PI;
  if (dlon < -M_PI) dlon += 2*M_PI;

  dlat = (*lat2-*lat1) * M_PI/180.0;

  arg = sin(theta1)*sin(theta2)*cos(dlon) + cos(theta1)*cos(theta2);
  arg = dmax(-1.0,dmin(1.0,arg)); 
  *delta = acos(arg) * 180.0/M_PI;

  if (fabs(dlat) > 1.0e-10)
    {
      *az = atan(tan(sin(theta2)*dlon)/sin(dlat)) * 180.0/M_PI;
      if (dlat < 0) *az += 180.0;
    }
  else 
    {
      if (dlon > 0) *az = 90.0; 
      else *az = 270.0; 
    } 
  if (*az < 0) *az += 360.0;
}

void cart2sphere(double x,double  y,double  z,double *phi_rad,double *lat_rad,double *r)
{
  double d;

  *r = sqrt(x * x + y * y + z * z);
  d = z / *r;
  *lat_rad = asin(d);
  if (x > 0)
    *phi_rad = atan(y / x);
  if (x < 0 && y >= 0)
    *phi_rad = atan(y / x) + M_PI;
  if (x < 0 && y < 0)
    *phi_rad = atan(y / x) - M_PI;
  if (fabs(x) < 1.0e-10 && y >= 0)
    *phi_rad = M_PI / 2;
  if (fabs(x) < 1.0e-10 && y < 0)
    *phi_rad = M_PI / -2;
  return;
}


void sphere2cart(double phi_rad,double  lat_rad,double  r,double  *x,double  *y,double  *z)
{
  *x = r * cos(phi_rad) * cos(lat_rad);
  *y = r * sin(phi_rad) * cos(lat_rad);
  *z = r * sin(lat_rad);
  return;
}


void sphere2cart_vec(double phi_rad,double lat_rad,double vlat,double vphi,double vminusr,
                     double *vx,double *vy,double *vz)
{
  double sp,st,cp,ct;

  sp = sin(phi_rad);
  cp = cos(phi_rad);
  st = sin(0.5*M_PI-lat_rad);
  ct = cos(0.5*M_PI-lat_rad);
  *vx = -sp*vphi - ct*cp*vlat - st*cp*vminusr;
  *vy =  cp*vphi - ct*sp*vlat - st*sp*vminusr;
  *vz =               st*vlat -    ct*vminusr;
} /* sphere2cart_vec */

void cart2sphere_vec(double phi_rad,double lat_rad,double vx,double vy,double vz,
                     double *vlat,double *vphi,double *vminusr)
{
  double sp,st,cp,ct;

  sp = sin(phi_rad);
  cp = cos(phi_rad);
  st = sin(0.5*M_PI-lat_rad);
  ct = cos(0.5*M_PI-lat_rad);

  *vlat =    -ct*cp*vx - ct*sp*vy + st*vz;
  *vphi =    -sp*vx + cp*vy;
  *vminusr = -st*cp*vx - st*sp*vy - ct*vz;

} /* cart2sphere_vec */


void IDENTITY_matrix(double A[3][3])
{
  A[0][0] = 1.0;
  A[1][0] = 0.0;
  A[2][0] = 0.0;

  A[0][1] = 0.0;
  A[1][1] = 1.0;
  A[2][1] = 0.0;

  A[0][2] = 0.0;
  A[1][2] = 0.0;
  A[2][2] = 1.0;
} /* IDENTITY_matrix */

void GEO2GSE_matrix(double fday,double A[3][3])
{
  double mjd,td;
  double theta,ct,st;
  double lambda,cl,sl;
  double eps,ce,se,time;

  mjd = floor(fday);
  td = (fday-mjd)*24.0;
  theta = (100.0 + 0.986*mjd + 15.04*td) * M_PI / 180.0; /* this is not the colatitude! */
  ct = cos(theta);
  st = sin(theta);
  time = 0.986*mjd + 0.04107*td;
  lambda = (279.97 + time
            + 1.915*sin((357.04+time)*M_PI/180.0)) * M_PI/180.0;
  cl = cos(lambda);
  sl = sin(lambda);
  eps = 23.44  * M_PI / 180.0;
  ce = cos(eps);
  se = sin(eps);

  A[0][0] = ct*cl + st*ce*sl;
  A[1][0] = -ct*sl + st*ce*cl;
  A[2][0] = -st*se;

  A[0][1] = -st*cl + ct*ce*sl;
  A[1][1] = st*sl + ct*ce*cl;
  A[2][1] = -ct*se;

  A[0][2] = se*sl;
  A[1][2] = se*cl;
  A[2][2] = ce;

} /* GEO2GSE_matrix */


void GSE2GEO_matrix(double fday,double A[3][3])
{
  double mjd,td;
  double theta,ct,st;
  double lambda,cl,sl;
  double eps,ce,se,time;

  mjd = floor(fday);
  td = (fday-mjd)*24.0;
  theta = (100.0 + 0.986*mjd + 15.04*td) * M_PI / 180.0; /* this is not the colatitude! */
  ct = cos(theta);
  st = sin(theta);
  time = 0.986*mjd + 0.04107*td;
    lambda = (279.97 + time
            + 1.915*sin((357.04+time)*M_PI/180.0)) * M_PI/180.0;
  cl = cos(lambda);
  sl = sin(lambda);
  eps = 23.44  * M_PI / 180.0;
  ce = cos(eps);
  se = sin(eps);

  A[0][0] = ct*cl + st*ce*sl;
  A[0][1] = -ct*sl + st*ce*cl;
  A[0][2] = -st*se;

  A[1][0] = -st*cl + ct*ce*sl;
  A[1][1] = st*sl + ct*ce*cl;
  A[1][2] = -ct*se;

  A[2][0] = se*sl;
  A[2][1] = se*cl;
  A[2][2] = ce;
} /* GSE2GEO_matrix */

#define LAMBDA (atan(H11/G11))
#define COSL cos(LAMBDA)
#define SINL sin(LAMBDA)
#define PHI ((M_PI/2.0) - asin((G11*COSL+H11*SINL)/G10))
#define XG (cos(PHI)*COSL)
#define YG (cos(PHI)*SINL)
#define ZG (sin(PHI))

#define LAMBDA_NP (atan(H11_NP/G11_NP))
#define COSL_NP cos(LAMBDA_NP)
#define SINL_NP sin(LAMBDA_NP)
#define PHI_NP ((M_PI/2.0) - asin((G11_NP*COSL_NP+H11_NP*SINL_NP)/G10_NP))
#define XG_NP (cos(PHI_NP)*COSL_NP)
#define YG_NP (cos(PHI_NP)*SINL_NP)
#define ZG_NP (sin(PHI_NP))

#define LAMBDA_SP (atan(H11_SP/G11_SP))
#define COSL_SP cos(LAMBDA_SP)
#define SINL_SP sin(LAMBDA_SP)
#define PHI_SP ((M_PI/2.0) - asin((G11_SP*COSL_SP+H11_SP*SINL_SP)/G10_SP))
#define XG_SP (cos(PHI_SP)*COSL_SP)
#define YG_SP (cos(PHI_SP)*SINL_SP)
#define ZG_SP (sin(PHI_SP))

void GEO2GSM_matrix(double fday,double B[3][3])
{ 
  double ye,ze,ct,st,theta;
  double A[3][3];

/*    printf("LAMBDA = %f, PHI = %f\n",LAMBDA*180/M_PI,PHI*180/M_PI); */

  GEO2GSE_matrix(fday,A);
  ye = A[1][0]*XG + A[1][1]*YG + A[1][2]*ZG;
  ze = A[2][0]*XG + A[2][1]*YG + A[2][2]*ZG;
  theta = atan(ye/ze);
  ct = cos(theta);
  st = sin(theta);
  B[0][0] = A[0][0];
  B[0][1] = A[0][1];
  B[0][2] = A[0][2];
  B[1][0] = ct*A[1][0] - st*A[2][0];
  B[1][1] = ct*A[1][1] - st*A[2][1];
  B[1][2] = ct*A[1][2] - st*A[2][2];
  B[2][0] = st*A[1][0] + ct*A[2][0];
  B[2][1] = st*A[1][1] + ct*A[2][1];
  B[2][2] = st*A[1][2] + ct*A[2][2];
} /* GEO2GSM_matrix */

void GSM2GEO_matrix(double fday,double B[3][3])
{ 
  double ye,ze,ct,st,theta;
  double A[3][3];

  GEO2GSE_matrix(fday,A);  /* A is now T2 T1^t */
  ye = A[1][0]*XG + A[1][1]*YG + A[1][2]*ZG;
  ze = A[2][0]*XG + A[2][1]*YG + A[2][2]*ZG;
  theta = atan(ye/ze);
  ct = cos(theta);
  st = sin(theta);

  B[0][0] = A[0][0];
  B[1][0] = A[0][1];
  B[2][0] = A[0][2];
  B[0][1] = ct*A[1][0] - st*A[2][0];
  B[1][1] = ct*A[1][1] - st*A[2][1];
  B[2][1] = ct*A[1][2] - st*A[2][2];
  B[0][2] = st*A[1][0] + ct*A[2][0];
  B[1][2] = st*A[1][1] + ct*A[2][1];
  B[2][2] = st*A[1][2] + ct*A[2][2];
} /* GSM2GEO_matrix */

void GEO2SM_matrix(double fday,double B[3][3],double xg,double yg,double zg)
{ 
  double xe,ye,ze,cm,sm,mu;
  double A[3][3];

/*    printf("LAMBDA = %f, PHI = %f\n",LAMBDA*180/M_PI,PHI*180/M_PI); */

  GEO2GSE_matrix(fday,A);
  xe = A[0][0]*xg + A[0][1]*yg + A[0][2]*zg;
  ye = A[1][0]*xg + A[1][1]*yg + A[1][2]*zg;
  ze = A[2][0]*xg + A[2][1]*yg + A[2][2]*zg;
  mu = atan(xe/sqrt(ye*ye+ze*ze));
  cm = cos(mu);
  sm = sin(mu);
  GEO2GSM_matrix(fday,A);
  B[0][0] = cm*A[0][0] - sm*A[2][0];
  B[0][1] = cm*A[0][1] - sm*A[2][1];
  B[0][2] = cm*A[0][2] - sm*A[2][2];
  B[1][0] = A[1][0];
  B[1][1] = A[1][1];
  B[1][2] = A[1][2];
  B[2][0] = sm*A[0][0] + cm*A[2][0];
  B[2][1] = sm*A[0][1] + cm*A[2][1];
  B[2][2] = sm*A[0][2] + cm*A[2][2];
} /* GEO2SM_matrix */

void SM2GEO_matrix(double fday,double B[3][3],double xg,double yg,double zg)
{ 
  double xe,ye,ze,cm,sm,mu;
  double A[3][3];

  GEO2GSE_matrix(fday,A);  /* A is now T2 T1^t */
  xe = A[0][0]*xg + A[0][1]*yg + A[0][2]*zg;
  ye = A[1][0]*xg + A[1][1]*yg + A[1][2]*zg;
  ze = A[2][0]*xg + A[2][1]*yg + A[2][2]*zg;
  mu = atan(xe/sqrt(ye*ye+ze*ze));
  cm = cos(mu);
  sm = sin(mu);
  GEO2GSM_matrix(fday,A);
  B[0][0] = cm*A[0][0] - sm*A[2][0];
  B[1][0] = cm*A[0][1] - sm*A[2][1];
  B[2][0] = cm*A[0][2] - sm*A[2][2];
  B[0][1] = A[1][0];
  B[1][1] = A[1][1];
  B[2][1] = A[1][2];
  B[0][2] = sm*A[0][0] + cm*A[2][0];
  B[1][2] = sm*A[0][1] + cm*A[2][1];
  B[2][2] = sm*A[0][2] + cm*A[2][2];

} /* SM2GEO_matrix */
void trans(int action,double fday,double lon_rad,double lat_rad,
                  double *lon_rad_s,double *lat_rad_s)
{
  double r,x,y,z,xs,ys,zs;
  double A[3][3];
  switch (action)
    {
    case GEO2GEO:  IDENTITY_matrix(A); break;
    case GEO2GSE:  GEO2GSE_matrix(fday,A); break;
    case GSE2GEO:  GSE2GEO_matrix(fday,A); break;
    case GEO2GSM:  GEO2GSM_matrix(fday,A); break;
    case GSM2GEO:  GSM2GEO_matrix(fday,A); break;
    case GEO2SM:   GEO2SM_matrix(fday,A,XG,YG,ZG); break;
    case SM2GEO:   SM2GEO_matrix(fday,A,XG,YG,ZG); break;
    case GEO2SM_NP:   GEO2SM_matrix(fday,A,XG_NP,YG_NP,ZG_NP); break;
    case SM2GEO_NP:   SM2GEO_matrix(fday,A,XG_NP,YG_NP,ZG_NP); break;
    case GEO2SM_SP:   GEO2SM_matrix(fday,A,XG_SP,YG_SP,ZG_SP); break;
    case SM2GEO_SP:   SM2GEO_matrix(fday,A,XG_SP,YG_SP,ZG_SP); break;
    }
  sphere2cart(lon_rad,lat_rad,1.0,&x,&y,&z);
  xs = A[0][0]*x + A[0][1]*y + A[0][2]*z;
  ys = A[1][0]*x + A[1][1]*y + A[1][2]*z;
  zs = A[2][0]*x + A[2][1]*y + A[2][2]*z;
  cart2sphere(xs,ys,zs,lon_rad_s,lat_rad_s,&r);
} /* trans */


void trans_vec(int action,double fday,double lon_rad,double lat_rad,
                      double vlat,double vphi,double vminusr,
                      double *lon_rad_s,double *lat_rad_s,double *vlat_s,double *vphi_s,double *vminusr_s)
{
  double r,x,y,z,xs,ys,zs,vx,vy,vz,vxs,vys,vzs;
  double A[3][3];
  switch (action)
    {
    case GEO2GEO:  IDENTITY_matrix(A); break;
    case GEO2GSE:  GEO2GSE_matrix(fday,A); break;
    case GSE2GEO:  GSE2GEO_matrix(fday,A); break;
    case GEO2GSM:  GEO2GSM_matrix(fday,A); break;
    case GSM2GEO:  GSM2GEO_matrix(fday,A); break;
    case GEO2SM:   GEO2SM_matrix(fday,A,XG,YG,ZG); break;
    case SM2GEO:   SM2GEO_matrix(fday,A,XG,YG,ZG); break;
    case GEO2SM_NP:   GEO2SM_matrix(fday,A,XG_NP,YG_NP,ZG_NP); break;
    case SM2GEO_NP:   SM2GEO_matrix(fday,A,XG_NP,YG_NP,ZG_NP); break;
    case GEO2SM_SP:   GEO2SM_matrix(fday,A,XG_SP,YG_SP,ZG_SP); break;
    case SM2GEO_SP:   SM2GEO_matrix(fday,A,XG_SP,YG_SP,ZG_SP); break;
    }

  sphere2cart(lon_rad,lat_rad,1.0,&x,&y,&z);
  sphere2cart_vec(lon_rad,lat_rad,vlat,vphi,vminusr,&vx,&vy,&vz);
  xs =  A[0][0]*x +  A[0][1]*y +  A[0][2]*z;
  ys =  A[1][0]*x +  A[1][1]*y +  A[1][2]*z;
  zs =  A[2][0]*x +  A[2][1]*y +  A[2][2]*z;
  vxs = A[0][0]*vx + A[0][1]*vy + A[0][2]*vz;
  vys = A[1][0]*vx + A[1][1]*vy + A[1][2]*vz;
  vzs = A[2][0]*vx + A[2][1]*vy + A[2][2]*vz;
  cart2sphere(xs,ys,zs,lon_rad_s,lat_rad_s,&r);
  cart2sphere_vec(*lon_rad_s,*lat_rad_s,vxs,vys,vzs,vlat_s,vphi_s,vminusr_s);

} /* trans_vec*/
