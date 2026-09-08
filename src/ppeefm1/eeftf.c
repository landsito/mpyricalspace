//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
/* The program calculate the day-time Equatorial Ionospheric Eastward Electric Field (EEF)
from the Interplanetary Soalr-wind Electric field data measured by  ACE satellite.
The input data is the 5 minute means of IEF (mV/m) data pre-processed by the OMNI group and available at
http://omniweb.gsfc.nasa.gov/form/omni_min.html . The program filters the input
data with a set of coefficients specified in the file TF.COF . These coefficents were estimated using 8 years of data from ACE satellite and JULIA RADAR
The output is predicted equatorial ionospheric eastward electric field (EEF, in mV/m) at 5 minutes interval with a time delay of 17 minutes.
Note that for ideal performance of the filter, one should feed IEF data atleast 1 hour prior to the desired start time of EEF.
Ref: Manoj, C., S. Maus, H. Lühr, and P. Alken (2008), Penetration characteristics of the interplanetary electric field to the daytime equatorial ionosphere, J. Geophys. Res., 113, A12310.
Questions ? mail manoj.c.nair@noaa.gov
Date March 9, 2009.
*/



int main()
{

	int i,j,k,na,nb,ndata=0;
	float a[50],b[50],*data, *eef, dummy;
	FILE *fp;
	float *eef_tf(float *, float *, float *, int , int , int );


	fp = fopen("TF.COF","r");
	if(!fp) {
		printf("Coef file, TF.COF not found\n");
		return 0;
		}
	else 	{

	fscanf(fp,"%d",&na);
	for(i=0;i<na;i++)
	fscanf(fp,"%f",&a[i]);

	fscanf(fp,"%d",&nb);
	for(i=0;i<nb;i++)
	fscanf(fp,"%f",&b[i]);

	fclose(fp);
		}
	
	fp = fopen("ief.dat","r");

	if (!fp)
	{
	perror ("Error opening data file, ief.dat");
	return 0;
	}
  	else
  	{
	while (!feof(fp)) {
	  fscanf(fp,"%f",&dummy);
	  ndata++;
	}
	fclose (fp);
	
  	}

	fp = fopen("ief.dat","r");

	if (!fp)
	{
	perror ("Error opening data file");
	return 0;
	}
  	else
  	{

	data = (float*) calloc (ndata,  sizeof(float));
	eef =  (float*) calloc (ndata,  sizeof(float));

	for(i=0;i<ndata;i++)
	fscanf(fp,"%f\n",&data[i]);
 
	fclose(fp);
	}
	

	eef =  eef_tf(b,a,data,na,nb,ndata);

	fp = fopen("eef_res.dat","w");
	if (!fp)
	{
	perror ("Error opening eef_res.dat for writing results\n");
	return 0;
	}
  	else
  	{
 	for(i=0;i<ndata;i++)
	fprintf(fp,"%f\n",eef[i]);
   
	fclose(fp);
	}

	return 0;
}
//---------------------------------------------------------------------------

float *eef_tf(float *b, float *a, float *ief, int na, int nb, int ndata)
{
// IIR implimentation 

	int    j,    k;
	float *eef,*dbuffer;

	dbuffer = (float*) calloc (nb,  sizeof(float));
	eef = (float*) calloc (ndata,  sizeof(float));
	na += 1;

	for( j = 0; j < ndata; j++ )
	{
		for( k = 0; k < (nb-1); k++ )
		{
			dbuffer[k] = dbuffer[k + 1];
		}
		dbuffer[nb-1] = 0.0;
		for( k = 0; k < nb; k++ )
		{
			dbuffer[k] += ief[j] * b[k];
		}
		for( k = 2; k < na; k++ )
		{
			dbuffer[k - 1] -= dbuffer[0] * a[k - 1];
		}
		eef[j] = dbuffer[0];
	}


		return eef;
}
