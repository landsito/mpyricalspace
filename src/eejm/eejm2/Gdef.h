#define GEO2GEO 0
#define GEO2GSE 1
#define GSE2GEO 2
#define GEO2GSM 3
#define GSM2GEO 4
#define GEO2SM 5
#define SM2GEO 6
#define GEO2SM_NP 7
#define SM2GEO_NP 8
#define GEO2SM_SP 9
#define SM2GEO_SP 10
/* IGRF-10 */
#define G10  (-29556.8)        /* for dipole latitude */
#define G11  ( -1671.8)  
#define H11     5080.0
#define MOM  (sqrt(G10*G10 + G11*G11 + H11*H11))  
#define COSTHETA0  ( -G10/MOM )
#define SINTHETA0  ( sin(acos(-G10/MOM)) )
#define PHI0       ( atan(H11/G11) ) /* here the minus sign cancels */

/* Hermann for Polar caps NP */
#define G10_NP  (-0.98823)
#define G11_NP  (-0.02393)
#define H11_NP    0.15109
#define MOM_NP  (sqrt(G10_NP*G10_NP + G11_NP*G11_NP + H11_NP*H11_NP))  
#define COSTHETA0_NP  ( -G10_NP/MOM_NP )
#define SINTHETA0_NP  ( sin(acos(-G10_NP/MOM_NP)) )
#define PHI0_NP       ( atan(H11_NP/G11_NP) ) /* here the minus sign cancels */
/* Hermann for Polar caps SP (corrected) */
/*  #define G10_SP  (-0.96269) */
/*  #define G11_SP  (-0.159055) */
/*  #define H11_SP    0.21892 */
/* Hermann for Polar caps SP (corrected) */
#define G10_SP  (-0.970296)
#define G11_SP  (-0.09453)
#define H11_SP   0.222691
#define MOM_SP  (sqrt(G10_SP*G10_SP + G11_SP*G11_SP + H11_SP*H11_SP))  
#define COSTHETA0_SP  ( -G10_SP/MOM_SP )
#define SINTHETA0_SP  ( sin(acos(-G10_SP/MOM_SP)) )
#define PHI0_SP       ( atan(H11_SP/G11_SP) ) /* here the minus sign cancels */






