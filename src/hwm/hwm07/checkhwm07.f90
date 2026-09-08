!
!  Test driver for HWM07 subroutines
!  The output of the program is given at the end of the file
!
!  AUTHOR
!    John Emmert, msishwmhelp@nrl.navy.mil
!
!  DATE
!    19 August 2008
!

!******************************************************************************

program checkhwm07

  implicit none
  INTEGER      :: IYD
  REAL(4)      :: SEC,ALT,GLAT,GLON,STL,F107A,F107,AP(2)
  REAL(4)      :: W(2), QW(2), DW(2)
  real(4)      :: mlt, mlat, kp, mmpwind, mzpwind
  real(4)      :: ut, apqt(2)
  integer      :: day
  real(4)      :: pershift
  integer      :: ialt,istl,ilat,ilon,iaptemp
  integer      :: imlat,imlt,ikp

! HEIGHT PROFILE
  day = 150
  iyd = 95000 + day
  ut = 12.0
  sec = ut * 3600.0
  glat = -45.0
  glon = -85.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 80.0
  apqt(2) = -1.0
  print *, 'HEIGHT PROFILE'
  print '(a5,i3, a5,f4.1, a7,f5.1, a7,f6.1, a6,f4.1, a5,f5.1)', &
            'DAY=',day, ', UT=',ut, ', GLAT=',glat,  &
            ', GLON=',glon, ', STL=',stl, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'ALT', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do ialt = 0, 400, 25
    alt = float(ialt)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.0,3(f12.3,f10.3))', alt, qw, dw, w
  end do
  print *
  print *

! LATITUDE PROFILE
  day = 305
  iyd = 95000 + day
  ut = 18.0
  sec = ut * 3600.0
  alt = 250.0
  glon = 30.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 48.0
  apqt(2) = -1.0
  print *, 'LATITUDE PROFILE'
  print '(a5,i3, a5,f4.1, a6,f5.1, a7,f6.1, a6,f4.1, a5,f5.1)', &
            'DAY=',day, ', UT=',ut, ', ALT=',alt,  &
            ', GLON=',glon, ', STL=',stl, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'GLAT', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do ilat = -90, 90, 10
    glat = float(ilat)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.1,3(f12.3,f10.3))', glat, qw, dw, w
  end do
  print *
  print *

! LOCAL TIME PROFILE
  day = 75
  iyd = 95000 + day
  alt = 125.0
  glat = 45.0
  glon = -70.0
  ap(2) = 30.0
  apqt(2) = -1.0
  print *, 'LOCAL TIME PROFILE'
  print '(a5,i3, a6,f5.1, a7,f5.1, a7,f6.1, a5,f5.1)', &
            'DAY=',day, ', ALT=',alt, ', GLAT=',glat,  &
            ', GLON=',glon, ', ap=',ap(2)
  print '(5x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a5,3(a12,a10))', 'STL', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do istl = 0,16
    stl = 1.5 * float(istl)
    sec = (stl - glon/15.0) * 3600.0
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f5.1,3(f12.3,f10.3))', stl, qw, dw, w
  end do
  print *
  print *

! LONGITUDE PROFILE
  day = 330
  iyd = 95000 + day
  ut = 6.0
  sec = ut * 3600.0
  alt = 40.0
  glat = -5.0
  ap(2) = 4.0
  apqt(2) = -1.0
  print *, 'LONGITUDE PROFILE'
  print '(a5,i3, a5,f4.1, a6,f5.1, a7,f5.1, a7,f6.1, a5,f5.1)', &
            'DAY=',day, ', UT=',ut, ', ALT=',alt, ', GLAT=',glat,  &
            ', GLON=',glon, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'GLON', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do ilon = -180, 180, 20
    glon = float(ilon)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.0,3(f12.3,f10.3))', glon, qw, dw, w
  end do
  print *
  print *

! DAY OF YEAR PROFILE
  ut = 21.0
  sec = ut * 3600.0
  alt = 200.0
  glat = -65.0
  glon = -135.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 15.0
  apqt(2) = -1.0
  print *, 'DAY OF YEAR PROFILE'
  print '(a4,f4.1, a6,f5.1, a7,f5.1, a7,f6.1, a6,f4.1, a5,f5.1)', &
            'UT=',ut, ', ALT=',alt, ', GLAT=',glat,  &
            ', GLON=',glon, ', STL=',stl, ', ap=',ap(2)
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'DAY', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do day = 0, 360, 20
    iyd = 95000 + day
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(i6,3(f12.3,f10.3))', day, qw, dw, w
  end do
  print *
  print *

! AP PROFILE
  day = 280
  iyd = 95000 + day
  ut = 21.0
  sec = ut * 3600.0
  alt = 350.0
  glat = 38.0
  glon = 125.0
  stl = pershift(ut + glon/15.0, (/0.0, 24.0/) )
  ap(2) = 48.0
  apqt(2) = -1.0
  print *, 'MAGNETIC ACTIVITY PROFILE'
  print '(a5,i3, a5,f4.1, a6,f5.1, a7,f5.1, a7,f6.1, a6,f4.1)', &
            'DAY=',day, ', UT=',ut, ', ALT=',alt,  &
            ', GLAT=',glat, ', GLON=',glon, ', STL=',stl
  print '(6x,3a22)', 'QUIET', 'DISTURBED', 'TOTAL'
  print '(a6,3(a12,a10))', 'ap', 'MER','ZON', 'MER','ZON', 'MER','ZON'
  do iaptemp = 0, 260, 20
    ap(2) = float(iaptemp)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,apqt,qw)
    call DWM07b_HWM_interface(iyd,sec,alt,glat,glon,ap,dw)
    call HWM07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
    print '(f6.1,3(f12.3,f10.3))', ap(2), qw, dw, w
  end do
  print *
  print *

! DWM: MLAT PROFILE
  kp = 6.0
  mlt = 3.0
  print *, 'DWM: MAGNETIC LATITUDE PROFILE'
  print '(a5,f4.1, a5,f3.1)', 'MLT=',mlt, ', Kp=',kp
  print '(a6,a12,a10)', 'MLAT', 'MAG MER', 'MAG ZON'
  do imlat = -90, 90, 10
    mlat = float(imlat)
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)
    print '(f6.1,f12.3,f10.3)', mlat, mmpwind, mzpwind
  end do
  print *
  print *

! DWM: MLT PROFILE
  kp = 6.0
  mlat = 45.0
  print *, 'DWM: MAGNETIC LOCAL TIME PROFILE'
  print '(a6,f5.1, a5,f3.1)', 'MLAT=',mlat, ', Kp=',kp
  print '(a6,a12,a10)', 'MLT', 'MAG MER', 'MAG ZON'
  do imlt = 0, 16
    mlt = float(imlt)*1.5
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)
    print '(f6.1,f12.3,f10.3)', mlt, mmpwind, mzpwind
  end do
  print *
  print *

! DWM: KP PROFILE
  mlat = -50.0
  mlt = 3.0
  print *, 'DWM: Kp PROFILE'
  print '(a6,f5.1, a6,f4.1)', 'MLAT=',mlat, ', MLT=',mlt
  print '(a6,a12,a10)', 'Kp', 'MAG MER', 'MAG ZON'
  do ikp = 0, 18
    kp = float(ikp)*0.5
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)
    print '(f6.1,f12.3,f10.3)', kp, mmpwind, mzpwind
  end do
  print *
  print *

end program checkhwm07

!******************************************************************************
!
!PERSHIFT
!JOHN EMMERT   9/12/03
!TRANSLATED TO FORTRAN-90 10/4/06. FORTRAN VERSION ONLY ALLOWS SCALAR INPUTS
!SHIFTS INPUT VALUES INTO A SPECIFIED PERIODIC INTERVAL
!
!CALLING SEQUENCE:   Result = PERSHIFT(x, range)
!
!ARGUMENTS
!      x:        The value to be shifted
!      perint:   2-element vector containing the start and end values
!                of the desired periodic interval.  The periodicity is
!                determined by the span of the range.
!
!ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
!      None

function pershift(x, perint)

  real(4), parameter :: tol=1e-4
  real(4)            :: x, perint(0:1)
  real(4)            :: a, span, offset, offset1, pershift

  pershift = x
  a = perint(0)
  span = perint(1) - perint(0)
  if (span .ne. 0) then
    offset = x-a
    offset1 = mod(offset,span)
    if (abs(offset1) .lt. tol) offset1 = 0
  endif
  pershift = a + offset1
  if ((offset .lt. 0) .and. (offset1 .ne. 0)) pershift = pershift + span

  return

end function pershift


!******************************************************************************
!  TEST OUTPUT
!******************************************************************************

! HEIGHT PROFILE
! DAY=150, UT=12.0, GLAT=-45.0, GLON= -85.0, STL= 6.3, ap= 80.0
!                       QUIET             DISTURBED                 TOTAL
!   ALT         MER       ZON         MER       ZON         MER       ZON
!    0.       0.484     7.859       0.000     0.000       0.484     7.859
!   25.       2.502    24.992       0.000     0.000       2.502    24.992
!   50.      -7.478    98.432       0.000     0.000      -7.478    98.432
!   75.      15.985    49.706       0.002    -0.001      15.987    49.705
!  100.     -38.402    45.775       0.298    -0.127     -38.104    45.648
!  125.      36.199   -27.612      22.279    -9.483      58.478   -37.095
!  150.      -3.955   -32.830      44.260   -18.838      40.304   -51.668
!  175.     -23.771   -38.211      44.556   -18.964      20.784   -57.175
!  200.     -11.783   -36.996      44.558   -18.965      32.775   -55.961
!  225.      -0.942   -38.801      44.558   -18.965      43.615   -57.766
!  250.       4.802   -42.668      44.558   -18.965      49.360   -61.633
!  275.       8.318   -45.691      44.558   -18.965      52.875   -64.656
!  300.      10.635   -47.684      44.558   -18.965      55.193   -66.649
!  325.      12.163   -48.998      44.558   -18.965      56.721   -67.963
!  350.      13.170   -49.864      44.558   -18.965      57.728   -68.830
!  375.      13.834   -50.435      44.558   -18.965      58.392   -69.401
!  400.      14.272   -50.812      44.558   -18.965      58.830   -69.777
! 
! 
! LATITUDE PROFILE
! DAY=305, UT=18.0, ALT=250.0, GLON=  30.0, STL=20.0, ap= 48.0
!                       QUIET             DISTURBED                 TOTAL
!  GLAT         MER       ZON         MER       ZON         MER       ZON
! -90.0     -46.412   110.719     -86.773    46.502    -133.186   157.221
! -80.0     -26.264    68.190    -133.904   -12.215    -160.169    55.974
! -70.0       9.717    28.536    -137.274  -112.417    -127.557   -83.881
! -60.0      32.894    11.525     -58.808  -158.941     -25.914  -147.415
! -50.0      32.972    20.006     -16.136  -111.817      16.836   -91.811
! -40.0      25.222    44.460     -10.932   -51.150      14.290    -6.690
! -30.0      26.510    67.704       7.971   -36.298      34.480    31.406
! -20.0      34.481    75.179      15.092   -34.427      49.574    40.753
! -10.0      35.540    68.603      11.120   -26.603      46.661    42.000
!   0.0      26.231    64.273       5.821   -16.657      32.051    47.616
!  10.0      18.178    72.452       0.659   -13.275      18.837    59.176
!  20.0      22.910    84.147      -5.329   -20.421      17.581    63.726
!  30.0      37.805    84.437     -12.740   -34.467      25.065    49.970
!  40.0      47.879    72.691     -15.446   -41.985      32.433    30.705
!  50.0      38.036    60.953      -6.688   -54.782      31.348     6.171
!  60.0       4.256    57.496      -7.001  -202.013      -2.745  -144.517
!  70.0     -42.170    63.602     -36.510   -52.515     -78.681    11.087
!  80.0     -79.970    83.380    -103.542    81.980    -183.511   165.360
!  90.0     -92.092   120.221     -75.838    90.679    -167.930   210.901
! 
! 
! LOCAL TIME PROFILE
! DAY= 75, ALT=125.0, GLAT= 45.0, GLON= -70.0, ap= 30.0
!                      QUIET             DISTURBED                 TOTAL
!  STL         MER       ZON         MER       ZON         MER       ZON
!  0.0     -25.523   -25.139      -8.475   -20.072     -33.998   -45.211
!  1.5     -29.545   -53.409      -7.022     2.427     -36.567   -50.981
!  3.0     -21.723   -65.269      -5.600    20.202     -27.324   -45.067
!  4.5      -8.261   -47.142      -8.034    21.516     -16.295   -25.627
!  6.0      -0.534   -15.152     -13.708     9.066     -14.242    -6.086
!  7.5       1.132     2.634     -16.795    -1.433     -15.663     1.200
!  9.0       8.180    -2.976     -16.107    -4.816      -7.926    -7.792
! 10.5      27.106   -14.648     -13.941    -4.061      13.165   -18.709
! 12.0      46.642   -11.904     -11.480    -1.964      35.162   -13.868
! 13.5      46.855     3.871      -8.243     0.139      38.612     4.010
! 15.0      21.394    12.923      -4.710    -0.930      16.684    11.993
! 16.5     -13.202     3.319      -3.931   -13.190     -17.134    -9.871
! 18.0     -33.181   -13.810      -8.103   -40.895     -41.284   -54.704
! 19.5     -31.037   -19.100     -13.157   -64.857     -44.194   -83.957
! 21.0     -20.447   -10.681     -12.756   -60.282     -33.203   -70.963
! 22.5     -18.238    -7.532     -10.061   -41.102     -28.299   -48.634
! 24.0     -25.523   -25.139      -8.473   -20.030     -33.996   -45.169
! 
! 
! LONGITUDE PROFILE
! DAY=330, UT= 6.0, ALT= 40.0, GLAT= -5.0, GLON= -70.0, ap=  4.0
!                       QUIET             DISTURBED                 TOTAL
!  GLON         MER       ZON         MER       ZON         MER       ZON
! -180.       0.125   -15.666       0.000     0.000       0.125   -15.666
! -160.       0.669   -18.202       0.000     0.000       0.669   -18.202
! -140.       0.871   -20.365       0.000     0.000       0.871   -20.365
! -120.       0.679   -21.247       0.000     0.000       0.679   -21.247
! -100.       0.227   -20.376       0.000     0.000       0.227   -20.376
!  -80.      -0.229   -17.950       0.000     0.000      -0.229   -17.950
!  -60.      -0.440   -14.761       0.000     0.000      -0.440   -14.761
!  -40.      -0.284   -11.873       0.000     0.000      -0.284   -11.873
!  -20.       0.175   -10.172       0.000     0.000       0.175   -10.172
!    0.       0.716   -10.009       0.000     0.000       0.716   -10.009
!   20.       1.063   -11.088       0.000     0.000       1.063   -11.088
!   40.       1.020   -12.651       0.000     0.000       1.020   -12.651
!   60.       0.564   -13.864       0.000     0.000       0.564   -13.864
!   80.      -0.137   -14.216       0.000     0.000      -0.137   -14.216
!  100.      -0.798   -13.754       0.000     0.000      -0.798   -13.754
!  120.      -1.146   -13.036       0.000     0.000      -1.146   -13.036
!  140.      -1.042   -12.827       0.000     0.000      -1.042   -12.827
!  160.      -0.543   -13.690       0.000     0.000      -0.543   -13.690
!  180.       0.125   -15.666       0.000     0.000       0.125   -15.666
! 
! 
! DAY OF YEAR PROFILE
! UT=21.0, ALT=200.0, GLAT=-65.0, GLON=-135.0, STL=12.0, ap= 15.0
!                       QUIET             DISTURBED                 TOTAL
!   DAY         MER       ZON         MER       ZON         MER       ZON
!     0      -1.953   -50.220       4.826     9.273       2.873   -40.947
!    20      -3.100   -45.982       4.827     9.270       1.727   -36.712
!    40      -8.144   -44.454       4.828     9.264      -3.316   -35.191
!    60     -16.767   -43.867       4.828     9.262     -11.939   -34.605
!    80     -28.168   -41.056       4.823     9.288     -23.345   -31.769
!   100     -41.057   -33.308       4.805     9.366     -36.252   -23.943
!   120     -53.790   -20.060       4.778     9.495     -49.012   -10.565
!   140     -64.620    -3.635       4.753     9.628     -59.867     5.993
!   160     -72.017    11.376       4.740     9.708     -67.278    21.084
!   180     -74.968    19.873       4.740     9.704     -70.228    29.577
!   200     -73.161    18.359       4.755     9.618     -68.406    27.977
!   220     -67.005     6.490       4.781     9.482     -62.224    15.973
!   240     -57.496   -12.646       4.807     9.356     -52.688    -3.290
!   260     -45.970   -33.666       4.824     9.284     -41.146   -24.382
!   280     -33.850   -50.949       4.829     9.262     -29.022   -41.687
!   300     -22.446   -60.734       4.828     9.264     -17.618   -51.470
!   320     -12.849   -62.349       4.827     9.270      -8.022   -53.079
!   340      -5.917   -58.075       4.826     9.273      -1.091   -48.802
!   360      -2.295   -51.752       4.826     9.273       2.531   -42.479
! 
! 
! MAGNETIC ACTIVITY PROFILE
! DAY=280, UT=21.0, ALT=350.0, GLAT= 38.0, GLON= 125.0, STL= 5.3
!                       QUIET             DISTURBED                 TOTAL
!    ap         MER       ZON         MER       ZON         MER       ZON
!   0.0      21.580   -43.920      -1.453    -5.318      20.127   -49.238
!  20.0      21.580   -43.920      -9.057    -2.548      12.523   -46.468
!  40.0      21.580   -43.920     -20.477   -17.007       1.103   -60.927
!  60.0      21.580   -43.920     -29.204   -30.416      -7.624   -74.336
!  80.0      21.580   -43.920     -34.047   -38.327     -12.467   -82.248
! 100.0      21.580   -43.920     -36.950   -43.408     -15.371   -87.328
! 120.0      21.580   -43.920     -38.558   -46.501     -16.979   -90.422
! 140.0      21.580   -43.920     -39.463   -48.508     -17.883   -92.428
! 160.0      21.580   -43.920     -39.902   -49.810     -18.323   -93.730
! 180.0      21.580   -43.920     -39.969   -50.507     -18.389   -94.427
! 200.0      21.580   -43.920     -39.767   -50.737     -18.188   -94.657
! 220.0      21.580   -43.920     -39.640   -50.729     -18.060   -94.649
! 240.0      21.580   -43.920     -39.640   -50.729     -18.060   -94.649
! 260.0      21.580   -43.920     -39.640   -50.729     -18.060   -94.649
! 
! 
! DWM: MAGNETIC LATITUDE PROFILE
! MLT= 3.0, Kp=6.0
!  MLAT     MAG MER   MAG ZON
! -90.0     157.613  -158.954
! -80.0     158.138  -141.865
! -70.0      56.293   -41.271
! -60.0      60.878    69.305
! -50.0     121.506    -5.004
! -40.0      40.882   -15.775
! -30.0      21.850   -25.517
! -20.0      11.640   -34.856
! -10.0       5.533   -40.228
!   0.0       5.327   -45.354
!  10.0       8.204   -50.906
!  20.0       9.422   -55.711
!  30.0       3.494   -56.017
!  40.0      -7.319   -49.619
!  50.0      -7.873    18.973
!  60.0     -33.112    80.795
!  70.0     -20.227    -8.189
!  80.0     -55.275   -64.783
!  90.0    -137.149  -170.526
! 
! 
! DWM: MAGNETIC LOCAL TIME PROFILE
! MLAT= 45.0, Kp=6.0
!   MLT     MAG MER   MAG ZON
!   0.0     -22.332  -129.270
!   1.5     -28.861   -78.026
!   3.0      -6.246   -33.077
!   4.5      -5.592   -33.656
!   6.0     -32.813   -39.796
!   7.5     -49.672   -34.709
!   9.0     -48.676   -26.754
!  10.5     -41.507   -20.200
!  12.0     -35.088   -15.157
!  13.5     -28.140   -12.072
!  15.0     -18.344   -13.870
!  16.5      -9.103   -27.605
!  18.0     -12.015   -57.829
!  19.5     -27.512   -86.677
!  21.0     -24.814  -103.292
!  22.5     -11.058  -127.315
!  24.0     -22.332  -129.270
! 
! 
! DWM: Kp PROFILE
! MLAT=-50.0, MLT= 3.0
!    Kp     MAG MER   MAG ZON
!   0.0      -8.066     5.787
!   0.5      -6.955     5.282
!   1.0      -3.365     3.849
!   1.5       3.202     1.647
!   2.0      13.413    -1.112
!   2.5      26.175    -3.461
!   3.0      39.805    -4.712
!   3.5      54.121    -5.081
!   4.0      68.787    -4.864
!   4.5      83.302    -4.428
!   5.0      97.009    -4.194
!   5.5     109.699    -4.429
!   6.0     121.506    -5.004
!   6.5     132.195    -5.750
!   7.0     141.548    -6.483
!   7.5     149.386    -7.019
!   8.0     155.577    -7.193
!   8.5     155.577    -7.193
!   9.0     155.577    -7.193
!
