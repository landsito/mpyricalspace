!
!  Test script for HL-TWiM
!  As a reference, the output of this test script is pasted at the end of this file.
!
!  AUTHORS:
!    Manbharat Singh Dhadly
!    John Emmert
!    Douglas Drob
!    ----------------------
!    Geospace Science and Technology Branch
!    Space Science Division
!    Naval Research Laboratory
!    4555 Overlook Ave.
!    Washington, DC 20375
!
!  Point of Contact
!   manbharat.dhadly.ctr@nrl.navy.mil
!
!   DATE: 05-June-2019
!******************************************
!  Model Call
!   hltwim(day, ut, glat, glon, kp, w, mw)
!
! Input arguments:
!        day - day of year
!        ut  - universal time
!        glat - geodetic latitude(deg)
!        glon - geodetic longitude(deg)
!        kp  - 3-hr geomagnetic activity index
! Output:
!        w(1) = geographic meridional wind (m/sec + northward)
!        w(2) = geograpgic zonal wind (m/sec + eastward)
!        mw(1) = geomagnetic meridional wind (m/sec + northward , QD coordinates)
!        mw(2) = geomagnetic zonal wind (m/sec + eastward , QD coordinates)
!
! If the magnetic latitude corresponding to the input geographic latitude is below 40 (|MLAT|<40), then 
! missing values (= -9999.0) are returned.
!
! To run HL-TWiM and writing output to a text file using gfortran or ifort, use the following from terminal:
!	    gfortran checkhltwim.f90 hltwim.f90 -o hltwim.output.exe
!	    ./hltwim.output.exe > result.txt
!or
!   	ifort checkhltwim.f90 hltwim.f90 -o hltwim.output.exe
!	    ./hltwim.output.exe > result.txt
!******************************************************************************

program checkhltwim

    implicit none

    integer :: day
    real(4) :: glat, glon, kp, ut
    real(4) :: w(2), mw(2)
    integer :: ilat, ilon, ikp, iday, iut
    real(4), external :: pershift

    !*************************
    ! Latitude profile
    !*************************
    day = 10
    ut = 21.0
    ! glat = 65.0
    glon = 215.0
    kp = 3.0
    print '(3a22)', 'Latitude Profile'
    print '(40x,3a22)', 'Geographic Wind', 'Geomagnetic Wind'
    print '(a6,a6,a6,3(a12,a10))', 'day', 'ut', 'kp', 'glat', 'glon', 'merid', 'zonal', 'merid', 'zonal'
    do ilat = -90, 90, 10
        glat = ilat
        call hltwim(day, ut, glat, glon, kp, w, mw)
        print '(I5,f8.1,f6.1,3(f12.3,f10.3))', day, ut, kp, glat, glon, w, mw
    end do
    print *
    print *


    !*************************
    ! Longitude profile
    !*************************
    day = 10
    ut = 21.0
    glat = 65.0
    ! glon = 215.0
    kp = 3.0
    print '(3a22)', 'Longitude Profile'
    print '(40x,3a22)', 'Geographic Wind', 'Geomagnetic Wind'
    print '(a6,a6,a6,3(a12,a10))', 'day', 'ut', 'kp', 'glat', 'glon', 'merid', 'zonal', 'merid', 'zonal'

    do ilon = 0, 360, 30
        glon = ilon
        call hltwim(day, ut, glat, glon, kp, w, mw)
        print '(I5,f8.1,f6.1,3(f12.3,f10.3))', day, ut, kp, glat, glon, w, mw
    end do
    print *
    print *


    !*************************
    ! Geomagnetic Activity Profile
    !*************************
    day = 10
    ut = 21.0
    glat = 65.0
    glon = 215.0
    !kp = 3.0
    print '(3a30)', 'Geomagnetic Activity Profile'
    print '(40x,3a22)', 'Geographic Wind', 'Geomagnetic Wind'
    print '(a6,a6,a6,3(a12,a10))', 'day', 'ut', 'kp', 'glat', 'glon', 'merid', 'zonal', 'merid', 'zonal'

    do ikp = 1, 7, 1
        kp = ikp
        call hltwim(day, ut, glat, glon, kp, w, mw)
        print '(I5,f8.1,f6.1,3(f12.3,f10.3))', day, ut, kp, glat, glon, w, mw
    end do
    print *
    print *

    !*************************
    ! Day of Year Profile
    !*************************
    ! day = 10
    ut = 21.0
    glat = 65.0
    glon = 215.0
    kp = 3.0
    print '(3a22)', 'Day of Year Profile'
    print '(40x,3a22)', 'Geographic Wind', 'Geomagnetic Wind'
    print '(a6,a6,a6,3(a12,a10))', 'day', 'ut', 'kp', 'glat', 'glon', 'merid', 'zonal', 'merid', 'zonal'

    do iday = 1, 365, 20
        day = iday
        call hltwim(day, ut, glat, glon, kp, w, mw)
        print '(I5,f8.1,f6.1,3(f12.3,f10.3))', day, ut, kp, glat, glon, w, mw
    end do
    print *
    print *

    !*************************
    ! Time Profile
    !*************************
    day = 10
    !ut = 21.0
    glat = 65.0
    glon = 215.0
    kp = 3.0
    print '(3a22)', 'UT Profile'
    print '(40x,3a22)', 'Geographic Wind', 'Geomagnetic Wind'
    print '(a6,a6,a6,3(a12,a10))', 'day', 'ut', 'kp', 'glat', 'glon', 'merid', 'zonal', 'merid', 'zonal'

    do iut = 0, 23, 1
        ut = iut
        call hltwim(day, ut, glat, glon, kp, w, mw)
        print '(I5,f8.1,f6.1,3(f12.3,f10.3))', day, ut, kp, glat, glon, w, mw
    end do
    print *
    print *

end program checkhltwim

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

            real(4), parameter :: tol = 1e-4
            real(4) :: x, perint(0:1)
            real(4) :: a, span, offset, offset1, pershift

            pershift = x
            a = perint(0)
            span = perint(1) - perint(0)
            if (span .ne. 0) then
                offset = x - a
                offset1 = mod(offset, span)
                if (abs(offset1) .lt. tol) offset1 = 0
            endif
            pershift = a + offset1
            if ((offset .lt. 0) .and. (offset1 .ne. 0)) pershift = pershift + span

            return

        end function pershift


        !******************************************************************************
        !  TEST OUTPUT: Default
        !******************************************************************************
        !      Latitude Profile
        !                                               Geographic Wind      Geomagnetic Wind
        !   day    ut    kp        glat      glon       merid     zonal       merid     zonal
        !   10    21.0   3.0     -90.000   215.000     180.839   -59.398     -65.540  -180.982
        !   10    21.0   3.0     -80.000   215.000     125.124  -129.992     -49.205  -173.275
        !   10    21.0   3.0     -70.000   215.000      30.165   -92.140     -32.143   -90.376
        !   10    21.0   3.0     -60.000   215.000     -15.729   -28.685     -24.929   -17.367
        !   10    21.0   3.0     -50.000   215.000     -17.181   -28.392     -23.898   -19.081
        !   10    21.0   3.0     -40.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0     -30.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0     -20.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0     -10.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0       0.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0      10.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0      20.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0      30.000   215.000   -9999.000 -9999.000   -9999.000 -9999.000
        !   10    21.0   3.0      40.000   215.000      73.187   -19.431      67.478   -36.659
        !   10    21.0   3.0      50.000   215.000      67.763   -48.344      51.721   -61.487
        !   10    21.0   3.0      60.000   215.000      43.899   -61.858      18.577   -66.610
        !   10    21.0   3.0      70.000   215.000      56.949   -46.422      31.016   -60.011
        !   10    21.0   3.0      80.000   215.000      65.542   -56.344       6.447   -79.246
        !   10    21.0   3.0      90.000   215.000      -1.567   -23.885     -19.550    14.101
        !
        !
        !     Longitude Profile
        !                                               Geographic Wind      Geomagnetic Wind
        !   day    ut    kp        glat      glon       merid     zonal       merid     zonal
        !   10    21.0   3.0      65.000     0.000     -34.654   -30.684     -30.548   -37.978
        !   10    21.0   3.0      65.000    30.000     -40.240   -15.356     -40.057   -25.327
        !   10    21.0   3.0      65.000    60.000     -45.699    21.443     -43.615    14.853
        !   10    21.0   3.0      65.000    90.000     -52.254    43.607     -48.098    43.271
        !   10    21.0   3.0      65.000   120.000     -44.859     2.404     -41.602     6.164
        !   10    21.0   3.0      65.000   150.000     -18.287   -47.500     -20.087   -41.338
        !   10    21.0   3.0      65.000   180.000      12.575   -56.846       0.073   -53.371
        !   10    21.0   3.0      65.000   210.000      39.349   -53.049      16.478   -57.603
        !   10    21.0   3.0      65.000   240.000      83.064   -58.433      48.921   -76.699
        !   10    21.0   3.0      65.000   270.000     107.171   -62.908      89.184   -75.382
        !   10    21.0   3.0      65.000   300.000      27.501   -46.094      43.948   -40.042
        !   10    21.0   3.0      65.000   330.000     -24.989   -46.337      -6.589   -46.716
        !   10    21.0   3.0      65.000   360.000     -34.654   -30.684     -30.548   -37.978
        !
        !
        !  Geomagnetic Activity Profile
        !                                               Geographic Wind      Geomagnetic Wind
        !   day    ut    kp        glat      glon       merid     zonal       merid     zonal
        !   10    21.0   1.0      65.000   215.000      67.733   -40.348      49.055   -56.186
        !   10    21.0   2.0      65.000   215.000      55.968   -45.070      35.040   -56.372
        !   10    21.0   3.0      65.000   215.000      45.492   -53.303      20.675   -59.953
        !   10    21.0   4.0      65.000   215.000      45.440   -65.411      14.951   -70.202
        !   10    21.0   5.0      65.000   215.000      55.810   -81.394      17.869   -87.118
        !   10    21.0   6.0      65.000   215.000      66.198   -93.341      22.695  -100.618
        !   10    21.0   7.0      65.000   215.000      67.497   -94.834      23.299  -102.305
        !
        !
        !   Day of Year Profile
        !                                               Geographic Wind      Geomagnetic Wind
        !   day    ut    kp        glat      glon       merid     zonal       merid     zonal
        !    1    21.0   3.0      65.000   215.000      41.760   -53.502      16.837   -58.911
        !   21    21.0   3.0      65.000   215.000      53.126   -57.016      26.594   -65.578
        !   41    21.0   3.0      65.000   215.000      72.336   -73.096      38.335   -85.444
        !   61    21.0   3.0      65.000   215.000      91.009   -95.840      46.415  -110.787
        !   81    21.0   3.0      65.000   215.000     100.700  -117.648      45.924  -132.421
        !  101    21.0   3.0      65.000   215.000      96.622  -132.267      34.986  -143.493
        !  121    21.0   3.0      65.000   215.000      79.449  -137.002      15.540  -141.935
        !  141    21.0   3.0      65.000   215.000      55.209  -133.512      -7.144  -131.111
        !  161    21.0   3.0      65.000   215.000      33.478  -127.155     -25.968  -118.671
        !  181    21.0   3.0      65.000   215.000      23.313  -124.385     -34.868  -113.024
        !  201    21.0   3.0      65.000   215.000      28.515  -128.845     -31.739  -118.494
        !  221    21.0   3.0      65.000   215.000      45.547  -138.640     -19.239  -132.325
        !  241    21.0   3.0      65.000   215.000      65.771  -147.471      -3.085  -146.373
        !  261    21.0   3.0      65.000   215.000      80.180  -148.826      10.735  -152.198
        !  281    21.0   3.0      65.000   215.000      83.192  -139.238      18.247  -145.046
        !  301    21.0   3.0      65.000   215.000      74.486  -119.338      18.834  -125.349
        !  321    21.0   3.0      65.000   215.000      59.096   -93.766      15.371   -98.674
        !  341    21.0   3.0      65.000   215.000      45.287   -69.913      12.689   -73.969
        !  361    21.0   3.0      65.000   215.000      40.688   -55.279      14.929   -60.070
        !
        !
        !            UT Profile
        !                                               Geographic Wind      Geomagnetic Wind
        !   day    ut    kp        glat      glon       merid     zonal       merid     zonal
        !   10     0.0   3.0      65.000   215.000      81.972   -63.134      52.667   -80.125
        !   10     1.0   3.0      65.000   215.000      82.288   -46.548      60.753   -66.165
        !   10     2.0   3.0      65.000   215.000      75.627   -24.950      64.185   -45.692
        !   10     3.0   3.0      65.000   215.000      63.290    -4.518      61.378   -24.365
        !   10     4.0   3.0      65.000   215.000      48.697     4.045      50.748   -12.371
        !   10     5.0   3.0      65.000   215.000      38.418   -11.433      33.187   -22.158
        !   10     6.0   3.0      65.000   215.000      34.836   -46.330      13.250   -50.584
        !   10     7.0   3.0      65.000   215.000      30.148   -77.139      -5.883   -75.184
        !   10     8.0   3.0      65.000   215.000      16.697   -88.447     -24.674   -80.407
        !   10     9.0   3.0      65.000   215.000      -6.666   -80.158     -44.230   -65.800
        !   10    10.0   3.0      65.000   215.000     -36.808   -58.376     -64.267   -37.552
        !   10    11.0   3.0      65.000   215.000     -67.058   -31.895     -82.213    -5.284
        !   10    12.0   3.0      65.000   215.000     -88.838   -10.155     -93.881    20.214
        !   10    13.0   3.0      65.000   215.000     -96.460     1.411     -96.110    32.494
        !   10    14.0   3.0      65.000   215.000     -89.529     3.101     -88.367    31.677
        !   10    15.0   3.0      65.000   215.000     -71.103    -1.635     -72.099    21.683
        !   10    16.0   3.0      65.000   215.000     -46.108    -8.682     -50.324     7.599
        !   10    17.0   3.0      65.000   215.000     -20.667   -15.024     -27.770    -6.033
        !   10    18.0   3.0      65.000   215.000       0.492   -20.606      -9.158   -17.630
        !   10    19.0   3.0      65.000   215.000      16.800   -27.952       3.763   -29.150
        !   10    20.0   3.0      65.000   215.000      30.897   -39.142      12.665   -43.211
        !   10    21.0   3.0      65.000   215.000      45.492   -53.303      20.675   -59.953
        !   10    22.0   3.0      65.000   215.000      60.895   -65.799      30.275   -75.546
        !   10    23.0   3.0      65.000   215.000      74.294   -70.187      41.661   -83.613

