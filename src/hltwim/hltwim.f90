!!!
!!!  High Latitude Thermospheric Wind Model - HL-TWiM
!!!
!!!  AUTHORS
!!!    Manbharat Singh Dhadly
!!!    John Emmert
!!!    Douglas Drob
!!!    ----------------------
!!!    Geospace Science and Technology Branch
!!!    Space Science Division
!!!    Naval Research Laboratory
!!!    4555 Overlook Ave.
!!!    Washington, DC 20375
!!!
!!!  Point of Contact
!!!   manbharat.dhadly.ctr@nrl.navy.mil
!!!
!!!   DATE
!!!    05-June-2019
!!!-------------------------------------
!!!  Model Call
!!!     hltwim(day, ut, glat, glon, kp, w, mw)
!!!
!!!  Input arguments:
!!!        day - day of year
!!!        ut  - universal time
!!!        kp  - 3-hr geomagnetic activity index
!!!        glat - geodetic latitude(deg)
!!!        glon - geodetic longitude(deg)
!!!
!!!  Output:
!!!        w(1) = geographic meridional wind (m/sec + northward)
!!!        w(2) = geograpgic zonal wind (m/sec + eastward)
!!!        mw(1) = geomagnetic meridional wind (m/sec + northward , QD coordinates)
!!!        mw(2) = geomagnetic zonal wind (m/sec + eastward , QD coordinates)
!!!-------------------------------------
!!! If the magnetic latitude corresponding to the input geographic latitude is below 40 (|MLAT|<40), then 
!!! missing values (= -9999.0) are returned.
!!!
!!! To run HL-TWiM and writing output to a text file using gfortran or ifort, use the following from terminal:
!!!	        gfortran checkhltwim.f90 hltwim.f90 -o hltwim.output.exe
!!!	        ./hltwim.output.exe > result.txt
!!! or
!!!	        ifort checkhltwim.f90 hltwim.f90 -o hltwim.output.exe
!!!	        ./hltwim.output.exe > result.txt

!####################################################################################
! Model Modules
!####################################################################################

module Base

    integer(4)           :: nmaxhltwim = 0        ! maximum degree of the model
    integer(4)           :: mmaxhltwim = 0        ! maximum order of the model     
    integer(4)           :: nmaxgeo = 0           ! maximum degree of coordinate coversion
    integer(4)           :: mmaxgeo = 0           ! maximum order of coordinate coversion
    real(8),allocatable  :: gpbar(:,:),gvbar(:,:),gwbar(:,:) ! alfs for geo coordinates
    real(8),allocatable  :: spbar(:,:),svbar(:,:),swbar(:,:) ! alfs MLT calculation
    real(8)              :: glatalf = -1.d32
    real(4)              :: xlat=0.0
    logical              :: baseinit = .true.

end module base

module base1

    implicit none

    integer(4)                 :: nterm             ! Number of terms in the model
    integer(4)                 :: nmax,mmax         ! Max latitudinal and local time garmonic (degree and order)
    integer(4)                 :: nvshterm          ! # of VSH basis functions
    integer(4)                 :: termarr(0:3,0:2175)    ! 3 x nterm index of coupled terms
    integer(4)                 :: termarr1(0:3,0:2175)    ! 3 x nterm index of coupled terms
    integer(4)                 :: termarr2(0:3,0:2175)    ! 3 x nterm index of coupled terms
    real(4)                    :: coeffx(0:2175)          ! Model coefficients
    real(4)                    :: coeffx1(0:2175)          ! Model coefficients
    real(4)                    :: coeffx2(0:2175)          ! Model coefficients
    real(4),allocatable        :: vshterms(:,:)     ! VSH basis values
    real(4),allocatable        :: termval(:,:)      ! Term values to which coefficients are applied
    real(8),allocatable        :: dpbar(:,:)        ! Associated lengendre fns
    real(8),allocatable        :: dvbar(:,:)        ! Associated lengendre fns
    real(8),allocatable        :: dwbar(:,:)        ! Associated lengendre fns
    real(8),allocatable        :: mltterms(:,:)     ! MLT Fourier terms
    real(8),allocatable        :: doyterms(:,:)     ! MLT Fourier terms
    real(8),allocatable        :: mlonterms(:,:)    ! MLT Fourier terms
    real(8), parameter         :: pi=3.1415926535897932
    real(8), parameter         :: dtor=pi/180.d0
    logical                    :: xinit = .true.
    character(128), parameter  :: xdefault1 = 'hltwimnorth.dat' !Fit coefficient file - North
    character(128), parameter  :: xdefault2 = 'hltwimsouth.dat' !Fit coefficient file - South
    
end module base1

! ################################################################################
! Portable utility to compute vector spherical harmonical harmonic basis functions
! ################################################################################

module alf

    implicit none

    integer(4)              :: nmax0,mmax0
    ! static normalizational coeffiecents
    real(8), allocatable    :: anm(:,:),bnm(:,:),dnm(:,:)
    real(8), allocatable    :: cm(:),en(:)
    real(8), allocatable    :: marr(:),narr(:)

contains

    ! -------------------------------------------------------------
    ! routine to compute vector spherical harmonic basis functions
    ! -------------------------------------------------------------

    subroutine alfbasis(nmax,mmax,theta,P,V,W)

        implicit none

        integer(4), intent(in)  :: nmax, mmax
        real(8), intent(in)     :: theta
        real(8), intent(out)    :: P(0:nmax,0:mmax)
        real(8), intent(out)    :: V(0:nmax,0:mmax)
        real(8), intent(out)    :: W(0:nmax,0:mmax)
        integer(8)              :: n, m
        real(8)                 :: x, y
        real(8), parameter      :: p00 = 0.70710678118654746d0

        ! P/V/W are intent(out) -> undefined on entry; the recursions below do not
        ! touch every element (e.g. V(0,0)), and gd2qd() reads V(0,0)*normadj(0),
        ! so an undefined NaN there poisons the result. Zero them first. L. Navarro
        P = 0.d0
        V = 0.d0
        W = 0.d0

        P(0,0) = p00
        x = dcos(theta)
        y = dsin(theta)
        do m = 1, mmax
            W(m,m) = cm(m) * P(m-1,m-1)
            P(m,m) = y * en(m) * W(m,m)
            do n = m+1, nmax
                W(n,m) = anm(n,m) * x * W(n-1,m) - bnm(n,m) * W(n-2,m)
                P(n,m) = y * en(n) * W(n,m)
                V(n,m) = narr(n) * x * W(n,m) - dnm(n,m) * W(n-1,m)
                W(n-2,m) = marr(m) * W(n-2,m)
            enddo
            W(nmax-1,m) = marr(m) * W(nmax-1,m)
            W(nmax,m) = marr(m) * W(nmax,m)
            V(m,m) = x * W(m,m)
        enddo
        P(1,0) = anm(1,0) * x * P(0,0)
        V(1,0) = -P(1,1)
        do n = 2, nmax
            P(n,0) = anm(n,0) * x * P(n-1,0) - bnm(n,0) * P(n-2,0)
            V(n,0) = -P(n,1)
        enddo        
        return

    end subroutine alfbasis

    ! -----------------------------------------------------
    ! routine to compute static normalization coeffiecents
    ! -----------------------------------------------------

    subroutine initalf(nmaxin,mmaxin)

        implicit none

        integer(4), intent(in) :: nmaxin, mmaxin
        integer(8)             :: n, m   ! 64 bits to avoid overflow for (m,n) > 60

        nmax0 = nmaxin
        mmax0 = mmaxin
        if (allocated(anm)) deallocate(anm, bnm, cm, dnm, en, marr, narr)
        allocate( anm(0:nmax0, 0:mmax0) )
        allocate( bnm(0:nmax0, 0:mmax0) )
        allocate( cm(0:mmax0) )
        allocate( dnm(0:nmax0, 0:mmax0) )
        allocate( en(0:nmax0) )
        allocate( marr(0:mmax0) )
        allocate( narr(0:nmax0) )
        do n = 1, nmax0
            narr(n) = dble(n)
            en(n)    = dsqrt(dble(n*(n+1)))
            anm(n,0) = dsqrt( dble((2*n-1)*(2*n+1)) ) / narr(n)
            bnm(n,0) = dsqrt( dble((2*n+1)*(n-1)*(n-1)) / dble(2*n-3) ) / narr(n)
        end do
        do m = 1, mmax0
            marr(m) = dble(m)
            cm(m)    = dsqrt(dble(2*m+1)/dble(2*m*m*(m+1)))
            do n = m+1, nmax0
                anm(n,m) = dsqrt( dble((2*n-1)*(2*n+1)*(n-1)) / dble((n-m)*(n+m)*(n+1)) )
                bnm(n,m) = dsqrt( dble((2*n+1)*(n+m-1)*(n-m-1)*(n-2)*(n-1)) &
                    / dble((n-m)*(n+m)*(2*n-3)*n*(n+1)) )
                dnm(n,m) = dsqrt( dble((n-m)*(n+m)*(2*n+1)*(n-1)) / dble((2*n-1)*(n+1)) )
            end do
        enddo
		
        return

    end subroutine initalf

end module alf

! ###################################################################################################################################################
!                         HL-TWiM Model Functions
! ###################################################################################################################################################

subroutine inithltwim1(nmaxout,mmaxout)

    use base
    use base1
    implicit none

    integer(4),intent(out)     :: nmaxout, mmaxout

    call findandopen(xdefault1,23)
    read(23) nterm, mmax, nmax, termarr1, coeffx1
    close(23)
    call findandopen(xdefault2,23)
    read(23) nterm, mmax, nmax, termarr2, coeffx2
    close(23)

    if (allocated(termval)) deallocate(termval,dpbar,dvbar,dwbar,mltterms,vshterms, doyterms,mlonterms)
    nvshterm = ( ((nmax+1)*(nmax+2) - (nmax-mmax)*(nmax-mmax+1))/2 - 1 ) * 4 - 2*nmax
    !nvshterm = 127.0 for nxm = 10x3
    allocate(termval(0:1, 0:nterm-1))
    allocate(dpbar(0:nmax,0:mmax),dvbar(0:nmax,0:mmax),dwbar(0:nmax,0:mmax))
    allocate(mltterms(0:mmax,0:1))
    allocate(doyterms(0:1,0:4))
    allocate(mlonterms(0:1,0:2))
    allocate(vshterms(0:1, 0:nvshterm-1))
    ! HLTWIMx only fills the entries it uses; zero the rest so an unfilled slot
    ! reached via termarr(...) can't poison the dot_product with NaN. L. Navarro
    dpbar = 0
    dvbar = 0
    dwbar = 0
    termval = 0
    mltterms = 0
    doyterms = 0
    mlonterms = 0
    vshterms = 0
    nmaxout = nmax
    mmaxout = mmax
    xinit = .false.

    return

end subroutine inithltwim1

subroutine inithltwim2()

    use base
    use base1
    use alf,only:initalf
    implicit none

    integer(4)           :: nmax0, mmax0

    if (xinit) call inithltwim1(nmaxhltwim, mmaxhltwim)
    nmaxgeo = 8 !  ! maximum degree of coordinate coversion in gd2qd.dat
    mmaxgeo = 4 !  ! maximum order of coordinate coversion in gd2qd.dat
    ! initalf must cover BOTH the wind model (nmax x mmax, from hltwim*.dat) and the
    ! gd2qd coordinate transform (nmaxgeo x mmaxgeo) -- gd2qd calls alfbasis with
    ! mmaxgeo=4 while the wind model has mmax=3, so sizing for `mmax` alone left
    ! alfbasis reading anm/bnm/dnm/cm one past the bound (nondeterministic NaN).
    ! L. Navarro
    nmax0 = max(nmax, nmaxgeo)
    mmax0 = max(mmax, mmaxgeo)
    call initalf(nmax0,mmax0)
    if (allocated(gpbar)) deallocate(gpbar,gvbar,gwbar)
    allocate(gpbar(0:nmaxgeo,0:mmaxgeo))
    allocate(gvbar(0:nmaxgeo,0:mmaxgeo))
    allocate(gwbar(0:nmaxgeo,0:mmaxgeo))
    gpbar = 0
    gvbar = 0
    gwbar = 0
    if (allocated(spbar)) deallocate(spbar,svbar,swbar)
    allocate(spbar(0:nmaxgeo,0:mmaxgeo))
    allocate(svbar(0:nmaxgeo,0:mmaxgeo))
    allocate(swbar(0:nmaxgeo,0:mmaxgeo))
    spbar = 0
    svbar = 0
    swbar = 0
    baseinit = .false.

    return

end subroutine inithltwim2


subroutine HLTWIM(DAY,UT,GLAT,GLON,KP,W,MW)

    use base
    use base1
    implicit none

    INTEGER,intent(in)      :: DAY
    REAL(4),intent(in)      :: UT,GLAT,GLON
    REAL(4),intent(in)      :: KP
    REAL(4),intent(out)     :: W(2), MW(2)
    real(4), save           :: mlat, mlon, mlt
    real(4)                 :: mmwind, mzwind
    real(4), save           :: f1e, f1n, f2e, f2n
    real(4), save           :: glatlast=1.0e16, glonlast=1.0e16
    real(4), save           :: daylast=1.0e16, utlast=1.0e16, kplast=1.0e16
    real(4), external       ::  mltcalc
	
    xlat=glat
    IF (DAY .le. 0) ERROR STOP "Hello!! day<=0 not allowed"
    IF (DAY .gt. 366) ERROR STOP "Hello!! day>366 not allowed"
    IF (UT .lt. 0) ERROR STOP "Hello!! UT<0 not allowed"    
    IF (UT .ge. 24) ERROR STOP "Hello!! UT>24 not allowed"
    IF (KP .lt. 0) ERROR STOP "Hello!! Kp< 0 not allowed"
    IF (KP .gt. 10) ERROR STOP "Hello!! Kp> 10 not allowed"
    IF (GLAT.GT.90) ERROR STOP "Hello!! GLAT>90 not allowed"
    IF (GLAT.lt.-90) ERROR STOP "Hello!! GLAT<-90 not allowed"

     if (baseinit) call inithltwim2()
    !CONVERT GEO LAT/LON TO QD LAT/LON
    if ((glat .ne. glatlast) .or. (glon .ne. glonlast)) then
        call gd2qd(glat,glon,mlat,mlon,f1e,f1n,f2e,f2n)
    endif
    !COMPUTE QD MAGNETIC LOCAL TIME (LOW-PRECISION)
    if ((day .ne. daylast) .or. (ut .ne. utlast) .or. &
        (glat .ne. glatlast) .or. (glon .ne. glonlast)) then
        mlt = mltcalc(mlat,mlon,day,ut)
    endif

    !Replace calcuated winds with missing value when |mlat|<40
    IF ((mlat .gt. - 40).and.(mlat .lt. 40)) then
        MW(1) = -9999.0
        MW(2) = -9999.0
        W(1) = -9999.0
        W(2) = -9999.0
    else
        !RETRIEVE HLTWIM geomagnetic winds
        call HLTWIMx(day, mlt, mlat, kp, mlon, mmwind, mzwind)
        MW(1) = mmwind
        MW(2) = mzwind
        !CONVERT GEOMAGNETIC WINDS TO GEOGRAPHIC COORDINATES
        W(1) = f2n * mmwind + f1n * mzwind
        W(2) = f2e * mmwind + f1e * mzwind
    endif

    glatlast = glat
    glonlast = glon
    daylast = day
    utlast = ut
    kplast = kp

    return

end subroutine HLTWIM

subroutine HLTWIMx(day, mlt, mlat, kp, mlon, mmwind, mzwind)

    use base
    use base1
    use alf,only:alfbasis
    implicit none
    
    INTEGER,intent(in)        :: day
    real(4),intent(in)        :: mlt       !Magnetic local time (hours)
    real(4),intent(in)        :: mlat      !Magnetic latitude (degrees)
    real(4),intent(in)        :: kp        !3-hour Kp
    real(4),intent(in)        :: mlon      !Magnetic Longitude (degree)
    real(4),intent(out)       :: mmwind   !Meridional wind (+north, QD coordinates)
    real(4),intent(out)       :: mzwind   !Zonal Wind (+east, QD coordinates)
    ! Local variables
    integer(4)                :: iterm, ivshterm, n, m, idoy, imlon
    real(4)                   :: termvaltemp(0:1)
    real(4),save              :: kpterms(0:2)
    real(4)                   :: kptermsxx(0:1,0:2)
    real(4),save              :: mltlast=1.e16, mlatlast=1.e16, kplast=1.e16, doylast=1.e16, mlonlast=1.e16
    real(8)                   :: theta, phi, mphi, xx, dxx, mlonx, mlonxx

    !LOAD MODEL PARAMETERS IF NECESSARY
    if (xinit) call inithltwim1(nmaxhltwim, mmaxhltwim)
	
    !COMPUTE LATITUDE PART OF VSH TERMS
    if (mlat .ne. mlatlast) then
        theta = (90.d0 - dble(mlat))*dtor
        call alfbasis(nmax,mmax,theta,dpbar,dvbar,dwbar)
    endif
    !COMPUTE MLT PART OF VSH TERMS
    if (mlt .ne. mltlast) then
        phi = dble(mlt)*dtor*15.d0
        do m = 0, mmax
            mphi = dble(m)*phi
            mltterms(m,0) = dcos(mphi)
            mltterms(m,1) = dsin(mphi)
        enddo
    endif
    !COMPUTE VSH TERMS
    if ((mlat .ne. mlatlast) .or. (mlt .ne. mltlast)) then
        ivshterm = 0
        do n = 1, nmax
!            print *,'n=',n
!            print *,'m=',m
            vshterms(0,ivshterm)   = -sngl(dvbar(n,0)*mltterms(0,0))
            vshterms(0,ivshterm+1) =  sngl(dwbar(n,0)*mltterms(0,0))
            vshterms(1,ivshterm)   = -vshterms(0,ivshterm+1)
            vshterms(1,ivshterm+1) =  vshterms(0,ivshterm)
            ivshterm = ivshterm + 2
            do m = 1, mmax
                if (m .gt. n) cycle
                vshterms(0,ivshterm)   = -sngl(dvbar(n,m)*mltterms(m,0))
                vshterms(0,ivshterm+1) =  sngl(dvbar(n,m)*mltterms(m,1))
                vshterms(0,ivshterm+2) =  sngl(dwbar(n,m)*mltterms(m,1))
                vshterms(0,ivshterm+3) =  sngl(dwbar(n,m)*mltterms(m,0))
                vshterms(1,ivshterm)   = -vshterms(0,ivshterm+2)
                vshterms(1,ivshterm+1) = -vshterms(0,ivshterm+3)
                vshterms(1,ivshterm+2) =  vshterms(0,ivshterm)
                vshterms(1,ivshterm+3) =  vshterms(0,ivshterm+1)
                ivshterm = ivshterm + 4
            enddo
        enddo
    endif
    !COMPUTE DOY TERMS
    if (day .ne. doylast) then
        xx = dble(day)*dtor*360.d0/366.d0
        idoy=0
        dxx = dble(0)*xx
        doyterms(0,0) = 1.0
        idoy=1
        dxx = dble(idoy)*xx
        doyterms(0,1) = dcos(dxx)
        doyterms(0,2)=dsin(dxx)  
        idoy=2
        dxx = dble(idoy)*xx
        doyterms(0,3) = dcos(dxx)
        doyterms(0,4)=dsin(dxx)
    endif
    doyterms(1,:) = doyterms(0,:)    
    !COMPUTE MLON TERMS
    if (mlon .ne. mlonlast) then
        mlonx = dble(mlon)*dtor
        imlon=0
        mlonxx = dble(imlon)*mlonx
        mlonterms(0,0) = 1.0
        imlon=1
        mlonxx = dble(imlon)*mlonx
        mlonterms(0,1) = dcos(mlonxx)
        mlonterms(0,2) = dsin(mlonxx)
    endif
    mlonterms(1,:) = mlonterms(0,:)
        !COMPUTE KP TERMS
    if (kp .ne. kplast) then
        call kpspl3(kp, kpterms)
    endif 
    kptermsxx(0,0) = kpterms(0)
    kptermsxx(0,1) = kpterms(1)
    kptermsxx(0,2) = kpterms(2)
    kptermsxx(1,:)=kptermsxx(0,:)
        !GENERATE COUPLED TERMS

    if (xlat.ge.0) termarr=termarr1
    if (xlat.le.0) termarr=termarr2
    if (xlat.ge.0) coeffx=coeffx1
    if (xlat.le.0) coeffx=coeffx2

    do iterm = 0, nterm-1
        termvaltemp = (/1.0, 1.0/)
        if (termarr(0,iterm) .ne. 999) termvaltemp = termvaltemp * vshterms(0:1,termarr(0,iterm))
        if (termarr(1,iterm) .ne. 999) termvaltemp = termvaltemp * doyterms(0:1,termarr(1,iterm))  
        if (termarr(2,iterm) .ne. 999) termvaltemp = termvaltemp * kptermsxx(0:1, termarr(2,iterm))
        if (termarr(3,iterm) .ne. 999) termvaltemp = termvaltemp * mlonterms(0:1,termarr(3,iterm))
        termval(0:1,iterm) = termvaltemp(0:1)
    enddo

    !APPLY COEFFICIENTS AND CALCAULTE QD WINDS
    mmwind = dot_product(coeffx, termval(0,0:nterm-1))
    mzwind = dot_product(coeffx, termval(1,0:nterm-1))
	
    mlatlast = mlat
    mltlast = mlt
    kplast = kp
    doylast = day
    mlonlast = mlon

    return

end subroutine HLTWIMx

!!=================================================================================
!!                           Convert Ap to Kp
!!=================================================================================
!
!function ap2kp(ap0)
!
!    real(4), parameter :: apgrid(0:27) = (/0.,2.,3.,4.,5.,6.,7.,9.,12.,15.,18., &
!        22.,27.,32.,39.,48.,56.,67.,80.,94., &
!        111.,132.,154.,179.,207.,236.,300.,400./)
!    real(4), parameter :: kpgrid(0:27) = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11., &
!        12.,13.,14.,15.,16.,17.,18.,19.,20.,21., &
!        22.,23.,24.,25.,26.,27./) / 3.0
!    real(4)            :: ap0, ap, ap2kp
!    integer(4)         :: i
!
!
!    ap = ap0
!    if (ap .lt. 0) ap = 0
!    if (ap .gt. 400) ap = 400
!
!    i = 1
!    do while (ap .gt. apgrid(i))
!        i = i + 1
!    end do
!    if (ap .eq. apgrid(i)) then
!        ap2kp = kpgrid(i)
!    else
!        ap2kp = kpgrid(i-1) + (ap - apgrid(i-1)) / (3.0 * (apgrid(i) - apgrid(i-1)))
!    end if
!
!    return
!
!end function ap2kp

! ########################################################################
!     Geographic <=> Geomagnetic Coordinate Transformations
!
!  Converts geodetic coordinates to Quasi-Dipole coordinates (Richmond, J. Geomag.
!  Geoelec., 1995, p. 191), using a spherical harmonic representation.
!
! ########################################################################

module gd2qdc

    implicit none

    integer(4)               :: nterm, nmax, mmax  !Spherical harmonic expansion parameters

    real(8), allocatable     :: coeff(:,:)         !Coefficients for spherical harmonic expansion
    real(8), allocatable     :: xcoeff(:)          !Coefficients for x coordinate
    real(8), allocatable     :: ycoeff(:)          !Coefficients for y coordinate
    real(8), allocatable     :: zcoeff(:)          !Coefficients for z coordinate
    real(8), allocatable     :: sh(:)              !Array to hold spherical harmonic fuctions
    real(8), allocatable     :: shgradtheta(:)     !Array to hold spherical harmonic gradients
    real(8), allocatable     :: shgradphi(:)       !Array to hold spherical harmonic gradients
    real(8), allocatable     :: normadj(:)         !Adjustment to VSH normalization factor
    real(4)                  :: epoch, alt
    real(8), parameter       :: pi = 3.1415926535897932d0
    real(8), parameter       :: dtor = pi/180.0d0
    real(8), parameter       :: sineps = 0.39781868d0
    logical                  :: gd2qdinit = .true.

contains

    subroutine initgd2qd()

        use base
        implicit none

        character(128), parameter   :: datafile='gd2qd.dat'
        integer(4)                  :: iterm, n
        integer(4)                  :: j

        call findandopen(datafile,23)
        read(23) nmax, mmax, nterm, epoch, alt
        if (allocated(coeff)) then
            deallocate(coeff,xcoeff,ycoeff,zcoeff,sh,shgradtheta,shgradphi,normadj)
        endif
        allocate( coeff(0:nterm-1, 0:2) )
        read(23) coeff
        close(23)

        allocate( xcoeff(0:nterm-1) )
        allocate( ycoeff(0:nterm-1) )
        allocate( zcoeff(0:nterm-1) )
        allocate( sh(0:nterm-1) )
        allocate( shgradtheta(0:nterm-1) )
        allocate( shgradphi(0:nterm-1) )
        allocate( normadj(0:nmax) )

        ! gd2qd() only fills the entries it uses; zero the rest so the trailing
        ! (uninitialised) elements don't poison dot_product(sh, *coeff) with NaN.
        ! L. Navarro
        sh = 0.d0
        shgradtheta = 0.d0
        shgradphi = 0.d0

        do iterm = 0, nterm-1
            xcoeff(iterm) = coeff(iterm,0)
            ycoeff(iterm) = coeff(iterm,1)
            zcoeff(iterm) = coeff(iterm,2)
        enddo

        do n = 0, nmax
            normadj(n) = dsqrt(dble(n*(n+1)))
        end do

        gd2qdinit = .false.

        return

    end subroutine initgd2qd

end module gd2qdc

subroutine gd2qd(glatin,glon,qlat,qlon,f1e,f1n,f2e,f2n)

    use base
    use gd2qdc
    use alf

    implicit none

    real(4), intent(in)         :: glatin, glon
    real(4), intent(out)        :: qlat, qlon
    real(4), intent(out)        :: f1e, f1n, f2e, f2n
    integer(4)               :: n, m, i
    real(8)                  :: glat, theta, phi
    real(8)                  :: mphi, cosmphi, sinmphi
    real(8)                  :: x, y, z
    real(8)                  :: cosqlat, cosqlon, sinqlon
    real(8)                  :: xgradtheta, ygradtheta, zgradtheta
    real(8)                  :: xgradphi, ygradphi, zgradphi
    real(8)                  :: qlonrad

    if (gd2qdinit) call initgd2qd()

    glat = dble(glatin)
    if (glat .ne. glatalf) then
        theta = (90.d0 - glat) * dtor
        call alfbasis(nmax,mmax,theta,gpbar,gvbar,gwbar)
        glatalf = glat
    endif
    phi = dble(glon) * dtor
	
    i = 0
    do n = 0, nmax
        sh(i) = gpbar(n,0)
        shgradtheta(i) =  gvbar(n,0) * normadj(n)
        shgradphi(i) = 0
        i = i + 1
    enddo
    do m = 1, mmax
        mphi = dble(m) * phi
        cosmphi = dcos(mphi)
        sinmphi = dsin(mphi)
        do n = m, nmax
           ! print *, 'gpbar(n,m)=',gpbar(n,m)
            sh(i)   = gpbar(n,m) * cosmphi
            sh(i+1) = gpbar(n,m) * sinmphi
            shgradtheta(i)   =  gvbar(n,m) * normadj(n) * cosmphi
            shgradtheta(i+1) =  gvbar(n,m) * normadj(n) * sinmphi
            shgradphi(i)     = -gwbar(n,m) * normadj(n) * sinmphi
            shgradphi(i+1)   =  gwbar(n,m) * normadj(n) * cosmphi
            i = i + 2
        enddo
    enddo

    x = dot_product(sh, xcoeff)
    y = dot_product(sh, ycoeff)
    z = dot_product(sh, zcoeff)

    qlonrad = datan2(y,x)
    cosqlon = dcos(qlonrad)
    sinqlon = dsin(qlonrad)
    cosqlat = x*cosqlon + y*sinqlon

    qlat = sngl(datan2(z,cosqlat) / dtor)
    qlon = sngl(qlonrad / dtor)

    xgradtheta = dot_product(shgradtheta, xcoeff)
    ygradtheta = dot_product(shgradtheta, ycoeff)
    zgradtheta = dot_product(shgradtheta, zcoeff)

    xgradphi = dot_product(shgradphi, xcoeff)
    ygradphi = dot_product(shgradphi, ycoeff)
    zgradphi = dot_product(shgradphi, zcoeff)

    f1e = sngl(-zgradtheta*cosqlat + (xgradtheta*cosqlon + ygradtheta*sinqlon)*z )
    f1n = sngl(-zgradphi*cosqlat   + (xgradphi*cosqlon   + ygradphi*sinqlon)*z )
    f2e = sngl( ygradtheta*cosqlon - xgradtheta*sinqlon )
    f2n = sngl( ygradphi*cosqlon   - xgradphi*sinqlon )

    return

end subroutine gd2qd

!==================================================================================
!                  (Function) Calculate Magnetic Local Time
!==================================================================================

function mltcalc(qlat,qlon,day,ut)

    use base
    use gd2qdc
    use alf

    implicit none

    real(4), intent(in)      :: qlat, qlon, ut
    real(4)                  :: mltcalc
    Integer(4), intent(in)   :: day
    integer(4)               :: n, m, i
    real(8)                  :: asunglat, asunglon, asunqlon
    real(8)                  :: glat, theta, phi
    real(8)                  :: mphi, cosmphi, sinmphi
    real(8)                  :: x, y
    real(8)                  :: cosqlat, cosqlon, sinqlon
    real(8)                  :: qlonrad

    if (gd2qdinit) call initgd2qd()

    !COMPUTE GEOGRAPHIC COORDINATES OF ANTI-SUNWARD DIRECTION (LOW PRECISION)
    asunglat = -asin(sin((dble(day)+dble(ut)/24.0d0-80.0d0)*dtor) * sineps) / dtor
    asunglon = -ut * 15.d0

    !COMPUTE MAGNETIC COORDINATES OF ANTI-SUNWARD DIRECTION
    theta = (90.d0 - asunglat) * dtor
    call alfbasis(nmax,mmax,theta,spbar,svbar,swbar)
    phi = asunglon * dtor
    i = 0
    do n = 0, nmax
        sh(i) = spbar(n,0)
        i = i + 1
    enddo
    do m = 1, mmax
        mphi = dble(m) * phi
        cosmphi = dcos(mphi)
        sinmphi = dsin(mphi)
        do n = m, nmax
            sh(i)   = spbar(n,m) * cosmphi
            sh(i+1) = spbar(n,m) * sinmphi
            i = i + 2
        enddo
    enddo
    x = dot_product(sh, xcoeff)
    y = dot_product(sh, ycoeff)
    asunqlon = sngl(datan2(y,x) / dtor)

    !COMPUTE MLT
    mltcalc = (qlon - asunqlon) / 15.0

    return

end function mltcalc

!================================================================================
!                           Cubic Spline interpolation of Kp
!================================================================================

subroutine kpspl3(kp, kpterms)

    implicit none

    real(4), intent(in)       :: kp
    real(4), intent(out)      :: kpterms(0:2)
    integer(4)                :: i, j
    real(4)                   :: x, kpspl(0:6)
    real(4), parameter        :: node(0:7)=(/-10., -8., 0., 2., 5., 6.5, 18., 20./)

    x = max(kp, 0.0)
    x = min(x,  6.5)
    kpterms(0:2) = 0.0
    do i = 0, 6
        kpspl(i) = 0.0
        if ((x .ge. node(i)) .and. (x .lt. node(i+1))) kpspl(i) = 1.0
    enddo
    do j = 2,3
        do i = 0, 8-j-1
            kpspl(i) = kpspl(i)   * (x - node(i))   / (node(i+j-1) - node(i)) &
                + kpspl(i+1) * (node(i+j) - x) / (node(i+j)   - node(i+1))
        enddo
    enddo
    kpterms(0) = kpspl(0) + kpspl(1)
    kpterms(1) = kpspl(2)
    kpterms(2) = kpspl(3) + kpspl(4)

    return

end subroutine kpspl3



! ========================================================================
! Utility to find and open the supporting data files
! ========================================================================

subroutine findandopen(datafile,unitid)

    implicit none

    character(128)      :: datafile
    integer             :: unitid
    character(128)      :: hwmpath
    logical             :: havefile
    integer             :: i

    i = index(datafile,'bin')
    if (i .eq. 0) then
        inquire(file=trim(datafile),exist=havefile)
        if (havefile) open(unit=unitid,file=trim(datafile),status='old',form='unformatted')
        if (.not. havefile) then
            call get_environment_variable('HWMPATH',hwmpath)
            inquire(file=trim(hwmpath)//'/'//trim(datafile),exist=havefile)
            if (havefile) open(unit=unitid, &
                file=trim(hwmpath)//'/'//trim(datafile),status='old',form='unformatted')
        endif
        if (.not. havefile) then
            inquire(file='../Meta/'//trim(datafile),exist=havefile)
            if (havefile) open(unit=unitid, &
                file='../Meta/'//trim(datafile),status='old',form='unformatted')
        endif
    else
        inquire(file=trim(datafile),exist=havefile)
        if (havefile) open(unit=unitid,file=trim(datafile),status='old',access='stream')
        if (.not. havefile) then
            call get_environment_variable('HWMPATH',hwmpath)
            inquire(file=trim(hwmpath)//'/'//trim(datafile),exist=havefile)
            if (havefile) open(unit=unitid, &
                file=trim(hwmpath)//'/'//trim(datafile),status='old',access='stream')
        endif
        if (.not. havefile) then
            inquire(file='../Meta/'//trim(datafile),exist=havefile)
            if (havefile) open(unit=unitid, &
                file='../Meta/'//trim(datafile),status='old',access='stream')
        endif
    endif

    if (havefile) then
        return
    else
        print *,"Can not find file ",trim(datafile)
        stop
    endif

end subroutine findandopen


!!!-------------------------------------------------------------------------------
!!!  Added vectorised driver hltwim_batch.  L. Navarro (luis.navarro@colorado.edu)
!!!
!!!  Loops the scalar hltwim() over n points in Fortran so python makes one call
!!!  instead of one per point.  Out-of-range inputs (which hltwim answers with
!!!  ERROR STOP) and the model's |mlat| < 40 missing value are both returned as
!!!  xnan here.
!!!
!!!  uvmw(i,1) = geographic zonal      wind (+ east)
!!!  uvmw(i,2) = geographic meridional wind (+ north)
!!!  uvmw(i,3) = quasi-dipole zonal      wind
!!!  uvmw(i,4) = quasi-dipole meridional wind
!!!  uvmw(i,5) = quasi-dipole (magnetic) latitude [deg]
!!!  uvmw(i,6) = magnetic local time [hours]
!!!-------------------------------------------------------------------------------
subroutine hltwim_batch(iday, ut, glat, glon, kp, xnan, uvmw, n)
!f2py integer,      intent(in)  :: iday
!f2py real(kind=4), intent(in)  :: ut
!f2py real(kind=4), intent(in)  :: glat
!f2py real(kind=4), intent(in)  :: glon
!f2py real(kind=4), intent(in)  :: kp
!f2py real(kind=4), intent(in)  :: xnan
!f2py real(kind=4), intent(out) :: uvmw
!f2py integer,      intent(hide), depend(iday) :: n = len(iday)
    implicit none
    integer      :: n, i
    integer      :: iday(n)
    real(kind=4) :: ut(n), glat(n), glon(n), kp(n), xnan
    real(kind=4) :: uvmw(n, 6)
    real(kind=4) :: w(2), mw(2), gln
    real(kind=4) :: qlat, qlon, f1e, f1n, f2e, f2n
    real(kind=4), external :: mltcalc
    external hltwim, gd2qd

    do i = 1, n
        if (iday(i) .lt. 1 .or. iday(i) .gt. 366 .or. ut(i) .lt. 0.0 .or. &
            ut(i) .ge. 24.0 .or. kp(i) .lt. 0.0 .or. kp(i) .gt. 10.0 .or. &
            abs(glat(i)) .gt. 90.0) then
            uvmw(i, :) = xnan
            cycle
        end if
        gln = mod(glon(i), 360.0)
        call hltwim(iday(i), ut(i), glat(i), gln, kp(i), w, mw)
        call gd2qd(glat(i), gln, qlat, qlon, f1e, f1n, f2e, f2n)
        uvmw(i, 5) = qlat
        uvmw(i, 6) = mltcalc(qlat, qlon, iday(i), ut(i))
        if (w(1) .lt. -9000.0) then
            uvmw(i, 1:4) = xnan               ! |magnetic latitude| < 40 deg
        else
            uvmw(i, 1) = w(2)
            uvmw(i, 2) = w(1)
            uvmw(i, 3) = mw(2)
            uvmw(i, 4) = mw(1)
        end if
    end do
    return
end subroutine hltwim_batch
