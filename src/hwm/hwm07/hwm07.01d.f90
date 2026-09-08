!
!  Horizontal Wind Model 08 (HWM08)
!  Version HWM071308E_DWM07B104i.01
!  See readme.txt file for detailed release notes.
!
!  AUTHORS
!    Douglas P. Drob, John T. Emmert, msishwmhelp@nrl.navy.mil
!
!  DATE
!    6 April 2009
!
!  REFERENCES
!    Drob, D. P, J. T. Emmert, G. Crowley, J. M. Picone, G. G. Shepherd, W. Skinner, 
!      Paul Hayes, R. J. Niciejewski, M. Larsen, C.Y. She, J. W. Meriwether, G. Hernandez, 
!      M. J. Jarvis, D. P. Sipler, C. A. Tepley, M. S. O'Brien, J. R. Bowman, Q. Wu, 
!      Y. Murayama, S. Kawamura, I.M. Reid, and R.A. Vincent (2008), An Empirical Model 
!      of the Earth's Horizontal Wind Fields: HWM07, J. Geophy. Res., doi:10.1029/2008JA013668.
!    Emmert, J. T., D. P. Drob, G. G. Shepherd, G. Hernandez, M. J. Jarvis, J. W. 
!      Meriwether, R. J. Niciejewski, D. P. Sipler, and C. A. Tepley (2008),
!      DWM07 global empirical model of upper thermospheric storm-induced 
!      disturbance winds, J. Geophys Res., 113, doi:10.1029/2008JA013541.
!
!================================================================================================
! Input arguments:
!        iyd - year and day as yyddd
!        sec - ut(sec)
!        alt - altitude(km)
!        glat - geodetic latitude(deg)
!        glon - geodetic longitude(deg)
!        stl - not used
!        f107a - not used
!        f107 - not used
!        ap - two element array with
!             ap(1) = not used
!             ap(2) = current 3hr ap index
!
! Output argument:
!        w(1) = meridional wind (m/sec + northward)
!        w(2) = zonal wind (m/sec + eastward)
!
!================================================================================================


subroutine hwm07(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)

    implicit none
    integer(4),intent(in)   :: iyd
    real(4),intent(in)      :: sec,alt,glat,glon,stl,f107a,f107
    real(4),intent(in)      :: ap(2)
    real(4),intent(out)     :: w(2)

    real(4)                 :: qw(2),dw(2)

    call hwmqt(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,qw)
    
    if (ap(2) .ge. 0.0) then
      call dwm07b_hwm_interface(iyd,sec,alt,glat,glon,ap,dw)
      w = qw + dw
    else
      w = qw
    endif
    
    return
    
end subroutine HWM07 

!******************************************************************************
!******************************************************************************
!                                 DATA MODULES
!******************************************************************************
!******************************************************************************
!==============================================================================
!                           Data Module HWMQTMODULE
!  Common data module for the quiet-time model definition. These parameters set by 
!  calling the subroutine loadmodel().
!   
!  Used by: HWMQT, LOADMODEL, LOADHWMQT
!=============================================================================

module hwmqtmodule

    implicit none

    integer(4)                 :: nbf              ! Count of basis terms per model level   
    integer(4)                 :: maxs             ! s seasonal
    integer(4)                 :: maxm             ! m stationary
    integer(4)                 :: maxl             ! l migrating
    integer(4)                 :: maxn             ! n latitude
    
    integer(4)                 :: p                ! B-splines order, p=4 cubic, p=3 quadratic 
    integer(4)                 :: nlev             ! e.g. Number of B-spline nodes
    integer(4)                 :: nnode            ! nlev + p

    real(8)                    :: alttns           ! Transition 1
    real(8)                    :: altsym           ! Transition 2
    real(8)                    :: altiso           ! Constant Limit

    integer(4),allocatable     :: nb(:)            ! total number of basis functions @ level
    integer(4),allocatable     :: order(:,:)       ! spectral content @ level
    real(8),allocatable        :: vnode(:)         ! Vertical Altitude Nodes
    real(8),allocatable        :: mparm(:,:)       ! Model Parameters
    real(8),allocatable        :: tparm(:,:)       ! Model Parameters

    ! Global store for quasi-static model space parameters
    ! These will change internally depending on the input parameters

    real(8)                    :: gprevious(1:5) = 1.0d32
    integer(4)                 :: gpriornb = 0
         
    real(8),allocatable        :: gfs(:,:),gfm(:,:),gfl(:,:)
    real(8),allocatable        :: gbz(:),gbm(:)
    
    real(8),allocatable        :: gzwght(:)
    integer(4)                 :: glev

    ! Miscellaneous flags and indicies
    
    integer(4)                 :: maxo
    integer(4)                 :: cseason = 0
    integer(4)                 :: cwave = 0
    integer(4)                 :: ctide = 0
    
    logical                    :: content(5) = .true.          ! Season/Waves/Tides
    logical                    :: component(0:1) = .true.      ! Compute zonal/meridional
    
    ! Initialization flags and information
    
    logical                    :: loadflag = .true.
    
end module hwmqtmodule


!=========================================================================================
!                           Data Module DWMMODULE
!  Common data module for the disturbance wind model definition. These parameters are set
!  by the calling the subroutine loadmodel().
!   
!  Used by: LOADDWM
!=========================================================================================

module dwmmodule

    implicit none

    integer(4)             :: nterm             ! Number of terms in the model
    integer(4)             :: nmax              ! Max latitudinal degree
    integer(4)             :: mmax              ! Max order of MLT variation
    integer(4)             :: nvshterm          ! # of VSH basis functions

    integer(4),allocatable :: termarr(:,:)      ! 3 x nterm index of coupled terms
    real(4),allocatable    :: coeff(:)          ! Model coefficients
    real(4),allocatable    :: vshterms(:,:)     ! VSH basis values
    real(4),allocatable    :: termval(:,:)      ! Term values to which coefficients are applied
    real(8),allocatable    :: dpbar(:,:)        ! Associated lengendre fns
    real(8),allocatable    :: dvbar(:,:)
    real(8),allocatable    :: dwbar(:,:)
    real(8),allocatable    :: mltterms(:,:)     ! MLT Fourier terms
    real(4)                :: twidth            ! Transition width of high-lat mask

    real(8), parameter      :: pi=3.1415926535897932
    real(8), parameter      :: dtor=pi/180.d0

    logical                 :: loadflag = .true.

end module dwmmodule


!=========================================================================================
!                        Data Module ALFGEOMODULE
!  Contains the associated Legendre functions and the latitude portion of the VSH
!  functions of geographic longitude and latitude. These values are used for the 
!  quiet-time component of the model and to compute magnetic coordinates and magnetic
!  local timefor the disturbance wind component. These parameters are set by calling the
!  subroutine loadmodel(). 
!   
!  Used by: LOADMODEL, LOADGD2QD, HWMUPDATE
!=========================================================================================

module alfgeomodule

    implicit none

    integer(4)              :: nmax_hwm         !Maximum degree in HWMQT
    integer(4)              :: omax_hwm         !Maximum order in HWMQT
    integer(4)              :: nmax_qd = 0      !Maximum degree of Quasi-Dipole coversion
    integer(4)              :: mmax_qd = 0      !Maximum order of Quasi-Dipole coversion
    integer(4)              :: nmax_geo         !Maximum of nmax_hwm, nmax_qd
    integer(4)              :: mmax_geo         !Maximum of omax_hwm, nmax_qd

    real(8)                 :: glat_alf = 1.d32  !Latitude of ALF arrays
    real(8),allocatable     :: gpbar(:,:),gvbar(:,:),gwbar(:,:)
    real(8),allocatable     :: spbar(:,:),svbar(:,:),swbar(:,:) !Used for MLT calculation
    
end module alfgeomodule


!=========================================================================================
!                        Data Module alfbasismodule
!  Contains the recursion coefficients needed for the computation of associated Legendre
!  functions and the latitude portion of VSH functions. These parameters are set by 
!  calling the subroutine loadmodel(). 
!   
!  Used by: ALFBASISINIT, ALFBASIS
!=========================================================================================

module alfbasismodule

    implicit none

    integer(4)           :: nmax0
    integer(4)           :: mmax0
    real(8), allocatable :: anm(:,:)
    real(8), allocatable :: cm(:)
    real(8), allocatable :: bnm(:,:)
    real(8), allocatable :: dnm(:,:)
    real(8), allocatable :: en(:)
    real(8), allocatable :: marr(:)
    real(8), allocatable :: narr(:)

end module alfbasismodule


!=========================================================================================
!                        Data Module gd2qdmodule
!  Contains the parameters for the conversion of geodetic to Quasi-Dipole coordinates via
!  a spherical harmonic expansion. These parameters are set by the calling the subroutine 
!  loadmodel(). 
!   
!  Used by: LOADGD2QD, GD2QD, MLTCALC
!=========================================================================================

module gd2qdmodule

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

    logical                  :: loadflag = .true.

end module gd2qdmodule

!***************************************************************************************************
!***************************************************************************************************
!                                 LOADING ROUTINES
!***************************************************************************************************
!***************************************************************************************************

!=========================================================================================
!                           Subroutine LOADMODEL
!  Calls all the individual initialization routines for the quiet-time model, the 
!  disturbance wind model, associated Legendre function calculation, and QD coordinate
!  calculation. Also allocates arrays in ALFGEOMODULE, which are used to store associated
!  Legendre function values for use by HWMUPDATE, GD2QD, AND MLTCALC
!   
!  Called by: HWMQT, DWM07B
!  Calls: LOADGD2QD, LOADHWMQT, LOADDWM, ALFBASISINIT
!=========================================================================================

subroutine loadmodel

    use alfgeomodule
    use hwmqtmodule
    
    implicit none
    
    integer(4)                  :: nmax_dwm, mmax_dwm
    integer(4)                  :: nmax0, mmax0

    call loadgd2qd
    call loadhwmqt
    call loaddwm(nmax_dwm, mmax_dwm)
    
    nmax_geo = max(nmax_hwm, nmax_qd)
    mmax_geo = max(omax_hwm, mmax_qd)
    nmax0 = max(nmax_geo, nmax_dwm)
    mmax0 = max(mmax_geo, mmax_dwm)
    call alfbasisinit(nmax0,mmax0)

    if (allocated(gpbar)) deallocate(gpbar,gvbar,gwbar)
    allocate(gpbar(0:nmax_geo,0:mmax_geo))
    allocate(gvbar(0:nmax_geo,0:mmax_geo))
    allocate(gwbar(0:nmax_geo,0:mmax_geo))
    gpbar = 0
    gvbar = 0
    gwbar = 0
    
    if (allocated(spbar)) deallocate(spbar,svbar,swbar)
    allocate(spbar(0:nmax_geo,0:mmax_geo))
    allocate(svbar(0:nmax_geo,0:mmax_geo))
    allocate(swbar(0:nmax_geo,0:mmax_geo))
    spbar = 0
    svbar = 0
    swbar = 0

    return

end subroutine loadmodel


!=========================================================================================
!                           Subroutine LOADHWMQT
!  Reads the parameters for the quiet-time model from the datafile. Allocates the arrays
!  in HWMQTMODULE.
!   
!  Called by: LOADMODEL
!  Calls: None
!=========================================================================================

subroutine loadhwmqt

    use hwmqtmodule
    use alfgeomodule
    implicit none

    integer(4)                     :: i,j
    integer(4)                     :: ncomp
    character(128), parameter   :: datafile = 'hwm071308e.dat'

    if (allocated(vnode)) then
        deallocate(order,nb,vnode,mparm)
        deallocate(gfs,gfm,gfl,gzwght,gbz,gbm)
    endif

    open(unit=23,file=trim(datafile),form='unformatted')
    read(23) nbf,maxs,maxm,maxl,maxn,ncomp
    read(23) nlev,p
    nnode = nlev + p
    allocate(nb(0:nnode))
    allocate(order(ncomp,0:nnode))
    allocate(vnode(0:nnode))
    read(23) vnode
    vnode(3) = 0.0
    allocate(mparm(nbf,0:nlev))
    mparm = 0.0d0
    do i = 0,nlev-p+1-2
        read(23) order(1:ncomp,i)
        read(23) nb(i)
        read(23) mparm(1:nbf,i)
    enddo
    close(23)

    ! Calculate the parity relationship permutations

    allocate(tparm(nbf,0:nlev))
    do i = 0,nlev-p+1-2
        call parity(order(:,i),nb(i),mparm(:,i),tparm(:,i))
    enddo

    ! Set transition levels
    
    alttns = vnode(nlev-2)
    altsym = vnode(nlev-1)
    altiso = vnode(nlev)

    ! Allocate the global store of quasi-static parameters
    
    maxo = max(maxs,maxm,maxl)
    omax_hwm = maxo
    nmax_hwm = maxn
    
    allocate(gfs(0:maxs,2),gfm(0:maxm,2),gfl(0:maxl,2))
    allocate(gbz(nbf),gbm(nbf))
    allocate(gzwght(0:p))
    
    gbz = 0.0d0
    gbm = 0.0d0
  
    ! Signal that the model has been initalized
    
    loadflag = .false.
    
    ! Signal a reset of the input variable comparison flags
    
    return

end subroutine loadhwmqt


!=========================================================================================
!                           Subroutine LOADDWM
!  Reads the parameters for the disturbance wind model from the datafile. Allocates the
!  arrays in DWMMODULE.
!   
!  Called by: LOADMODEL
!  Calls: None
!=========================================================================================

subroutine loaddwm(nmaxout,mmaxout)

    use dwmmodule

    implicit none

    integer(4),intent(out)     :: nmaxout, mmaxout
    character(128), parameter  :: datafile = 'dwm07b_104i.dat'

    open(unit=23,file=trim(datafile),form='unformatted')
    if (allocated(termarr)) deallocate(termarr,coeff)
    read(23) nterm, mmax, nmax
    allocate(termarr(0:2, 0:nterm-1))
    read(23) termarr
    allocate(coeff(0:nterm-1))
    read(23) coeff
    read(23) twidth
    close(23)
    
    nvshterm = ( ((nmax+1)*(nmax+2) - (nmax-mmax)*(nmax-mmax+1))/2 - 1 ) * 4 - 2*nmax
    allocate(termval(0:1, 0:nterm-1))
    allocate(dpbar(0:nmax,0:mmax),dvbar(0:nmax,0:mmax),dwbar(0:nmax,0:mmax))
    allocate(mltterms(0:mmax,0:1))
    allocate(vshterms(0:1, 0:nvshterm-1))
    dpbar = 0
    dvbar = 0
    dwbar = 0

    nmaxout = nmax
    mmaxout = mmax
    loadflag = .false.

    return

end subroutine loaddwm


!=========================================================================================
!                           Subroutine LOADGD2QD
!  Reads the parameters for conversion of geodetic to Quasi-Dipole Coordinates. Allocates
!  the arrays in GD2QDMODULE
!   
!  Called by: LOADMODEL
!  Calls: None
!=========================================================================================

subroutine loadgd2qd

    use gd2qdmodule
    use alfgeomodule

    implicit none

    character(128), parameter   :: datafile='gd2qd.dat'
    integer(4)                  :: iterm, n
    integer(4)                  :: j

    open(unit=23,file=trim(datafile),status='old',form='unformatted')
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

    do iterm = 0, nterm-1
      xcoeff(iterm) = coeff(iterm,0)
      ycoeff(iterm) = coeff(iterm,1)
      zcoeff(iterm) = coeff(iterm,2)
    enddo

    do n = 0, nmax
      normadj(n) = dsqrt(dble(n*(n+1)))
    end do

    nmax_qd = nmax
    mmax_qd = mmax
    
    loadflag = .false.
   
    return

end subroutine loadgd2qd


!======================================================================================
!                           Subroutine ALFBASISINIT
!  Computes the recursion coefficients for the computation of associated Legendre
!  functions. Allocates the arrays in ALFBASISMODULE.
!   
!  Called by: LOADMODEL
!  Calls: None
!======================================================================================

subroutine alfbasisinit(nmax0in,mmax0in)

    use alfbasismodule

    implicit none

    integer(4), intent(in) :: nmax0in, mmax0in
    integer(8)             :: n, m   ! define as integer(8) for (m,n) > 60 to avoid overflow

    nmax0 = nmax0in
    mmax0 = mmax0in

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

end subroutine alfbasisinit



!***************************************************************************************************
!***************************************************************************************************
!                               QUIET-TIME MODEL ROUTINES
!***************************************************************************************************
!***************************************************************************************************

!=========================================================================================
!                           Subroutine HWMQT
!  Evaluates the quiet-time model.
!   
!  Called by: HWM07
!  Calls: LOADMODEL, HWMUPDDATE
!=========================================================================================

subroutine HWMQT(IYD,SEC,ALT,GLAT,GLON,STL,F107A,F107,AP,W)
    
    use hwmqtmodule
    implicit none

    INTEGER,intent(in)      :: IYD
    REAL(4),intent(in)      :: SEC,ALT,GLAT,GLON,STL,F107A,F107
    REAL(4),intent(in)      :: AP(2)
    REAL(4),intent(out)     :: W(2)

    real(8)                 :: input(1:5)
    real(8)                 :: u,v
    
    input(1) = dble(mod(IYD,1000))
    input(2) = dble(sec)
    input(3) = dble(glon)
    input(4) = dble(glat)
    input(5) = dble(alt)

    if (loadflag) then
        call loadmodel()
    endif
    
    call HWMupdate(input,gprevious,gfs,gfl,gfm,gbz,gbm,gzwght,glev,gpriornb,u,v)
              
    w(1) = sngl(v)
    w(2) = sngl(u)

    return

end subroutine HWMQT


!=========================================================================================
!                           Subroutine HWMUPDATE
!  Computes the quiet-time model terms, updating only the terms that have not changed
!  since the last call. Applies the model parameters to the terms to get the winds.
!   
!  Programming Note: This subroutine is only OPENMP/THREAD SAFE when no calls to
!  LOADMODEL are made.
!
!  Called by: HWM07QT
!  Calls: ALFBASIS, VERTWGHT
!=========================================================================================

subroutine HWMupdate(input,previous,fs,fl,fm,bz,bm,zwght,lev,priornb,u,v)

    use hwmqtmodule
    use alfgeomodule
    implicit none
    
    real(8),intent(in)              :: input(5) ! jday,utsec,glon,glat,alt
    real(8),intent(inout)           :: previous(5)
    real(8),intent(inout)           :: fs(0:maxs,2)
    real(8),intent(inout)           :: fm(0:maxm,2)
    real(8),intent(inout)           :: fl(0:maxl,2)
    
    real(8),intent(inout)           :: bz(nbf)
    real(8),intent(inout)           :: bm(nbf)

    real(8),intent(inout)           :: zwght(0:p)
    integer,intent(inout)           :: lev
    integer,intent(inout)           :: priornb   ! the previous number of basis function used last p
   
    real(8),intent(out)             :: u,v
 
    ! Local variables

    real(8)                   :: cs,ss,cm,sm,cl,sl
    real(8)                   :: cmcs,smcs,cmss,smss
    real(8)                   :: clcs,slcs,clss,slss	
    real(8)                   :: AA,BB,CC,DD
    real(8)                   :: vb,wb
    
    integer(4)                :: b,c,d,m,n,s,l
    
    integer(4)                :: amaxs,amaxn
    integer(4)                :: pmaxm,pmaxs,pmaxn
    integer(4)                :: tmaxl,tmaxs,tmaxn
 
    logical                   :: refresh(5)
       
    real(8),parameter         :: twoPi = 2.0d0*3.1415926535897932384626433832795d0
    real(8),parameter         :: deg2rad = twoPi/360.0d0
    
! ====================================================================
! Update VSH model terms based on any change in the input parameters
! ====================================================================
    
    refresh(1:5) = .false.

    ! Seasonal variations

    if (input(1) .ne. previous(1)) then
        AA = input(1)*twoPi/365.25d0
        do s = 0,MAXS
            BB = dble(s)*AA
            fs(s,1) = dcos(BB)
            fs(s,2) = dsin(BB)
        enddo
        refresh(1:5) = .true.
        previous(1) = input(1)
    endif

    ! Hourly time changes, tidal variations

    if (input(2) .ne. previous(2) .or. input(3) .ne. previous(3)) then
        AA = mod(input(2)/3600.d0 + input(3)/15.d0 + 48.d0,24.d0)
        BB = AA*twoPi/24.d0
        do l = 0,MAXL
            CC = dble(l)*BB
            fl(l,1) = dcos(CC)
            fl(l,2) = dsin(CC)
        enddo
        refresh(3) = .true.
        previous(2) = input(2)
    endif

    ! Longitudinal variations, planetary waves
    
    if (input(3) .ne. previous(3)) then
        AA = input(3)*deg2rad
        do m = 0,MAXM
            BB = dble(m)*AA
            fm(m,1) = dcos(BB)
            fm(m,2) = dsin(BB)
        enddo
        refresh(2) = .true.
        previous(3) = input(3)
    endif

    ! Latitude

    if (input(4) .ne. glat_alf) then
        AA = (90.0d0 - input(4))*deg2rad        ! theta = colatitude in radians
        call alfbasis(nmax_geo,mmax_geo,AA,gpbar,gvbar,gwbar)
        refresh(1:4) = .true.
        glat_alf = input(4)
        previous(4) = input(4)
    endif

    ! Altitude
    
    if (input(5) .ne. previous(5)) then
        call vertwght(input(5),zwght,lev)
        previous(5) = input(5)
    endif

    ! ====================================================================
    ! Calculate the VSH functions
    ! ====================================================================
       
    u = 0.0d0
    v = 0.0d0

    do b = 0,p

        if (zwght(b) .eq. 0.d0) cycle
        
        d = b + lev

        if (priornb .ne. nb(d)) refresh(1:5) = .true. ! recalculate basis functions
        priornb = nb(d)

        if (.not. any(refresh)) then
            c = nb(d)
            if (component(0)) u = u + zwght(b)*dot_product(bz(1:c),mparm(1:c,d))
            if (component(1)) v = v + zwght(b)*dot_product(bz(1:c),tparm(1:c,d))
            cycle
        endif
    
        amaxs = order(1,d)
        amaxn = order(2,d)
        pmaxm = order(3,d)
        pmaxs = order(4,d)
        pmaxn = order(5,d)
        tmaxl = order(6,d)
        tmaxs = order(7,d)
        tmaxn = order(8,d)

        c = 1

        ! ------------- Seasonal - Zonal average (m = 0) ----------------

        if (refresh(1) .and. content(1)) then  
            do n = 1,amaxn                   ! s = 0
                bz(c)   = -0.5d0*gvbar(n,0)  ! Cr     A    
                bz(c+1) =  0.0d0             ! Br     B
                c = c + 2
            enddo             
            do s = 1,amaxs                   ! Seasonal variations
                cs = fs(s,1)
                ss = fs(s,2)
                do n = s,amaxn
                    vb = gvbar(n,s)
                    wb = gwbar(n,s)
                    bz(c) = -vb*cs     ! Cr     A
                    bz(c+1) =  vb*ss   ! Ci     B
                    bz(c+2) = -wb*ss   ! Br     C
                    bz(c+3) = -wb*cs   ! Bi     D
                    c = c + 4
                enddo
            enddo
            cseason = c
        else
            c = cseason
        endif
            
        ! ---------------- Stationary planetary waves --------------------

        if (refresh(2) .and. content(2)) then
            do m = 1,pmaxm
               cm = fm(m,1)
               sm = fm(m,2)
               do n = m,pmaxn           ! s = 0
                    vb = gvbar(n,m)
                    wb = gwbar(n,m)                                          
                    bz(c) =   -vb*cm    ! Cr * (cm) * -vb   A
                    bz(c+1) =  vb*sm    ! Ci * (sm) *  vb   B
                    bz(c+2) = -wb*sm	! Br * (sm) * -wb   C
                    bz(c+3) = -wb*cm	! Bi * (cm) * -wb   D
                    c = c + 4
               enddo
               do s = 1,pmaxs
                  cs = fs(s,1)
                  ss = fs(s,2)                  
                  do n = m,pmaxn
                     vb = gvbar(n,m)
                     wb = gwbar(n,m)
                     bz(c) =   -vb*cm*cs	! Crc * (cmcs) * -vb   A
                     bz(c+1) =  vb*sm*cs    ! Cic * (smcs) *  vb   B
                     bz(c+2) = -wb*sm*cs	! Brc * (smcs) * -wb   C
                     bz(c+3) = -wb*cm*cs	! Bic * (cmcs) * -wb   D
                     bz(c+4) = -vb*cm*ss	! Crs * (cmss) * -vb   E
                     bz(c+5) =  vb*sm*ss    ! Cis * (smss) *  vb   F
                     bz(c+6) = -wb*sm*ss	! Brs * (smss) * -wb   G
                     bz(c+7) = -wb*cm*ss	! Bis * (cmss) * -wb   H
                     c = c + 8
                  enddo
               enddo
               cwave = c      
            enddo        
        else
            c = cwave  
        endif

        ! ---------------- Migrating Solar Tides ---------------------

        if (refresh(3) .and. content(3)) then
             do l = 1,tmaxl
               cl = fl(l,1)
               sl = fl(l,2)
               do n = l,tmaxn           ! s = 0
                    vb = gvbar(n,l)
                    wb = gwbar(n,l)                                          
                    bz(c) =   -vb*cl    ! Cr * (cl) * -vb
                    bz(c+1) =  vb*sl    ! Ci * (sl) *  vb
                    bz(c+2) = -wb*sl	! Br * (sl) * -wb
                    bz(c+3) = -wb*cl	! Bi * (cl) * -wb
                    c = c + 4
               enddo
               do s = 1,tmaxs
                  cs = fs(s,1)
                  ss = fs(s,2)
                  do n = l,tmaxn
                     vb = gvbar(n,l)
                     wb = gwbar(n,l)
                     bz(c) =   -vb*cl*cs	! Crc * (clcs) * -vb
                     bz(c+1) =  vb*sl*cs    ! Cic * (slcs) *  vb
                     bz(c+2) = -wb*sl*cs	! Brc * (slcs) * -wb
                     bz(c+3) = -wb*cl*cs	! Bic * (clcs) * -wb
                     bz(c+4) = -vb*cl*ss	! Crs * (clss) * -vb
                     bz(c+5) =  vb*sl*ss    ! Cis * (slss) *  vb
                     bz(c+6) = -wb*sl*ss	! Brs * (slss) * -wb
                     bz(c+7) = -wb*cl*ss	! Bis * (clss) * -wb
                     c = c + 8
                  enddo
               enddo
               ctide = c
            enddo
        else
            c = ctide  
        endif
            
        ! ---------------- Non-Migrating Solar Tides ------------------
        
        ! TBD
            
        c = c - 1
        
        ! ====================================================================
        ! Calculate the wind components 
        ! ====================================================================

        if (component(0)) u = u + zwght(b)*dot_product(bz(1:c),mparm(1:c,d))
        if (component(1)) v = v + zwght(b)*dot_product(bz(1:c),tparm(1:c,d))

    enddo

    return 

end subroutine HWMupdate


!=========================================================================================
!                           Subroutine VERTWGHT
!  Updates the vertical basis functions of the quiet-time model.
!
!  Called by: HWMUPDATE
!  Calls: None.
!=========================================================================================

subroutine vertwght(alt,wght,iz)

    use hwmqtmodule
    implicit none
    
    real(8),intent(in)  :: alt
    real(8),intent(out) :: wght(4)
    integer(4),intent(out) :: iz
    
    real(8)             :: we(0:4)

    real(8)             :: e1(0:4) = &
        (/1.d0, 0.428251121076233d0,0.192825112107623d0,0.484304932735426d0,0.0d0/)
    real(8)             :: e2(0:4) = &
        (/0.d0, 0.571748878923767d0,0.807174887892377d0,-0.484304932735426d0,1.0d0/)

    real(8),parameter   :: H = 60.0d0
        
    iz = findspan(nnode-p-1_4,p,alt,vnode) - p

    iz = min(iz,26)
    
    wght(1) = bspline(p,nnode,vnode,iz,alt)
    wght(2) = bspline(p,nnode,vnode,iz+1_4,alt)
    if (iz .le. 25) then
        wght(3) = bspline(p,nnode,vnode,iz+2_4,alt)
        wght(4) = bspline(p,nnode,vnode,iz+3_4,alt)
        return
    endif
        if (alt .gt. 250.0d0) then
            we(0) = 0.0d0
            we(1) = 0.0d0
            we(2) = 0.0d0
            we(3) = exp(-(alt - 250.0d0)/H) 
            we(4) = 1.0d0
        else
            we(0) = bspline(p,nnode,vnode,iz+2_4,alt)
            we(1) = bspline(p,nnode,vnode,iz+3_4,alt)
            we(2) = bspline(p,nnode,vnode,iz+4_4,alt)
            we(3) = 0.0d0
            we(4) = 0.0d0
        endif
        wght(3) = dot_product(we,e1)
        wght(4) = dot_product(we,e2)
        
    return
    
contains

    function bspline(p,m,V,i,u)

        implicit none
        
        real(8)     :: bspline
        integer(4)  :: p,m
        real(8)     :: V(0:m)
        integer(4)  :: i
        real(8)     :: u
        
        real(8)     :: N(0:p+1)
        real(8)     :: Vleft,Vright
        real(8)     :: saved,temp
        integer(4)  :: j,k
                
        if ((i .eq. 0) .and. (u .eq. V(0))) then
            bspline = 1.d0
            return
        endif
        
        if ((i .eq. (m-p-1)) .and. (u .eq. V(m))) then
            bspline = 1.d0
            return
        endif

        if (u .lt. V(i) .or. u .ge. V(i+p+1)) then
            bspline = 0.d0
            return
        endif
        
        N = 0.0d0
        do j = 0,p
            if (u .ge. V(i+j) .and. u .lt. V(i+j+1)) then
                N(j) = 1.0d0
            else
                N(j) = 0.0d0
            endif
        enddo
        
        do k = 1,p
            if (N(0) .eq. 0.d0) then
                saved = 0.d0
            else
                saved = ((u - V(i))*N(0))/(V(i+k) - V(i))
            endif
            do j = 0,p-k
                Vleft = V(i+j+1)
                Vright = V(i+j+k+1)
                if (N(j+1) .eq. 0.d0) then
                    N(j) = saved
                    saved = 0.d0
                else
                    temp = N(j+1)/(Vright - Vleft)
                    N(j) = saved + (Vright - u)*temp
                    saved = (u - Vleft)*temp
                endif
            enddo
        enddo
        
        bspline = N(0)

        return

    end function bspline

    ! =====================================================
    ! Function to locate the knot span
    ! =====================================================

    integer(4) function findspan(n,p,u,V)

        implicit none
        
        integer(4),intent(in)   :: n,p
        real(8),intent(in)      :: u
        real(8),intent(in)      :: V(0:n+1)
        integer(4)              :: low,mid,high
        
        if (u .ge. V(n+1)) then
            findspan = n
            return
        endif
        
        low = p
        high = n+1
        mid = (low + high)/2

        do while (u .lt. V(mid) .or. u .ge. V(mid + 1))
            if (u .lt. V(mid)) then
                high = mid
            else
                low = mid
            endif
            mid = (low + high)/2
        end do
    
        findspan = mid
        return

    end function findspan
    
end subroutine vertwght


!=========================================================================================
!                           Subroutine ALFBASIS
!  Concurrently computes the latitude-dependent portions of scalar and vector spherical
!  harmonics.
!
!  Programming Note: This subroutine is only OPENMP/THREAD SAFE when no calls to
!  LOADMODEL are made.
!
!  Called by: HWMUPDATE, DWM07B
!  Calls: None.
!=========================================================================================

subroutine alfbasis(nmax,mmax,theta,P,V,W)

    use alfbasismodule

    implicit none

    integer(4), intent(in)  :: nmax, mmax
    real(8), intent(in)     :: theta
    real(8), intent(out)    :: P(0:nmax,0:mmax)
    real(8), intent(out)    :: V(0:nmax,0:mmax)
    real(8), intent(out)    :: W(0:nmax,0:mmax)

    integer(8)              :: n, m
    real(8)                 :: x, y
    real(8), parameter      :: p00=0.70710678118654746d0

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



!***************************************************************************************************
!***************************************************************************************************
!                             DISTURBANCE WIND MODEL ROUTINES
!***************************************************************************************************
!***************************************************************************************************

!=========================================================================================
!                           Subroutine DWM07B_HWM_INTERFACE
!  Using HWM inputs, computes Quasi-dipole latitude and local time, and Kp. Retrieves DWM
!  results for these conditions, converts to geographic directions, and applies artificial
!  height profile.
!
!  Called by: HWM07
!  Calls: AP2KP, GD2QD, MLTCALC, DWM07B
!=========================================================================================

subroutine dwm07b_hwm_interface(IYD,SEC,ALT,GLAT,GLON,AP,DW)

    implicit none

    INTEGER,intent(in)      :: IYD
    REAL(4),intent(in)      :: SEC,ALT,GLAT,GLON
    REAL(4),intent(in)      :: AP(2)
    REAL(4),intent(out)     :: DW(2)

    real(4), save           :: day, ut, mlat, mlon, mlt, kp
    real(4)                 :: mmpwind, mzpwind
    real(4), save           :: f1e, f1n, f2e, f2n
    real(4), save           :: glatlast=1.0e16, glonlast=1.0e16
    real(4), save           :: daylast=1.0e16, utlast=1.0e16, aplast=1.0e16
    real(4), parameter      :: talt=125.0, twidth=5.0

    real(4), external       :: ap2kp, mltcalc

    !CONVERT AP TO KP
    if (ap(2) .ne. aplast) then
      kp = ap2kp(ap(2))
    endif

    !CONVERT GEO LAT/LON TO QD LAT/LON
    if ((glat .ne. glatlast) .or. (glon .ne. glonlast)) then
      call gd2qd(glat,glon,mlat,mlon,f1e,f1n,f2e,f2n)
    endif

    !COMPUTE QD MAGNETIC LOCAL TIME (LOW-PRECISION)
    day = real(mod(iyd,1000))
    ut = sec / 3600.0
    if ((day .ne. daylast) .or. (ut .ne. utlast) .or. &
        (glat .ne. glatlast) .or. (glon .ne. glonlast)) then
      mlt = mltcalc(mlat,mlon,day,ut)
    endif

    !RETRIEVE DWM WINDS
    call dwm07b(mlt, mlat, kp, mmpwind, mzpwind)

    !CONVERT TO GEOGRAPHIC COORDINATES
    dw(1) = f2n*mmpwind + f1n*mzpwind
    dw(2) = f2e*mmpwind + f1e*mzpwind

    !APPLY HEIGHT PROFILE
    dw = dw / (1 + exp(-(alt - talt)/twidth))

    glatlast = glat
    glonlast = glon
    daylast = day
    utlast = ut
    aplast = ap(2)
    
    return

end subroutine dwm07b_hwm_interface


!=========================================================================================
!                           Subroutine AP2KP
!  Converts ap values to Kp values, via linear interpolation on the lookup table.
!
!  Called by: DWM07B_HWM_INTERFACE
!  Calls: None.
!=========================================================================================

function ap2kp(ap0)

  real(4), parameter :: apgrid(0:27) = (/0.,2.,3.,4.,5.,6.,7.,9.,12.,15.,18., &
                                         22.,27.,32.,39.,48.,56.,67.,80.,94., &
                                       111.,132.,154.,179.,207.,236.,300.,400./)
  real(4), parameter :: kpgrid(0:27) = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11., &
                                         12.,13.,14.,15.,16.,17.,18.,19.,20.,21., &
                                         22.,23.,24.,25.,26.,27./) / 3.0
  real(4)            :: ap0, ap, ap2kp
  integer(4)         :: i


  ap = ap0
  if (ap .lt. 0) ap = 0
  if (ap .gt. 400) ap = 400

  i = 1
  do while (ap .gt. apgrid(i))
    i = i + 1
  end do
  if (ap .eq. apgrid(i)) then
    ap2kp = kpgrid(i)
  else
    ap2kp = kpgrid(i-1) + (ap - apgrid(i-1)) / (3.0 * (apgrid(i) - apgrid(i-1)))
  end if

  return

end function ap2kp


!=========================================================================================
!                           Subroutine GD2QD
!  Converts geodetic coordinates to Quasi-Dipole coordinates (Richmond, J. Geomag. 
!  Geoelec., 1995, p. 191), using a spherical harmonic representation.
!
!  Programming Note: This subroutine is only OPENMP/THREAD SAFE when no calls to
!  LOADMODEL are made.
!
!  Called by: DWM07B_HWM_INTERFACE
!  Calls: ALFBASIS
!=========================================================================================

subroutine gd2qd(glatin,glon,qlat,qlon,f1e,f1n,f2e,f2n)

    use gd2qdmodule
    use alfgeomodule

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
            
    if (loadflag) call loadmodel
    
    glat = dble(glatin)
    if (glat .ne. glat_alf) then
      theta = (90.d0 - glat) * dtor
      call alfbasis(nmax,mmax,theta,gpbar,gvbar,gwbar)
      glat_alf = glat
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


!=========================================================================================
!                           Function MLTCALC
!  Low-precision calculation of magnetic local time.
!
!  Programming Note: This subroutine is only OPENMP/THREAD SAFE when no calls to
!  LOADMODEL are made.
!
!  Called by: DWM07B_HWM_INTERFACE
!  Calls: ALFBASIS
!=========================================================================================

function mltcalc(qlat,qlon,day,ut)

    use gd2qdmodule
    use alfgeomodule

    implicit none

    real(4), intent(in)      :: qlat, qlon, day, ut
    real(4)                  :: mltcalc

    integer(4)               :: n, m, i
    real(8)                  :: asunglat, asunglon, asunqlon
    real(8)                  :: glat, theta, phi
    real(8)                  :: mphi, cosmphi, sinmphi
    real(8)                  :: x, y
    real(8)                  :: cosqlat, cosqlon, sinqlon
    real(8)                  :: qlonrad
            
    if (loadflag) call loadmodel
    
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


!=========================================================================================
!                           Subroutine DWM07B
!  Evaluates the disturbance wind model.
!
!  Programming Note: This subroutine is only OPENMP/THREAD SAFE when no calls to
!  LOADMODEL are made.
!
!  Called by: DWM07B_HWM_INTERFACE
!  Calls: ALFBASIS, DWM_KPSPL3
!=========================================================================================

subroutine dwm07b(mlt, mlat, kp, mmpwind, mzpwind)

    use dwmmodule
    implicit none

    real(4),intent(in)        :: mlt            !Magnetic local time (hours)
    real(4),intent(in)        :: mlat           !Magnetic latitude (degrees)
    real(4),intent(in)        :: kp             !3-hour Kp

    real(4),intent(out)       :: mmpwind        !Mer. disturbance wind (+north, QD coordinates)
    real(4),intent(out)       :: mzpwind        !Zon. disturbance wind (+east, QD coordinates)

    ! Local variables
    integer(4)                :: iterm, ivshterm, n, m
    real(4)                   :: termvaltemp(0:1)
    real(4),save              :: kpterms(0:2)
    real(4)                   :: latwgtterm
    real(4),save              :: mltlast=1.e16, mlatlast=1.e16, kplast=1.e16
    real(8)                   :: theta, phi, mphi

    real(4),external          :: dwm_latwgt2

    
    !LOAD MODEL PARAMETERS IF NECESSARY
    if (loadflag) call loadmodel

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

    !COMPUTE KP TERMS
    if (kp .ne. kplast) then
      call dwm_kpspl3(kp, kpterms)
    endif

    !COMPUTE LATITUDINAL WEIGHTING TERM
    latwgtterm = dwm_latwgt2(mlat, mlt, kp, twidth)

    !GENERATE COUPLED TERMS
    do iterm = 0, nterm-1
      termvaltemp = (/1.0, 1.0/)
      if (termarr(0,iterm) .ne. 999) termvaltemp = termvaltemp * vshterms(0:1,termarr(0,iterm))
      if (termarr(1,iterm) .ne. 999) termvaltemp = termvaltemp * kpterms(termarr(1,iterm))
      if (termarr(2,iterm) .ne. 999) termvaltemp = termvaltemp * latwgtterm
      termval(0:1,iterm) = termvaltemp(0:1)
    enddo
  
    !APPLY COEFFICIENTS
    mmpwind = dot_product(coeff, termval(0,0:nterm-1))
    mzpwind = dot_product(coeff, termval(1,0:nterm-1))

    mlatlast = mlat
    mltlast = mlt
    kplast = kp

    return

end subroutine dwm07b


!=========================================================================================
!                           Subroutine DWM_KPSPL3
!  Computes a quadratic spline basis for the Kp dependence. Covers the interval [0,8],
!  with 2 interior nodes at Kp = 2 and 5. The basis is constrained to have zero slope at
!  the bounds.
!
!  Called by: DWM07B
!  Calls: None.
!=========================================================================================

subroutine dwm_kpspl3(kp, kpterms)

    implicit none

    real(4), intent(in)       :: kp
    real(4), intent(out)      :: kpterms(0:2)

    integer(4)                :: i, j
    real(4)                   :: x, kpspl(0:6)    
    real(4), parameter        :: node(0:7)=(/-10., -8., 0., 2., 5., 8., 18., 20./)

    x = max(kp, 0.0)
    x = min(x,  8.0)

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

end subroutine dwm_kpspl3


!=========================================================================================
!                           Function DWM_LATWGT2
!  Computes a latitude dependent weighting function that goes to Zero at low latitudes and
!  one at high latitudes. The transition Is an exponential s-curve, with the transition
!  latitude determined by a MLT/Kp model.
!
!  Called by: DWM07B
!  Calls: None.
!=========================================================================================

function dwm_latwgt2(mlat, mlt, kp0, twidth)

    implicit none

    real(4)                   :: dwm_latwgt2
    real(4)                   :: mlat, mlt, kp0, kp, twidth
    real(4)                   :: mltrad, sinmlt, cosmlt, tlat

    real(4), parameter :: coeff(0:5) = (/ 65.7633,  -4.60256,  -3.53915,  &
                                         -1.99971,  -0.752193,  0.972388 /)
    real(4), parameter :: pi=3.141592653590
    real(4), parameter :: dtor=pi/180.d0


    mltrad = mlt * 15.0 * dtor
    sinmlt = sin(mltrad)
    cosmlt = cos(mltrad)
    kp = max(kp0, 0.0)
    kp = min(kp,  8.0)
    tlat = coeff(0) + coeff(1)*cosmlt + coeff(2)*sinmlt +   &
           kp*(coeff(3) + coeff(4)*cosmlt + coeff(5)*sinmlt)
    dwm_latwgt2 = 1.0 / ( 1 + exp(-(abs(mlat)-tlat)/twidth) )

    return

end function dwm_latwgt2

!=========================================================================================
!
!=========================================================================================


subroutine parity(order,nb,mparm,tparm)

    implicit none

    integer(4),intent(in)     :: order(8)
    integer(4),intent(in)     :: nb
    real(8),intent(in)        :: mparm(nb)
    real(8),intent(out)       :: tparm(nb)

    integer(4)                :: c,m,n,s,l
    
    integer(4)                :: amaxs,amaxn
    integer(4)                :: pmaxm,pmaxs,pmaxn
    integer(4)                :: tmaxl,tmaxs,tmaxn

    amaxs = order(1)
    amaxn = order(2)
    pmaxm = order(3)
    pmaxs = order(4)
    pmaxn = order(5)
    tmaxl = order(6)
    tmaxs = order(7)
    tmaxn = order(8)

    c = 1

    ! ------------- Seasonal - Zonal average (m = 0) ----------------

    do n = 1,AMAXN
        tparm(c) = mparm(c+1)
        tparm(c+1) = -mparm(c)
        c = c + 2
    enddo                
    do s = 1,AMAXS
        do n = s,AMAXN
            tparm(c) = mparm(c+2)
            tparm(c+1) = mparm(c+3)
            tparm(c+2) = -mparm(c)
            tparm(c+3) = -mparm(c+1)
            c = c + 4
        enddo
    enddo
    
    ! ---------------- Stationary planetary waves --------------------

    do m = 1,pmaxm   
        do n = m,pmaxn
            tparm(c) = mparm(c+2)
            tparm(c+1) = mparm(c+3)
            tparm(c+2) = -mparm(c)
            tparm(c+3) = -mparm(c+1)    
            c = c + 4
        enddo
        do s = 1,pmaxs
            do n = m,pmaxn            
                tparm(c) = mparm(c+2)
                tparm(c+1) = mparm(c+3)
                tparm(c+2) = -mparm(c)
                tparm(c+3) = -mparm(c+1)
                tparm(c+4) = mparm(c+6)
                tparm(c+5) = mparm(c+7)
                tparm(c+6) = -mparm(c+4)
                tparm(c+7) = -mparm(c+5)              
                c = c + 8          
            enddo  
        enddo

    enddo

    ! ---------------- Migrating Solar Tides ---------------------

    do l = 1,tmaxl
        do n = l,tmaxn
            tparm(c) = mparm(c+2)
            tparm(c+1) = mparm(c+3)
            tparm(c+2) = -mparm(c)
            tparm(c+3) = -mparm(c+1)    
            c = c + 4  
        enddo   
        do s = 1,tmaxs
            do n = l,tmaxn
                tparm(c) = mparm(c+2)
                tparm(c+1) = mparm(c+3)
                tparm(c+2) = -mparm(c)
                tparm(c+3) = -mparm(c+1)
                tparm(c+4) = mparm(c+6)
                tparm(c+5) = mparm(c+7)
                tparm(c+6) = -mparm(c+4)
                tparm(c+7) = -mparm(c+5)              
                c = c + 8          
            enddo
        enddo
    enddo
    
    return

end subroutine parity
