!  Replacement modules for the four TIE-GCM 2.0 modules that heelis.F imports (params_module, cons_module,
!  magfield_module, aurora_module), reduced to what the Heelis potential needs.
!
!  This software is part of the NCAR TIE-GCM.  Use is governed by the
!  Open Source Academic Research License Agreement contained in the file
!  tiegcmlicense.txt.
!
!  Derived from TIE-GCM 2.0 (params.F, cons.F, magfield.F and, above all, aurora.F: heelis_cons below is the part of
!  sub aurora_cons that feeds heelis.F, with its comments).  Written/changed by L. Navarro on 2026-09-25 for
!  mpyricalspace.  What is different from TIE-GCM:
!    - the dimensions and the magnetic grid are not needed (the driver passes the points directly), so
!      params_module / cons_module / magfield_module only declare what heelis.F needs to compile;
!    - heelis_cons takes ctpoten and by as arguments instead of reading TIE-GCM's input_module, always applies the
!      Heelis limits on By, and leaves out the AMIE branch;
!    - hvariant, heelis_prm, paper_cons and the variables they set are new: the optional latitude function of
!      Heelis et al. (1982) (see gpaper in heelis.F) and overrides of the pattern constants.
!  Needs -fdefault-real-8, as TIE-GCM does.

module params_module
    implicit none
    ! magnetic grid of TIE-GCM (params.F); heelis.F declares its own 2-d array phihm(nmlonp1,nmlat) and the
    ! routine potm with them, which the driver does not use
    integer, parameter :: nmlat = 97, nmlon = 80, nmlonp1 = nmlon + 1
    ! nlon is TIE-GCM's number of GEOGRAPHIC longitudes. flwv32 only uses it to pick the element of dlat that decides
    ! the hemisphere, dlat(max0(1,nlon/2)). The driver calls flwv32 with the points of one hemisphere only, so any
    ! element does; nlon = 2 selects the first one and keeps the index valid for a call of any length.
    integer, parameter :: nlon = 2
end module params_module

module cons_module
    use params_module, only: nmlat, nmlonp1
    implicit none
    real, parameter :: pi_dyn = 3.14159265358979312       ! pi used in dynamo calculations (cons.F)
    real :: ylatm(nmlat) = 0., ylonm(nmlonp1) = 0.         ! magnetic grid lats/lons: only potm (not used here) reads them
end module cons_module

module magfield_module
    implicit none
    real :: sunlons(1) = 0.                                ! sun's longitude in dipole coordinates: only potm reads it
end module magfield_module

module aurora_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none

    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: dtr = pi/180.                       ! degrees to radians (pi/180)
    real, parameter :: h2deg = 15.                         ! convert from hours to degrees
    integer, parameter :: isouth = 1, inorth = 2

    ! The following parameters are used by heelis.F (dimension 2 is for south, north hemispheres):
    real :: &
        theta0(2),   &  ! convection reversal boundary in radians
        offa(2),     &  ! offset of oval towards 0 MLT relative to magnetic pole (rad)
        dskofa(2),   &  ! offset of oval in radians towards 18 MLT (f(By))
        phid(2),     &  ! dayside convection entrance in MLT (subtract 12h since 0=noon) converted to radians (f(By))
        phin(2),     &  ! night convection entrance in MLT (subtract 12h since 0=noon) converted to radians (f(By))
        offc(2),     &  ! offset of convection towards 0 MLT relative to magnetic pole (rad)
        dskofc(2),   &  ! offset of convection in radians towards 18 MLT (f(By))
        psim(2),     &  ! maximum potential in the morning cell (V)
        psie(2),     &  ! minimum potential in the evening cell (V)
        pcen(2),     &  ! potential at the center (offc,dskofc) of the convection pattern (V, f(By))
        phidp0(2),   &  ! angle curvature of convection on plus side of dayside entrance (rad)
        phidm0(2),   &  ! angle curvature of convection on minus side of dayside entrance (rad)
        phinp0(2),   &  ! angle curvature of convection on plus side of nightside entrance (rad)
        phinm0(2),   &  ! angle curvature of convection on minus side of nightside entrance (rad)
        rr1(2)          ! exponential fall-off of convection from convection radius (r1 of Heelis et al. 1982 if hvariant=1)

    ! NEW (not in TIE-GCM): the variant of the latitude function and the constants of the paper's G(theta).
    integer :: hvariant = 0        ! 0 = TIE-GCM expressions (default), 1 = latitude function of Heelis et al. (1982)
    real :: r2, thetac             ! index of the polar-cap function and its phase angle [rad] (paper: 2 and 14 deg)
    real :: dth1, dth2             ! theta1 = theta0 + dth1 (equatorward), theta2 = theta0 - dth2 (poleward), the edges
                                   ! of the elliptical region that smooths the flow reversal [rad]
    real :: ga1(2), gb1(2), ga2(2), gb2(2)      ! A1, B1, A2, B2 of the paper, set by paper_cons

    ! Overrides of the constants, set from python: NaN = keep the value computed here. In order: theta0 [deg],
    ! psim, psie, pcen [kV], phid, phin [MLT h], phidp0, phidm0, phinp0, phinm0 [deg], offc, dskofc [deg], r1,
    ! r2, thetac, dth1, dth2 [deg]
    integer, parameter :: nprm = 17
    real :: ovr(nprm) = 0.

contains

    subroutine heelis_prm(prm, np)
        ! store the overrides; entries beyond np, or NaN, are not overridden
        integer, intent(in) :: np
        real, intent(in) :: prm(np)
        integer :: k
        ovr = ieee_value_nan()
        do k = 1, min(np, nprm)
            ovr(k) = prm(k)
        end do
    end subroutine heelis_prm

    real function ieee_value_nan()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        ieee_value_nan = ieee_value(1., ieee_quiet_nan)
    end function ieee_value_nan

    subroutine heelis_cons(ctpoten, byimf, variant, ok)
        ! Set the parameters of the Heelis potential from the cross-cap potential and the IMF By. This is the part of
        ! TIE-GCM's sub aurora_cons (aurora.F, "called from sub advance once per time-step") that feeds heelis.F,
        ! kept as it is there, comments included, except for the changes listed at the top of this file.
        real, intent(in) :: ctpoten   ! cross-cap potential (kV)      (e.g., 45.)
        real, intent(in) :: byimf     ! BY component of IMF (nT)      (e.g., 0.)
        integer, intent(in) :: variant
        logical, intent(out) :: ok    ! .false. if the paper's latitude function has no solution for these constants
        real :: byloc                 ! local by
        integer :: n

        !
        ! Add limits to byimf if use the Heelis convection pattern,this is to have
        ! asymmetric dawn and dusk convection cells and By effect. --Wenbin Wang 12/02/2008
        !
        byloc = byimf ! init local from original namelist input
        if (byloc > 7.0) byloc = 7.0
        if (byloc < -11.0) byloc = -11.0

        theta0(isouth) = (-3.80+8.48*(ctpoten**0.1875))*dtr
        theta0(inorth) = theta0(isouth)
        offa(isouth) = 1.0*dtr
        offa(inorth) = 1.0*dtr
        dskofa(isouth) = 0.
        dskofa(inorth) = 0.
        !
        ! The following parameters (offc through rr1) are used only in heelis
        !   potential calculation for the dynamo (see heelis.F)
        !
        !  tiegcm1.9
        !     offc(isouth) = 1.*dtr
        !     offc(inorth) = 1.*dtr
        !     dskofc(isouth) = 0.
        !     dskofc(inorth) = 0.
        !     phid(isouth) = 0.
        !     phid(inorth) = 0.
        !     phin(isouth) = 180.*dtr
        !     phin(inorth) = 180.*dtr
        !     psim(:) =  0.50 * ctpoten * 1000.
        !     psie(:) = -0.50 * ctpoten * 1000.
        !     pcen(isouth) = 0.
        !     pcen(inorth) = 0.
        !     phidp0(:) = 90.*dtr
        !     phidm0(:) = 90.*dtr
        !     phinp0(:) = 90.*dtr
        !     phinm0(:) = 90.*dtr
        !     rr1(:) = -2.6
        !
        ! tiegcm original with assymmetry and By effect
        !
        offc(isouth) = 1.1*dtr
        offc(inorth) = 1.1*dtr
        dskofc(isouth) = (-0.08 + 0.15*byloc)*dtr
        dskofc(inorth) = (-0.08 - 0.15*byloc)*dtr
        phid(isouth) = (9.39 + 0.21*byloc - 12.) * h2deg * dtr
        phid(inorth) = (9.39 - 0.21*byloc - 12.) * h2deg * dtr
        phin(isouth) = (23.50 + 0.15*byloc - 12.) * h2deg * dtr
        phin(inorth) = (23.50 - 0.15*byloc - 12.) * h2deg * dtr
        psim(:) =  0.44 * ctpoten * 1000.
        psie(:) = -0.56 * ctpoten * 1000.
        pcen(isouth) = (-0.168 + 0.027*byloc) * ctpoten * 1000.
        pcen(inorth) = (-0.168 - 0.027*byloc) * ctpoten * 1000.
        phidp0(:) = 90.*dtr
        phidm0(:) = 90.*dtr
        phinp0(:) = 90.*dtr
        phinm0(:) = 90.*dtr
        rr1(:) = -2.6

        ! ---- new: variant and overrides (not in TIE-GCM) -------------------------------------------------------
        hvariant = variant
        if (variant == 1) then              ! typical values of Heelis et al. (1982), "r1 = -4, r2 = 2, thetac = 14 deg"
            rr1(:) = -4.
            r2 = 2.
            thetac = 14.*dtr
            dth1 = 1.*dtr                   ! the paper's example has a 2 deg wide region; how it is split about
            dth2 = 1.*dtr                   ! theta0 is not stated: +-1 deg
        end if
        if (.not. ieee_is_nan(ovr(1)))  theta0(:) = ovr(1)*dtr
        if (.not. ieee_is_nan(ovr(2)))  psim(:) = ovr(2)*1000.
        if (.not. ieee_is_nan(ovr(3)))  psie(:) = ovr(3)*1000.
        if (.not. ieee_is_nan(ovr(4)))  pcen(:) = ovr(4)*1000.
        if (.not. ieee_is_nan(ovr(5)))  phid(:) = (ovr(5) - 12.)*h2deg*dtr
        if (.not. ieee_is_nan(ovr(6)))  phin(:) = (ovr(6) - 12.)*h2deg*dtr
        if (.not. ieee_is_nan(ovr(7)))  phidp0(:) = ovr(7)*dtr
        if (.not. ieee_is_nan(ovr(8)))  phidm0(:) = ovr(8)*dtr
        if (.not. ieee_is_nan(ovr(9)))  phinp0(:) = ovr(9)*dtr
        if (.not. ieee_is_nan(ovr(10))) phinm0(:) = ovr(10)*dtr
        if (.not. ieee_is_nan(ovr(11))) offc(:) = ovr(11)*dtr
        if (.not. ieee_is_nan(ovr(12))) dskofc(:) = ovr(12)*dtr
        if (.not. ieee_is_nan(ovr(13))) rr1(:) = ovr(13)
        if (.not. ieee_is_nan(ovr(14))) r2 = ovr(14)
        if (.not. ieee_is_nan(ovr(15))) thetac = ovr(15)*dtr
        if (.not. ieee_is_nan(ovr(16))) dth1 = ovr(16)*dtr
        if (.not. ieee_is_nan(ovr(17))) dth2 = ovr(17)*dtr
        do n = 1, 2                         ! flwv32 takes asin(dskofc/sqrt(offc**2+dskofc**2)): keep it defined
            if (offc(n)**2 + dskofc(n)**2 == 0.) offc(n) = 1.e-9
        end do

        ok = .true.
        if (variant == 1) call paper_cons(ok)
    end subroutine heelis_cons

    subroutine paper_cons(ok)
        ! A1, B1, A2, B2 of G(theta) (Heelis et al. 1982, eq. 3), from the requirement that G and its derivative are
        ! continuous at theta1 and theta2 and that G(theta0) = 1.
        !   outer: A1*(sin t/sin th0)**r1 meets the ellipse sqrt(1-(t-th0)**2/B1) at th1:
        !          B1 = d1*(d1 - tan(th1)/r1),  d1 = th1-th0,   A1 = sqrt(1-d1**2/B1)/(sin th1/sin th0)**r1
        !   inner: A2*H(t), H = (sin(t+thc)/sin th0)**r2 - (sin thc/sin th0)**r2, meets the ellipse at th2:
        !          B2 = d2*(d2 - 1/R),  d2 = th2-th0,  R = H'(th2)/H(th2),   A2 = sqrt(1-d2**2/B2)/H(th2)
        logical, intent(out) :: ok
        real :: th0, th1, th2, d1, d2, hh, dh, rr, c0
        integer :: n

        ok = .true.
        do n = 1, 2
            th0 = theta0(n); th1 = th0 + dth1; th2 = th0 - dth2
            d1 = th1 - th0; d2 = th2 - th0
            if (.not. (dth1 > 0. .and. dth2 > 0. .and. th2 > 0. .and. th1 < 0.5*pi .and. thetac >= 0. &
                       .and. rr1(n) < 0. .and. r2 > 0.)) then
                ok = .false.; return
            end if
            ! outer side
            gb1(n) = d1*(d1 - tan(th1)/rr1(n))
            if (gb1(n) <= d1**2) then
                ok = .false.; return
            end if
            ga1(n) = sqrt(1. - d1**2/gb1(n))/(sin(th1)/sin(th0))**rr1(n)
            ! inner side
            c0 = (sin(thetac)/sin(th0))**r2
            hh = (sin(th2 + thetac)/sin(th0))**r2 - c0
            dh = r2*(sin(th2 + thetac)/sin(th0))**(r2 - 1.)*cos(th2 + thetac)/sin(th0)
            if (hh <= 0. .or. dh <= 0.) then
                ok = .false.; return
            end if
            rr = dh/hh
            gb2(n) = d2*(d2 - 1./rr)
            if (gb2(n) <= d2**2) then
                ok = .false.; return
            end if
            ga2(n) = sqrt(1. - d2**2/gb2(n))/hh
        end do
    end subroutine paper_cons

end module aurora_module
