!  Vectorised HWM93 driver: loops GWS5 in Fortran.  Added by L. Navarro.
!  local solar time is derived from sec + longitude (as the check drivers do);
!  storm-time winds are already inside GWS5 via the ap array.

subroutine gws5_batch(iyd, sec, alt, glat, glong, f107a, f107, ap, xnan, uvdd, n)
!f2py integer, intent(in)  :: iyd
!f2py real,    intent(in)  :: sec
!f2py real,    intent(in)  :: alt
!f2py real,    intent(in)  :: glat
!f2py real,    intent(in)  :: glong
!f2py real,    intent(in)  :: f107a
!f2py real,    intent(in)  :: f107
!f2py real,    intent(in)  :: ap
!f2py real,    intent(in)  :: xnan
!f2py real,    intent(out) :: uvdd
!f2py integer, intent(hide), depend(iyd) :: n = len(iyd)
    implicit none
    integer :: n, i
    integer :: iyd(n)
    real :: sec(n), alt(n), glat(n), glong(n), f107a(n), f107(n), ap(n), xnan
    real :: uvdd(n, 4)
    real :: apv(2), ww(2), stl
    external gws5

    do i = 1, n
        apv(1) = ap(i)
        apv(2) = ap(i)
        stl = mod(sec(i) / 3600.0 + glong(i) / 15.0, 24.0)
        call gws5(iyd(i), sec(i), alt(i), glat(i), mod(glong(i), 360.0), &
                  stl, f107a(i), f107(i), apv, ww)
        uvdd(i, 1) = ww(2)      ! u  (GWS5 returns w = [meridional, zonal])
        uvdd(i, 2) = ww(1)      ! v
        uvdd(i, 3) = xnan       ! HWM93 has no separate disturbance model
        uvdd(i, 4) = xnan
    end do
    return
end subroutine gws5_batch
