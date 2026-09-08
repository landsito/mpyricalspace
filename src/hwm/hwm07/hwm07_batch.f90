!  Vectorised HWM07 driver: quiet wind (hwmqt) + disturbance (dwm07b), looping in
!  Fortran.  Added by L. Navarro.  Calls hwmqt directly (not hwm07, which would
!  auto-add DWM), so the disturbance is applied exactly once, here.
!  local solar time is derived from sec + longitude.

subroutine hwm07_batch(iyd, sec, alt, glat, glon, f107a, f107, ap, xnan, uvdd, n)
!f2py integer, intent(in)  :: iyd
!f2py real,    intent(in)  :: sec
!f2py real,    intent(in)  :: alt
!f2py real,    intent(in)  :: glat
!f2py real,    intent(in)  :: glon
!f2py real,    intent(in)  :: f107a
!f2py real,    intent(in)  :: f107
!f2py real,    intent(in)  :: ap
!f2py real,    intent(in)  :: xnan
!f2py real,    intent(out) :: uvdd
!f2py integer, intent(hide), depend(iyd) :: n = len(iyd)
    implicit none
    integer :: n, i
    integer :: iyd(n)
    real :: sec(n), alt(n), glat(n), glon(n), f107a(n), f107(n), ap(n), xnan
    real :: uvdd(n, 4)
    real :: apv(2), qw(2), dw(2), stl
    external hwmqt, dwm07b_hwm_interface

    do i = 1, n
        apv(1) = xnan
        apv(2) = ap(i)
        stl = mod(sec(i) / 3600.0 + glon(i) / 15.0, 24.0)
        call hwmqt(iyd(i), sec(i), alt(i), glat(i), mod(glon(i), 360.0), &
                   stl, f107a(i), f107(i), apv, qw)
        if (apv(2) .ge. 0.0) then
            call dwm07b_hwm_interface(iyd(i), sec(i), alt(i), glat(i), mod(glon(i), 360.0), apv, dw)
        else
            dw = 0.0
        end if
        uvdd(i, 1) = qw(2) + dw(2)      ! u = zonal (total)
        uvdd(i, 2) = qw(1) + dw(1)      ! v = meridional (total)
        uvdd(i, 3) = dw(2)             ! du = disturbance zonal
        uvdd(i, 4) = dw(1)             ! dv = disturbance meridional
    end do
    return
end subroutine hwm07_batch
