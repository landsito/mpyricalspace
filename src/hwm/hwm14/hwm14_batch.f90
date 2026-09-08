!  Vectorised HWM14 driver: total wind (hwm14) + disturbance wind (dwm07),
!  looping in Fortran so python makes one call instead of one per point.
!  Added by L. Navarro.

subroutine hwm14_batch(iyd, sec, alt, glat, glon, ap, xnan, uvdd, n)
!f2py integer,      intent(in)  :: iyd
!f2py real(kind=4), intent(in)  :: sec
!f2py real(kind=4), intent(in)  :: alt
!f2py real(kind=4), intent(in)  :: glat
!f2py real(kind=4), intent(in)  :: glon
!f2py real(kind=4), intent(in)  :: ap
!f2py real(kind=4), intent(in)  :: xnan
!f2py real(kind=4), intent(out) :: uvdd
!f2py integer,      intent(hide), depend(iyd) :: n = len(iyd)
    implicit none
    integer :: n, i
    integer :: iyd(n)
    real(kind=4) :: sec(n), alt(n), glat(n), glon(n), ap(n), xnan
    real(kind=4) :: uvdd(n,4)
    real(kind=4) :: w(2), d(2), apv(2)
    external hwm14, dwm07

    apv(1) = xnan
    do i = 1, n
        apv(2) = ap(i)
        call hwm14(iyd(i), sec(i), alt(i), glat(i), mod(glon(i), 360.0), &
                   xnan, xnan, xnan, apv, w)
        call dwm07(iyd(i), sec(i), alt(i), glat(i), mod(glon(i), 360.0), apv, d)
        uvdd(i,1) = w(2)      ! u  (hwm14 returns w = [meridional, zonal])
        uvdd(i,2) = w(1)      ! v
        uvdd(i,3) = d(2)      ! du
        uvdd(i,4) = d(1)      ! dv
    end do
    return
end subroutine hwm14_batch
