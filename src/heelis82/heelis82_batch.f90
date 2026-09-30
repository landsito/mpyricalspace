!  Vectorised driver for the Heelis convection potential as TIE-GCM 2.0 computes it: electric potential at n
!  points, looping in Fortran so python makes one call instead of one per point.  Written by L. Navarro.
!
!  Every row i is an independent point with its own drivers:
!     cp    cross-cap potential [kV]; NaN -> derived from kp (TIE-GCM's empirical relation, kp_module)
!     kp    Kp index (0..9), used only where cp is NaN
!     by    IMF By [nT] (GSM); limited to [-11, 7] nT as TIE-GCM does for the Heelis model
!     mlat  magnetic latitude [deg], signed (negative = southern hemisphere)
!     mlt   magnetic local time [h]
!  (|mlat| may exceed 90 by up to 1e-6 deg, a rounding error of the caller's grid; it is clipped to 90.)
!  variant: 0 = TIE-GCM expressions, 1 = latitude function of Heelis et al. (1982) (see heelis.F)
!  prm:     overrides of the pattern constants, in the order listed in aurora_module (NaN = keep); may be empty
!
!  Outputs: epot [kV]; cpused [kV], the cross-cap potential each row used (NaN where the row could not be computed).
!  Poleward of 30 deg (|mlat| > 30, TIE-GCM's own limit in sub potm) epot is computed; at or below it TIE-GCM sets the
!  potential to zero, and so does this driver.  NaN where any input is NaN, Kp is outside 0..9 (and cp is NaN) or
!  cp <= 0.
!  ier is a bit mask: 1 = some row had cp <= 0 or Kp outside 0..9 (those rows are NaN)
!                     2 = the constants of the paper's latitude function have no solution (variant 1): NaN
!
!  flwv32 is called once per run of consecutive rows that share cp, by and hemisphere (a whole grid at one time step
!  is one or two runs), because it takes the hemisphere from its first point.

subroutine heelis82_batch(cp, kp, by, mlat, mlt, variant, prm, xnan, epot, cpused, ier, n, np)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use kp_module, only: ctpoten_from_kp
    use aurora_module, only: heelis_cons, heelis_prm, dtr, h2deg
    use heelis_module, only: flwv32
    implicit none
    integer :: n, np, variant, ier
    real(kind=8) :: cp(n), kp(n), by(n), mlat(n), mlt(n), prm(np)
    real(kind=8) :: xnan, epot(n), cpused(n)

    real(kind=8), allocatable :: ce(:), dlat(:), dlon(:), ratio(:), pot(:)
    integer, allocatable :: idx(:), iflag(:)
    logical, allocatable :: ok(:)
    logical :: north, good
    real(kind=8) :: c
    integer :: i, j, k, m

    ier = 0
    epot = xnan
    cpused = xnan
    if (allocated(ce)) deallocate (ce, dlat, dlon, ratio, pot, idx, iflag, ok)     ! -fno-automatic keeps locals
    allocate (ce(n), dlat(n), dlon(n), ratio(n), pot(n + 1), idx(n), iflag(n), ok(n))
    call heelis_prm(prm, np)

    ! the cross-cap potential of every row: given, or from Kp
    do i = 1, n
        ok(i) = .false.
        ce(i) = xnan
        if (ieee_is_nan(by(i)) .or. ieee_is_nan(mlat(i)) .or. ieee_is_nan(mlt(i))) cycle
        if (abs(mlat(i)) > 90.000001_8) cycle           ! a rounding error above 90 (np.arange) is clipped below
        if (.not. ieee_is_nan(cp(i))) then
            c = cp(i)
        else if (.not. ieee_is_nan(kp(i))) then
            c = ctpoten_from_kp(kp(i))
        else
            cycle                                       ! no driver: missing data, not an error
        end if
        if (ieee_is_nan(c) .or. c <= 0._8) then
            ier = ior(ier, 1)
            cycle
        end if
        ce(i) = c
        ok(i) = .true.
    end do
    cpused = ce

    i = 1
    do while (i <= n)
        if (.not. ok(i)) then
            i = i + 1
            cycle
        end if
        north = mlat(i) >= 0._8
        m = 0
        j = i
        do while (j <= n)
            if (.not. ok(j)) exit
            if (by(j) /= by(i) .or. ce(j) /= ce(i) .or. ((mlat(j) >= 0._8) .neqv. north)) exit
            if (abs(mlat(j)) > 30._8) then              ! TIE-GCM (sub potm): potential is zero at |lat| <= 30 deg
                m = m + 1
                idx(m) = j
            else
                epot(j) = 0._8
            end if
            j = j + 1
        end do
        if (m > 0) then
            call heelis_cons(real(ce(i)), real(by(i)), variant, good)
            if (.not. good) then
                ier = ior(ier, 2)
            else
                do k = 1, m
                    dlat(k) = max(-90._8, min(90._8, mlat(idx(k))))*dtr
                    dlon(k) = (mlt(idx(k)) - 12._8)*h2deg*dtr       ! longitude from the sun: 0 at noon, + eastward
                    ratio(k) = 1._8
                    iflag(k) = 1                                    ! must be reset for every call, as sub potm does
                end do
                call flwv32(dlat(1:m), dlon(1:m), ratio(1:m), iflag(1:m), m, pot(1:m + 1), 1)
                do k = 1, m
                    epot(idx(k)) = pot(k)*1.e-3_8                   ! V -> kV
                end do
            end if
        end if
        i = j
    end do

    deallocate (ce, dlat, dlon, ratio, pot, idx, iflag, ok)
end subroutine heelis82_batch
