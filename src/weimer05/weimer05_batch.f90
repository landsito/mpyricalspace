!  Vectorised Weimer 2005 driver: electric potential (epotval) + field-aligned current (mpfac) at n points,
!  looping in Fortran so python makes one call instead of one per point.
!  Added by L. Navarro.
!
!  Every row i is an independent point with its own drivers:
!     by, bz  IMF (GSM) [nT]      tilt  dipole tilt [deg]      swvel  solar wind speed [km/s]
!     swden   proton density [cm-3]                            mlat, mlt  magnetic latitude [deg] / local time [h]
!  The model keeps its set-up (coefficients, cap boundary, basis-function table) in module state, so consecutive
!  rows with identical drivers -- a whole grid at one time step -- share a single set-up.  Put all the points of a
!  time step next to each other; a time series is then just the rows of each step in order.
!
!  Outputs: epot [kV], fac [uA/m2]; NaN outside the auroral cap, for the southern hemisphere (the model is
!  defined on the northern cap only) and wherever any input is NaN.
!  ier: 0 ok | 1 a data file was not found in `path` (nothing computed) | 2 some time step had a cap boundary
!       outside the range the basis-function tables support (those rows are NaN).
!  Must be compiled with -fdefault-real-8 (as the model itself requires); the real(kind=8) declarations below
!  make the compiler reject a mismatch instead of silently mixing precisions.

subroutine weimer05_batch(by, bz, tilt, swvel, swden, mlat, mlt, path, xnan, epot, fac, ier, n)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use read_data, only: verbose, th0s
    use w05sc, only: setmodel, epotval, mpfac, bndyfitr, mxtablesize
    implicit none
    integer :: n, ier
    real(kind=8) :: by(n), bz(n), tilt(n), swvel(n), swden(n), mlat(n), mlt(n)
    real(kind=8) :: xnan, epot(n), fac(n)
    character(len=*) :: path

    real(kind=8), parameter :: fill = 1.e36_8      ! what epotval / mpfac return outside the cap
    character(len=:), allocatable :: dir
    character(len=14), parameter :: files(4) = [character(len=14) :: 'W05scEpot.dat', 'W05scBpot.dat', &
                                                                    'SCHAtable.dat', 'W05scBndy.dat']
    character(len=4), parameter :: models(2) = [character(len=4) :: 'epot', 'bpot']
    real(kind=8) :: cur(5), now(5)
    logical :: have, ok, exists
    integer :: i, k

    verbose = .false.
    ier = 0
    epot = xnan
    fac  = xnan

    ! the model reads file_path//name, so the directory must end in a separator
    dir = trim(path)
    if (len(dir) == 0) then
        dir = './'
    else if (dir(len(dir):len(dir)) /= '/') then
        dir = dir // '/'
    end if
    do k = 1, size(files)
        inquire(file=dir // trim(files(k)), exist=exists)
        if (.not. exists) then      ! a failed open() inside the model would end the whole python process
            ier = 1
            return
        end if
    end do

    ! pass 1 = electric potential ('epot' coefficients), pass 2 = magnetic potential -> field-aligned current
    do k = 1, 2
        have = .false.
        ok   = .false.
        do i = 1, n
            now = (/ by(i), bz(i), tilt(i), swvel(i), swden(i) /)
            if (any(ieee_is_nan(now)) .or. ieee_is_nan(mlat(i)) .or. ieee_is_nan(mlt(i))) cycle

            if (.not. have .or. any(now /= cur)) then       ! new drivers: (re)do the set-up; otherwise reuse it
                call setmodel(now(1), now(2), now(3), now(4), now(5), dir, models(k))
                cur = now
                have = .true.
                ! the basis-function table needs 3*nint(cap) <= mxtablesize points and a cap size the SCHA table
                ! interpolates (it returns 0 below its 2nd entry); the model itself would `stop` the process
                ok = (bndyfitr >= th0s(2)) .and. (3 * nint(bndyfitr) <= mxtablesize)
            end if
            if (.not. ok) then
                ier = 2
                cycle
            end if

            if (k == 1) then
                call epotval(mlat(i), mlt(i), fill, epot(i))
                if (epot(i) == fill) epot(i) = xnan
            else
                call mpfac(mlat(i), mlt(i), fill, fac(i))
                if (fac(i) == fill) fac(i) = xnan
            end if
        end do
    end do
end subroutine weimer05_batch
