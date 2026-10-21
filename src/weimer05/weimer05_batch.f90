!  Vectorised Weimer 2005 driver: electric potential (epotval), magnetic potential (bpotval, below) and
!  field-aligned current (mpfac) at n points,
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
!  Outputs: epot [kV], fac [uA/m2], bpot = the magnetic-potential series sum(P * bsphc), in the model's own
!  units (the same series mpfac differentiates; the python wrapper scales it); NaN outside the auroral cap, for the southern hemisphere (the model is
!  defined on the northern cap only) and wherever any input is NaN.
!  ier: 0 ok | 1 a data file was not found in `path` (nothing computed) | 2 some time step had a cap boundary
!       outside the range the basis-function tables support (those rows are NaN).
!  Must be compiled with -fdefault-real-8 (as the model itself requires); the real(kind=8) declarations below
!  make the compiler reject a mismatch instead of silently mixing precisions.

subroutine weimer05_batch(by, bz, tilt, swvel, swden, mlat, mlt, path, xnan, epot, fac, bpot, ier, n)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use read_data, only: verbose, th0s, csize, ms, ab
    use w05sc, only: setmodel, epotval, mpfac, bndyfitr, mxtablesize, checkinputs, scplm, bsphc
    implicit none
    integer :: n, ier
    real(kind=8) :: by(n), bz(n), tilt(n), swvel(n), swden(n), mlat(n), mlt(n)
    real(kind=8) :: xnan, epot(n), fac(n), bpot(n)
    character(len=*) :: path

    real(kind=8), parameter :: fill = 1.e36_8      ! what epotval / mpfac return outside the cap
    character(len=:), allocatable :: dir
    character(len=14), parameter :: files(4) = [character(len=14) :: 'W05scEpot.dat', 'W05scBpot.dat', &
                                                                    'SCHAtable.dat', 'W05scBndy.dat']
    character(len=4), parameter :: models(2) = [character(len=4) :: 'epot', 'bpot']
    real(kind=8) :: cur(5), now(5)
    logical :: have, ok, exists
    integer :: i, k, m

    verbose = .false.
    ier = 0
    epot = xnan
    fac  = xnan
    bpot = xnan

    ! the model reads file_path//name, so the directory must end in a separator.
    ! Explicit allocate: with -std=legacy / -fno-automatic (meson.build) an assignment does not
    ! reallocate dir, so `dir = dir // '/'` dropped the '/' and dir kept the length of a previous call
    m = len_trim(path)
    if (allocated(dir)) deallocate(dir)
    if (m == 0) then
        allocate(character(len=2) :: dir)
        dir = './'
    else if (path(m:m) == '/') then
        allocate(character(len=m) :: dir)
        dir = path(1:m)
    else
        allocate(character(len=m+1) :: dir)
        dir = path(1:m) // '/'
    end if
    do k = 1, size(files)
        inquire(file=dir // trim(files(k)), exist=exists)
        if (.not. exists) then      ! a failed open() inside the model would end the whole python process
            ier = 1
            return
        end if
    end do

    ! pass 1 = electric potential ('epot' coefficients), pass 2 = magnetic potential and, from it, the
    ! field-aligned current ('bpot' coefficients)
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
                call bpotval(mlat(i), mlt(i), bpot(i))
            end if
        end do
    end do

contains

    ! The magnetic potential: epotval's series with the magnetic coefficients (bsphc, set by setmodel('bpot'))
    ! in place of the electric ones -- the model has no routine for it; mpfac sums the same series weighted by
    ! n(n+1), truncated at degree 10. All degrees are kept here, as in the potential itself. NaN off the cap.
    subroutine bpotval(lat, mlt, z)
        real(kind=8), intent(in)  :: lat, mlt
        real(kind=8), intent(out) :: z
        real(kind=8) :: phir, colat, nlm, plm, cospm(2), sinpm(2)
        integer :: inside, j, m
        logical :: skip

        call checkinputs(lat, mlt, inside, phir, colat)
        if (inside == 0) then
            z = xnan
            return
        end if
        cospm = cos([phir, 2. * phir])
        sinpm = sin([phir, 2. * phir])
        z = 0.
        skip = .false.
        do j = 1, csize
            if (skip) then                  ! the sine coefficient of the previous (m > 0) term
                skip = .false.
                cycle
            end if
            if (ab(j) /= 1) cycle
            m = ms(j)
            plm = scplm(j, colat, nlm)
            if (m == 0) then
                z = z + plm * bsphc(j)
            else
                z = z + plm * (bsphc(j) * cospm(m) + bsphc(j + 1) * sinpm(m))
                skip = .true.
            end if
        end do
    end subroutine bpotval
end subroutine weimer05_batch
