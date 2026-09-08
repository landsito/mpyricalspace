!  Vectorised IRI driver: loops iri_sub over n points in Fortran.
!  jf(50) and the oarr(100) template are constant across the grid (they come from
!  a single get_options() call on the python side), so python passes them once.
!  f107 / f107a are per-point external F10.7: a negative entry means "not
!  supplied, let IRI read apf107.dat", a non-negative entry overwrites oarr(41) /
!  oarr(46) for that point (the matching jf switch must be turned off python-side).
!  Rz12 / IG12 are always left to IRI's own ig_rz.dat.  Added by L. Navarro.

subroutine iri_batch(jf, oarr_in, iyyyy, mmdd, dhour, alati, along, hei, &
                     f107, f107a, outv, oarrv, n)
!f2py logical, intent(in)  :: jf
!f2py real,    intent(in)  :: oarr_in
!f2py integer, intent(in)  :: iyyyy
!f2py integer, intent(in)  :: mmdd
!f2py real,    intent(in)  :: dhour
!f2py real,    intent(in)  :: alati
!f2py real,    intent(in)  :: along
!f2py real,    intent(in)  :: hei
!f2py real,    intent(in)  :: f107
!f2py real,    intent(in)  :: f107a
!f2py real,    intent(out) :: outv
!f2py real,    intent(out) :: oarrv
!f2py integer, intent(hide), depend(iyyyy) :: n = len(iyyyy)
    implicit none
    integer :: n, i, k
    logical :: jf(50)
    real    :: oarr_in(100)
    integer :: iyyyy(n), mmdd(n)
    real    :: dhour(n), alati(n), along(n), hei(n)
    real    :: f107(n), f107a(n)
    real    :: outv(n, 20), oarrv(n, 100)
    real    :: outf(20, 1000), oarr(100)
    external iri_sub

    do i = 1, n
        do k = 1, 100
            oarr(k) = oarr_in(k)
        end do
        if (f107(i)  .ge. 0.0) oarr(41) = f107(i)
        if (f107a(i) .ge. 0.0) oarr(46) = f107a(i)
        call iri_sub(jf, 0, alati(i), mod(along(i), 360.0), iyyyy(i), mmdd(i), &
                     dhour(i), hei(i), hei(i) + 1.0, 1.0, outf, oarr)
        do k = 1, 20
            outv(i, k) = outf(k, 1)
        end do
        do k = 1, 100
            oarrv(i, k) = oarr(k)
        end do
    end do
    return
end subroutine iri_batch
