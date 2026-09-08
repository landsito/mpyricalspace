!!!
!!!  ROCSAT-1 Vertical Drift Model
!!!
!!! Point of contact
!!!    B. G. Fejer (bela.fejer@usu.edu)
!!!
!!! References:
!!!   Fejer, B. G., J. W. Jensen, and S.-Y. Su (2008),
!!!   Quiet time equatorial F region vertical plasma drift model
!!!   derived from ROCSAT-1 observations, J. Geophys. Res., 113, A05304
!!!   doi:10.1029/2007JA012801.
!!!
!!! History:
!!!   Added driver subroutine getverticaldrift. L. Navarro (luis.navarro@colorado.edu)
!!!
!!!================================================================================
!!! Input arguments:
!!!        f107 - daily solar flux
!!!        idoy - day of the year [integer]
!!!        ttl - local time [hours]
!!!        gglon - geographic longitude [radians]
!!!
!!! Output argument:
!!!        viv - vertical plasma drift [m/s]
!!!
!!! How to use it:
!!!    1) initialize drift arrays (fill in missing values)
!!!        CALL vfjmodelrocstart()
!!!    2) initialize data array for solar activity and day of year
!!!        CALL vfjmodelrocinit( f107, idoy )
!!!    3) retrieve vertical drift parsing hour and geographic longitude
!!!        CALL vfjmodelroc(ttl,gglon,viv)
!!!
!!! How to compile it (using example attached):
!!!    gfortran -std=legacy -ffixed-line-length-none check.F
!!!             vfjmodroc1.F -undefined dynamic_lookup -o check.x
!!!================================================================================

module vfjmodroc

    parameter( MSOL= 11, MLON= 24, MLT= 58, MSN= 4 )
    parameter( MLON1= MLON+1, MLT1=MLT+1 )

    dimension sfl(MSOL), glon(MLON1), tl(MLT1), seas(MSN)
    dimension fjrocvz(MLT1,MLON1,MSN,MSOL)

    integer :: jsfl= 1, jseas= 1

    data sfl /   100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200.  /
    data seas /   59., 120., 243., 304. /
    data glon / -180.,-165.,-150.,-135.,-120.,-105., -90., -75., -60., &
                 -45., -30., -15.,   0.,  15.,  30.,  45.,  60.,  75., &
                  90., 105., 120., 135., 150., 165., 180.  /
    data tl /   0.00,  0.50,  1.00,  1.50,  2.00,  2.50,  3.00,  3.50,  4.00, &
                4.50,  5.00,  5.50,  6.00,  6.50,  7.00,  7.50,  8.00,  8.50, &
                9.00,  9.50, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, &
               13.50, 14.00, 14.50, 15.00, 15.50, 16.00, 16.50, 17.00, 17.25, &
               17.50, 17.75, 18.00, 18.25, 18.50, 18.75, 19.00, 19.25, 19.50, &
               19.75, 20.00, 20.25, 20.50, 20.75, 21.00, 21.25, 21.50, 21.75, &
               22.00, 22.50, 23.00, 23.50, 24.00 /

    INCLUDE "fjrocdata.h"

    save

    contains

    subroutine vfjmodelroc( ttl, gglon, viv )
!  Inputs: ttl - local time [hour]
!          gglon - geographic longitude [radians]
!  Output: viv - vertical plasma drift [m/s]
    parameter( RD =   57.29577951308232 )

    xgglon= gglon*RD
    if( xgglon.gt.180. )  xgglon= xgglon - 360.
    xtl= ttl
    if( xtl.lt.0.  )  xtl= xtl+24.
    if( xtl.gt.24. )  xtl= xtl-24.
    call fjlin22dex( viv, xtl, xgglon, tl, glon, fjrocvz(1,1,jseas,jsfl), MLT1, MLON1 )
    return
    end subroutine vfjmodelroc

    subroutine vfjmodelrocstart
!  initialize drift arrays (fill in missing values)
    do 40 is= 1,MSOL
        do 40 isn= 1,MSN
            do 40 il= 1,MLON1
                do 30 it= 1,MLT1
                    if( fjrocvz(it,il,isn,is).lt.-900. )  then
                        itm= it-1
                        if( itm.lt.1 )  itm= itm+MLT
                        fjrocvz(it,il,isn,is)= fjrocvz(itm,il,isn,is)
                    endif
30  continue
40  continue
    return
    end subroutine vfjmodelrocstart

    subroutine vfjmodelrocinit( f107, idoy )
!  initialize data array for solar activity and day of year
    doy= idoy
    if( 0.lt.doy .and. doy.le.seas(1) )  then
        jseas= 1
    else if( seas(1).lt.doy .and. doy.le.seas(2) )  then
        jseas= 2
    else if( seas(2).lt.doy .and. doy.le.seas(3) )  then
        jseas= 3
    else if( seas(3).lt.doy .and. doy.le.seas(4) )  then
        jseas= 4
    else
        jseas= 1
    endif

    if( f107.lt.sfl(2) )  then
        jsfl= 1
    else if( f107.ge.sfl(MSOL)  )  then
        jsfl= MSOL
    else
        do 10 i= 1,MSOL-1
            if( sfl(i).le.f107 .and. f107.lt.sfl(i+1) )  then
                jsfl= i
                return
            endif
10      continue
    endif
    return
    end subroutine vfjmodelrocinit

    subroutine fjlin22dex( z, x, y, xd, yd, d, nx, ny )
! linear interp. in two dimensions  z = d(x,y)
    dimension d(nx,ny), xd(nx), yd(ny)

    call fjlocate( kx, x, xd, nx )
    kx= max( 1, min( kx, nx-1 ) )
    hx= xd(kx+1)-xd(kx)
    ax= (xd(kx+1)-x)/hx
    bx= (x-xd(kx))/hx
    call fjlocate( ky, y, yd, ny )
    ky= max( 1, min( ky, ny-1 ) )
    hy= yd(ky+1)-yd(ky)
    ay= (yd(ky+1)-y)/hy
    by= (y-yd(ky))/hy

    z=  ax*ay*d(kx,ky)   + bx*ay*d(kx+1,ky) &
      + ax*by*d(kx,ky+1) + bx*by*d(kx+1,ky+1)
    return
    end subroutine fjlin22dex

    subroutine fjlocate( kl, x, xa, na )
!  locates nearest item in array xa to item x
!   (items xa(kl) and xa(kl+1) bracket x if x is in the range of xa)
    dimension xa(na)

    if( na.le.1 )  then
        kl= 0
        return
    endif

    klo=1
    if( xa(2).gt.xa(1) )  then
        if( x.lt.xa(1) )  then
            kl= 0
            return
        else if( x.ge.xa(na) )  then
            kl= na
            return
        else if( xa(klo).le.x .and. x.lt.xa(klo+1) )  then
            kl= klo
            return
        endif
        klo=1
        khi=na
1       if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(xa(k).gt.x)  then
                khi=k
            else
                klo=k
            endif
            goto 1
        endif
        kl= klo

    else
!      xa is in descending order
        if( x.gt.xa(1) )  then
            kl= 0
            return
        else if( x.lt.xa(na) )  then
            kl= na
            return
        else if( xa(klo).ge.x .and. x.gt.xa(klo+1) )  then
            kl= klo
            return
        endif
        klo=1
        khi=na
2       if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(xa(k).ge.x)  then
                klo=k
            else
                khi=k
            endif
            goto 2
        endif
        kl= klo
    endif
    return
    end subroutine fjlocate

end module vfjmodroc



subroutine getverticaldrift(GL,IYD,F107,SLT,VD)
!f2py real,intent(in) :: GL
!f2py integer,intent(in) :: IYD
!f2py real,intent(in) :: F107
!f2py real,intent(in) :: SLT
!f2py real,intent(out) :: VD
    use vfjmodroc, only: vfjmodelrocstart, vfjmodelrocinit, vfjmodelroc
    implicit none
    INTEGER      :: IYD
    REAL*4       :: GL,F107,SLT,VD,GL1

    GL1=GL*0.017453292519943295

    CALL vfjmodelrocstart()
    CALL vfjmodelrocinit( F107, IYD )
    CALL vfjmodelroc(SLT,GL1,VD)
    RETURN
end subroutine getverticaldrift



subroutine getverticaldrift_batch(GL,IYD,F107,SLT,VD,N)
!  vectorised getverticaldrift: loops in Fortran, and runs the (f107/doy-independent)
!  gap-fill vfjmodelrocstart() only once instead of once per point.
!f2py real*4,    intent(in)  :: GL
!f2py integer,   intent(in)  :: IYD
!f2py real*4,    intent(in)  :: F107
!f2py real*4,    intent(in)  :: SLT
!f2py real*4,    intent(out) :: VD
!f2py integer,   intent(hide), depend(GL) :: N = len(GL)
    use vfjmodroc, only: vfjmodelrocstart, vfjmodelrocinit, vfjmodelroc
    implicit none
    integer :: N, I
    integer :: IYD(N)
    real*4  :: GL(N), F107(N), SLT(N), VD(N)
    real*4  :: RD
    parameter( RD = 0.017453292519943295 )

    call vfjmodelrocstart()
    do I = 1, N
        call vfjmodelrocinit( F107(I), IYD(I) )
        call vfjmodelroc( SLT(I), GL(I)*RD, VD(I) )
    end do
    return
end subroutine getverticaldrift_batch
