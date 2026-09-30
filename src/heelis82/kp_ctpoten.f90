!  Cross-cap potential from Kp: function ctpoten_from_kp of TIE-GCM 2.0 (util.F), which sub getgpi (gpi.F) calls
!  when the user provides no ctpoten.  The formula and its comments are TIE-GCM's.
!
!  This software is part of the NCAR TIE-GCM.  Use is governed by the
!  Open Source Academic Research License Agreement contained in the file
!  tiegcmlicense.txt.
!
!  Changed by L. Navarro on 2026-09-25 for mpyricalspace: a Kp outside 0..9 makes TIE-GCM print a message and stop
!  the run (call shutdown); here the function returns NaN instead, so that one bad value does not end python.
!  Needs -fdefault-real-8, as TIE-GCM does.

module kp_module
    implicit none
contains

    real function ctpoten_from_kp(kp)
        !
        ! Calculate cross-tail potential from Kp.
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real, intent(in) :: kp
        !
        ctpoten_from_kp = 0.
        if (kp < 0. .or. kp > 9.) then
            ! TIE-GCM:  write(6,"('>>> ctpoten_from_kp: Bad kp=',e12.4)") kp
            !           call shutdown('ctpoten_from_kp')
            ctpoten_from_kp = ieee_value(kp, ieee_quiet_nan)
            return
        endif
        !
        ! Formula for potential:
        ! modified by LQIAN, 2007
        ! formula given by Wenbin based on data fitting
        !
        ctpoten_from_kp = 15. + 15.*kp + 0.8*kp**2
    end function ctpoten_from_kp

end module kp_module
