     PROGRAM MAIN
     PARAMETER (MAX=366*24*4)
     IMPLICIT REAL*8 (A-H,O-Z)
     CHARACTER*1 c1
     CHARACTER*5 c5
     CHARACTER*20 FLname
     INTEGER FLAG,jF,jL,yymmdd(1:MAX),YMD(1:MAX),iYD,iP
     REAL*8 UT(1:MAX),AE(1:MAX),time,long,param(2)
     REAL*8 F107(1:MAX),Vd_mean,Vd_storm,PromptVd,DynamoVd
!
     do i=1,max
        YYMMDD(i)=0.0D0
        YMD(i)=0.0D0
        AE(i)=0.0D0
        UT(i)=0.0D0
        F107(i)=0.0D0
     end do
!
     WRITE(*,*) 'Input AE data file name: ae_15min.dat'
     FLname='ae_15min.dat'
!
     OPEN(20,file=FLname)
     READ(20,*,END=30)
     DO I=1,MAX
          READ(20,*,END=30) yymmdd(I),UT(I),AE(I)
          jL=I
     ENDDO
30   CLOSE(20)
!
     IF (ABS(UT(jL)-UT(jL-1)).GT.0.5D0.AND.ABS(UT(2)-UT(1)).GT.0.5D0) THEN
          FLAG=0   !1h resolution;
          jF=49
          PRINT *,'1H resolution'
     ELSE
          FLAG=-1  !15 min resolution;
          jF=193
          PRINT *,'15M resolution'
     ENDIF
! Read F107;
     WRITE(*,*) 'Input F107 data file name: f107.dat'
     FLname='f107.dat'
     OPEN(30,file=FLname)
     DO I=1,MAX
          READ(30,35,END=40) iYEAR,c1,MONTH,c1,iDAY,c5,F107(I)
          YMD(I)=iYEAR*10000+MONTH*100+iDAY
     ENDDO
35   FORMAT(I4,A1,I2,A1,I2,A7,F5.1)
40   CLOSE(30)
!
     WRITE(*,*) 'INPUT GLON: -76.859012'
!     READ(*,*) long      ! Geographic Longitude
     long=-76.859012
!
     OPEN(10,file='drift.dat')
     write(10,*) '#yymmdd(iP)UT(iP)   Vd_quiet  PromptVd  DynamoVd  Vd_storm    VdTot'
     DO iP=jF,jL
          time=UT(iP)+long/15.0D0
          IF (time.GT.24.0D0)  time=time-24.0D0
! Year, Month and Day;
          iYEAR=yymmdd(iP)/10000
          MONTH=(yymmdd(iP)-iYEAR*10000)/100
          iDAY=yymmdd(iP)-(yymmdd(iP)/100)*100

          CALL MODA(0,iYEAR,MONTH,iDAY,iYD)

          param(1)=iYD   ! Day of year (e.g.,January 1st: param(1)=1.)
          param(2)=0.0   ! F10.7cm solar flux
          ij=1
45        IF (yymmdd(iP).EQ.YMD(ij)) THEN
             param(2)=F107(ij)  !F10.7cm solar flux;
             GOTO 50
          ELSE
             ij=ij+1
             IF (ij.LT.MAX) GOTO 45
             STOP 'F107 data error!'
          END IF
! Vertical drifts during quiet time;
50        CALL vdrift_quiet_model(time,long,param,Vd_mean)
! Vertical drifts during storm times;
          CALL vdrift_storm_model(FLAG,iP,AE,time,PromptVd,DynamoVd,Vd_storm)
          VdTot=Vd_mean+Vd_storm
          write(10,200) yymmdd(iP),UT(iP),Vd_mean,PromptVd,DynamoVd,Vd_storm,VdTot
          GOTO 60
     END DO
60   CLOSE(10)
100  format(5x,a3,6x,a11,2x,a11,6x,a3,8x,a7)
200  format(1X,I8,F8.3,5F10.3)
300  format(1X,I8,F8.3,5F10.3)
     END

