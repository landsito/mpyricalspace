*************************************************************************
C      TEST DRIVER
      PROGRAM TESTER
      USE vfjmodroc
      
      PARAMETER (K=48)
      DIMENSION VD(K),HOUR(K)
      DIMENSION GL(1)
      
      
      GL(1)=-76.859012*0.017453292519943295
      F107=150.0
      IDOY=201
      
      DO 10, J=1,K
       HOUR(J)=J*0.5
   10 CONTINUE
     
      CALL vfjmodelrocstart()
      CALL vfjmodelrocinit( F107, IDOY )

      DO 11, I=1,K
         CALL vfjmodelroc(HOUR(I),GL(1),VD(I))
         WRITE(6,100) HOUR(I),VD(I)
   11 CONTINUE
  100 FORMAT (1X,2F10.2)
      STOP
      END PROGRAM TESTER
