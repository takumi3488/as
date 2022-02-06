*#RUN * :L=TSPLIB L=LIB/ASL7 F=/LAGA2/D1,R(03) F=(09)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JB(2,3),KKM(3,7000)
      CHARACTER*2 CN
      DIMENSION KB(3)
c     CALL FPARAM(1,132)
      READ(1,*)
      READ(1,*) IL,NGEN,INV
      CALL TSPACE(IL)
C     CALL TSPHDS
C     CALL TSOPDS
      DO 1 I=1,NGEN
      READ(1,*) JA,((JB(J,K),J=1,2),K=1,3)
      CALL TSGENR(JA,JB)
    1 CONTINUE
      CALL TSPGRP(INV)
      CALL TSPGDS
      READ(1,*) A,B,C
      READ(1,'(3E23.20)') CA,CB,CC
      CALL TSLATC(A,B,C,CA,CB,CC0)
      READ(5,*) NX,NY,NZ
      CALL TSKPG8(NX,NY,NZ,KKM,ICC,NK)
      DO 10 K=1,NK
      KX=KKM(1,K)
      KY=KKM(2,K)
      KZ=KKM(3,K)
      CALL TSNMKP(KX,KY,KZ,ICC,CN)
      WRITE(6,600) KX,KY,KZ,ICC,CN
      WRITE(9,600) KX,KY,KZ,ICC,CN
  600 FORMAT(4I4,2X,A2)
      KB(1)=KX
      KB(2)=KY
      KB(3)=KZ
      CALL TSIREP(KB,ICC,0)
      CALL TSPKDS
   10 CONTINUE
      END
C
C SUBROUTINE TSKPG8 ====*====3====*====4====*====5====*====6====*====7
C
C   KPOINTS IN THE FIRST B.Z. ARE GENERATED IN KKM.
C   ICC IS COMMON DENOMINATOR AND NK IS NUMBER OF KPOINT.
C   AXIS OF GAMMA TO ZONE BOUNDARY CENTER IN THE X-DIRECTION
C        IS DIVIDED BY NX.
C   AXIS OF GAMMA TO ZONE BOUNDARY CENTER IN THE Y-DIRECTION
C        IS DIVIDED BY NY.
C   AXIS OF GAMMA TO ZONE BOUNDARY CENTER IN THE Z-DIRECTION
C        IS DIVIDED BY NZ.
C   ZONE BOUNDARY CENTER MEANS, HERE, MID-POINT BETWEEN GAMMA AND
C   THE NEAREST REXIPROCAL LATTICE POINT IN THE RESPECTIVE DIRECTION.
C
C                 1988.10.18 :  A. YANASE
C    MODIFIED FRON KPGEN 2022.01.09 : H.HARIMA
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSKPG8(NX,NY,NZ,KKM,ICC,NK)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
C
      DIMENSION KKM(3,*)
      IF(NX*NY*NZ.LE.0) THEN
       WRITE(6,900) NX,NY,NZ
  900 FORMAT(' ALL OF NX,NY,NZ SHOULD BE POSITVE',
     &       ' HOWEVER NX=',I3,' NY=',I3,' NZ=',I3,' THEN STOP HERE')
       STOP
      END IF
      NTLATC=-1
      CALL ADDINV(NTLATC)
      CALL ZZZY53(NX,NY,KK)
      LCM=NX*(NY/KK)
      CALL ZZZY53(LCM,NZ,KK)
      LCM=LCM*(NZ/KK)
C     WRITE(6,*) NX,NY,NZ,LCM
      IF(IL.EQ.-1) THEN
         IF(MOD(LCM,3).NE.0) THEN
           KDX=3*(LCM/NX)
           KDY=3*(LCM/NY)
           KDZ=3*(LCM/NZ)
           ICC=2*LCM
         ELSE
           KDX=LCM/NX
           KDY=LCM/NY
           KDZ=LCM/NZ
           ICC=2*(LCM/3)
         END IF
      ELSE IF(IL.EQ.0.OR.IL.EQ.1) THEN
           KDX=LCM/NX
           KDY=LCM/NY
           KDZ=LCM/NZ
           ICC=LCM*2
      ELSE IF(IL.EQ.2.OR.IL.EQ.3) THEN
           KDX=LCM/NX
           KDY=LCM/NY
           KDZ=LCM/NZ
           ICC=LCM
      ELSE IF(IL.EQ.4) THEN
           KDX=2*(LCM/NX)
           KDY=2*(LCM/NY)
           KDZ=LCM/NZ
           ICC=LCM*2
           IF(MOD(KDZ,2).EQ.0) THEN
               KDZ=KDZ/2
               KDY=KDY/2
               KDX=KDX/2
               ICC=ICC/2
           END IF
      END IF
      WRITE(6,901)
  901 FORMAT(/' GENERATED KPOINT'
     &       /'    NO  KX, KY, KZ')
      NK=0
      DO 1 IZ=0,NZ
         KKZ=IZ*KDZ
      DO 2 IY=0,NY
         IF(NTLATC.NE.1.OR.IY*IZ.EQ.0) THEN
               JEY=1
         ELSE
               JEY=2
         END IF
         KY=IY*KDY
      DO 3 IX=0,NX
         IF(NTLATC.GE.4.OR.NTLATC.LE.0.OR.IX.EQ.0) THEN
              JEX=1
         ELSE IF(IZ.EQ.0.AND.IY.EQ.0) THEN
              JEX=1
         ELSE
              JEX=2
         END IF
         KX=IX*KDX
       WRITE(6,*) KX,KY,KZ,ICC
      DO 12 JY=1,JEY
         KKY=KY*(3-2*JY)
      DO 13 JX=1,JEX
         KKX=KX*(3-2*JX)
         CALL TSKFBZ(KKX,KKY,KKZ,ICC,IND)
         WRITE(6,601) KKX,KKY,KKZ,ICC,IND
  601    FORMAT(' B.Z. TEST',4I6)
         IF(IND.EQ.0) GO TO 13
         CALL KALRST(KKX,KKY,KKZ,ICC,KKM,NK,IND)
         WRITE(6,602) KKX,KKY,KKZ,ICC,IND
  602    FORMAT(' TEST IF ALREADY REGISTERED',4I5)
         IF(IND.EQ.1) THEN
            NK=NK+1
            KKM(1,NK)=KKX
            KKM(2,NK)=KKY
            KKM(3,NK)=KKZ
            WRITE(6,910) NK,KKX,KKY,KKZ,ICC,IND
  910 FORMAT(1H ,I5,3I4,'/',I4,' IND=',I2)
         END IF
c        KKX=-KKX;KKY=-KKY;KKZ=-KKZ
c        CALL TSKFBZ(KKX,KKY,KKZ,ICC,IND)
C        WRITE(6,601) KKX,KKY,KKZ,ICC,IND
C 601    FORMAT(' B.Z. TEST',4I6)
c        IF(IND.EQ.0) GO TO 13
c        CALL KALRST(KKX,KKY,KKZ,ICC,KKM,NK,IND)
c        WRITE(6,602) KKX,KKY,KKZ,ICC,IND
c 602    FORMAT(' TEST IF ALREADY REGISTERED',4I5)
c        IF(IND.EQ.1) THEN
c           NK=NK+1
c           KKM(1,NK)=KKX
c           KKM(2,NK)=KKY
c           KKM(3,NK)=KKZ
C           WRITE(6,900) NK,KKX,KKY,KKZ,ICC,IND
C 900 FORMAT(1H ,I5,3I4,'/',I4,' IND=',I2)
c        END IF
   13 CONTINUE
   12 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
      CALL REMINV
      RETURN
      END
