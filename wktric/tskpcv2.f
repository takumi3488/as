*#RUN * :L=TSPLIB L=LIB/ASL7 F=/LAGA2/D1,R(03) F=(09)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JB(2,3),KKM(3,1000),E(1000)
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
      J=1
   10 READ(2,200,END=99) (KKM(I,J),I=1,3),IIC,E(J)
c     write(6,*) (KKM(I,J),I=1,3),IIC
  200 FORMAT(4X,3I4,2X,I3,F10.6)
      J=J+1
      GOTO 10
   99 J=J-1
      WRITE(6,*) 'HAS BEEN READ ',J,' k-points'
      NK=J
C
      READ(3,*)
      J=1
   20 READ(3,*,END=98)KX,KY,KZ,IC,KXO,KYO,KZO,ICO
c     write(6,*) " IC,IIC=",IC,IIC
      IF(IC.NE.IIC) STOP
      CALL EQUIVK(KX,KY,KZ,IC,KKM,NK,IND)
      print *," IND=",IND,KX,KY,KZ,IC
      if(IND.EQ.0) THEN
      KX=-KX
      KY=-KY
      KZ=-KZ
      CALL EQUIVK(KX,KY,KZ,IC,KKM,NK,IND)
      print *," IND=",IND,KX,KY,KZ,IC
      endif
      IF(IND.EQ.0) STOP
      WRITE(6,300)IND,KX,KY,KZ,IC,E(IND),(KKM(I,IND),I=1,3),IIC
  300 FORMAT(I4,3I3,2X,I3,F10.6,3I3,2X,I3)
      WRITE(7,200)KXO,KYO,KZO,ICO,E(IND)
      J=J+1
      GOTO 20
   98 WRITE(6,*) 'FINISHED FOR ',J-1,' NEW POINTS'
      STOP
      END
C SUBROUTINE EQUIVK ===*====3====*====4====*====5====*====6====*====7
C
C    IF STAR OF (KX,KY,KZ)/IC IS ALREADY IN KKM, IND=THE INDEX
C    (KX,KY,KZ)/IC IS NEW                      , IND=0
C               1988.10.18  AKIRA YANASE
C               1996.11.19  HISATOMO HARIMA FROM KALRST
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE EQUIVK(KX,KY,KZ,IC,KKM,NK,IND)
        IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CW
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &  ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG2/,/SPG3/,/STK/
      DIMENSION KKM(3,*),KKB(3),KA(3),KC(3)
      KKB(1)=KB(1)
      KKB(2)=KB(2)
      KKB(3)=KB(3)
      ICBB=ICB
      KB(1)=KX
      KB(2)=KY
      KB(3)=KZ
      ICB=IC
C---------------------------------------------
C   CALCULATION OF STAR
      CALL ZZZY38
C---------------------------------------------
      IF(NK.EQ.0) GO TO 11
      DO 1 IS=1,NS
      DO 2 J=1,NK
         JN=J
         INDC=0
      DO 3 K=1,3
         KA(K)=KS(K,IS)-KKM(K,J)
         IF(MOD(IABS(KA(K)),ICB).NE.0) GO TO 2
         KA(K)=KA(K)/ICB
         KC(K)=IABS(KA(K))
         IF(KA(K).NE.0) INDC=1
    3 CONTINUE
      IF(INDC.EQ.1) THEN
         IF(IL.EQ.-1) GO TO 4
         IF(IL.EQ.0) GO TO 5
         IF(IL.EQ.1) GO TO 5
         IF(IL.EQ.2) GO TO 6
         IF(IL.EQ.3) GO TO 7
         IF(IL.EQ.4) GO TO 8
    4    IIWW=-KA(1)+KA(2)+KA(3)
         IF(MOD(IABS(IIWW),3).NE.0) GO TO 2
         GO TO 5
    6    IF(MOD(KC(1),2).NE.MOD(KC(2),2)) GO TO 2
         IF(MOD(KC(2),2).NE.MOD(KC(3),2)) GO TO 2
         GO TO 5
    7    IF(MOD(KC(1)+KC(2)+KC(3),2).NE.0) GO TO 2
         GO TO 5
    8    IF(MOD(KC(1)+KC(2),2).NE.0) GO TO 2
         GO TO 5
       ELSE
         GO TO 5
       END IF
    2 CONTINUE
    1 CONTINUE
   11    IND=0
         GO TO 13
    5    IND=JN
   13 KB(1)=KKB(1)
      KB(2)=KKB(2)
      KB(3)=KKB(3)
      ICB=ICBB
      RETURN
      END
