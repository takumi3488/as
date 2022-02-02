*#RUN *:FLTCHK OPT=3 INLINE=2 IAP L=C62479/TSPLIB;
*#      F=/CAVO3/D1,R(01) F=/CAVO3/ELSO,R(02) ;
*#      F=(07) F=(08)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'prmdos.f'
c     PARAMETER (LMNKP1=800,LMNEIG=62)
c     PARAMETER (LMCOMP=5)
c     SAVE /EIGEN0/
      COMMON/KPDATA/ NX,NY,NZ,KINTM,WIDTH
      CHARACTER*4 NCCHR
      COMMON/LCM/COMP(LMNEIG,LMNKP1,LMCOMP,2),
     &          IPEN(5,2,LMCOMP),NKCOMP(LMCOMP),NCOMP,NCCHR(LMCOMP)
C
      WRITE(6,*) ' Input NX,NY,NZ(for full BZ division)'
      READ(5,*)NX,NY,NZ
      WRITE(6,*) NX,NY,NZ
      READ(5,*)
      READ(5,*) NDI,NDE
      READ(5,*)
      READ(5,*) NSPIN
      READ(5,*)
      READ(5,*) NCOMP
      IF(NCOMP.LT.0.OR.NCOMP.GT.LMCOMP) STOP ' NCOMP ERROR'
      IF(NCOMP.NE.0)THEN
       DO 10 ICOMP=1,NCOMP
       READ(5,100) NKCOMP(ICOMP),NCCHR(ICOMP),
     &    (IPEN(NPEN,1,ICOMP),IPEN(NPEN,2,ICOMP),NPEN=1,NKCOMP(ICOMP))
  100  FORMAT(I2,1X,A4,10I2)
       WRITE(6,*) NKCOMP(ICOMP),(IPEN(NPEN,1,ICOMP),IPEN(NPEN,2,ICOMP),
     &           NPEN=1,NKCOMP(ICOMP)),NCCHR(ICOMP)
       IF(NKCOMP(ICOMP).GT.5) STOP ' NKCOMP ERROR'
   10  CONTINUE
      ENDIF
C
      CALL CRYRED(1)
      CALL SPACEG
      CALL KPGEN
      CALL ENRRED(2,NKP1)
      CALL ENERDS(7)
      CALL BANDDS
      CALL TDOSC0(NSPIN,NKP1,NDI,NDE,8)
      STOP
      END
      SUBROUTINE CRYRED(IFL1)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'prmdos.f'
C
c     PARAMETER (LMNATM=6,LMNKAT=3)
      COMMON/STR1  / NAT,NKAT,NAK(LMNKAT),KAT(2,LMNKAT),NAME(LMNKAT)
      COMMON/STR2  / ZNUCL(LMNKAT), SNUCL(LMNKAT)
      COMMON/STR3  / RATOM(3,LMNATM), SATOM(LMNKAT)
      COMMON/LATIN / A,B,C,CA,CB,CC
      COMMON/TSIN  / IL,NGEN,INV,IDU,JA(3),JB(2,3,3)
      COMMON/NEIGHB/ NR(LMNKAT,LMNKAT), NREG(LMNKAT,LMNKAT,100),
     &               REG(LMNKAT,LMNKAT,100)
      COMMON/VOL0  / ALATT,VCELL,VIN(LMNKAT),VOUT,WUNIT,CR(LMNKAT)
      COMMON/ELE0  / VAL,TOE,NSPIN
      COMMON/MAD0  / TMAD(LMNKAT,LMNKAT)
      COMMON/ATOM0 / NCS(LMNKAT), NNLZ2(35,LMNKAT), WWNL2(35,LMNKAT,2)
      DIMENSION      ANA(LMNKAT)
      INTEGER KION(LMNKAT)
      CHARACTER*40 CM
      CHARACTER*8  DA
C
      READ(IFL1,99) CM,DA
   99 FORMAT(1H ,A40,A8)
      WRITE(6,98) IFL1,CM,DA
   98 FORMAT(1H ,'WE USE THE DATA ON FILE(',I2,'):  ',A40,A8)
      READ(IFL1,102) IL,NGEN,INV,IDU
  102 FORMAT(14I5)
      WRITE(6,*) ' IL=',IL,'   NGEN=',NGEN,'   INV=',INV,'   IDU=',IDU
      DO 50 I=1,NGEN
      READ(IFL1,102) JA(I),((JB(I1,I2,I),I1=1,2),I2=1,3)
   50 WRITE(6,*) ' JA=',JA(I),' WITH (',JB(1,1,I),'/',JB(2,1,I),
     &                            ',  ',JB(1,2,I),'/',JB(2,2,I),
     &                            ',  ',JB(1,3,I),'/',JB(2,3,I),'  )'
      READ(IFL1,110) A,B,C
  110 FORMAT(3D23.16)
      WRITE(6,200) A,B,C
  200 FORMAT(1H ,'A,B,C=',3D24.16)
      READ(IFL1,110) CA,CB,CC
      WRITE(6,210) CA,CB,CC
  210 FORMAT(1H ,'CA,CB,CC=',3D24.16)
C     CALL TSLATC(A,B,C,CA,CB,CC)
C
      READ(IFL1,120) NAT,NKAT,(KION(I),I=1,NKAT)
  120 FORMAT(14I5)
      WRITE(6,220) NAT,NKAT,(KION(I),I=1,NKAT)
  220 FORMAT(1H ,'NAT,NKAT,(KION(I),I=1,NKAT)',15I5)
      IF(NKAT.GT.LMNKAT) THEN
        WRITE(6,*) ' ========== STOP DUE TO NKAT>LNNKAT IN SUB.CRYRED.'
        STOP
      END IF
      IF(NAT.GT.LMNATM) THEN
        WRITE(6,*) ' ========== STOP DUE TO NAT>LMNATM IN SUB.CRYRED.'
        STOP
      END IF
C
      DO 10 I=1,NAT
      READ(IFL1,110) (RATOM(J,I),J=1,3)
      WRITE(6,230) I,(RATOM(J,I),J=1,3)
  230 FORMAT(1H ,'I, (RATOM(J,I),J=1,3)',I5,3F12.7)
   10 CONTINUE
C
      J=0
      DO 20 I=1,NKAT
      JJ=J+1
      KAT(1,I)=JJ
      KAT(2,I)=KION(I)
      J=KION(I)
      NAK(I)=KAT(2,I)-KAT(1,I)+1
      ANA(I)=NAK(I)
   20 CONTINUE
C---*
C---*
      READ(IFL1,105) ((NR(I,J),J=1,NKAT),I=1,NKAT)
  105 FORMAT(14I5)
C     WRITE(6,204) ((NR(I,J),J=1,NKAT),I=1,NKAT)
C 204 FORMAT(1H ,'((NR(I,J),J=1,NKAT),I=1,NKAT)',20I5)
      DO 2 I=1,NKAT
      DO 2 J=1,NKAT
        IF(NR(I,J).LE.0) THEN
          READ(IFL1,105)
          READ(IFL1,101)
  101     FORMAT(3D23.16)
C         WRITE(6,205) I,J
          WRITE(6,*) ' ***** WARNING ***** NR(I,J) <= 0'
        ELSE IF(NR(I,J).LE.100) THEN
          READ(IFL1,105) (NREG(I,J,K),K=1,NR(I,J))
          READ(IFL1,101) (REG(I,J,K),K=1,NR(I,J))
C         WRITE(6,205) I,J, (NREG(I,J,K),K=1,NR(I,J))
C 205     FORMAT(1H ,2I5,7I15 / (11X,7I15) )
C         WRITE(6,206) (REG(I,J,K),K=1,NR(I,J))
C 206     FORMAT(1H ,10X,7F15.7)
        ELSE
          WRITE(6,*) ' ===== STOP ST SUB.CRYRED.  NR(I,J)=',NR(I,J)
          STOP ' ==== STOP IN SUB.CRYRED. (NR(I,J)>100) ==='
        END IF
    2 CONTINUE
      DO 3 I=1,NKAT
      READ(IFL1,101) (TMAD(I,J),J=1,NKAT)
      WRITE(6,207) I, (TMAD(I,J),J=1,NKAT)
  207 FORMAT(1H ,'TMAD(I,J),I=',I2,':  ',7F15.7 / (18X,7F15.7))
    3 CONTINUE
C---*                                    LATTICE CONTSTANTS
      READ(IFL1,101) VP,AS
      WRITE(6,208) VP,AS
  208 FORMAT(1H ,'VP,AS=',2F15.7)
C     ALATT=A/0.52917706D0
      ALATT=A
      VCELL=VP*(ALATT**3)
C                    ..... VOLUME OF PRIMITIVE CELL
      PAI  =4.D0*DATAN(1.D0)
      WUNIT=(2.0*PAI/ALATT)*AS
C                    ..... UNIT OF K-SPACE (LENGTH OF K=(1,0,0))
C---*                              ATOMIC SPECIES
      DO 12 I=1,NKAT
      READ(IFL1,103) NAME(I),ZNUCL(I),NCS(I)
  103 FORMAT(A2,D23.16,I5)
      READ(IFL1,102) (NNLZ2(J,I),J=1,NCS(I))
      READ(IFL1,101) (WWNL2(J,I,1),J=1,NCS(I))
   12 READ(IFL1,101) (WWNL2(J,I,2),J=1,NCS(I))
      VAL=0.0
      TOE=0.0
      DO 13 I=1,NKAT
      DO 14 J=1,NCS(I)
      TOE=WWNL2(J,I,1)*ANA(I)+TOE
      TOE=WWNL2(J,I,2)*ANA(I)+TOE
      IF(MOD(NNLZ2(J,I),10).LT.5) GO TO 14
      VAL=WWNL2(J,I,1)*ANA(I)+VAL
      VAL=WWNL2(J,I,2)*ANA(I)+VAL
   14 CONTINUE
   13 CONTINUE
C     NTOE=TOE+0.5
C     NVAL=VAL+0.5
      WRITE(6,*) ' ====== ELECTRON NUMBER ===='
      WRITE(6,*) ' TOTAL=',TOE,'    VALENCE=',VAL
      WRITE(6,*) ' ==========================='
C
      DO 60 I=1,NKAT
        CALL NUCLEI(ZNUCL(I),AMASS,SNUCL(I))
        WRITE(6,*) ' IKAT=',I,'  Z,MASS,RADIUS=',ZNUCL(I),AMASS,SNUCL(I)
   60 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE SPACEG
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmdos.f'
c     PARAMETER (LMNATM=6,LMNKAT=3,LMMESH=380)
      COMMON/STR1  / NAT,NKAT,NAK(LMNKAT),KAT(2,LMNKAT),NAME(LMNKAT)
      COMMON/STR3  / RATOM(3,LMNATM), SATOM(LMNKAT)
      COMMON/TSIN  / ILT,NGEN,INV,IDU,JA(3),JB(2,3,3)
      COMMON/LATIN / AI,BI,CI,CAI,CBI,CCI
C     COMMON/SPG2  / IL,NG,IG(48),JV(2,3,48)
      INTEGER KS(11,20)
C
      CALL TSPACE(ILT)
C     CALL TSOPDS
      DO 10 I=1,NGEN
   10 CALL TSGENR(JA(I),JB(1,1,I))
      CALL TSPGRP(INV)
      CALL TSLATC(AI,BI,CI,CAI,CBI,CCI)
      CALL LATT01(6)
      CALL TSPGDS
      CALL TSCRST(RATOM,KAT,NKAT,NAT,KS)
C
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TDOSC0(NSPIN,NKP1,NDI,NDE,JD)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'prmdos.f'
c     PARAMETER (LMNEIG=62,LMNKP1=800,LMNKP0=20000,NEMAX=10001)
c     PARAMETER (LMCOMP=5,LNCOMP=LMCOMP+1)
      SAVE /KPDATA/,/EIGEN0/,/LCM/
      COMMON/KPDATA/ NX,NY,NZ,KINTM,WIDTH
      COMMON/EIGEN0/ EIG(LMNEIG,LMNKP1,2),NEIG,ISO
      CHARACTER*4 NCCHR
      COMMON/LCM/COMP(LMNEIG,LMNKP1,LMCOMP,2),
     &          IPEN(5,2,LMCOMP),NKCOMP(LMCOMP),NCOMP,NCCHR(LMCOMP)
      DIMENSION EA(LMNKP0),WT(LMCOMP,LMNKP0)
      DIMENSION D(LNCOMP,NEMAX),S(LNCOMP,NEMAX)
      DIMENSION DD(LNCOMP,NEMAX,2),SS(LNCOMP,NEMAX,2)
      IXMAX=NX
      IYMAX=NY
      IZMAX=NZ
      JS=1
c     NCOMP=5
      EMIN=EIG(1,1,1)
      EMAX=EIG(1,1,1)
      WRITE(6,*) ' 2050 NSPIN,NKP1,NEIG,NCOMP=',NSPIN,NKP1,NEIG,NCOMP
      DO 10 IS=1,NSPIN
      DO 10 IKP1=1,NKP1
      DO 10 I=1,NEIG
        IF(EMIN.GT.EIG(I,IKP1,IS)) EMIN=EIG(I,IKP1,IS)
        IF(EMAX.LT.EIG(I,IKP1,IS)) EMAX=EIG(I,IKP1,IS)
c     WRITE(6,*) ' I,IKP,EMIN,EMAX=',I,IKP1,EMIN,EMAX
   10 CONTINUE
c     WRITE(6,*) ' 2122 EMIN,EMAX=',EMIN,EMAX
      IF(EMIN.GE.0.D0) THEN
        EMIN=DBLE(INT(EMIN*10))/10
      ELSE
        EMIN=DBLE(INT(EMIN*10))/10-0.1D0
      END IF
      IF(EMAX.GE.0.D0) THEN
        EMAX=DBLE(INT(EMAX*10))/10+0.1D0
      ELSE
        EMAX=DBLE(INT(EMAX*10))/10
      END IF
      IE=INT((EMAX-EMIN)*100)/100 + 1
      DE=IE*0.0001D0
      NE=INT((EMAX-EMIN)/DE+1.D-3)+1
      IF(NE.GT.NEMAX) STOP ' === STOP IN SUB.TDOSC0. (NE>101) ==='
      IF(NSPIN.EQ.1.AND.ISO.EQ.0) THEN
        WEI=2.D0
      ELSE
        WEI=1.D0
      END IF
      WRITE(6,*)' 2155 NE,NEIG,EMIN,EMAX=',NE,NEIG,EMIN,EMAX
      DO 20 IS=1,NSPIN
      DO 20 I=1,NE
      DO 20,ICOMP=1,NCOMP+1
        DD(ICOMP,I,IS)=0.D0
        SS(ICOMP,I,IS)=0.D0
   20 CONTINUE
C
      NDOSI=1
      NDOSE=NEIG
      IF(NDI.GT.NDOSI)NDOSI=NDI
      IF(NDE.LT.NDOSE)NDOSE=NDE
      WRITE(6,*) ' *** START CAL. FROM ',NDOSI,' TO',NDOSE,' ***'
C
      DO 30 IS=1,NSPIN
c     DO 30 IB=1,NEIG
      DO 30 IB=NDOSI,NDOSE
      WRITE(6,*) ' ***** DOS FOR ',IB,'*****'
C       CALL TDOSC1(NE,EMIN,DE,IB,IS,1,D,S)
        CALL ENSORT(EA,WT,IB,IXMAX,IYMAX,IZMAX,IS,JS)
C     WRITE(6,*) NE,EMIN,DE,NCOMP,IXMAX,IYMAX,IZMAX,JS,(WT(II,1),II=1,5)
        CALL PDOS(EA,WT,NE,EMIN,DE,NCOMP,IXMAX,IYMAX,IZMAX,JS,D,S)
        DO 40 I=1,NE
        DO 150,ICOMP=1,NCOMP+1
          DD(ICOMP,I,IS)=DD(ICOMP,I,IS)+D(ICOMP,I)*WEI
          SS(ICOMP,I,IS)=SS(ICOMP,I,IS)+S(ICOMP,I)*WEI
  150    CONTINUE
   40   CONTINUE
   30 CONTINUE
      WRITE(JD,100)(NCCHR(ICOMP),ICOMP=1,NCOMP)
      IFLG1=0
      IFLG0=1
      DO 51,IS=1,NSPIN
      DO 50 I=1,NE
      IFLG00=IFLG0
      IFLG0=IFLG1
      IFLG1=1
      IF(DD(1,I,1).LT.1.D-4.AND.DD(1,I,2).LT.1.D-4) IFLG1=0
      IF(IFLG0+IFLG1.EQ.0) GOTO 50
C
        DO 1000 ICOMP=1,NCOMP+1
         IF(DD(ICOMP,I,IS).GT.9999.999)THEN
          DD(ICOMP,I,IS)=9999.999
         ENDIF
 1000   CONTINUE
C
      IF(IFLG00+IFLG0.EQ.0) THEN
        WRITE(JD,110) EMIN+DE*(I-2),
     & (DD(ICOMP,I-1,IS),SS(ICOMP,I-1,IS),ICOMP=1,NCOMP+1)
      ENDIF
        WRITE(JD,110) EMIN+DE*(I-1),
     & (DD(ICOMP,I,IS),SS(ICOMP,I,IS),ICOMP=1,NCOMP+1)
  100 FORMAT(1H ,1X,'ENERGY(Ry)',2X,'D.O.S.',6X,'N.O.S.',2X,
     &10(A4,' DOS   NOS '))
  110 FORMAT(1H ,'E=',F7.4,' D=',F9.3,' S=',F7.3,1X,10(F8.3,F7.3))
   50 CONTINUE
   51 CONTINUE
      RETURN
      END
C SUBROUTINE NUCLEI ====*====3====*====4====*====5====*====6====*====7
C
C      TO GET MASS NUMBER AND RADIUS OF THE NUCLEUS WITH CHARGE Z.
C
C      1983.8.4. :  N. HAMADA
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE NUCLEI(Z,A,R)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C          REF) E. FERMI, 'NUCLEAR PHYSICS',
C               THE UNIVERSITY OF CHICAGO PRESS (1950,CHICAGO).
C
C          Z = A / (P + Q * A**(2/3) )
C
C          Z :  NUCLEAR CHARGE
C          A :  MASS NUMBER
C
      DATA P /1.98D0/, Q /0.015D0/, C2 /0.666666666666666D0/
C
C          R = S * A**(1/3)   :     RADIUS OF THE NUCLEUS.
C
      DATA S /0.208D-4/, C1 /0.333333333333333D0/
C
      A1=2.0D0*Z
   10 CONTINUE
      A2=Z*(P+Q*(A1**C2))
      A0=DABS(A2-A1)
      A1=A2
      IF(A0.GT.0.1D-5) GO TO 10
      A=A1
C
      R = S * (A**C1)
C
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE KPGEN
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmdos.f'
c     PARAMETER (LMNKP1=800,LMNKP0=20000)
      COMMON/KPDATA/ NX,NY,NZ,KINTM,WIDTH
      COMMON/SPG2  / IL,NG,IG(48),JV(2,3,48)
      COMMON/KP0   / KP0(4,LMNKP0),IKP0(LMNKP0),NKP0
      COMMON/KP1   / KP1(4,LMNKP1),NKST1(LMNKP1),NKP1,NKP2
C
      CALL KPGEN1(IL,NG,IG,NX,NY,NZ,LMNKP0,NKP0,KP0,IKP0,
     &            LMNKP1,NKP1,KP1)
      CALL NKSTAR(IDU,NKP1,KP1,NKST1,NKP2)
c     WRITE(6,*) ' === K POINTS === (NX,NY,NZ)=(',NX,NY,NZ,')',
c    &    ' ===== NKP1=',NKP1,'   NKP2=',NKP2
      WRITE(6,100) NX,NY,NZ,NKP1,NKP2
  100 FORMAT(' === K POINTS === (NX,NY,NZ)=(',3I3,')',
     &       ' ===== NKP1=',I5,'   NKP2=',I5)
      WRITE(6,120) ((KP1(J,JKP1),J=1,4),NKST1(JKP1),JKP1=1,NKP1)
  120 FORMAT(1H ,4('  (',3I5,' ) /',I5,' (',I3,')'))
C     WRITE(6,*) '======= NKP0=',NKP0
C     WRITE(6,120) ((KP0(J,JKP0),J=1,4),IKP0(JKP0),JKP0=1,NKP0)
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE KPGEN1(IL,NG,IG,NX,NY,NZ,LMNKP0,NKP0,KP0,IKP0,
     &                  LMNKP1,NKP1,KP1)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IG(NG)
      INTEGER KP0(4,LMNKP0),IKP0(LMNKP0),KP1(4,LMNKP1)
c
      INTEGER KP2(3,LMNKP1)
c
C
      CALL KPGEN0(IL,NX,NY,NZ,LMNKP0,NKP0,KP0)
c
      CALL ZZZY53(NY*NZ,NZ*NX,NGC1)
      CALL ZZZY53(NX*NY,NGC1,NGCB)
      do 21 I0=1,NKP0
      do 21 J=1,4
      if(MOD(KP0(J,I0),NGCB).NE.0) then
       write(6,*) '****WARN in KPGEN1 ***'
      endif
      KP0(J,I0)=KP0(J,I0)/NGCB
   21 continue
      CALL ADDINV(NTLATC)
C
      DO 20 J=1,4
   20 KP1(J,1)=KP0(J,1)
      IKP0(1)=1
      NKP1=1
      DO 22 I0=2,NKP0
C       WRITE(6,*) ' IL=',IL,'   I0=',I0
c       DO 24 I1=1,NKP1
c       I11=I1
c       CALL EQUIVK(IL,NG,IG,KP1(1,I1),KP0(1,I0),IND)
c
        DO 24 I1=1,NKP1
        DO 24 J=1,3
        KP2(J,I1)=KP1(J,I1)
   24 CONTINUE
      CALL KPFIND(KP0(1,I0),KP0(2,I0),KP0(3,I0),KP0(4,I0),KP2,NKP1,IND)
        I11=IND
c     WRITE(6,*) 'IND=',IND
c     WRITE(6,*) ' KP1=',(KP1(J,I11),J=1,4),'   KP0=',(KP0(J,I0),J=1,4)
        IF(IND.NE.0) GO TO 28
c
c     WRITE(6,*) 'IND=',IND
c     WRITE(6,*) ' KP1=',(KP1(J,I11),J=1,4),'   KP0=',(KP0(J,I0),J=1,4)
c       IF(IND.EQ.0) GO TO 28
c  24   CONTINUE
          NKP1=NKP1+1
          IF(NKP1.GT.LMNKP1) THEN
            WRITE(6,*) ' NKP1,LMNKP1=',NKP1,LMNKP1
            STOP ' === STOP IN SUB.KPGEN1. (NKP1) ==='
          END IF
          DO 26 J=1,4
   26     KP1(J,NKP1)=KP0(J,I0)
          IKP0(I0)=NKP1
C         WRITE(6,*) '      I0=',I0,'   IKP0=',IKP0(I0)
        GO TO 22
   28   CONTINUE
          IKP0(I0)=I11
   22 CONTINUE
      CALL REMINV
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE KPGEN0(IL,NX,NY,NZ,LMNKP, NKP,KP)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INTEGER KP(4,LMNKP)
C
      IF(IL.LT.-1) THEN
        STOP ' === STOP IN SUB.KPGEN0. (IL.LT.-1) ==='
      ELSE IF(IL.EQ.-1) THEN
        NXX=NX*3
        NYY=NY*3
        NZZ=NZ*3
      ELSE IF(IL.LE.1) THEN
        NXX=NX
        NYY=NY
        NZZ=NZ
      ELSE IF(IL.LE.3) THEN
        NXX=NX*2
        NYY=NY*2
        NZZ=NZ*2
      ELSE IF(IL.EQ.4) THEN
        NXX=NX*2
        NYY=NY*2
        NZZ=NZ
      ELSE
        STOP ' === STOP IN SUB.KPGEN0. (IL.GT.4) ==='
      END IF
      ND=NX*NY*NZ
      NKP=0
      DO 10 IX=0,NXX
      DO 10 IY=0,NYY
      DO 10 IZ=0,NZZ
        NKP=NKP+1
        IF(NKP.GT.LMNKP) THEN
          WRITE(6,*) 'IX,IY,IZ=',IX,IY,IZ,'   NKP,LMNKP=',NKP,LMNKP
          STOP ' === STOP IN SUB.KPGEN0. (NKP>LMNKP) ==='
        END IF
        KP(1,NKP)=IX*NY*NZ
        KP(2,NKP)=NX*IY*NZ
        KP(3,NKP)=NX*NY*IZ
        KP(4,NKP)=ND
   10 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE EQUIVK(IL,NG,IG,K1,K2,IND)
C
      INTEGER IG(NG)
      INTEGER K1(4),K2(4),KK1(3),K(4)
C
      K(4)=K1(4)*K2(4)
      DO 10 M=1,NG
        CALL ROTATK(IL,IG(M),K1,KK1)
        DO 20 J=1,3
   20   K(J)=KK1(J)*K2(4)-K2(J)*K1(4)
C
        CALL LATPK(IL,K,IND)
C
        IF(IND.EQ.0) RETURN
   10 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE LATPK(IL,K,IND)
C
C     IL :  LATTICE TYPE.
C     (K(1),K(2),K(3))/K(4) IS LATTICE POINT OR NOT? (IN K SPACE).
C        IND=0 :  LATTECE POINT
C            1 :  NOT
C
      DIMENSION K(4),KK(3)
C
      I=K(4)
      DO 10 J=1,3
        IF(MOD(K(J),I).NE.0) THEN
          IND=1
          RETURN
        END IF
        KK(J)=K(J)/I
   10 CONTINUE
C
      IND=1
      IF(IL.EQ.-1) THEN
        IF(MOD(-KK(1)+KK(2)+KK(3),3).EQ.0) IND=0
      ELSE IF(IL.EQ.0) THEN
        IND=0
      ELSE IF(IL.EQ.1) THEN
        IND=0
      ELSE IF(IL.EQ.2) THEN
        KK1=IABS(KK(1))
        KK2=IABS(KK(2))
        KK3=IABS(KK(3))
        IF(MOD(KK1,2).EQ.0 .AND.
     &     MOD(KK2,2).EQ.0 .AND.
     &     MOD(KK3,2).EQ.0      ) THEN
           IND=0
        ELSE IF(MOD(KK1,2).EQ.1 .AND.
     &          MOD(KK2,2).EQ.1 .AND.
     &          MOD(KK3,2).EQ.1      ) THEN
           IND=0
        END IF
      ELSE IF(IL.EQ.3) THEN
        IF(MOD(KK(1)+KK(2)+KK(3),2).EQ.0) IND=0
      ELSE IF(IL.EQ.4) THEN
        IF(MOD(KK(1)+KK(2),2).EQ.0) IND=0
      ELSE
        STOP ' === STOP IN SUB.LATPK. (IL) ==='
      END IF
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ROTATK(IL,JA,KM,KK)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IT(3,48),IH(3,24),KK(3),KM(3)
      DATA IT/
     &  4,5,6, 4,1,0, 2,5,0, 2,1,6, 6,4,5, 0,4,1,
     &  0,2,5, 6,2,1, 5,6,4, 5,0,2, 1,6,2, 1,0,4,
     &  5,4,0, 1,2,0, 6,1,4, 2,6,5, 0,1,2, 2,0,1,
     &  4,0,5, 6,5,2, 1,4,6, 4,6,1, 0,5,4, 5,2,6,
     &  2,1,0, 2,5,6, 4,1,6, 4,5,0, 0,2,1, 6,2,5,
     &  6,4,1, 0,4,5, 1,0,2, 1,6,4, 5,0,4, 5,6,2,
     &  1,2,6, 5,4,6, 0,5,2, 4,0,1, 6,5,4, 4,6,5,
     &  2,6,1, 0,1,4, 5,2,0, 2,0,5, 6,1,2, 1,4,0/
      DATA IH/
     &  4,5,6, 1,3,6, 7,4,6, 2,1,6, 5,7,6, 3,2,6,
     &  2,3,0, 3,1,0, 1,2,0, 4,7,0, 7,5,0, 5,4,0,
     &  2,1,0, 5,7,0, 3,2,0, 4,5,0, 1,3,0, 7,4,0,
     &  4,7,6, 7,5,6, 5,4,6, 2,3,6, 3,1,6, 1,2,6/
      DO 1 K=1,3
        MA=IT(K,JA)-3
        IF(IL.LE.0) MA=IH(K,JA)-3
        IF(MA.EQ.0) GO TO 2
        IF(MA.EQ.4) GO TO 3
        MB=IABS(MA)
        KK(K)=KM(MB)
        IF(MA.LT.0) KK(K)=-KK(K)
        GO TO 1
    2   KK(K)=KM(1)+KM(2)
        GO TO 1
    3   KK(K)=-(KM(1)+KM(2))
    1 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE NKSTAR(IDU,NK,K,NKST,NK2)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER K(4,1),NKST(1)
      COMPLEX*16 CW
      COMMON/STK   / KS(3,48),JS(48),NS,ICBB,CW(48,12)
C
      NK2=0
      DO 10 I=1,NK
        CALL TSIREP(K(1,I),K(4,I),IDU)
C       CALL TSPKDS
        CALL ZZZY38
        NKST(I)=NS
        NK2=NK2+NS
   10 CONTINUE
      RETURN
      END
      SUBROUTINE ENRRED(JF,KKMAX)
C**********************************************************************
C************ LIBRARY HEADER ******************************************
C     ***READ ROUTINE FROM APW-OUTPUT FORMAT FILE                   ***
C     *** JF:LAPW-OUTPUT FILE                                       ***
C     *** NKAT:NUMBER OF KINDS OF ATOMS                             ***
C     ***                                    BY H.HARIMA    SEP.1986***
C     ***                         MODIFIED   BY H.HARIMA    APR.1990***
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'prmdos.f'
c     PARAMETER (LMNKP1=800,LMNEIG=62)
c     PARAMETER (LMNKAT=3,LCMAX=4)
c     PARAMETER (LMCOMP=5)
      COMMON/STR1/NAT,NKAT,NAK(LMNKAT),KAT(2,LMNKAT),NAME(LMNKAT)
      DIMENSION MRR(LMNKP1),MRN(LMNKP1)
      DIMENSION NEIG(LMNKP1)
      SAVE /KP1/,/EIGEN0/,/LCM/
      COMMON/KP1   / KP1(4,LMNKP1),NKST1(LMNKP1),NKP1,NKP2
      COMMON/EIGEN0/EIG(LMNEIG,LMNKP1,2),JMAX,ISO
      CHARACTER*4 NCCHR
      COMMON/LCM/COMP(LMNEIG,LMNKP1,LMCOMP,2),
     &          IPEN(5,2,LMCOMP),NKCOMP(LMCOMP),NCOMP,NCCHR(LMCOMP)
      DIMENSION KV(3,LMNKP1),IC(LMNKP1),IUD(LMNKP1)
      DIMENSION KKV(3,LMNKP1),ICC(LMNKP1)
      DIMENSION CO(LMNEIG,LMNKAT,LCMAX),CC(LMNEIG,LMCOMP)
      DIMENSION ENR(LMNEIG)
      SAVE KMAX
      KMAX=KKMAX
      K=0
      MXIC=1
      DO 10 KK=1,NKP1
      KKV(1,KK)=KP1(1,KK)
      KKV(2,KK)=KP1(2,KK)
      KKV(3,KK)=KP1(3,KK)
      ICC(KK)=KP1(4,KK)
      NEIG(KK)=0
C     WRITE(6,*) (KKV(J,KK),J=1,3),ICC(KK)
   10 CONTINUE
  100 CONTINUE
      READ(JF,150,END=110) KX,KY,KZ,IIC,IU,MR,MN,ND,NGIV
c     IF(MOD(ND,2).NE.0)STOP ' ND=ODD'
      MWEI=ND
c     MWEI=ND/2
      WRITE(6,150) KX,KY,KZ,IIC,IU,MR,MN,ND,NGIV
c 150 FORMAT(3X,3X,3I3,1X,2I3,2X,A2,I2,I3,3X,I3)
  150 FORMAT(2X,3X,3I3,1X,2I3,2X,A2,I2,I3,3X,I3)
      IF(NGIV.EQ.0) GOTO 100
      READ(JF,*) (ENR(J),J=1,NGIV)
      READ(JF,*)(((CO(NL,NK,L),L=1,4),NK=1,NKAT),NL=1,NGIV)
      DO 11 NL=1,NGIV
      OUT=1.0
      DO 310 NK=1,NKAT
      DO 310 L=1,4
      OUT=OUT-CO(NL,NK,L)
  310 CONTINUE
      DO 330 I=1,NCOMP
      CC(NL,I)=0.D0
      DO 340,NKC=1,NKCOMP(I)
      IF(IPEN(NKC,1,I).NE.0)THEN
      CC(NL,I)=CC(NL,I)+CO(NL,IPEN(NKC,1,I),IPEN(NKC,2,I))
      ELSE
      CC(NL,I)=CC(NL,I)+OUT
      ENDIF
  340 CONTINUE
  330 CONTINUE
   11 CONTINUE
  200 FORMAT(8F8.4)
C
      IF(MOD(KX*ICC(1),IIC).NE.0)GOTO 100
      KKX=KX*ICC(1)/IIC
      IF(MOD(KY*ICC(1),IIC).NE.0)GOTO 100
      KKY=KY*ICC(1)/IIC
      IF(MOD(KZ*ICC(1),IIC).NE.0)GOTO 100
      KKZ=KZ*ICC(1)/IIC
      CALL KPFIND(KKX,KKY,KKZ,ICC(1),KKV,NKP1,KK)
      IF(KK.EQ.0)THEN
C     WRITE(6,*) KX,KY,KZ,IIC
C     WRITE(6,*) KKX,KKY,KKZ,ICC(1)
      GOTO 100
c     STOP '  AT 8260 UNKNOWN POINT'
      ENDIF
C
C     REGISTRATION OLD POINT
C
      JMAX=NEIG(KK)
      IF(JMAX.EQ.0)GOTO 20
      DO 30 I=1,NGIV
      J1=JMAX
      DO 45 JJ=1,J1
      J=J1-JJ+1
C     DLTAE=ABS(ENR(I)-EIG(J,KK,IU))
C     IF(DLTAE.LT.0.00001) GOTO 30
   45 CONTINUE
      DO 40 JJ=1,J1
      J=J1-JJ+1
      IF(ENR(I).GT.EIG(J,KK,IU))GOTO 50
      EIG(J+MWEI,KK,IU)=EIG(J,KK,IU)
      DO 350 L=1,NCOMP
      COMP(J+MWEI,KK,L,IU)=COMP(J,KK,L,IU)
  350 CONTINUE
   40 CONTINUE
      J=0
   50 DO 51 MW=1,MWEI
      EIG(J+MW,KK,IU)=ENR(I)
      DO 351,L=1,NCOMP
      COMP(J+MW,KK,L,IU)=CC(I,L)
  351 CONTINUE
   51 CONTINUE
      JMAX=JMAX+MWEI
   30 CONTINUE
      NEIG(KK)=JMAX
      WRITE(6,*)' NEIG=',NEIG(KK)
      GOTO 100
C
C     REGISTRATION NEW POINT
C
   20 CONTINUE
      K=K+1
      MRR(KK)=MR
      MRN(KK)=MN
C     WRITE(6,*)'K=',K,KX,KY,KZ,IIC,MWEI,NGIV
      KV(1,KK)=KX
      KV(2,KK)=KY
      KV(3,KK)=KZ
      IC(KK)=IIC
C     BKV(1,KK)=DBLE(KX)/DBLE(IIC)
C     BKV(2,KK)=DBLE(KY)/DBLE(IIC)
C     BKV(3,KK)=DBLE(KZ)/DBLE(IIC)
      IF(IIC.GT.MXIC) MXIC=IIC
      DO 80 N=1,NGIV
      DO 81 MW=1,MWEI
      NN=MWEI*(N-1)
      EIG(NN+MW,KK,IU)=ENR(N)
      DO 360 L=1,NCOMP
      COMP(NN+MW,KK,L,IU)=CC(N,L)
  360 CONTINUE
   81 CONTINUE
   80 CONTINUE
      NEIG(KK)=NGIV*MWEI
      WRITE(6,*)' NEIG=',NEIG(KK)
      IUD(KK)=IU
      GOTO 100
  110 CONTINUE
      KMAX=K
      WRITE(6,*) KMAX,NKP1
C      IF(KMAX.NE.NKP1)STOP ' AT 8825 IN ENRRED'
      JMAX=NEIG(1)
      DO 461 K=1,KMAX
      IF(JMAX.GT.NEIG(K))JMAX=NEIG(K)
  461 CONTINUE
      WRITE(6,*)'8830  JMAX=',JMAX
      KKMAX=KMAX
      RETURN
C
      ENTRY ENERDS(KF)
      DO 400 I=1,KMAX
      WRITE(KF,600) MRR(I),MRN(I),(KV(J,I),J=1,3),IC(I),NEIG(I)
  600 FORMAT(3X,A2,I2,'(',3I3,')/',I3,'     N.O.ENERGY=',I3)
      WRITE(KF,601) (EIG(J,I,1),J=1,NEIG(I))
  601 FORMAT(8F9.5)
  400 CONTINUE
      RETURN
      ENTRY BANDDS
      DO 500 J=1,JMAX
      EMIN=EIG(J,KMAX,1)
      EMAX=EIG(J,1,1)
      DO 510 I=1,KMAX
      IF(EMIN.GT.EIG(J,I,1))EMIN=EIG(J,I,1)
      IF(EMAX.LT.EIG(J,I,1))EMAX=EIG(J,I,1)
  510 CONTINUE
      WRITE(6,700)J,EMIN,EMAX
  700 FORMAT(3X,'==',I3,'TH BAND EMIN=',F7.4,' EMAX=',F7.4)
  500 CONTINUE
      RETURN
      END
C SUBROUTINE KPFIND ===*====3====*====4====*====5====*====6====*====7
C
C    IF STAR OF (KX,KY,KZ)/IC IS ALREADY IN KKM, IND=(INDEX OF KKM)
C    (KX,KY,KZ)/IC IS NEW                      , IND=0
C               1988.10.18  AKIRA YANASE
C               1988.12.20  H. HARIMA (KALRST YORI KAITEI)
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE KPFIND(KX,KY,KZ,IC,KKM,NK,IND)
        IMPLICIT REAL*8(A-H,O-Z)
      SAVE SPG2,SPG3,STK
      COMPLEX*16 CW
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &  ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      DIMENSION KKM(3,1),KKB(3),KA(3),KC(3)
      CALL ADDINV(NTLATC)
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
         JJ=J
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
    5    IND=JJ
   13 KB(1)=KKB(1)
      KB(2)=KKB(2)
      KB(3)=KKB(3)
      ICB=ICBB
      CALL REMINV
      RETURN
      END
      SUBROUTINE PDOS(EA,WT,NE,E0,DE,NW,IXMAX,IYMAX,IZMAX,JS,D,S)
C     ***  DENSITY OF STATES IS OBTAINED FOR A SINGLE BAND  ***
C     ***  BY LINEAR INTERPOLATION IN TETRAHEDORONS         ***
C     ***  EA  ENERGIES AT MESH POINTS
C     ***  WT  WEIGHTS AT MESH POINTS
C     ***  NE  NUMBER OF ENERGY POINTS                      ***
C     ***  E0  STARTING ENERGY                              ***
C     ***  DE  STEP OF ENERGY                               ***
C     ***  NW  NUMBER OF COMPONENTS TO BE CALCULATED
C     ***  IXMAX  NUMBER OF MESH POINTS IN B.Z.             ***
C     ***         IN X, Y, Z-DIRECTION FOR JS=3
C     ***         IN X,Y-DIRECTION FOR JS=2
C     ***         IN X-DIRECTION FOR JS=1
C     ***  IYMAX  NUMBER OF MESH POINTS IN B.Z.
C     ***         IN Y-DIRECTION FOR JS=1
C     ***  IZMAX  NUMBER OF MESH POINTS IN B.Z.             ***
C     ***         IN Z-DIRECTION FOR JS=1 OR 2
C     ***  JS=1  ORTHORHOMBIC
C     ***  JS=2  TETRAGONAL
C     ***  JS=3  CUBIC SYMMETRY
C     ***  D DENSITY OF STATES
C     ***  S NUMBER OF STATES
C     *** BY A. YANASE  JAN.1983                            ***
C     ***    MODIFIED FROM TDOSB
      IMPLICIT  REAL*8(A-H,O-Z)
      INCLUDE 'prmdos.f'
c     PARAMETER (LMCOMP=5,LNCOMP=LMCOMP+1)
      DIMENSION D(LNCOMP,*),S(LNCOMP,*),WT(LMCOMP,*),WET(LMCOMP)
      DIMENSION EA(*),EB(4),EC(8),ET(4),IQMAT(6,2),PWT(LMCOMP,8)
      DIMENSION ECUB(2,2,2),PWCUB(LMCOMP,2,2,2)
      EQUIVALENCE(EC(1),ECUB(1,1,1))
      EQUIVALENCE(PWT(1,1),PWCUB(1,1,1,1))
      DATA IQMAT/2,2,5,3,3,5, 4,6,6,4,7,7/
C     WRITE(6,*)' DE=',DE,E0,NE
      DO 250 I=1,NE
      DO 250 K=1,NW+1
      D(K,I)=0.
  250 S(K,I)=0.
      EMAX=E0
      EMIN=E0+DBLE(NE-1)*DE
c     WRITE(6,*)' *PDOS',E0,DE,NE,EMIN,EMAX
      DO 21 JX=1,IXMAX
      LX=JX-1
      JYM=JX
      IF(JS.EQ.1) JYM=IYMAX
      DO 22 JY=1,JYM
      LY=JY-1
      JZM=JY
      IF(JS.NE.3) JZM=IZMAX
      DO 23 JZ=1,JZM
      LZ=JZ-1
      IF(JS.EQ.3) NII=(LX*(LX+1)*(LX+2))/6+(LY*(LY+1))/2+LZ+1
      IF(JS.EQ.2) NII=((LX*(LX+1))/2+LY)*IZMAX+LZ+1
      IF(JS.EQ.1) NII=(LX*IYMAX+LY)*IZMAX+LZ+1
      IF(EA(NII).GT.EMAX) EMAX=EA(NII)
      IF(EA(NII).LT.EMIN) EMIN=EA(NII)
C     WRITE(6,*)JX,JY,JZ,NII,EA(NII)
   23 CONTINUE
   22 CONTINUE
   21 CONTINUE
      WRITE(6,*)EMIN,EMAX
      NSTA=(EMIN-E0)/DE
      NSTA=NSTA-2
      IF(NSTA.LT.1) NSTA=1
C     **  INTEGRATION OVER B.Z. STARTS    ***
      NTT=0
      K=0
      JXM=IXMAX-1
      DO 353 IX=1,JXM
      LX=IX-1
      JYM=IX
      IF(JS.EQ.1) JYM=IYMAX-1
      DO 352 IY=1,JYM
      LY=IY-1
      JZM=IY
      IF(JS.NE.3) JZM=IZMAX-1
      DO 351 IZ=1,JZM
      LZ=IZ-1
C     ***  ENERGIES AT CUBE CORNERS  ***
      IF(JS.EQ.3) NI=(LX*(LX+1)*(LX+2))/6+(LY*(LY+1))/2+LZ+1
      IF(JS.EQ.2) NI=IZMAX*((LX*(LX+1))/2+LY)+LZ+1
      IF(JS.EQ.1) NI=IZMAX*(IYMAX*LX+LY)+LZ+1
      NN=NI
      NP=0
      DO 321 KX=1,2
      DO 322 KY=1,2
      DO 323 KZ=1,2
      NP=NP+1
      ECUB(KX,KY,KZ)=EA(NN)
      DO 330 K=1,NW
  330 PWCUB(K,KX,KY,KZ)=WT(K,NN)
      NN=NN+1
  323 CONTINUE
      IF(JS.EQ.3) NN=NN-2+LY+KY
      IF(JS.EQ.2) NN=NN-2+IZMAX
      IF(JS.EQ.1) NN=NN-2+IZMAX
  322 CONTINUE
      IF(JS.EQ.3) NN=NI+(IX*(IX+1))/2
      IF(JS.EQ.2) NN=NI+IZMAX*(LX+1)
      IF(JS.EQ.1) NN=NI+IZMAX*IYMAX
  321 CONTINUE
C     *** SIX TETRAHEDORONS  ***
      ET(1)=EC(1)
      ET(4)=EC(8)
      IT=1
  325 CONTINUE
      DO 331 K=1,NW
  331 WET(K)=PWT(K,1)+PWT(K,8)
      DO 324 IP=1,2
      IQ=IQMAT(IT,IP)
      ET(IP+1)=EC(IQ)
      DO 332 K=1,NW
  332 WET(K)=WET(K)+PWT(K,IQ)
  324 CONTINUE
      DO 333 K=1,NW
  333 WET(K)=WET(K)/4.0
      NTT=NTT+1
      DO 340 M=1,4
      EB(M)=ET(M)
  340 CONTINUE
C     ***  EB(1).LE.EB(2).LE.EB(3).LE.EB(4)  ***
      DO 342 M=2,4
      J=M
  341 J=J-1
      IF(EB(J).LE.EB(J+1))GOTO 342
      C=EB(J)
      EB(J)=EB(J+1)
      EB(J+1)=C
      IF(J.NE.1)GOTO 341
  342 CONTINUE
C     ***  JEPSON AND ANDERSON FORMULA  ***
      E1=EB(1)
      E2=EB(2)
      E3=EB(3)
      E4=EB(4)
      ISTA=(E1-E0)/DE
      ISTA=ISTA+2
      IF(ISTA.LT.1) ISTA=1
      IEND=(E4-E0)/DE
      IEND=IEND+1
      IF(IEND.GT.NE) IEND=NE
      IF(IEND.LT.0) IEND=0
      IF(ISTA.GT.IEND) GO TO 513
      D21=E2-E1
      D31=E3-E1
      D41=E4-E1
      D32=E3-E2
      D42=E4-E2
      D43=E4-E3
      DEL=E4+E3-E2-E1
      EM=(E4*E3-E2*E1)/DEL
      A21=1.0/DEL
      A22=(D31*D41+D32*D42)/(DEL*DEL)
      A20=0.0
      IF(D32.GT.1.0E-4) A20=DEL/(D31*D41*D32*D42)
      A1=0.0
      IF(D21.GT.1.0E-4) A1=1.0/(D21*D31*D41)
      A3=0.0
      IF(D43.GT.1.0E-4) A3=1.0/(D41*D42*D43)
      DO 514 I=ISTA,IEND
      E=E0+DE*FLOAT(I-1)
      IF(E.LE.E2) GO TO 512
      IF(E.GE.E3) GO TO 511
      DM=E-EM
      WD=DM*DM*A20
      WS=WD*DM
      S(1,I)=S(1,I)+A22+3.0*DM*A21-WS
      D(1,I)=D(1,I)+A21-WD
      DO 334 K=1,NW
      S(K+1,I)=S(K+1,I)+(A22+3.0*DM*A21-WS)*WET(K)
  334 D(K+1,I)=D(K+1,I)+(A21-WD)*WET(K)
      GO TO 514
  511 D4=E-E4
      WD=D4*D4*A3
      WS=WD*D4
      S(1,I)=S(1,I)+1.0+WS
      D(1,I)=D(1,I)+WD
      DO 335 K=1,NW
      S(K+1,I)=S(K+1,I)+(1.0+WS)*WET(K)
  335 D(K+1,I)=D(K+1,I)+WD*WET(K)
      GO TO 514
  512 D1=E-E1
      WD=D1*D1*A1
      WS=WD*D1
      S(1,I)=S(1,I)+WS
      D(1,I)=D(1,I)+WD
      DO 336 K=1,NW
      S(K+1,I)=S(K+1,I)+WS*WET(K)
  336 D(K+1,I)=D(K+1,I)+WD*WET(K)
  514 CONTINUE
  513 IIEND=IEND+1
      IF(IIEND.GT.NE) GO TO 510
      DO 515 I=IIEND,NE
      S(1,I)=S(1,I)+1.0
      DO 337 K=1,NW
  337 S(K+1,I)=S(K+1,I)+WET(K)
  515 CONTINUE
  510 CONTINUE
C
      IT=IT+1
      IF(IT.EQ.7) GO TO 351
      IF(JS.EQ.2) GO TO 326
      IF(JS.EQ.1) GO TO 325
      IF((IT.EQ.2).AND.(IZ.EQ.IY)) IT=IT+2
      IF((IT.EQ.4).AND.(IY.EQ.IX)) GO TO 351
      IF((IT.EQ.6).AND.(IZ.EQ.IY)) GO TO 351
      GO TO 325
  326 IF((IT.EQ.4).AND.(IY.EQ.IX)) GO TO 351
      GO TO 325
C
  351 CONTINUE
  352 CONTINUE
  353 CONTINUE
      FACTS=1.0/FLOAT(NTT)
      FACT=3.0*FACTS
      DO 530 I=1,NE
      DO 530 K=1,NW+1
      D(K,I)=FACT*D(K,I)
      S(K,I)=S(K,I)*FACTS
  530 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ENSORT(EA,WT,IB,IXMAX,IYMAX,IZMAX,ISPIN,JS)
C
      IMPLICIT  REAL*8(A-H,O-Z)
C
      INCLUDE 'prmdos.f'
c     PARAMETER (LMNEIG=62,LMNKP1=800,LMNKP0=20000)
c     PARAMETER (LMCOMP=5,LNCOMP=LMCOMP+1)
      DIMENSION EA(LMNKP0),WT(LMCOMP,LMNKP0)
      SAVE /KPDATA/,/SPG2/,/KP0/,/KP1/,/EIGEN0/,/LCM/
      COMMON/KPDATA/ NX,NY,NZ,KINTM,WIDTH
      COMMON/SPG2  / IL,NG,IG(48),JV(2,3,48)
      COMMON/KP0   / KP0(4,LMNKP0),IKP0(LMNKP0),NKP0
      COMMON/KP1   / KP1(4,LMNKP1),NKST1(LMNKP1),NKP1,NKP2
      COMMON/EIGEN0/ EIG(LMNEIG,LMNKP1,2),NEIG,ISO
      CHARACTER*4 NCCHR
      COMMON/LCM/COMP(LMNEIG,LMNKP1,LMCOMP,2),
     &          IPEN(5,2,LMCOMP),NKCOMP(LMCOMP),NCOMP,NCCHR(LMCOMP)
      DIMENSION KB(3)
      IF(IB.GT.NEIG) STOP ' === STOP IN SUB.ENSORT. (IB>NEIG) ==='
      IF(IL.LT.-1) THEN
        STOP ' === STOP IN SUB.ENSORT. (IL.LT.-1) ==='
      ELSE IF(IL.EQ.-1) THEN
        NXX=NX*3
        NYY=NY*3
        NZZ=NZ*3
      ELSE IF(IL.LE.1) THEN
        NXX=NX
        NYY=NY
        NZZ=NZ
      ELSE IF(IL.LE.3) THEN
        NXX=NX*2
        NYY=NY*2
        NZZ=NZ*2
      ELSE IF(IL.EQ.4) THEN
        NXX=NX*2
        NYY=NY*2
        NZZ=NZ
      ELSE
        STOP ' === STOP IN SUB.ENSORT. (IL.GT.4) ==='
      END IF
      IXMAX=NXX+1
      IYMAX=NYY+1
      IZMAX=NZZ+1
C
C     ***  PREPARE ENERGIES AT THREE DIMENSIONAL MESH POINTS  ***
      ND=NX*NY*NZ
      CALL ZZZY53(NY*NZ,NZ*NX,NGC1)
      CALL ZZZY53(NX*NY,NGC1,NGCB)
      ND=ND/NGCB
      IKP00=0
      DO 21 IX=0,NXX
      KB(1)=IX*NY*NZ/NGCB
      NYM=IX
      IF(JS.EQ.1) NYM=NYY
      DO 22 IY=0,NYM
      KB(2)=NX*IY*NZ/NGCB
      NZM=IY
      IF(JS.NE.3) NZM=NZZ
      DO 23 IZ=0,NZM
      KB(3)=NX*NY*IZ/NGCB
      IKP00=IKP00+1
      IF(KB(1).NE.KP0(1,IKP00) .OR.
     &   KB(2).NE.KP0(2,IKP00) .OR.
     &   KB(3).NE.KP0(3,IKP00) .OR.
     &      ND.NE.KP0(4,IKP00)     ) THEN
        WRITE(6,*) 'IKP00=',IKP00,' KB=',(KB(I),I=1,3),' / ',ND
     &            ,'    KP0=',(KP0(I,IKP00),I=1,3),' / ', KP0(4,IKP00)
        STOP ' === STOP IN SUB.ENSORT. (KB.NE.KP0) ==='
      END IF
      LX=IX
      LY=IY
      LZ=IZ
      IF(JS.EQ.3) NII=(LX*(LX+1)*(LX+2))/6+(LY*(LY+1))/2+LZ+1
      IF(JS.EQ.2) NII=((LX*(LX+1))/2+LY)*IZMAX+LZ+1
      IF(JS.EQ.1) NII=(LX*IYMAX+LY)*IZMAX+LZ+1
      EA(NII)=EIG(IB,IKP0(IKP00),ISPIN)
C     WRITE(6,*)NII,IKP00,IKP0(IKP00),EA(NII)
      DO 140,ICOMP=1,NCOMP
      WT(ICOMP,NII)=COMP(IB,IKP0(IKP00),ICOMP,ISPIN)
  140 CONTINUE
C     WRITE(6,*) ' IKP00=',IKP00,' IKP0=',IKP0(IKP00),' IB=',IB,
C    &           '   E=',EA(NII)
   23 CONTINUE
   22 CONTINUE
   21 CONTINUE
      IF(NKP0.NE.IKP00) STOP ' === STOP IN SUB.ENSORT. (NKP0) ==='
      RETURN
      END
