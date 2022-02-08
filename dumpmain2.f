*#RUN *:OPT=3 INLINE=2 CORE=400 F=/LAGA2FRM/D1,R(01) ;
*#    F=/LAGA2FRM/BAND18,R(03) F=(07)
      IMPLICIT REAL*8 (A-H,O-Z)
      JF=3
      KF=7
      IF=1
      READ(5,*)JS
      READ(IF,*)
      READ(IF,*)IL,NG
      DO 10 N=1,NG
   10 READ(IF,*)
      READ(IF,*)AA,BB,CC
      READ(IF,*)
      READ(IF,*)NATM,NKAT
      AM=1.20
      READ(5,*)KXM,KYM,KZM
      E0=0.2D0
      NE=601
      DE=1.D-3
      IBM=28
      MX=9
      MY=5
      MZ=5
C     CALL SPLDOS(JF,IL,JS,NKAT,AM,KXM,KYM,KZM,E0,NE,DE,IBM,MX,MY,MZ,KF)
      CALL FERMIS(JF,IL,JS,NKAT,AM,KXM,KYM,KZM)
      STOP 'NORMAL END'
      END
      SUBROUTINE ENDUMP(KXM,KYM,KZM,JS,IL)
      DIMENSION IINV(3,8)
      DATA IINV/1,1,1, 1,1,-1, 1,-1,1, 1,-1,-1,
     &         -1,1,1, -1,1,-1, -1,-1,1, -1,-1,-1/
      DIMENSION EA(3000000)
      REAL*8 X(3),TRYFUN
      DO 40 ISXX=1,8
      ISX=IINV(1,ISXX)
      ISY=IINV(2,ISXX)
      ISZ=IINV(3,ISXX)
      N=0
      MIC=1
      IF(IL.EQ.1)MIC=2
      DO 10 IX=1,KXM
      IMY=IX
      IF(JS.EQ.1) IMY=KYM
      DO 20 IY=1,IMY
      IMZ=IY
      IF(JS.NE.3) IMZ=KZM
      DO 30 IZ=1,IMZ
      N=N+1
      X(1)=0.D0
      X(2)=0.D0
      X(3)=0.D0
      IF(KXM.NE.1)X(1)=REAL(IX-1)*ISX/(KXM-1)
      IF(IMY.NE.1)X(2)=REAL(IY-1)*ISY/(KYM-1)
      IF(IMZ.NE.1)X(3)=REAL(IZ-1)*ISZ/(KZM-1)/MIC
      EA(N)=TRYFUN(X,1)
      WRITE(8,"(3F6.2,2X,F10.7)")(X(I),I=1,3),EA(N)
   30 CONTINUE
   20 CONTINUE
   10 CONTINUE
      NMAX=N
      WRITE(7,'(6I5,I7)') ISX,ISY,ISZ,KXM,KYM,KZM,NMAX
      do 50 NN=1,NMAX
      WRITE(7,*) EA(NN)
   50 continue
   40 CONTINUE
      RETURN
      END
C
      SUBROUTINE FERMIS(JF,IL,JS,NKAT,AM,KXM,KYM,KZM)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NBMAX=91999,KPMAX=48999)
      PARAMETER (NEMAX=100,IBMAX=2001)
      EXTERNAL TRYFUN
      DIMENSION BKV(3,KPMAX),ENR(KPMAX,NEMAX),NGIV(KPMAX)
      DIMENSION EIG(KPMAX),BK(3)
      DATA C1/1.D-1/,C2/1.D-2/
      REAL*4 EF,DR(3),XC(3)
      DATA DR/0.0,1.0,0.0/
      DATA XC/0.0,1.0,0.0/
      CALL CLOCKM(IT0)
c     CALL FPARAM(1,132)
      WRITE(6,*)'  BEFORE ENRRED'
      CALL ENRRED(JF,NKAT,BKV,ENR,NGIV,NKP)
      WRITE(6,*)'  AFTER ENRRED'
C     KF=7
C     CALL ENERDS(KF)
      CALL CLOCKM(IT1)
      CALL SPLIN(BKV,NKP,AM,KXM,KYM,KZM,IL,JS,C1,C2)
      CALL CLOCKM(IT2)
      WRITE(6,701) (IT2-IT1)/1000.0,(IT2-IT0)/1000.0
  701 FORMAT(3X,'TIME FOR SPLIN=',F10.3,
     & ' SEC.  ALL CPU TIME=',F10.3,' SEC')
C
      ID=0
C
      DO 1201 IB=1,1
      WRITE(6,1605) IB
 1605 FORMAT(' ***   DENSITY OF STATES FOR ',I3,'TH BAND ***')
      DO 1211 J=1,NKP
      EIG(J)=ENR(J,IB)
 1211 CONTINUE
      CALL CLOCKM(IT1)
      CALL AMGIV(EIG)
      CALL CLOCKM(IT2)
      WRITE(6,702) (IT2-IT1)/1000.0,(IT2-IT0)/1000.0
  702 FORMAT(3X,'TIME FOR AMGIV=',F10.3,
     & ' SEC.  ALL CPU TIME=',F10.3,' SEC')
      CALL CLOCKM(IT1)
C     EF=0.4934
C     CALL FRMCUT(EF,DR,XC)
      CALL ENRCHK(BKV,ENR,IB,NKP)
      CALL ENDUMP(61,61,61,JS,IL)
      CALL CLOCKM(IT2)
      WRITE(6,703) (IT2-IT1)/1000.0,(IT2-IT0)/1000.0
  703 FORMAT(3X,'TIME FOR CHECK=',F10.3,
     & ' SEC.  ALL CPU TIME=',F10.3,' SEC')
 1201 CONTINUE
      RETURN
      END
C
      SUBROUTINE ENRCHK(BKV,ENR,IB,NKP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KPMAX=48999,NEMAX=100)
      DIMENSION BKV(3,KPMAX),ENR(KPMAX,NEMAX)
      DO 10 IK=1,NKP
      EE=TRYFUN(BKV(1,IK),1)
      DE=ENR(IK,IB)-EE
      WRITE(6,'(3F8.4,3F10.5)') (BKV(J,IK),J=1,3),ENR(IK,IB),EE,DE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE SPLDOS(JF,IL,JS,NKAT,AM,KXM,KYM,KZM,E0,NE,DE,IBM,
     & MX,MY,MZ,KF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NBMAX=91999,KPMAX=48999)
      PARAMETER (NEMAX=100,IBMAX=2001)
      EXTERNAL TRYFUN
      DIMENSION BKV(3,KPMAX),ENR(KPMAX,NEMAX),NGIV(KPMAX)
      DIMENSION EIG(KPMAX),BK(3)
      DIMENSION D(IBMAX),S(IBMAX),DD(IBMAX),SS(IBMAX)
      DATA C1/1.D-1/,C2/1.D-2/
      CALL CLOCKM(IT0)
c     CALL FPARAM(1,132)
      WRITE(6,*)'  BEFORE ENRRED'
      CALL ENRRED(JF,NKAT,BKV,ENR,NGIV,NKP)
      WRITE(6,*)'  AFTER ENRRED'
C     KF=7
C     CALL ENERDS(KF)
      CALL CLOCKM(IT1)
      CALL SPLIN(BKV,NKP,AM,KXM,KYM,KZM,IL,JS,C1,C2)
      CALL CLOCKM(IT2)
      WRITE(6,701) (IT2-IT1)/1000.0,(IT2-IT0)/1000.0
  701 FORMAT(3X,'TIME FOR SPLIN=',F10.3,
     & ' SEC.  ALL CPU TIME=',F10.3,' SEC')
C
      ID=0
C
      DO 1101 I=1,NE
      SS(I)=0.D0
 1101 DD(I)=0.D0
      DO 1201 IB=1,1
      WRITE(6,1605) IB
 1605 FORMAT(' ***   DENSITY OF STATES FOR ',I3,'TH BAND ***')
      DO 1211 J=1,NKP
      EIG(J)=ENR(J,IB)
 1211 CONTINUE
      CALL CLOCKM(IT1)
      CALL AMGIV(EIG)
      CALL CLOCKM(IT2)
      WRITE(6,702) (IT2-IT1)/1000.0,(IT2-IT0)/1000.0
  702 FORMAT(3X,'TIME FOR AMGIV=',F10.3,
     & ' SEC.  ALL CPU TIME=',F10.3,' SEC')
      CALL CLOCKM(IT1)
      CALL ENRCHK(BKV,ENR,IB,NKP)
      CALL TDOSA(TRYFUN,NE,E0,DE,MX,MY,MZ,ID,JS,D,S)
C
      DO 1031 I=1,NE
      DD(I)=DD(I)+D(I)
 1031 SS(I)=SS(I)+S(I)
 1201 CONTINUE
C
      WRITE(6,1602) MX-1,MY-1,MZ-1
 1602 FORMAT(' *** MESH MX=',I3,'  MY=',I3,'  MZ=',I3,'  ***')
      DO 1001 I=1,NE
      E=E0+DE*(I-1)
C     IF(DD(I).GT.1.D-4) WRITE(6,1601) E,DD(I),SS(I)
      IF(DD(I).GT.1.D-4) WRITE(KF,1601) E,DD(I),SS(I)
 1601 FORMAT(F8.4,F9.4,F8.4)
 1001 CONTINUE
      CALL CLOCKM(IT2)
      WRITE(6,703) (IT2-IT1)/1000.0,(IT2-IT0)/1000.0
  703 FORMAT(3X,'TIME FOR CHECK=',F10.3,
     & ' SEC.  ALL CPU TIME=',F10.3,' SEC')
      RETURN
      END
      SUBROUTINE SPLIN(BKV,NKP,AM,KXM,KYM,KZM,IL,JS,C1,C2)
C**********************************************************************
C************ LIBRARY HEADER ******************************************
C     ***PREPARATION OF INTERPOLATION FOR 3-D FUNCTION              ***
C     *** BKV:MESH POINTS                                           ***
C     *** NKP:NUMBER OF MESH POINTS (LESS THAN KPMAX)               ***
C     *** AM:AM*NKP IS MINIMUM OF NUMBER OF BASIC FUNCTIONS         ***
C     ***                                 (LESS THAN NBMAX)         ***
C     *** KXM:MINIMUM MESH OF BASIC FUNCTIONS  X DIRECTION          ***
C     *** KYM:MINIMUM MESH OF BASIC FUNCTIONS  Y DIRECTION          ***
C     *** KZM:MINIMUM MESH OF BASIC FUNCTIONS  Z DIRECTION          ***
C     *** IL=1:SIMPLE LATTICE                                       ***
C     *** IL=2:FACE CENTERED LATTICE                                ***
C     *** IL=3:BODY CENTERED LATTICE                                ***
C     *** IL=4:BASE CENTERED LATTICE                                ***
C     *** JS=1:ORTHORHOMBIC                                         ***
C     *** JS=2:TETRAGONAL                                           ***
C     *** JS=3:CUBIC SYMMETRY                                       ***
C     *** C1,C2:PARAMETERS TO GIVE C                                ***
C     ***                                    BY H.HARIMA    SEP.1986***
C     ***                            KAITEI  BY H.HARIMA    APR.1989***
C**********************************************************************
      PARAMETER (NBMAX=91999,KPMAX=48999)
      IMPLICIT REAL*8 (A-H,O-Z)
      complex(kind(0d0)) FMXN,EI,ANN,C
      COMMON/BCOF/MRBC(48*NBMAX,3),INM(48*NBMAX),NRM(NBMAX),MRIC
      COMMON/LTIC/IIL,JJS,NNKP,NNBM,AAM,IWT,KKXM,KKYM,KKZM
      DIMENSION BKV(3,KPMAX)
      DIMENSION FMXN(NBMAX,KPMAX),ALOU(NBMAX),A(KPMAX,KPMAX)
      COMMON/MATC/C(NBMAX,KPMAX)
      DIMENSION IW(2*KPMAX)
      EI=cmplx(0.0,1.0,kind(0d0))
      IIL=IL
      JJS=JS
      NNKP=NKP
      AAM=AM
      KKXM=KXM
      KKYM=KYM
      KKZM=KZM
C
C     GIVE BASIC FUNCTIONS TO COMMON/BCOF/
C
      CALL GEBCOE
      CALL MRBCDS(1)
C     JF=1
C     CALL MRBCDS(JF)
      WRITE(6,102) C1,C2
  102 FORMAT(1X,'*** PARAMETER C1=',D12.5,' C2=',D12.5,'***')
      NBM=NNBM
      NBF=NBM/IWT
      DO 10 I=1,NBF
      J0=IWT*(I-1)+1
      JE=J0+IWT-1
      DO 20 NX=1,NKP
      BKX=BKV(1,NX)
      BKY=BKV(2,NX)
      BKZ=BKV(3,NX)
      FAC=2.0*3.14159265358D0
      FAC=FAC/4.0
      WA=0.D0
      WB=0.D0
      DO 25 J=J0,JE
      WWA=BKX*MRBC(J,1)+BKY*MRBC(J,2)+BKZ*MRBC(J,3)
      WA=WA+COS(FAC*DBLE(WWA)/MRIC)
      WB=WB+SIN(FAC*DBLE(WWA)/MRIC)
   25 CONTINUE
c     print *,I,NX,WA,WB,BKX,BKY,BKZ,(MRBC(J0,I3),I3=1,3)
c     IF(DABS(WB).GT.1.D-6) STOP 'IMAGINARY PART'
      FMXN(I,NX)=cmplx(WA/DBLE(IWT),WB/DBLE(IWT),kind(0d0))
   20 CONTINUE
      ALOU(I)=1.0/(1.D0+C1*NRM(I)+C2*NRM(I)**2)
   10 CONTINUE
      DO 30 N1=1,NKP
      DO 30 N2=1,N1
      ANN=cmplx(0.0,0.0)
      DO 40 I=1,NBF
      ANN=ANN+ALOU(I)*FMXN(I,N1)*conjg(FMXN(I,N2))
   40 continue
      if(aimag(ANN).GT.1d-6) then
        print *,ANN
        stop "IMAGINARY PART"
      end if
      A(N1,N2)=real(ANN)
      A(N2,N1)=A(N1,N2)
   30 CONTINUE
C
C     WRITE(6,*) 'MATRIX SIZE =',NKP,'X',NKP
C     DO 510 N1=1,NKP
C     WRITE(6,*) 'N1=',N1
C     WRITE(6,511) (A(N1,N2),N2=1,NKP)
C 510 CONTINUE
C
      E=1.0D-15
      CALL TINVDD(NKP,A,DT,E,KPMAX,IW,INDER)
      IF(INDER.NE.0) STOP 'AFTER TINBDD BECAUSE INDER.NE.0'
C
C     WRITE(6,*) 'INDER=',INDER
C     DO 520 N1=1,NKP
C     WRITE(6,511) (A(N1,N2),N2=1,NKP)
C 511 FORMAT(3X,6F10.5)
C 520 CONTINUE
C
      DO 50 I=1,NBF
      DO 50 N1=1,NKP
      ANN=cmplx(0.0,0.0)
      DO 60 N2=1,NKP
      ANN=ANN+A(N2,N1)*FMXN(I,N2)
   60 CONTINUE
      C(I,N1)=ALOU(I)*ANN
   50 CONTINUE
      RETURN
      END
      SUBROUTINE AMGIV(ENR)
C**********************************************************************
C     ***GIVE COEFFICIENT AM                                        ***
C     *** ENR:ENERGIES AT MESH POINTS                               ***
C     ***                                    BY H.HARIMA    SEP.1986***
C***********************************************************************
      PARAMETER (NBMAX=91999,KPMAX=48999)
      IMPLICIT REAL*8 (A-H,O-Z)
      complex(kind(0d0)) C,AM,WA
      COMMON/LTIC/IL,JS,NKP,NBF,AAM,IWT,KXM,KYM,KZM
      COMMON/MATC/C(NBMAX,KPMAX)
      COMMON/COEF/AM(NBMAX)
      DIMENSION ENR(KPMAX)
      NB=NBF/IWT
      DO 10 I=1,NB
      WA=(0.0,0.0)
      DO 20 N=1,NKP
      WA=WA+C(I,N)*ENR(N)
   20 CONTINUE
      AM(I)=WA
   10 CONTINUE
      RETURN
      END
      REAL*8 FUNCTION TRYFUN(BK,ID)
C**********************************************************************
C     ***INTERPOLATED FUNCTION                                      ***
C     *** BK:POINT TO BE CALCULATED                                 ***
C     *** ID:SUB-SUFFIX OF FUNCTION                                 ***
C     ***                                    BY H.HARIMA    SEP.1986***
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      complex(kind(0d0)) ANN,AM
      PARAMETER (NBMAX=91999,KPMAX=48999)
      COMMON/BCOF/MRBC(48*NBMAX,3),INM(48*NBMAX),NRM(NBMAX),MRIC
      COMMON/LTIC/IL,JS,NKP,NBF,AAM,IWT,KXM,KYM,KZM
      COMMON/COEF/AM(NBMAX)
      DIMENSION BK(3)
      FAC=2.0*3.14159265358D0
      FAC=FAC/4.0
      ANN=cmplx(0.0,0.0)
      DO 20 I=1,NBF
      WB=FAC*(BK(1)*MRBC(I,1)+BK(2)*MRBC(I,2)+BK(3)*MRBC(I,3))/MRIC
      ANN=ANN+AM(INM(I))*cmplx(COS(WB),SIN(WB),kind(0d0))
   20 CONTINUE
      TRYFUN=real(ANN/DBLE(IWT))
      RETURN
      END
c     SUBROUTINE CLOCKM(T)
c     INTEGER T
c     REAL*8 TIME
c     CALL CLOCK(TIME)
c     T=1000.*TIME
C     CALL CPTIME(T)
c     RETURN
c     END
      SUBROUTINE ENRRED(JF,NKAT,BKV,ENR,NGIV,NKP)
C**********************************************************************
C************ LIBRARY HEADER ******************************************
C     *** READ ROUTINE FROM SHORT E(K) FILE                         ***
C     *** JF:E(K) FILE                                              ***
C     *** NKAT:NUMBER OF KINDS OF ATOMS                             ***
C     *** BKV:MESH POINTS                                           ***
C     *** ENR:ENERGIES OF MESH POINTS                               ***
C     *** NGIV:NUMBER OF ENERGIES OF MESH POINTS                    ***
C     *** NKP:NUMBER OF MESH POINTS                                 ***
C     ***                                    BY H.HARIMA    SEP.1986***
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KPMAX=48999,NEMAX=100)
      DIMENSION BKV(3,KPMAX),ENR(KPMAX,NEMAX),NGIV(KPMAX)
      DIMENSION KV(3,KPMAX),IC(KPMAX)
      DIMENSION EIG(NEMAX)
      K=0
  100 CONTINUE
      NDIG=1
      NEIG=1
      READ(JF,*,END=110) KX,KY,KZ,IIC,EIG(NEIG)
  150 FORMAT(1X,I3,2X,3I2,2X,2I3,2X,A4,2I3)
C 150 FORMAT(I3,2X,3I2,1X,2I3,2X,A4,2I3)
      K=K+1
C     WRITE(6,*)'K=',K,KX,KY,KZ,IIC,NDIG,NEIG
      KV(1,K)=KX
      KV(2,K)=KY
      KV(3,K)=KZ
      IC(K)=IIC
      BKV(1,K)=DBLE(KX)/DBLE(IIC)
      BKV(2,K)=DBLE(KY)/DBLE(IIC)
      BKV(3,K)=DBLE(KZ)/DBLE(IIC)
      J=0
      DO 80 N=1,NEIG
      DO 80 I=1,NDIG
      J=J+1
      ENR(K,J)=EIG(N)
   80 CONTINUE
      NGIV(K)=J
      GOTO 100
  110 CONTINUE
      NKP=K
      RETURN
C
c     ENTRY ENERDS(KF)
c     DO 400 I=1,NKP
c     WRITE(KF,600)(KV(J,I),J=1,3),IC(I),NGIV(I)
c 600 FORMAT(3X,'(',3I3,')/',I3,'     N.O.ENERGY=',I3)
c     WRITE(KF,601) (ENR(I,J),J=1,NGIV(I))
c 601 FORMAT(8F8.4)
c 400 CONTINUE
c     RETURN
      END
*#RUN *:FORM NLNO OPT=3 IAP
C     IMPLICIT REAL*8 (A-H,O-Z)
C     COMMON/LTIC/IL,JS,KPM,NBM,AM,IWT,KXM,KYM,KZM
C     IL=4
C     JS=1
C     KPM=2
C     AM=1.1D0
C     KXM=2
C     KYM=2
C     KZM=2
C     WRITE(6,*) ' ENTER IL,JS,KPM,AM,KXM,KYM,KZM'
C     WRITE(6,*)  IL,JS,KPM,AM,KXM,KYM,KZM
C     READ(5,*,END=99) IL,JS,KPM,AM,KXM,KYM,KZM
C  99 CALL GEBCOE
C     CALL MRBCDS(1)
C     STOP
C     END
      SUBROUTINE GEBCOE
C**********************************************************************
C     ***GENERATION OF BASIC FUNCTIONS                              ***
C     *** CALL FROM SPLIN                                           ***
C     ***                                    BY H.HARIMA    SEP.1986***
C     ***                          KAITEI    BY H.HARIMA    APR.1989***
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NBMAX=91999,KPMAX=48999)
      COMMON/LTIC/IL,JS,KPM,NBM,AM,IWT,KXM,KYM,KZM
      DIMENSION IROT(3,6)
      DATA IROT/1,2,3, 2,1,3, 1,3,2,
     &          3,1,2, 2,3,1, 3,2,1/
      DIMENSION IINV(3,8)
      DATA IINV/1,1,1, 1,1,-1, 1,-1,1, 1,-1,-1,
     &         -1,1,1, -1,1,-1, -1,-1,1, -1,-1,-1/
      DIMENSION IWIT(3)
      DATA IWIT/8,16,48/
      COMMON/BCOF/MRBC(48*NBMAX,3),INM(48*NBMAX),NRM(NBMAX),MRIC
      DIMENSION MRB(20*NBMAX,3),NORM(20*NBMAX)
      IF(AM.LT.1.0) STOP 'AM SHOULD BE GREATER THAN 1.0'
C
c     IWT=IWIT(JS)
      IWT=1
      KXMM=KXM
      KYMM=KYM
      KZMM=KZM
      IVOL=1
      IF(IL.EQ.2)IVOL=2
      IF(IL.EQ.3)IVOL=4
      IF(IL.EQ.4)IVOL=2
      IF(JS.EQ.3) THEN
   33 IF(AM*KPM*IVOL.GT.(KXM+1)*(KXM+2)*(KXM+3)/6*8)THEN
      KXM=KXM+1
      GOTO 33
      END IF
      END IF
      IF(JS.EQ.2) THEN
   32 IF(AM*KPM*IVOL.GT.(KXM+1)*(KXM+1)*(KZM+2)/2*8)THEN
      KXM=KXM+1
      KZM=KZM+1
      GOTO 32
      END IF
      END IF
      IF(JS.EQ.1) THEN
   31 IF(AM*KPM*IVOL.GT.(KXM+1)*(KYM+1)*(KZM+1)*8)THEN
      KXM=KXM+1
      KYM=KYM+1
      KZM=KZM+1
      GOTO 31
      END IF
      END IF
c     IF(IL.EQ.-1)STOP 'NOT SUPPORT IL=-1'
      IF(IL.EQ.0)STOP 'NOT SUPPORT IL=0'
C
      MRIC=2
      IF(IL.EQ.1)MRIC=1
      WRITE(6,220) IL,JS
  220 FORMAT(1X,'***  IL=',I2,'     JS=',I2,25X,'***')
      WRITE(6,240) KPM
  240 FORMAT(1X,'***     NUMBER OF GIVEN MESH POINTS     =',I4,'***')
      WRITE(6,230) KXMM,AM,KXM
  230 FORMAT(1X,'***  GIVEN KXM=',I3,'  GIVEN AM=',F4.2,'  KXM=',I3,
     &3X,'***')
      IF(JS.EQ.1) THEN
      WRITE(6,231) KYMM,AM,KYM
  231 FORMAT(1X,'***  GIVEN KYM=',I3,'  GIVEN AM=',F4.2,'  KYM=',I3,
     &3X,'***')
      ENDIF
      IF(JS.LE.2)THEN
      WRITE(6,232) KZMM,AM,KZM
  232 FORMAT(1X,'***  GIVEN KZM=',I3,'  GIVEN AM=',F4.2,'  KZM=',I3,
     &3X,'***')
      ENDIF
C
      IK=1
      DO 41 IX=-KXM,KXM
      KY=IX
      IF(JS.EQ.1) KY=KYM
      DO 43 IY=-KY,KY
      KZ=IY
      IF(JS.NE.3) KZ=KZM
      IDZ=1
      IF(IL.EQ.4)IDZ=2
      DO 45 IZ=-KZ,KZ,IDZ
      IF(IL.EQ.0) GO TO 460
      IF(IL.EQ.1) GO TO 460
      IF(IL.EQ.2) GO TO 463
      IF(IL.EQ.3) GO TO 462
      IF(IL.EQ.4) GO TO 464
      IF(IL.EQ.-1) GO TO 461
      IF(MOD(IABS(IX+IY),2).NE.0) GO TO 45
      GO TO 460
  461 ISUM=-IX+IY+IZ
      IF(MOD(IABS(ISUM),3).NE.0) GO TO 45
      GO TO 460
  464 IF(MOD(IABS(IX+IY),2).NE.0) GO TO 45
      GO TO 460
  463 IF(MOD(IABS(IX+IY+IZ),2).NE.0) GO TO 45
      GO TO 460
  462 IF(MOD(IABS(IX),2).NE.MOD(IABS(IY),2)) GO TO 45
      IF(MOD(IABS(IY),2).NE.MOD(IABS(IZ),2)) GO TO 45
  460 CONTINUE
C
C     REGISTRATION
C
      MRB(IK,1)=IX
      MRB(IK,2)=IY
      MRB(IK,3)=IZ
      NORM(IK)=IX*IX+IY*IY+IZ*IZ
      IK=IK+1
      IF(IK.GT.20*NBMAX) GO TO 53
   45 CONTINUE
   43 CONTINUE
   41 CONTINUE
C
      IK=IK-1
      WRITE(6,210) IK
  210 FORMAT(1X,'***  NUMBER OF GENERATED BASIC FUNCTIONS=',I4,'***')
      NIM=1
      NRMX=1
      IF(JS.EQ.3) NRMX=6
      IF(JS.EQ.2) NRMX=2
      LK=1
      DO 20 JK=1,IK
      DO 10 NR=1,NRMX
      DO 10 NI=1,NIM
      MRBC(LK,1)=IINV(1,NI)*MRB(JK,IROT(1,NR))
      MRBC(LK,2)=IINV(2,NI)*MRB(JK,IROT(2,NR))
      MRBC(LK,3)=IINV(3,NI)*MRB(JK,IROT(3,NR))
      INM(LK)=JK
      LK=LK+1
   10 CONTINUE
      NRM(JK)=NORM(JK)
   20 CONTINUE
      NBM=LK-1
      RETURN
C
C     WRITE ROUTINE:JF=1 4-PARTS;JF=OTERS ALL
C
      ENTRY MRBCDS(JF)
      DO 110 I=1,NBM/IWT
      WRITE(9,700) I,NRM(I)
  700 FORMAT(3X,I3,' NORM=',I5)
      IJ0=IWT*(I-1)+1
      IJE=IWT*I
      IF(JF.EQ.1)IJE=IJ0+3
      WRITE(9,710) ((MRBC(IJ,K),K=1,3),IJ=IJ0,IJE)
  710 FORMAT(4(2X,' (',3I4,')'))
  110 CONTINUE
      RETURN
   53 WRITE(9,660)
  660 FORMAT(' STOP AT 53 IN GEBCOE')
      STOP
      END
CTINVDD       0002F ; COMPUTER CENTER TOHOKU UNIVERSITY
C  ********** LIBRARY HEADER ********************
C00  F1/CTU0002F/TINVDD;COMPUTER CENTER TOHOKU UNIVERSITY
C10  TINVDD;MATRIX INVERSION,SWEEP OUT    (DOUBLE PRECISION)
C20  SUBROUTINE
C30  710430;SAKATA MASATO
C40
C50
C60  ACOS SERIES 77 NEAC SYSTEM 700;ACOS-6,R3.1
C61  ACOS SERIES 77 NEAC SYSTEM 900;ACOS-6,R5.2(BIN,HEX)
C70  LP=6
C80  CM= 579W
C90
C  **********************************************
      SUBROUTINE    TINVDD(N,A,DT,E,NN,IW,INDER)
C     INVERSION OF DOUBLE PRECISION MATRIX  A(NXN)
C     SWEEP OUT, COMPLETE POSITIONING
      DOUBLE PRECISION     A,DT,E,AA,AZ,WORK,EPS,PIV
      DIMENSION     A(NN,N)  ,IW(1)
      INDER=0
      IF(N-1)    910,930,101
  101 IF(N.GT.NN)    GO TO  900
      EPS=0.0D0
      DT=1.0D0
      DO  100     K=1,N
      PIV=0.0D0
      DO  110       I=K,N
      DO  110       J=K,N
      IF(DABS(A(I,J)).LE.DABS(PIV))  GO TO  110
      IPIV=I
      JPIV=J
      PIV=A(I,J)
  110 CONTINUE
      DT=DT*PIV
      IF(DABS(PIV).LE.EPS)   GO TO  920
      IF(K.EQ.1)   EPS=DABS(PIV)*E
      IF(IPIV.EQ.K)      GO TO 130
      DT=-DT
      DO 120   J=1,N
      WORK=A(IPIV,J)
      A(IPIV,J)=A(K,J)
  120 A(K,J)=WORK
  130 IF(JPIV.EQ.K)      GO TO  150
      DT=-DT
      DO 140   I=1,N
      WORK=A(I,JPIV)
      A(I,JPIV)=A(I,K)
  140 A(I,K)=WORK
  150 IW(2*K-1)=IPIV
      IW(2*K)=JPIV
      AA=1.0D0/PIV
      DO 210   J=1,N
  210 A(K,J)=A(K,J)*AA
      DO 220  I=1,N
      IF(I.EQ.K)    GO TO  220
      AZ=A(I,K)
      IF(AZ.EQ.0.0D0)    GO TO  220
      DO 230   J=1,N
  230 A(I,J)=A(I,J)-A(K,J)*AZ
      A(I,K)=-AA*AZ
  220 CONTINUE
  100 A(K,K)=AA
      DO  400 KK=2,N
      K=N+1-KK
      IJ=IW(2*K)
      IF(IJ.EQ.K)   GO TO  420
      DO 410   J=1,N
      WORK=A(IJ,J)
      A(IJ,J)=A(K,J)
  410 A(K,J)=WORK
  420 IJ=IW(2*K-1)
      IF(IJ.EQ.K)   GO TO  400
      DO 430   I=1,N
      WORK=A(I,IJ)
      A(I,IJ)=A(I,K)
  430 A(I,K)=WORK
  400 CONTINUE
      RETURN
  910 INDER=-1
      WRITE(6,691)  N
  691 FORMAT(1H0,8H *** N = ,I5,5X,30HN SHOULD BE POSITIVE IN TINVDD)
      RETURN
  900 INDER=-1
      WRITE(6,690)    N,NN
  690 FORMAT(1H0,8H *** N =  ,I5,5X,4HNN =,I5/
     & 55H    N SHOULD BE LESS THAN OR EQUAL TO NN IN TINVDD       )
      RETURN
  920 DT=0.0D0
      INDER=N-K+1
      NNN=K-1
      WRITE(6,692)  NNN
  692 FORMAT(1H0,48H*** GIVEN MATRIX A TO TINVDD IS ILL CONDITIONED
     &,24HOR SINGULAR WITH RANK = ,I4,27H.   RETURN WITH NO FURTHER
     &,11HCALCULATION       )
      RETURN
  930 DT=A(1,1)
      K=1
      IF(DT.EQ.0.0D0)   GO TO  920
      A(1,1)=1.0D0/A(1,1)
      RETURN
      END
CTDOSA        0268F ; COMPUTER CENTER TOHOKU UNIVERSITY
C  ********** LIBRARY HEADER ********************
C00  Y9/CTU0268F/TDOSA ;COMPUTER CENTER TOHOKU UNIVERSITY
C10  TDOSA ;DENSITY OF STATES-A
C20  SUBROUTINE
C30  790723;YANASE AKIRA,HASEGAWA AKIRA
C40
C50
C60  ACOS SERIES 77 NEAC SYSTEM 900;ACOS-6,R7.1
C70
C80  CM=7440W
C90
C  **********************************************
      SUBROUTINE TDOSA(ENR,NE,E0,DE,IXMAX,IYMAX,IZMAX,IB,JS,D,S)
C     ***  DENSITY OF STATES IS OBTAINED FOR A SINGLE BAND  ***
C     ***  BY LINEAR INTERPOLATION IN TETRAHEDORONS         ***
C     ***  NE  NUMBER OF ENERGY POINTS                      ***
C     ***  E0  STARTING ENERGY                              ***
C     ***  DE  STEP OF ENERGY                               ***
C     ***  IXMAX  NUMBER OF MESH POINTS IN B.Z.             ***
C     ***         IN X, Y, Z-DIRECTION FOR JS=3
C     ***         IN X,Y-DIRECTION FOR JS=2
C     ***         IN X-DIRECTION FOR JS=1
C     ***  IYMAX  NUMBER OF MESH POINTS IN B.Z.
C     ***         IN Y-DIRECTION FOR JS=1
C     ***  IZMAX  NUMBER OF MESH POINTS IN B.Z.             ***
C     ***         IN Z-DIRECTION FOR JS=1 OR 2
C     ***  IB  BAND NUMBER                                  ***
C     ***  JS=1  ORTHORHOMBIC
C     ***  JS=2  TETRAGONAL
C     ***  JS=3  CUBIC SYMMETRY
C     ***  D DENSITY OF STATES
C     ***  S NUMBER OF STATES
C     *** BY A. YANASE  APR.1977                            ***
C     ***    MODIFIED IN FEB. 1979
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KB
      DIMENSION D(1),S(1)
      DIMENSION EA(6545),EB(4),EC(8),ECUB(2,2,2),ET(4),IQMAT(6,2)
      DIMENSION KB(3)
      EQUIVALENCE(EC(1),ECUB(1,1,1))
      DATA IQMAT/2,2,5,3,3,5, 4,6,6,4,7,7/
      DO 250 I=1,NE
      D(I)=0.D0
  250 S(I)=0.D0
C     ***  PREPARE ENERGIES AT THREE DIMENSIONAL MESH POINTS  ***
      PX=0.5D0/DBLE(IXMAX-1)
      PY=PX
      IF(JS.EQ.1) PY=0.5D0/DBLE(IYMAX-1)
      PZ=PX
      IF(JS.NE.3) PZ=0.5D0/DBLE(IZMAX-1)
      EMAX=E0
      EMIN=E0+DBLE(NE-1)*DE
      DO 21 JX=1,IXMAX
      LX=JX-1
      KB(1)=PX*DBLE(LX)
      JYM=JX
      IF(JS.EQ.1) JYM=IYMAX
      DO 22 JY=1,JYM
      LY=JY-1
      KB(2)=PY*DBLE(LY)
      JZM=JY
      IF(JS.NE.3) JZM=IZMAX
      DO 23 JZ=1,JZM
      LZ=JZ-1
      KB(3)=PZ*DBLE(LZ)
      IF(JS.EQ.3) NII=(LX*(LX+1)*(LX+2))/6+(LY*(LY+1))/2+LZ+1
      IF(JS.EQ.2) NII=((LX*(LX+1))/2+LY)*IZMAX+LZ+1
      IF(JS.EQ.1) NII=(LX*IYMAX+LY)*IZMAX+LZ+1
      EA(NII)=ENR(KB,IB)
      IF(EA(NII).GT.EMAX) EMAX=EA(NII)
      IF(EA(NII).LT.EMIN) EMIN=EA(NII)
   23 CONTINUE
   22 CONTINUE
   21 CONTINUE
      NSTA=(EMIN-E0)/DE
      NSTA=NSTA-2
      NEND=(EMAX-E0)/DE
      NEND=NEND+3
      IF(NSTA.LT.1) NSTA=1
      IF(NEND.GT.NE) NEND=NE
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
      DO 321 KX=1,2
      DO 322 KY=1,2
      DO 323 KZ=1,2
      ECUB(KX,KY,KZ)=EA(NN)
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
      DO 324 IP=1,2
      IQ=IQMAT(IT,IP)
      ET(IP+1)=EC(IQ)
  324 CONTINUE
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
      IF(ISTA.GT.IEND) GO TO 513
      D21=E2-E1
      D31=E3-E1
      D41=E4-E1
      D32=E3-E2
      D42=E4-E2
      D43=E4-E3
      DEL=E4+E3-E2-E1
      EM=(E4*E3-E2*E1)/DEL
      A21=1.D0/DEL
      A22=(D31*D41+D32*D42)/(DEL*DEL)
      A20=0.D0
      IF(D32.GT.1.0D-10) A20=DEL/(D31*D41*D32*D42)
      A1=0.D0
      IF(D21.GT.1.0D-10) A1=1.D0/(D21*D31*D41)
      A3=0.D0
      IF(D43.GT.1.0D-10) A3=1.D0/(D41*D42*D43)
      DO 514 I=ISTA,IEND
      E=E0+DE*DBLE(I-1)
      IF(E.LE.E2) GO TO 512
      IF(E.GE.E3) GO TO 511
      DM=E-EM
      WD=DM*DM*A20
      WS=WD*DM
      S(I)=S(I)+A22+3.0*DM*A21-WS
      D(I)=D(I)+A21-WD
      GO TO 514
  511 D4=E-E4
      WD=D4*D4*A3
      WS=WD*D4
      S(I)=S(I)+1.D0+WS
      D(I)=D(I)+WD
      GO TO 514
  512 D1=E-E1
      WD=D1*D1*A1
      WS=WD*D1
      S(I)=S(I)+WS
      D(I)=D(I)+WD
  514 CONTINUE
  513 IIEND=IEND+1
      IF(IIEND.GT.NEND) GO TO 510
      DO 515 I=IIEND,NEND
  515 S(I)=S(I)+1.D0
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
      FACTS=1.D0/DBLE(NTT)
      FACT=3.0*FACTS
      DO 530 I=1,NE
      D(I)=FACT*D(I)
      S(I)=S(I)*FACTS
  530 CONTINUE
      NNEND=NEND+1
      IF(NNEND.GT.NE) GO TO 531
      DO 532 I=NNEND,NE
  532 S(I)=1.D0
  531 CONTINUE
      RETURN
      END
