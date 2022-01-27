C************************************************************************
C      BAND PLOT PROGRAM                                                *
C               PROGRAMMED BY H.HARIMA                                  *
C       SPIN FERRO TAIOU  1989/10/30                                    *
C               MODIFIED BY A.YANASE
C                         1994/02/10
C     01:D1 FILE                                                        *
C     02:ENERGY OUTPUT FILE                                             *
C     03:FORMAT FILE                                                    *
C     40:40-BAN FILE                                                    *
C************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SCL/WX,WY,EO,EM,XM,YM,IPR,ISO
      CHARACTER*28 D
C
      CALL TSPPRP
      WRITE(6,*) ' END OF TSPPRP'
C
      CALL STRUCT(NLCOMP,NSPIN,JMARK)
      WRITE(6,*) ' END OF STRUCT'
       write(6,*) 'jmark=',jmark
C
      CALL ENERED(NLCOMP)
      WRITE(6,*) ' END OF ENERED'
C
      IPRR=IPR
      JD=0
      IF(ISO.EQ.3) JD=1
      CALL AXENER(JD,IPRR,NSPIN)
      WRITE(6,*) ' END OF AXENR' 
C
c      CALL TEST(7)
      CALL DRFSTR(0)
      CALL CHIOP(8)
      CALL DEGREE
      CALL ORIGIN(50.0,40.0)
      write(6,*) ' END OF ORIGIN'
C
      CALL PLOTR
      write(6,*) ' END OF PLOTR'
C
      IF(NSPIN.NE.2) IUD=1
      IF(NSPIN.EQ.2) IUD=2
      CALL PLOTSC(IUD)
      IF(JMARK.EQ.1) CALL MARKPL(IUD)
      IF(NSPIN.EQ.3) THEN
         CALL SETLN2(2.0,1.0,0.0,1)
         CALL PLOTSC(2)
         IF(JMARK.EQ.1) CALL MARKPL(2)
      END IF
      write(6,*) ' END OF PLOTSC'
      CALL SETPEN(6)
      READ(3,*,END=10)EF
      GOTO 20
   10 WRITE(6,*)' FERMI LEVEL'
      READ(5,*,ERR=30)EF
   20 CALL SETPEN(5)
      CALL FERMIE(EF)
      write(6,*) ' END OF FERMIE'
   30 READ(3,'(A28)',err=40)D
      GOTO 50
   40 WRITE(6,*)' TITLE <28 1****+****2****+***'
      READ(5,'(A28)',END=60)D
   50 CALL SETPEN(4)
      CALL TITLE(D,28)
      WRITE(6,*) ' END OF TITLE'
   60 CALL DRFEND
      STOP
      END
      SUBROUTINE TITLE(D,K)
C***********************************************************************
C WRITE TITLE                                                          *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SCL/WX,WY,EO,EM,XM,YMM,IPR,ISO
      CHARACTER*28 D
C     REAL*4 YM
      YM=YMM
      CALL BTYPE(2.0,YM+5.0,7.0,D,0.0,K,1,1,-1.0)
      RETURN
      END
      SUBROUTINE FERMIE(EF)
C***********************************************************************
C WRITE FERMI ENERGY                                                   *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/SCL/WX,WY,EO,EM,XM,YM,IPR,ISO
C     REAL*4 XO,XE,YF,XMM,YMM
      YMM=YM
      XMM=XM
      YF=(EF-EO)*WY
C         write(6,*) naxm,'fermie','yf=',yf,ymm
      IF(YF.GT.YMM.OR.YF.LT.0.0) RETURN
      XO=0.0
      DO 10 N=1,NAXM
C         write(6,*) n,icon(n),naxm
      IF(ICON(N).EQ.1.OR.N.EQ.NAXM) THEN
         XE=WX*PS(2,N)
C          write(6,*) xo,yf,xe
         CALL LINE3(XO,YF,XE,YF,1.0,2.0,1)
         IF(N.NE.NAXM) XO=WX*PS(1,N+1)
      ENDIF
   10 CONTINUE
      CALL BTYPE(XMM+2.0,YF-2.5,5.0,'E$DF',0.0,4,1,1,-1.0)
      RETURN
      END
      SUBROUTINE TSPPRP
C***********************************************************************
C PREPERATION TO CALL TSPACE                                           *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SCL/WX,WY,EO,EM,XM,YM,IPR,ISO
      DIMENSION JB(2,3)
      READ(1,*)
      READ(1,*) IL,NGEN,INV
      CALL TSPACE(IL)
      DO 1 I=1,NGEN
      READ(1,*) JA,((JB(J,K),J=1,2),K=1,3)
      CALL TSGENR(JA,JB)
    1 CONTINUE
      CALL TSPGRP(INV)
      IF(IPR.GE.3) CALL TSPGDS
      READ(1,104) A,B,C
      READ(1,104) CA,CB,CC
  104 FORMAT(3D23.16)
      CALL TSLATC(A,B,C,CA,CB,CC)
C 102 FORMAT(14I5)
C     READ(1,102) NATOM,NKAT
C     NK=NKAT
C     CALL ENRRED(2,NK)
C     PRINT ,'  END OF ENRRED'
C     IF(IPR.GE.4)CALL ENERDS(6)
      RETURN
      END
      SUBROUTINE STRUCT(NLCOMP,NSPIN,JMARK)
C***********************************************************************
C SELECTION OF POINTS TO BE PLOTTED                                    *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)                
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX)
     &   ,MK(3,NAXMAX),NAXM
      COMMON/KKEND/KK1M(3,NAXMAX),KK2M(3,NAXMAX)
     &      ,ICC1M(NAXMAX),ICC2M(NAXMAX)
      COMMON/SCL/WX,WY,EO,EM,XM,YM,IPR,ISO
      DIMENSION KK1(3),KK2(3),KK3(3),MR(3),KTV(3),KG(3)
      CHARACTER*2 CCR1,CCR2,CCR3
      CHARACTER*4 CHARA4
      CHARACTER*10 MAGNET
C
      READ(3,*) MAGNET
      CHARA4=MAGNET
      IF(CHARA4.EQ.'NONM') ISO=1
      IF(CHARA4.EQ.'MAGN') ISO=2
      IF(CHARA4.EQ.'SPIN') ISO=3
      READ(3,*) NLCOMP,NSPIN
      IF((ISO.EQ.1.OR.ISO.EQ.3).AND.NSPIN.NE.0) THEN
         WRITE(6,*) ' FOR PARAMAGNETIC BANDS NSPIN MUST BE 0',NSPIN
         STOP
      END IF
      IF(ISO.EQ.2.AND.(NSPIN.LE.0.OR.NSPIN.GT.3)) THEN
         WRITE(6,*) ' FOR FERROMAGNETIC BANDS'
         WRITE(6,*) ' NSPIN=1 ONLY FOR UP   SPIN BANDS'
         WRITE(6,*) ' NSPIN=2 ONLY FOR DOWN SPIN BANDS'
         WRITE(6,*) ' NSPIN=3 FOR BOTH SPIN BANDS'
         STOP
      END IF
      READ(3,*) IPR,JMARK
      READ(3,*) EO,EM,YM,XM
      WY=YM/(EM-EO)
C
      PS(1,1)=0.0
      NK=0
      READ(3,*) NAXM
      IF(NAXM.GT.NAXMAX) THEN
         WRITE(6,*) ' NAXMAX=',NAXMAX,' NAXM=',NAXM
         STOP
      END IF
      DO 99 NAX=1,NAXM
      READ(3,*) (KK1(J),J=1,3),ICC1,(KK2(J),J=1,3),ICC2
      IF(IPR.GE.2) WRITE(6,600) (KK1(J),J=1,3),ICC1,(KK2(J),J=1,3),ICC2
  600 FORMAT( '(',3I4,')/',I4,'---','(',3I4,')/',I4)
C
      DO 75 J=1,3
      IIX=J
      IF(KK1(J)*ICC2.NE.KK2(J)*ICC1) GOTO 81
   75 CONTINUE
      write(6,*)  'KK1 AND KK2 ARE THE SAME POINT'
      STOP
   81 IX(NAX)=IIX
      DO 71 J=1,3
      KK1M(J,NAX)=KK1(J)
      KK2M(J,NAX)=KK2(J)
   71 CONTINUE
      ICC1M(NAX)=ICC1
      ICC2M(NAX)=ICC2
      KTV(1)=KK2(1)*ICC1-KK1(1)*ICC2
      KTV(2)=KK2(2)*ICC1-KK1(2)*ICC2
      KTV(3)=KK2(3)*ICC1-KK1(3)*ICC2
      ICC=ICC1*ICC2
      CALL ZZZY37(KTV(1),KTV(2),KTV(3),ICC,WW)
      RT(NAX)=WW
      BLK=0.0
      ICON(NAX)=0
      IF(NAX.NE.1) THEN
        CALL EQUIKK(KK1,ICC1,KK3,ICC3,KG,IND)
        IF(IND.EQ.0) THEN
            ICON(NAX-1)=1
            BLK=0.3*ABS(RT(1))
        ENDIF
        PS(1,NAX)=PS(2,NAX-1)+BLK
      ENDIF
      PS(2,NAX)=RT(NAX)+PS(1,NAX)
C      IF(KTV(IIX).LT.0.) RT(NAX)=-RT(NAX)
      IF(IPR.GE.2) WRITE(6,601) PS(1,NAX),PS(2,NAX),RT(NAX)
  601 FORMAT(3F10.5)
      DO 34 J=1,3
      KK3(J)=KK2(J)
   34 CONTINUE
      ICC3=ICC2
      KX=KK1(1)
      KY=KK1(2)
      KZ=KK1(3)
      ICC=ICC1
      CALL KPNAME(KX,KY,KZ,ICC,CCR1,KG)
      READ(CCR1,'(A2)')MR(1)
      MK(1,NAX)=MR(1)
      KX=KK2(1)
      KY=KK2(2)
      KZ=KK2(3)
      ICC=ICC2
      CALL KPNAME(KX,KY,KZ,ICC,CCR3,KG)
      READ(CCR3,'(A2)')MR(1)
      MK(3,NAX)=MR(1)
      KX=KK1(1)*ICC2+KK2(1)*ICC1
      KY=KK1(2)*ICC2+KK2(2)*ICC1
      KZ=KK1(3)*ICC2+KK2(3)*ICC1
      ICC=ICC1*ICC2*2
      CALL KPNAME(KX,KY,KZ,ICC,CCR2,KG)
      READ(CCR2,'(A2)')MR(1)
      MK(2,NAX)=MR(1)
      IF(IPR.GE.2) WRITE(6,602) CCR1,CCR2,CCR3
  602 FORMAT(3(3X,A2))
   99 CONTINUE
      WX=XM/PS(2,NAXM)
      RETURN
      END
      SUBROUTINE PLOTR
C***********************************************************************
C PLOT FRAMES AND SCALES                                               *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/SCL/WX,WY,EO,EM,XM,YMM,IPR,ISO
C     REAL*4 X,Y,XO,XE,YM,E
      WRITE(6,*) ' NAXM=',NAXM
      YM=YMM
      XO=0.
      DO 10 N=1,NAXM
      IF(IPR.GE.3) WRITE(6,600) N,PS(1,N),PS(2,N),ICON(N)
  600 FORMAT(I4,2F10.6,I2)
      IF(ICON(N).EQ.1.OR.N.EQ.NAXM) THEN
         XE=WX*PS(2,N)
         CALL LINE2(XO,0.0,2)
         CALL LINE2(XE,0.0,1)
         CALL LINE2(XE,YM,2)
         CALL LINE2(XO,YM,1)
         XO=WX*PS(1,N+1)
      ENDIF
   10 CONTINUE
      DO 301 I=1,NAXM
      JM=1
      IF(ICON(I).EQ.1.OR.I.EQ.NAXM) JM=2
      DO 302 J=1,JM
      X=PS(J,I)*WX
      CALL LINE2(X,0.0,2)
      CALL LINE2(X,YM,1)
  302 CONTINUE
  301 CONTINUE
c      WRITE(6,*) ' 301'
      CALL LINE2(-7.0,0.0,2)
      CALL LINE2(-10.0,0.0,1)
      CALL LINE2(-10.0,YM,1)
      CALL LINE2(-7.0,YM,1)
      ESDV=0.1
      IF(EM-EO.LT.0.1) ESDV=0.01
      E0A=(EO-0.00001)/ESDV
      IF(E0A.GE.0.0) I0=INT(E0A)+1
      IF(E0A.LT.0.0) I0=INT(E0A)
      EMA=(EM+0.00001)/ESDV
      IF(EMA.GE.0.0) IM=INT(EMA)
      IF(EMA.LT.0.0) IM=INT(EMA)-1
c      WRITE(6,*) IM,I0,ESDV,E0A,EMA
      DO 100 I=1,IM-I0+1
      E=(I0+I-1)*ESDV
      Y=(E-EO)*WY
      IF(ABS(Y).LT.1.0) GO TO 101
      IF(ABS(Y-YM).LT.1.0) GO TO 101
      CALL LINE2(-10.0,Y,2)
      CALL LINE2(-7.0,Y,1)
c      WRITE(6,*) I,E,Y
  101 IF(ESDV.GE.0.1)
     &CALL FTYPE2(-13.0,Y-2.5,5.0,E,0.0,4,1,2,1,-1.0)
      IF(ESDV.LT.0.1)
     &CALL FTYPE2(-13.0,Y-2.5,5.0,E,0.0,4,2,2,1,-1.0)
  100 CONTINUE
C      write(6,*) ' btype'
      CALL BTYPE(-35.0,YM/2.0,5.0,'Energy (Ry)',90.0,11,3,1,-1.0)
C      write(6,*) ' btype'
      DO 110 N=1,NAXM
      MMM=MK(1,N)
      X=WX*PS(1,N)
c      write(6,*) ' ctype',x
      CALL CTYPE(X,-10.0,5.0,MMM)
      MMM=MK(2,N)
      X=WX*(PS(1,N)+PS(2,N))/2.0
c      write(6,*) ' ctype',x
      CALL CTYPE(X,-7.5,4.0,MMM)
      IF(ICON(N).EQ.1.OR.N.EQ.NAXM) THEN
      MMM=MK(3,N)
      X=WX*PS(2,N)
c      write(6,*) ' ctype',x
      CALL CTYPE(X,-10.0,5.0,MMM)
      ENDIF
  110 CONTINUE
      RETURN
      END
      SUBROUTINE CTYPE(X,Y,H,M)
C***********************************************************************
C WRITE NAMES OF POINTS AND AXES                                       *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4 D
C     REAL*4 X,Y,H
      WRITE(D(1:2),'(A2)')M
      K=2
      IF(D(2:2).EQ.' ')K=1
      IF(D(1:2).EQ.'GM') THEN
      D(1:3)='$GC'
      K=3
      ELSEIF(D(1:2).EQ.'XD') THEN
      D(1:1)='X'
      K=1
      ELSEIF(D(1:2).EQ.'SM') THEN
      D(1:3)='$GS'
      K=3
      ELSEIF(D(1:2).EQ.'LD') THEN
      D(1:3)='$GL'
      K=3
      ELSEIF(D(1:2).EQ.'DT') THEN
      D(1:3)='$GD'
      K=3
      else if(D(1:2).eq.'TP') then
      D(1:2)="T'"
      else if(D(1:2).eq.'SP') then
      D(1:2)="S'"
      ELSEIF(D(1:2).EQ.'MX') THEN
      D(1:1)='M'
      K=1
      ELSEIF(D(1:2).EQ.'LX') THEN
      D(1:1)='L'
      K=1
      ENDIF
      CALL BTYPE(X,Y,H,D(1:K),0.0,K,3,1,-1.0)
      RETURN
      END
      SUBROUTINE PLOTSC(IUD)
C***********************************************************************
C PLOT BAND                                                            *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)
      PARAMETER (MZXKPT=500,MAXEIG=960)
      PARAMETER (MAXFKP=1300,MAXIRP=367)
      PARAMETER (N65=65)
      COMMON/AXE/EAX(N65,MAXEIG,MAXIRP),NKPT(MAXIRP)
     &      ,NLIN(MAXIRP),JRR(MAXIRP),JUD(MAXIRP),KPM
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/SCL/WX,WY,EO,EM,XM,YMM,IPR,ISO
C     REAL*4 YM,XC,XA,YA,XB,YB
      CHARACTER*4 DDD
      YM=YMM
      DO 110 K=1,KPM
C        write(6,*) iud,k,jud(k)
      IF(IUD.NE.JUD(K)) GO TO 110
      NP=NKPT(K)
      WRITE(DDD(1:2),'(A2)') MK(2,JRR(K))
      IF(IPR.GE.1) WRITE(6,600) K,NP,DDD,NLIN(K)
  600 FORMAT(2I5,1X,A2,I4)
      IF(NLIN(K).LE.0) GO TO 110
      DO 100 J=1,NLIN(K)
      IF(IPR.GE.4) WRITE(6,*) ' BNAD NO=',J
      XA=WX*PS(1,JRR(K))
      YA=(EAX(1,J,K)-EO)*WY
      IF(YA.GT.0.0.AND.YA.LT.YM) CALL LINE4(XA,YA,2)
C      IF(IPR.GE.2) WRITE(6,*) '2',XA,YA
      DO 180 I=2,NP
C
      XB=XA
      YB=YA
      XA=PS(1,JRR(K))+(RT(JRR(K))*(I-1))/(NP-1)
      XA=XA*WX
      YA=(EAX(I,J,K)-EO)*WY
C      write(6,*) xa,ya,eax(i,j,k)
      IF(YB.GT.YM) GO TO 291
      IF(YA.LT.YM) GO TO 295
      XC=XB+(XA-XB)*((YM-YB)/(YA-YB))
      CALL LINE4(XC,YM,1)
      GO TO 180
  291 IF(YA.GT.YM) GO TO 180
      XC=XB+(XA-XB)*((YB-YM)/(YB-YA))
      CALL LINE4(XC,YM,2)
      GO TO 299
  295 IF(YB.LT.0.0) GO TO 296
      IF(YA.GT.0.0) GO TO 299
      XC=XB+(XA-XB)*(YB/(YB-YA))
      CALL LINE4(XC,0.0,1)
      GO TO 180
  296 IF(YA.LT.0.0) GO TO 180
      XC=XB+(XA-XB)*(-YB/(YA-YB))
      CALL LINE4(XC,0.0,2)
  299 CALL LINE4(XA,YA,1)
C      WRITE(6,*) '1',XA,YA
  180 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
      SUBROUTINE ENERED(NLCOMP)
      PARAMETER (MAXKPT=3600,MAXEIG=960)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*2 MRR,MRRM
      COMMON/ENR/KXM(3,MAXKPT),ICM(MAXKPT),IUDM(MAXKPT)
     &   ,MRRM(MAXKPT),MRNM(MAXKPT),MWEIM(MAXKPT),NSTM(MAXKPT)
     &   ,NEIGM(MAXKPT),EIGM(MAXEIG,MAXKPT),NKPINT
      COMMON/MRK/NSIT(MAXKPT),MSIT(3,MAXKPT),XIMS(3,MAXKPT)
      DIMENSION KX(3),EIG(MAXEIG),KBB(3,10),KG(3,10)
      REWIND 2 
C
      N=0
   40 READ(2,150,END=41) KK,(KX(J),J=1,3),IC,IUD,MRR,MRN,MWEI,NST,NEIG
      print *," ENERGY READ KK=",KK
      IF(NEIG.GT.MAXEIG) THEN
         WRITE(6,*) ' MAXEIG=',MAXEIG,' NEIG=',NEIG
         STOP
      END IF
      IF(NEIG.GT.0) READ(2,*) (EIG(J),J=1,NEIG)
  150 FORMAT(I4,2X,3I3,1X,2I3,2X,A2,I2,3I3)
C  151 FORMAT(8F8.4)
      IF(NEIG.EQ.0) GO TO 40
      IF(NLCOMP.NE.0) THEN
         NFACT=(NLCOMP-1)/2+1
         DO 50 K=1,NEIG*NFACT
   50       READ(2,*)
      END IF
      N=N+1
      IF(N.GT.MAXKPT) THEN
         WRITE(6,*) ' MAXKPT=',MAXKPT,' NKPINT=',N
         STOP
      END IF
      CALL NEAREC(KX,IC,KBB,KG,NG)
      DO 10 J=1,3
         KXM(J,N)=KBB(J,1)
   10 CONTINUE
      ICM(N)=IC
C       write(6,*) n,(kxm(j,n),j=1,3),icm(n)
      IUDM(N)=IUD
      MRRM(N)=MRR
      MRNM(N)=MRN
      MWEIM(N)=MWEI
      NSTM(N)=NST
      NEIGM(N)=NEIG
      DO 20 J=1,NEIG
        EIGM(J,N)=EIG(J)
   20 CONTINUE
      GO TO 40
   41 CONTINUE
      NKPINT=N
      DO 1 I=1,NKPINT
         NSIT(I)=0
    1 CONTINUE
      RETURN
      END
      SUBROUTINE MARKPL(IUD)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXKPT=3600,MAXEIG=960)
      PARAMETER (NAXMAX=20)                
      CHARACTER*2 MRRM
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/ENR/KXM(3,MAXKPT),ICM(MAXKPT),IUDM(MAXKPT)
     &   ,MRRM(MAXKPT),MRNM(MAXKPT),MWEIM(MAXKPT),NSTM(MAXKPT)
     &   ,NEIGM(MAXKPT),EIGM(MAXEIG,MAXKPT),NKPINT
      COMMON/SCL/WX,WY,EO,EM,XMM,YMM,IPR,ISO
      COMMON/MRK/NSIT(MAXKPT),MSIT(3,MAXKPT),XIMS(3,MAXKPT)

C     REAL*4 XM,YM,XA,YA
      DIMENSION NMARK(12)
      DATA NMARK/1,2,3,4,5,6,7,8,9,11,12,13/
C      DATA NMARK/73, 74, 75, 76, 80, 63, 68, 94, 70, 69, 72, 71/
C
      XM=XMM
      YM=YMM
      DO 10 I=1,NKPINT
         IF(IUD.NE.IUDM(I)) GO TO 10
C       write(6,*) i,mrrm(i),mrnm(i),nsit(i)
         IF(NSIT(I).EQ.0) GO TO 10
         DO 11 N=1,NSIT(I)
            XA=(PS(1,MSIT(N,I))+RT(MSIT(N,I))*XIMS(N,I))*WX
            DO 12 KK=1,NEIGM(I)
               YA=(EIGM(KK,I)-EO)*WY
               IF(YA.LT.YM.AND.YA.GT.0.0) THEN
                  CALL MARK(XA,YA,4.0,NMARK(MRNM(I)))
               END IF
   12       CONTINUE
   11    CONTINUE
   10 CONTINUE
      RETURN
      END  
      SUBROUTINE AXENER(JDOUB,IPR,NSPIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)                
      PARAMETER (MAXKPT=3600,MAXEIG=960)
      PARAMETER (MAXFKP=1300,MAXIRP=367)
      PARAMETER (N65=65)
      CHARACTER*2 MRRM
      COMMON/AXE/EAX(N65,MAXEIG,MAXIRP),NKPT(MAXIRP)
     &      ,NLIN(MAXIRP),JRR(MAXIRP),JUD(MAXIRP),NAXEN
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/KKEND/KK1M(3,NAXMAX),KK2M(3,NAXMAX)
     &      ,ICC1M(NAXMAX),ICC2M(NAXMAX)
      COMMON/ENR/KXM(3,MAXKPT),ICM(MAXKPT),IUDM(MAXKPT)
     &   ,MRRM(MAXKPT),MRNM(MAXKPT),MWEIM(MAXKPT),NSTM(MAXKPT)
     &   ,NEIGM(MAXKPT),EIGM(MAXEIG,MAXKPT),NKPINT
      DIMENSION KB(3),JTR(12),IPA(12),ND(12)
      DIMENSION NRECPT(3,20)
      JD=JDOUB
      NAXEN=0
      DO 11 I=1,NAXM
         II=I
         IF(IPR.GE.2) WRITE(6,601) II
  601    FORMAT(' RECFND',I3)
         CALL RECFND(II,NRECPT,NRPT,IPR)
         IF(IPR.GE.2) WRITE(6,602) II,NRPT
  602    FORMAT(' KPFIEN',I3,' NUMBER OF RECPR. LAT. POINTS=',I3)
         CALL KPFIEN(II,NRECPT,NRPT,IPR)
         DO 21 K=1,3
            KB(K)=KK1M(K,I)*ICC2M(I)+KK2M(K,I)*ICC1M(I)
   21    CONTINUE
         IC=ICC1M(I)*ICC2M(I)*2
         CALL TSIREP(KB,IC,JD)
         CALL DGTRST(JDO,NR,NH,NSTR,ND,JTR,IPA)
         IF(IPR.GE.2) WRITE(6,603) KB,IC
  603    FORMAT(' CMPTBL (',3I3,')/',I4)
         CALL CMPTBL(JD,II,KB,IC,NR,NH,IPR)   
         DO 22 IR=1,NR
            IF(JTR(IR).EQ.0.AND.IPA(IR).LT.IR) GO TO 22
            IIR=IR
            NAXEN=NAXEN+1
            IF(NAXEN.GT.MAXIRP) THEN
               WRITE(6,*) ' MAXIRP=',MAXIRP,' NAXEN=',NAXEN
               STOP
            END IF
            IF(NSPIN.NE.2) IUD=1
            IF(NSPIN.EQ.2) IUD=2 
            CALL INTPOR(JD,II,KB,IC,NH,IIR,IUD,IPR)
            JUD(NAXEN)=IUD
            IF(NSPIN.EQ.3) THEN
               NAXEN=NAXEN+1
               IF(NAXEN.GT.MAXIRP) THEN
                  WRITE(6,*) ' MAXIRP=',MAXIRP,' NAXEN=',NAXEN
                  STOP
               END IF
               IUD=2
               CALL INTPOR(JD,II,KB,IC,NH,IIR,IUD,IPR)
               JUD(NAXEN)=IUD
            END IF
   22    CONTINUE
   11 CONTINUE
      RETURN
      END
      SUBROUTINE RECFND(II,NRECPT,NRPT,IPR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)                
C
      COMMON/KKEND/KK1M(3,NAXMAX),KK2M(3,NAXMAX)
     &      ,ICC1M(NAXMAX),ICC2M(NAXMAX)
C
      DIMENSION KB(3),NRECPT(3,20),KG(3,10),KBB(3,10)
      JJ=0
      DO 10 J=1,4
         IF(J.EQ.1) THEN
            DO 11 K=1,3
               KB(K)=2*KK1M(K,II)*ICC2M(II)-KK2M(K,II)*ICC1M(II)
   11       CONTINUE
            ICC=ICC2M(II)*ICC1M(II)
         ELSE IF(J.EQ.2) THEN
            DO 12 K=1,3
               KB(K)=KK1M(K,II)
   12       CONTINUE
            ICC=ICC1M(II)
         ELSE IF(J.EQ.3) THEN
            DO 13 K=1,3
               KB(K)=KK2M(K,II)
   13       CONTINUE
            ICC=ICC2M(II)
         ELSE IF(J.EQ.4) THEN
            DO 14 K=1,3
               KB(K)=2*KK2M(K,II)*ICC1M(II)-KK1M(K,II)*ICC2M(II)
   14       CONTINUE
            ICC=ICC1M(II)*ICC2M(II)
         END IF
         CALL NEAREC(KB,ICC,KBB,KG,NNG)
         DO 16 N=1,NNG
            IF(JJ.EQ.0) GO TO 17
            DO 18 IJ=1,JJ
               DO 19 K=1,3
                  IF(KG(K,N).NE.NRECPT(K,IJ)) GO TO 18
   19          CONTINUE
               GO TO 16
   18       CONTINUE
   17       JJ=JJ+1
            DO 15 K=1,3
               NRECPT(K,JJ)=KG(K,N)
   15       CONTINUE
            IF(IPR.GE.4) WRITE(6,601) II,J,JJ,(NRECPT(K,JJ),K=1,3)
  601       FORMAT(6I5)
   16    CONTINUE
   10 CONTINUE
      NRPT=JJ
      RETURN
      END
      SUBROUTINE KPFIEN(II,NRECPT,NRPT,IPR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)                
      PARAMETER (MAXKPT=3600,MAXEIG=960)
      PARAMETER (MAXFKP=1300,MAXIRP=367)
      CHARACTER*2 MRRM
C
      COMPLEX*16 CW
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
C
      COMMON/KKEND/KK1M(3,NAXMAX),KK2M(3,NAXMAX)
     &      ,ICC1M(NAXMAX),ICC2M(NAXMAX)
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/ENR/KXM(3,MAXKPT),ICM(MAXKPT),IUDM(MAXKPT)
     &   ,MRRM(MAXKPT),MRNM(MAXKPT),MWEIM(MAXKPT),NSTM(MAXKPT)
     &   ,NEIGM(MAXKPT),EIGM(MAXEIG,MAXKPT),NKPINT
      COMMON/AXP/IFKP(MAXFKP),JSKP(MAXFKP),KKB(3,MAXFKP)
     &   ,ICCB(MAXFKP),ICPTBL(12,MAXFKP),INDM(MAXFKP),NPOINT
C
      DIMENSION KB(3),KBC(3)
      DIMENSION NRECPT(3,20),JIND(MAXKPT),IMM(MAXKPT)
C
      CALL ADDINV(NTLATC)
      X1=DBLE(KK1M(1,II))/ICC1M(II)
      Y1=DBLE(KK1M(2,II))/ICC1M(II)
      Z1=DBLE(KK1M(3,II))/ICC1M(II)
      DX=DBLE(KK2M(1,II))/ICC2M(II)-X1
      DY=DBLE(KK2M(2,II))/ICC2M(II)-Y1
      DZ=DBLE(KK2M(3,II))/ICC2M(II)-Z1
C       write(6,*) x1,y1,z1
C       write(6,*) dx,dy,dz
      D6=1.0D-6
      NPOINT=0
      DO 10 I=1,NKPINT
         JIND(I)=0
   10 CONTINUE
      DO 1 I=1,NKPINT
         IF(JIND(I).NE.0) GO TO 1
         KB(1)=KXM(1,I)
         KB(2)=KXM(2,I)
         KB(3)=KXM(3,I)
         IC=ICM(I)
c          write(6,*) 'kb',kb,ic
         NJ=1
         IMM(1)=I
         IF(I.LE.NKPINT) THEN
            DO 4 JI=I+1,NKPINT
               IF(JIND(JI).NE.0) GO TO 4
               DO 5 K=1,3
                  IF(KB(K)*ICM(JI).NE.KXM(K,JI)*IC) GO TO 4
    5          CONTINUE
               NJ=NJ+1
               IMM(NJ)=JI
    4       CONTINUE
         END IF
C-----------------------------------------
C        STAR OF K IS OBTAINED
         CALL TSIREP(KB,IC,0)
         CALL ZZZY38
C------------------------------------------
         DO 2 J=1,NS
c          write(6,*) j,(ks(k,j),k=1,3),icbb
           DO 3 KGG=1,NRPT
              XK=DBLE(KS(1,J))/ICBB+DBLE(NRECPT(1,KGG))
              YK=DBLE(KS(2,J))/ICBB+DBLE(NRECPT(2,KGG))
              ZK=DBLE(KS(3,J))/ICBB+DBLE(NRECPT(3,KGG))
c           write(6,*) xk,yk,zk
              IF(DABS(DX).LE.D6.AND.DABS(XK-X1).GT.D6) GO TO 3 
              IF(DABS(DY).LE.D6.AND.DABS(YK-Y1).GT.D6) GO TO 3
              IF(DABS(DZ).LE.D6.AND.DABS(ZK-Z1).GT.D6) GO TO 3
              IF(DABS(DX).GT.D6) THEN
                 IF(DABS(DY).GT.D6) THEN
                    IF(DABS((XK-X1)/DX-(YK-Y1)/DY).GT.D6) GO TO 3
                 END IF
                 IF(DABS(DZ).GT.D6) THEN
                    IF(DABS((XK-X1)/DX-(ZK-Z1)/DZ).GT.D6) GO TO 3
                 END IF
              END IF
              IF(DABS(DY).GT.D6.AND.DABS(DZ).GT.D6) THEN
                 IF(DABS((YK-Y1)/DY-(ZK-Z1)/DZ).GT.D6) GO TO 3
              END IF 
              IF(IX(II).EQ.1.AND.
     &         ((XK-X1)/DX.LT.-1.0D0.OR.(XK-X1)/DX.GT.2.0D0)) GO TO 3
              IF(IX(II).EQ.2.AND.
     &         ((YK-Y1)/DY.LT.-1.0D0.OR.(YK-Y1)/DY.GT.2.0D0)) GO TO 3
              IF(IX(II).EQ.3.AND.
     &         ((ZK-Z1)/DZ.LT.-1.0D0.OR.(ZK-Z1)/DZ.GT.2.0D0)) GO TO 3
              KBC(1)=KS(1,J)+NRECPT(1,KGG)*ICBB
              KBC(2)=KS(2,J)+NRECPT(2,KGG)*ICBB
              KBC(3)=KS(3,J)+NRECPT(3,KGG)*ICBB
              IF(NPOINT.EQ.0) GO TO 20
              DO 21 N=1,NPOINT
                 DO 22 K=1,3
                    IF(KBC(K)*ICCB(N).NE.KKB(K,N)*ICBB) GO TO 21
   22            CONTINUE
                 IF(MRNM(IFKP(N)).EQ.MRNM(I)) GO TO 3
   21         CONTINUE
   20         DO 6 JI=1,NJ
                 JIND(IMM(JI))=I
                 NPOINT=NPOINT+1
                 IF(NPOINT.GT.MAXFKP) THEN
                    WRITE(6,*) ' MAXFKP=',MAXFKP,' NPOINT=',NPOINT
                    STOP
                 END IF
                 IFKP(NPOINT)=IMM(JI)
                 JSKP(NPOINT)=JS(J)
                 KKB(1,NPOINT)=KBC(1)
                 KKB(2,NPOINT)=KBC(2)
                 KKB(3,NPOINT)=KBC(3)
                 ICCB(NPOINT)=ICBB
                 IF(IPR.GE.4) THEN
                    N=NPOINT
                    WRITE(6,601) N,KBC,ICBB,IFKP(N)
     &                          ,KB,IC,MRRM(IFKP(N)),MRNM(IFKP(N))
  601               FORMAT(I3,'(',3I4,')/',I4,I3,
     &                     '(',3I3,')/',I3,2X,A2,I3)
                 END IF
    6         CONTINUE
c             write(6,*) xk,yk,zk
c             write(6,*) npoint,kbc,icbb,i,j,kgg
    3      CONTINUE
    2   CONTINUE
    1 CONTINUE
      CALL REMINV
      RETURN
      END
      SUBROUTINE CMPTBL(JD,II,KBB,ICC,NRT,NH,IPR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)    
      PARAMETER (MAXKPT=3600,MAXEIG=960)
      PARAMETER (MAXFKP=1300,MAXIRP=367)
      PARAMETER (N65=65)
      CHARACTER*2 MRRM
      COMMON/AXE/EAX(N65,MAXEIG,MAXIRP),NKPT(MAXIRP)
     &      ,NLIN(MAXIRP),JRR(MAXIRP),JUD(MAXIRP),NAXEN
      COMMON/AXP/IFKP(MAXFKP),JSKP(MAXFKP),KKB(3,MAXFKP)
     &   ,ICCB(MAXFKP),ICPTBL(12,MAXFKP),INDM(MAXFKP),NPOINT
      COMMON/ENR/KXM(3,MAXKPT),ICM(MAXKPT),IUDM(MAXKPT)
     &   ,MRRM(MAXKPT),MRNM(MAXKPT),MWEIM(MAXKPT),NSTM(MAXKPT)
     &   ,NEIGM(MAXKPT),EIGM(MAXEIG,MAXKPT),NKPINT
C
      DIMENSION KB(3),KBX(3),KBB(3)
      DIMENSION ITCR(12),ND(12),JTR(12),IPA(12)
      DIMENSION NDX(12),JTRX(12),IPAX(12)
      DIMENSION ICP(12,12),ICP0(12,12)
C
      JDOB=JD
      DO 11 I=1,NPOINT
         INDM(I)=0
   11 CONTINUE
      DO 1 I=1,NPOINT
         IF(INDM(I).NE.0) GO TO 1
         KB(1)=KKB(1,I)
         KB(2)=KKB(2,I)
         KB(3)=KKB(3,I)
         IC=ICCB(I)
         KBX(1)=KXM(1,IFKP(I))
         KBX(2)=KXM(2,IFKP(I))
         KBX(3)=KXM(3,IFKP(I))
         ICX=ICM(IFKP(I))
         CALL CORRES(KB,IC,KBX,ICX,JDOB,ITCR,NRR,IND)
         IF(IND.EQ.0) THEN
            write(6,*) ' kb=',kb,iv,' kbx=',kbx,icx
            WRITE(6,*) ' STOP AT CORRES IN CMPTBL'
            STOP
         END IF
         CALL TSIREP(KBX,ICX,JDOB)
         CALL DGTRST(JDO,NRX,NHX,NSTRX,NDX,JTRX,IPAX)
         CALL TSIREP(KB,IC,JDOB)
         CALL DGTRST(JDO,NRR,NHH,NSTR,ND,JTR,IPA)
         DO 4 IR1=1,NRX
            IF(JTR(IR1).EQ.0.AND.IR1.LT.IPAX(IR1).AND.
     &           ITCR(IR1).GT.IPA(ITCR(IR1))) THEN
                 IW1=ITCR(IR1)
                 ITCR(IR1)=ITCR(IPAX(IR1))
                 ITCR(IPAX(IR1))=IW1
            END IF
    4    CONTINUE
C          write(6,*) ' cor',kb,ic,kbx,icx,(itcr(L),L=1,nrr)
         IF(NH.NE.NHH) THEN
            CALL COMPAT(KBB,ICC,NRTT,KB,IC,NRR,JDOB,ICP,INDC)
C            write(6,*) kbb,icc,kb,ic
C            do 1000 i1=1,nrt
C            write(6,*) ' icp',(icp(L,i1),L=1,nrr)
C 1000       continue
            IF(INDC.NE.0) THEN
               WRITE(6,*) 'STOP AT COMPAT IN COMTBL'
               STOP
            END IF
            CALL CMPTRV(KBB,ICC,KB,IC,JDOB,ICP,ICP0)
C             do 1001 i1=1,nrt
C             write(6,*) ' icp0',(icp0(L,i1),L=1,nrr)
C 1001        continue
            DO 2 IIR=1,NRT
               ICPTBL(IIR,I)=ICP0(ITCR(MRNM(IFKP(I))),IIR)
    2       CONTINUE
            INDM(I)=1
            IF(IPR.GE.4) WRITE(6,601) I,INDM(I),(ICPTBL(K,I),K=1,NRT)
  601       FORMAT(14I3)
            IF(I.EQ.NPOINT) GO TO 1
            DO 20 J=I+1,NPOINT
               DO 22 K=1,3
                  IF(KB(K)*ICCB(J).NE.KKB(K,J)*IC) GO TO 1
   22          CONTINUE
               DO 23 IIR=1,NRT
                  ICPTBL(IIR,J)=ICP0(ITCR(MRNM(IFKP(J))),IIR)
   23          CONTINUE
               INDM(J)=1
               IF(IPR.GE.4)
     &            WRITE(6,601) J,INDM(J),(ICPTBL(K,J),K=1,NRT)
   20       CONTINUE
         ELSE
            DO 3 IIR=1,NRT
               ICPTBL(IIR,I)=0
    3       CONTINUE
            ICPTBL(ITCR(MRNM(IFKP(I))),I)=1
            INDM(I)=-1
            IF(IPR.GE.4)
     &            WRITE(6,601) I,INDM(I),(ICPTBL(K,I),K=1,NRT)
            IF(I.EQ.NPOINT) GO TO 1
            DO 30 J=I+1,NPOINT
               DO 32 K=1,3
                  IF(KB(K)*ICCB(J).NE.KKB(K,J)*IC) GO TO 1
   32          CONTINUE
               DO 33 IIR=1,NRT
                  ICPTBL(IIR,J)=0
   33          CONTINUE
               ICPTBL(ITCR(MRNM(IFKP(J))),J)=1
               INDM(J)=-1
               IF(IPR.GE.4)
     &            WRITE(6,601) J,INDM(J),(ICPTBL(K,J),K=1,NRT)
   30       CONTINUE
         END IF
    1 CONTINUE
      RETURN
      END
      SUBROUTINE INTPOR(JD,II,KBB,ICC,NH,IIR,IUD,IPR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NAXMAX=20)                
      PARAMETER (MAXKPT=3600,MAXEIG=960)
      PARAMETER (MAXFKP=1300,MAXIRP=367)
      PARAMETER (N65=65)
      CHARACTER*2 MRRM,MRR
      COMMON/AXE/EAX(N65,MAXEIG,MAXIRP),NKPT(MAXIRP)
     &      ,NLIN(MAXIRP),JRR(MAXIRP),JUD(MAXIRP),NAXEN
      COMMON/AXP/IFKP(MAXFKP),JSKP(MAXFKP),KKB(3,MAXFKP)
     &   ,ICCB(MAXFKP),ICPTBL(12,MAXFKP),INDM(MAXFKP),NPOINT
      COMMON/ISP/IX(NAXMAX),RT(NAXMAX),PS(2,NAXMAX),ICON(NAXMAX),
     &          MK(3,NAXMAX),NAXM
      COMMON/KKEND/KK1M(3,NAXMAX),KK2M(3,NAXMAX)
     &      ,ICC1M(NAXMAX),ICC2M(NAXMAX)
      COMMON/ENR/KXM(3,MAXKPT),ICM(MAXKPT),IUDM(MAXKPT)
     &   ,MRRM(MAXKPT),MRNM(MAXKPT),MWEIM(MAXKPT),NSTM(MAXKPT)
     &   ,NEIGM(MAXKPT),EIGM(MAXEIG,MAXKPT),NKPINT
      COMMON/MRK/NSIT(MAXKPT),MSIT(3,MAXKPT),XIMS(3,MAXKPT)
C
      DIMENSION KBB(3)
      DIMENSION XX(N65),YY(N65),X(MAXFKP),YM(MAXEIG,MAXFKP)
      DIMENSION SM(MAXFKP),JMAX(MAXFKP)
      DIMENSION IXM(MAXFKP),Y(MAXFKP),XXI(MAXFKP)
C
      D6=1.0D-6
      IF(IPR.GE.2) THEN
        KX=KBB(1)
        KY=KBB(2)
        KZ=KBB(3)
        CALL KPNAME(KX,KY,KZ,ICC,MRR,KG)
        WRITE(6,600) MRR,IIR
  600   FORMAT(' INTERPORATION START FOR ',A2,I2)
      END IF
      JRR(NAXEN)=II
      JDOB=JD
      WIDE=DBLE(KK2M(IX(II),II))/ICC2M(II)
     &    -DBLE(KK1M(IX(II),II))/ICC1M(II)
      NXI=0
      DO 1 I=1,NPOINT
         IF(IUD.NE.IUDM(IFKP(I))) GO TO 1
         XI=(DBLE(KKB(IX(II),I))/ICCB(I)
     &      -DBLE(KK1M(IX(II),II))/ICC1M(II))/WIDE
C       write(6,*) i,xi,nxi
         IF(IIR.NE.1) GO TO 13
         IF(XI+D6.LT.0.0) GO TO 13
         IF(XI-D6.GT.1.0D0) GO TO 13
         IF(XI+D6.LT.1.0D0.OR.II.EQ.NAXM.OR.ICON(II).EQ.1) THEN
            NSIT(IFKP(I))=NSIT(IFKP(I))+1
            IF(NSIT(IFKP(I)).GT.3) THEN
               WRITE(6,*) ' NSIT MUST BE LESS THAN 3',NSIT(IFKP(I))
               STOP
            END IF
            MSIT(NSIT(IFKP(I)),IFKP(I))=II
            XIMS(NSIT(IFKP(I)),IFKP(I))=XI
         END IF
C         write(6,*) i,ifkp(i),nsit(ifkp(i)),ii,naxm,xi,'13up'
   13    CONTINUE
         IF(ICPTBL(IIR,I).EQ.0) GO TO 1
         IF(NXI.EQ.0) GO TO 11
         DO 15 NI=1,NXI
            IF(DABS(XI-X(NI)).LT.1.0D-6) GO TO 12
   15    CONTINUE
   11    NXI=NXI+1
         JMAX(NXI)=0
         X(NXI)=XI
         NI=NXI
   12    continue
C   12    IF(INDM(I).EQ.1) THEN
            DO 5 L=1,NEIGM(IFKP(I))
               DO 2 K=1,ICPTBL(IIR,I)
                  IF(JMAX(NI).EQ.0) GO TO 51
                  DO 52 JJ=1,JMAX(NI)
                     J1=JMAX(NI)-JJ+1
                     IF(YM(J1,NI).LT.EIGM(L,IFKP(I))) GO TO 53
                     YM(J1+1,NI)=YM(J1,NI)
   52             CONTINUE
   51             J1=0
   53             YM(J1+1,NI)=EIGM(L,IFKP(I))
                  JMAX(NI)=JMAX(NI)+1
    2          CONTINUE
    5       CONTINUE
C         ELSE
C            DO 54 L=1,NEIGM(IFKP(I))
C               YM(L,NI)=EIGM(L,IFKP(I))
C   54       CONTINUE
C            JMAX(NI)=NEIGM(IFKP(I))
C         END IF
    1 CONTINUE
      DO 6 I=1,NXI
         IXM(I)=I
    6 CONTINUE
      XXI(1)=X(1)
      DO 61 I=2,NXI
         DO 62 J=1,I-1
            JJ=I-J
            IF(XXI(JJ).LT.X(I)) GO TO 63
            XXI(JJ+1)=XXI(JJ)
            IXM(JJ+1)=IXM(JJ)
   62    CONTINUE
         JJ=0
   63    JJ=JJ+1
         XXI(JJ)=X(I)
         IXM(JJ)=I
   61 CONTINUE
C        write(6,*) ' xxi',(xxi(i),i=1,nxi)
      NLIN(NAXEN)=0
      NKPT(NAXEN)=N65
      IF(NXI.LT.3) GO TO 73
      JJMAX=0
      DO 21 I=1,NXI
         IF(JMAX(I).GT.JJMAX) JJMAX=JMAX(I)
   21 CONTINUE
C        write(6,*) 'jjmax=',jjmax
C        write(6,*) (jmax(ixm(i)),i=1,nxi)
      DO 70 J=1,JJMAX 
         NNN=0
         NCO=0
         DO 71 I=1,NXI
            IF(JMAX(IXM(I)).LT.J) GO TO 71
            NNN=NNN+1
            X(NNN)=XXI(I)
            IF(XXI(I)+D6.GE.0.0.AND.XXI(I)-D6.LT.1.0D0) NCO=NCO+1
            Y(NNN)=YM(J,IXM(I))
   71    CONTINUE
         IF(NCO.LT.3) GO TO 70
         IF(IPR.GE.4) THEN
            WRITE(6,*) ' BAND NO=',J
            WRITE(6,666) (X(I),I=1,NNN)
            WRITE(6,666) (Y(I),I=1,NNN)
  666       FORMAT(10F7.4)
         END IF
         C1=4.0*((Y(2)-Y(1))/(X(2)-X(1)))
         AMU1=2.0
         CN=4.0*((Y(NNN)-Y(NNN-1))/(X(NNN)-X(NNN-1)))
         ALMN=2.0
         NLIN(NAXEN)=J
         INIT=0
         ICO=0
         DO 100 I=1,N65
            DD=DBLE(I-1)/(N65-1)
            IF(DD.LT.X(1)) GO TO 100
            IF(DD.GT.X(NNN)) GO TO 101
            IF(INIT.EQ.0) INIT=I
            ICO=ICO+1
            XX(ICO)=DD
  100    CONTINUE
  101    CONTINUE
         CALL S3N(X,Y,SM,XX,YY,NNN,ICO,C1,CN,AMU1,ALMN)
         DO 72 I=1,N65
            EAX(I,J,NAXEN)=99.0
   72    CONTINUE
         DO 75 I=1,ICO
            EAX(INIT+I-1,J,NAXEN)=YY(I)
   75    CONTINUE
   70 CONTINUE
   73 RETURN
      END
      SUBROUTINE S3N(X,Y,SM,XX,YY,N,NN,C1,CN,AMU1,ALMN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(99),Y(99),SM(99),XX(200),YY(200)
      DIMENSION H(99),ALM(99),AMU(99),C(99),P(99),Q(99),U(99)
      N1=N-1
      DO 110 I=2,N
      H(I)=X(I)-X(I-1)
  110 CONTINUE
      DO 120 I=2,N1
      ALM(I)=H(I+1)/(H(I)+H(I+1))
      AMU(I)=1.0-ALM(I)
  120 CONTINUE
      DO 130 I=2,N1
      C(I)=3.0*(ALM(I)*(Y(I)-Y(I-1))/H(I)+AMU(I)*(Y(I+1)-Y(I))/H(I+1))
  130 CONTINUE
      C(1)=C1
      C(N)=CN
      AMU(1)=AMU1
      ALM(N)=ALMN
      P(1)=2.0
      Q(1)=-AMU(1)/P(1)
      U(1)=C(1)/P(1)
      DO 140 K=2,N
      P(K)=ALM(K)*Q(K-1)+2.0
      Q(K)=-AMU(K)/P(K)
      U(K)=(C(K)-ALM(K)*U(K-1))/P(K)
  140 CONTINUE
      SM(N)=U(N)
      DO 150 K=1,N1
      K1=N1-K+1
      SM(K1)=Q(K1)*SM(K1+1)+U(K1)
  150 CONTINUE
      DO 160 I=1,NN
      XXI=XX(I)
      DO 170 K=2,N
      IF(XXI.GT.X(K)) GO TO 170
      J1=K
      GO TO 180
  170 CONTINUE
  180 J=J1-1
      SMJ=SM(J)
      SMJ1=SM(J1)
      YJ=Y(J)
      YJ1=Y(J1)
      HJ1=H(J1)
      XJ1=X(J1)-XXI
      XJ=XXI-X(J)
      HJ2=HJ1*HJ1
      HJ3=HJ2*HJ1
      YY(I)=SMJ*XJ1*XJ1*XJ/HJ2-SMJ1*XJ*XJ*XJ1/HJ2+YJ*XJ1*XJ1*(2.0*XJ+
     &      HJ1)/HJ3+YJ1*XJ*XJ*(2.0*XJ1+HJ1)/HJ3
  160 CONTINUE
      RETURN
      END
