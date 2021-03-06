      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NKMAX=40000)
      DIMENSION IOKP(NKMAX,3),IHKP(NKMAX,3)
      NX=36
      NY=36
      NZ=36
      IC=72
      NIX=6
      NIY=2
      NIZ=1
      WRITE(6,*) " NX,NY,NZ,IC=",NX,NY,NZ,IC
      READ(5,*,END=5) NX,NY,NZ,IC
    5 WRITE(6,*) " NIX,NIY,NIZ=",NIX,NIY,NIZ 
      READ(5,*,END=6)NIX,NIY,NIZ
C
    6 NI=1
      DO 10 KX=-NX,NX,NIX
      DO 10 KY=-NY,NY,NIY
      DO 10 KZ=-NZ,NZ,NIZ
      IOKP(NI,1)=KX
      IOKP(NI,2)=KY
      IOKP(NI,3)=KZ
      NI=NI+1
   10 CONTINUE
      NIMAX=NI-1
C
      DO 30 NI=1,NIMAX
c     IHKP(NI,1)=(IOKP(NI,1)-IOKP(NI,2))/2
      IHKP(NI,1)=IOKP(NI,1)-IOKP(NI,2)
      IHKP(NI,2)=IOKP(NI,2)
      IHKP(NI,3)=IOKP(NI,3)
   30 CONTINUE
C
      DO 20 NI=1,NIMAX
      WRITE(6,100) (IOKP(NI,I),I=1,3),IC,(IHKP(NI,I),I=1,3),IC
      WRITE(7,105) (IHKP(NI,I),I=1,3),IC,(IOKP(NI,I),I=1,3),IC
  100 FORMAT(2(2X,"(",3I4,")/",I4))
  105 FORMAT(3I4,I4,1x,3I4,I4)
   20 CONTINUE
      WRITE(6,110) NIMAX
  110 FORMAT(2X,I5," POINTS REGISTERED")
C
      STOP
      END
