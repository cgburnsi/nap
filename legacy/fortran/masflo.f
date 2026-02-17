	SUBROUTINE MASFLO(ISURF)                                          MAS   10
C                                                                     MAS   20
C	**************************************************************    MAS   30
C                                                                     MAS   40
C	THIS SUBROUTINE CALCULTES THE INITIAL-DATA OR SOLUTION SURFACE    MAS   50
C   MASS FLOW AND THRUST                                              MAS   60
C                                                                     MAS   70
C	**************************************************************    MAS   80
C                                                                     MAS   90
	COMMON /AV/ IAV,CAV,NST,SMP,LSS,CTA,XMU,XLA,RKMU,QUT(81,21),QVT(81MAS  100
   1,21),QPT(81,21)                                                   MAS  110
	COMMON /ONESID/ UD(4),VD(4),PD(4),ROD(4)                          MAS  120
	COMMON /SOLUTN/ U(8l,21,2),V(8l,21,2),P(8l,21,2),RO(81,2l,2)      MAS  130
	COMMON /CNTRLC/ LMAX,MMAX,NMAX,NPRINT,TCONV,FDT,GAMMA,RGAS,GAM1,GAMAS  140
   1M2,L1,L2,L3,Ml,M2,DX,DY,DT,N,N1,N3,NASM,IVEL,ICHAR,N1D,LJET,JFLAG,MAS  150
   2IERR,IUI,IUO,DXR,DYR,LD,MD,LMD1,LMD3,IB,RSTAR,RSTARS,NPLOT,G,PC,TCMAS  160
   3,LC,PLOW,ROLOW                                                    MAS  170
	COMMON /GEMTRYC/ NGEOM,XI,RI,XT,RT,XE,RE,RCI,RCT,ANGI,ANGE,XW(8l),MAS  180
   1YW(8l),XWI(8l),YWI(8l),NXNY(8l),NWPTS,IINT,IDIF,LT,NDIM           MAS  190
	COMMON /GCB/ NGCB,XICB,RICB,XTCB,RTCB,XECB,RECB,RCICB,RCTCB,ANGICBMAS  200
   1,ANGECB,XCB(8l),YCB(8l),XCBI(81),YCBI(8l),NXNYCB(8l),NCBPTS,IINTCBMAS  210
   2,IDIFCB,LECB                                                      MAS  220
	COMMON /BCC/ PT(21),TT(21),THETA(21),PE,MASSE,MASSI,MASST,THRUST,NMAS  230
   3STAG                                                              MAS  240
	REAL MN3,NXNY,MASSI,MASST,NXNYCB,MASSE                            MAS  250
C                                                                     MAS  260
    LC2=LC*LC                                                         MAS  270
    LDUM=LMAX-1                                                       MAS  280
    IF (LT,EQ,LMAX) LT=LMAX-1                                         MAS  290
    IF (JFLAG,NE,0) LDUM=LJET-1                                       MAS  300
    IF (ISURF,EQ,1,OR,N1D,EQ,0) GO TO 30                              MAS  310
C                                                                     MAS  320
C   CALCULATE THE MASS FLOW AND THRUST FOR THE 1-D INITIAL-DATA       MAS  330
C   SURFACE                                                           MAS  340
C                                                                     MAS  350
    IF (NDIM,EQ,1) GO TO 10                                           MAS  360
    AREAI=(YW(1)-YCB(1))/LC2                                          MAS  370
    AREAT=(YW(LT)-YCB(LT))/LC2                                        MAS  380
    AREAE=(YW(LDUM)-YCB(LDUM))/LC2                                    MAS  390
    GO TO 20                                                          MAS  400
10  AREAI=3.141593*(YW(1)**2-YCB(1)**2)/LC2                           MAS  410
    AREAT=3.141593*(YW(LT)**2-YCB(LT)**2)/LC2                         MAS  420
    AREAE=3.141593*(YW(LDUM)**2-YCB(LDUM)**2)/LC2                     MAS  430
    GO TO 20                                                          MAS  440
20  VMI=SQRT(U(1,1,1)**2+V(1,1,1)**2)                                 MAS  450
    VMT=SQRT(U(LT,1,1)**2+V(LT,1,1)**2)                               MAS  460
    VME=SQRT(U(LDUM,1,1)**2+V(LDUM,1,1)**2)                           MAS  470
    MASSI=RO(1,1,1)*VMI*AREAI*G                                       MAS  480
    MASST=RO(LT,1,1)*VMT*AREAT*G                                      MAS  490
    MASSE=RO(LDUM,1,1)*VME*AREAE*G                                    MAS  500
    THRUST=RO(LDUM,1,1)*U(LDUM,1,1)**2*AREAE                          MAS  510
C                                                                     MAS  520
C   CALCULATE THE MASS FLOW AND THRUST FOR THE 2-D INITIAL-DATA       MAS  530
C   SURFACE                                                           MAS  540
C                                                                     MAS  550
30  MASSI=0.0                                                         MAS  560
    MASST=0.0                                                         MAS  570
    MASSE=0.0                                                         MAS  580
    THRUST=0.0                                                        MAS  590
    DYI=DY*(YW(1)-YCB(1))                                             MAS  600
    DYT=DY*(YW(LT)-YCB(LT))                                           MAS  610
    DYE=DY*(YW(LDUM)-YCB(LDUM))                                       MAS  620
    ND=1                                                              MAS  630
    IF (ISURF,EQ,1) ND=N3                                             MAS  640
    DO 60 M=1,M1                                                      MAS  650
    RADI=(M-1)*DYI+YCB(1)                                             MAS  660
    RADT=(M-1)*DYT+YCB(LT)                                            MAS  670
    RADE=(M-1)*DYE+YCB(LDUM)                                          MAS  680
    IF (NDIM,EQ,1) GO TO 40                                           MAS  690
    AREAI=DYI/LC2                                                     MAS  700
    AREAT=DYT/LC2                                                     MAS  710
    AREAE=DYE/LC2                                                     MAS  720
    GO TO 50                                                          MAS  730
40  AREAI=3.141593*((RADI+DYI)**2-RADI**2)/LC2                        MAS  740
    AREAT=3.141593*((RADT+DYT)**2-RADT**2)/LC2                        MAS  750
    AREAE=3.141593*((RADE+DYE)**2-RADE**2)/LC2                        MAS  760
    ROUI=(RO(1,M,ND)*U(1,M,ND)+RO(1,M+1,ND)*U(1,M+1,ND))*0.5          MAS  770
    ROUT=(RO(LT,M,ND)*U(LT,M,ND)+RO(LT,M+1,ND)*U(LT,M+1,ND))*0.5      MAS  780
    ROUE=(RO(LDUM,M,ND)*U(LDUM,M,ND)+RO(LDUM,M+1,ND)*U(LDUM,M+1,ND))*0MAS  790
   1.5                                                                MAS  800
    ROUE2=(RO(LDUM,M,ND)*U(LDUM,M,ND)**2+RO(LDUM,M+1,ND)*U(LDUM,M+1,NDMAS  810
   1)**2)*0.5                                                         MAS  820
    MASSI=MASSI+ROUI*AREAI*G                                          MAS  830
    MASST=MASST+ROUT*AREAT*G                                          MAS  840
    MASSE=MASSE+ROUE*AREAE*G                                          MAS  850
    THRUST=THRUST+ROU2E*AREAE                                         MAS  860
60  CONTINUE                                                          MAS  870
    RETURN                                                            MAS  880
    END                                                               MAS  890
