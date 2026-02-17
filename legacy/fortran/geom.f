      SUBROUTINE GEOM                                                   GEO   10
C                                                                       GEO   20
C     **************************************************************    GEO   30
C                                                                       GEO   40
C	  THIS SUBROUTINE CALCULATES THE NOZZLE RADIUS AND OUTER NORMAL     GEO   50
C                                                                       GEO   60
C     **************************************************************    GEO   70
C                                                                       GEO   80
      COMMON /AV/ IAV,CAV,NST,SMP,LSS,CTA,XMU,XLA,RKMU,QUT(81,21),QVT(81GEO   90
     1,21),QPT(81,21)                                                   GEO  100
      COMMON /ONESID/ UD(4),VD(4),PD(4),ROD(4)                          GEO  110
      COMMON /SOLUTN/ U(81,21,2),V(81,21,2),P(81,21,2),RO(81,21,2)      GEO  120
      COMMON /CNTRLC/ LMAX,MMAX,NMAX,NPRINT,TCONV,FDT,GAMMA,RGAS,GAM1,GAGEO  130
     1M2,L1,L2,L3,M1,M2,DX,DY,DT,N,N1,N3,NASM,IVEL,ICHAR,N1D,LJET,JFLAG,GEO  140
     2IERR,IUI,IUO,DXR,DYR,LD,MD,LMD1,LMD3,IB,RSTAR,RSTARS,NPLOT,G,PC,TCGEO  150
     3,LC,PLOW,ROLOW                                                    GEO  160
      COMMON /GEMTRYC/ NGEOM,XI,RI,XT,RT,XE,RE,RCI,RCT,ANGI,ANGE,XW(81),GEO  170
     1YW(81),XWI(81),YWI(81),NXNY(81),NWPTS,IINT,IDIF,LT,NDIM           GEO  180
      COMMON /GCB/ NGCB,XICB,RICB,XTCB,RTCB,XECB,RECB,RCICB,RCTCB,ANGICBGEO  190
     2,ANGECB,XCB(81),YCB(81),XCBI(81),YCBI(81),NXNYCB(81),NCBPTS,IINTCBGEO  200
     3,IDIFCB,LECB                                                      GEO  210
      COMMON /BCC/ PT(21),TT(21),THETA(21),PE,MASSE,MASSI,MASST,THRUST,NGEO  220
     1STAG                                                              GEO  230
      REAL MN3,NXNY,MASSI,MASST,NXNYCB,MASSE                            GEO  240
C                                                                       GEO  250
      GO TO (10,30,120,170), NGEOM                                      GEO  260
C                                                                       GEO  270
C     CONSTANT AREA DUCT CASE                                           GEO  280
C                                                                       GEO  290
10    PRINT 230                                                         GEO  300
      IF (IUI,EQ,1) PRINT 250, XI,RI,XE                                 GEO  310
      IF (IUI,EQ,2) PRINT 260, XI,RI,XE                                 GEO  320
      LT=LMAX                                                           GEO  330
      DX=(XE-XI)/(LMAX-1)                                               GEO  340
      XT=XE                                                             GEO  350
      RT=RI                                                             GEO  360
      RE=RI                                                             GEO  370
      DO 20 L=1,LMAX                                                    GEO  380
      YW(L)=RI                                                          GEO  390
      NXNY(L)=0.0                                                       GEO  400
20    CONTINUE                                                          GEO  410
      IF (JFLAG,EQ,0) GO TO 210                                         GEO  420
C                                                                       GEO  430
      XWL=XI+(LJET-2)*DX                                                GEO  440
      IF (IUI,EQ,1) PRINT 370, XWL,LJET,LMAX                            GEO  450
      IF (IUI,EQ,2) PRINT 380, XWL,LJET,LMAX                            GEO  460 
      GO TO 210                                                         GEO  470
C                                                                       GEO  480
C     CIRCULAR-ARC, CONICAL NOZZLE CASE                                 GEO  490
C                                                                       GEO  500
30    PRINT 230                                                         GEO  510
      IF (RCI,EQ,0.0,OR,RCT,EQ,0.0) GO TO 200                           GEO  520 
      ANI=ANGI*3.141593/180.0                                           GEO  530
      ANE=ANGE*3.141593/180.0                                           GEO  540
      XTAN=XI+RCI*SIN(ANI)                                              GEO  550
      RTAN=RI+RCI*(COS(ANI)-1.0)                                        GEO  560
      RT1=RT-RCT*(COS(ANI)-1.0)                                         GEO  570
      XT1=XTAN+(RTAN-RT1)/TAN(ANI)                                      GEO  580
      IF (XT1,GE,XTAN) GO TO 40                                         GEO  590
      XT1=XTAN                                                          GEO  600
      RT1=RTAN                                                          GEO  610
40    XT=XT1+RCT*SIN(ANI)                                               GEO  620
      XT2=XT+RCT*SIN(ANE)                                               GEO  630
      RT2=RT+RCT*(1.0-COS(ANE))                                         GEO  640
      RE=RT2+(XE-XT2)*TAN(ANE)                                          GEO  650
      LT=1                                                              GEO  660
      DX=(XE-XI)/(LMAX-1)                                               GEO  670
      IF (IUI,EQ,1) PRINT 270, XI,RI,RT,XE,RCI,RCT,ANGI,ANGE,XT,RE      GEO  680
      IF (IUI,EQ,2) PRINT 280, XI,RI,RT,XE,RCI,RCT,ANGI,ANGE,XT,RE      GEO  690
      DO 110 L=1,LMAX                                                   GEO  700
      X=XI+(L-1)*DX                                                     GEO  710
      IF (X,GE,XI,AND,X,LE,XTAN) GO TO 50                               GEO  720
      IF (X,GT,XTAN,AND,X,LE,XT1) GO TO 60                              GEO  730
      IF (X,GT,XT1,AND,X,LE,XT) GO TO 70                                GEO  740
      IF (X,GT,XT,AND,X,LE,XT2) GO TO 80                                GEO  750
      IF (X,GT,XT2,AND,X,LE,XE) GO TO 90                                GEO  760
C                                                                       GEO  770
50    YW(L)=RI+RC*(COS(ASIN((X-XI)/RCI))-1.0)                           GEO  780
      NXNY(L)=(XI-XI)/(YW(L)-RI+RCI)                                    GEO  790
      GO TO 100                                                         GEO  800
C                                                                       GEO  810
60    YW(L)=RT1+(XT1-X)*TAN(ANI)                                        GEO  820
      NXNY(L)=TAN(ANI)                                                  GEO  830
      GO TO 100                                                         GEO  840
C                                                                       GEO  850

      
      
