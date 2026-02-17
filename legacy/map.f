      SUBROUTINE MAP(IP,L,M,AL,BE,DE,LD1,AL1,BE1,DE1)                   MAP   10
C                                                                       MAP   20
C     **************************************************************    MAP   30
C                                                                       MAP   40
C	  THIS SUBROUTINE CALCULATES THE MAPPING FUNCTIONS                  MAP   50
C                                                                       MAP   60
C     **************************************************************    MAP   70
C                                                                       MAP   80
      COMMON /AV/ IAV,CAV,NST,SMP,LSS,CTA,XMU,XLA,RKMU,QUT(81,21),QVT(81MAP   90
     1,21),QPT(81,21)                                                   MAP  100
      COMMON /ONESID/ UD(4),VD(4),PD(4),ROD(4)                          MAP  110
      COMMON /SOLUTN/ U(81,21,2),V(81,21,2),P(81,21,2),RO(81,21,2)      MAP  120
      COMMON /CNTRLC/ LMAX,MMAX,NMAX,NPRINT,TCONV,FTD,GAMMA,RGAS,GAM1,GAMAP  130
     1M2,L1,L2,L3,M1,M2,DX,DY,DT,N,N1,N3,NASM,IVEL,ICHAR,N1D,LJET,JFLAG,MAP  140
     2IERR,IUI,IUO,DXR,DYR,LD,MD,LMND,LMD3,IB,RSTAR,RSTARS,NPLOT,G,PC,TCMAP  150
     3,LC,PLOW,ROLOW                                                    MAP  160
      COMMON /GEMTRYC/ NGEOM,XI,RI,XT,RT,XE,RE,RCI,RCT,ANGI,ANGE,XW(81),MAP  170
     1YW(81),XWI(81),YWI(81),NXNY(81),NWPTS,IINT,IDIF,LT,NDIM           MAP  180
      COMMON /GCB/ NGCB,XICB,RICB,XTCB,RTCB,XECB,RECB,RCICB,RCTCB,ANGICBMAP  190
     2,ANGECB,XCB(81),YCB(81),XCBI(81),YCBI(81),NXNYCB(81),NCBPTS,IINTCBMAP  200
     3,IDIFCB,LECB                                                      MAP  210
      COMMON /BCC/ PT(21),TT(21),THETA(21),PE,MASSE,MASSI,MASST,THRUST,NMAP  220
     1STAG                                                              MAP  230
      REAL MN3,NXNY,MASSI,MASST,NXNYCB,MASSE                            MAP  240
C                                                                       MAP  250
      BE=1.0/(YW(L)-YCB(L))                                             MAP  260
      IF (IP,EQ,0) RETURN                                               MAP  270
      Y=(M-1)*DY                                                        MAP  280
      AL=BE*(NXNYCB(L)+Y*(NXNY(L)-NXNYCB(L)))                           MAP  290
      DE=-BE*Y*XWI(L)                                                   MAP  300
      IF (IP,EQ,1) RETURN                                               MAP  310
      BE1=1.0/(YW(LD1)-YCB(LD1))                                        MAP  320
      AL1=BE1*(NXNYCB(LD1)+Y*(NXNY(LD1)-NXNYCB(LD1)))                   MAP  330
      DE1=-BE1*Y*XWI(LD1)                                               MAP  340
      RETURN                                                            MAP  350
      END                                                               MAP  360