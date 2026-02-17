
MODULE DataStorage
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Double precision kind
  INTEGER, PARAMETER :: max1D = 81, max2D = 21

  REAL(dp) :: IAV, CAV, NST, SMP, LSS, CTA, XMU, XLA, RKMU
  REAL(dp), DIMENSION(max1D, max2D) :: QUT, QVT, QPT
  REAL(dp), DIMENSION(4) :: UD, VD, PD, ROD
  REAL(dp), DIMENSION(max1D, max2D, 2) :: U, V, P, RO
  INTEGER :: LMAX, MMAX, NMAX, NPRINT, NASM, IVEL, ICHAR, N1D, LJET, JFLAG
  REAL(dp) :: TCONV, FDT, GAMMA, RGAS, GAM1, GAM2, DX, DY, DT, DXR, DYR
  REAL(dp) :: RSTAR, RSTARS, G, PC, TC, PLOW, ROLOW
  REAL(dp), DIMENSION(max1D) :: XW, YW, XWI, YWI, NXNY
  REAL(dp), DIMENSION(max1D) :: XCB, YCB, XCBI, YCBI, NXNYCB
  INTEGER :: NGCB, LT, NDIM
  REAL(dp), DIMENSION(max2D) :: PT, TT, THETA
  REAL(dp) :: PE, MASSE, MASSI, MASST, THRUST, NSTAG

END MODULE DataStorage

SUBROUTINE OneDim
  USE DataStorage
  IMPLICIT NONE

  INTEGER :: L, M, ITER, NXCK
  REAL(dp) :: MN3, GRGAS, ACOEF, BCOEF, CCOEF, X, RAD, ARATIO, OMN3
  REAL(dp) :: ABM, ABMC, FM, FPM, DEM, DEMP, TEMP, Q, DN, DNS, SIGN

  ! Initialize variables
  MN3 = 0.01
  IF (N1D == -1 .OR. N1D > 2) MN3 = 2.0

  GRGAS = 1.0 / (RGAS * G)
  NXCK = 0
  ACOEF = 2.0 / (GAMMA + 1.0)
  BCOEF = (GAMMA - 1.0) / (GAMMA + 1.0)
  CCOEF = (GAMMA + 1.0) / (2.0 * (GAMMA - 1.0))

  ! Main computation loop
  DO L = 1, LMAX
    IF ((L == 1 .AND. N1D == -1) .OR. (L == 1 .AND. N1D > 2)) CYCLE

    X = XI + DX * (L - 1)
    IF (N1D < 0) CYCLE

    IF (NGCB /= 0) THEN
      RAD = YW(L) - YCB(L)
      ARATIO = RAD / RSTAR
    ELSE
      RSTAR = RT
      RSTARS = RT ** 2
    END IF

    ! Newton-Raphson iteration for MN3
    DO ITER = 1, 20
      ABM = ACOEF + BCOEF * MN3 ** 2
      ABMC = ABM ** CCOEF
      FM = ABMC / MN3 - ARATIO
      FPM = ABMC * (2.0 * BCOEF * CCOEF / ABM - 1.0 / MN3 ** 2)
      OMN3 = MN3
      MN3 = OMN3 - FM / FPM

      IF (ABS(MN3 - OMN3) / OMN3 <= 0.0005) EXIT
    END DO

    ! Fill in 2-D arrays
    DEM = 1.0 + GAM2 * MN3 ** 2
    DEMP = DEM ** GAM1
    DN = (NXNY(L) - NXNYCB(L)) / M1

    DO M = 1, MMAX
      P(L, M, 1) = PT(M) / DEMP
      TEMP = TT(M) / DEM
      RO(L, M, 1) = P(L, M, 1) * GRGAS / TEMP
      Q = MN3 * SQRT(GAMMA * P(L, M, 1) / RO(L, M, 1))
      DNS = DN ** 2

      IF (DNS == 0.0_DP) THEN
        U(L, M, 1) = Q
        V(L, M, 1) = 0.0_DP
      ELSE
        SIGN = -1.0
        IF (DN > 0.0) SIGN = 1.0
        U(L, M, 1) = Q / SQRT(1.0 + DNS)
        V(L, M, 1) = SIGN * Q / SQRT(1.0 + 1.0 / DNS)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE OneDim
