MODULE inter_module
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_L = 81, MAX_M = 21
  REAL :: GAMMA, GAM1, GAM2, PLOW, PC, RLOW, G, DXR, DYR, DT
  REAL :: U(MAX_L, MAX_M, 2), V(MAX_L, MAX_M, 2), P(MAX_L, MAX_M, 2), RO(MAX_L, MAX_M, 2)
  REAL :: QUT(MAX_L, MAX_M), QVT(MAX_L, MAX_M), QPT(MAX_L, MAX_M)
  REAL :: AL, BE, DE, AL1, BE1, DE1, LMD1, LMD2, LMD3, DX, DY, FDT
  INTEGER :: LMAX, MMAX, NMAX, NPRINT, N1, N2, N3, IVEL, ICHAR, N1D, LJET, JFLAG
CONTAINS

  SUBROUTINE INTER
    IMPLICIT NONE
    INTEGER :: L, M, MDUM, LMN1, LMN3, L1MN1, LM1N1
    REAL :: UB, VB, PB, ROB, ASB, DUDX, DPDX, DRODX, DVDY, URHS, RORHS, PRHS
    REAL :: ATERM, UVB, DUDY, DVDX, DUDX, DPDY, DRODY, VRHS
    LOGICAL :: IAV

    ATERM = 0.0
    IF (ICHAR .EQ. 2) RETURN

    MDUM = 1
    IF (NGCB .NE. 0) MDUM = 2

    ! Compute the tentative solution at T + DT
    DO L = 2, LMAX
      DO M = MDUM, M1
        CALL MAP(1, L, M, AL, BE, DE, LD1, AL1, BE1, DE1)
        LMD2 = LD * (M - 1)
        LMN1 = L + LMD2 + LMD1
        LMN3 = L + LMD2 + LMD3
        L1MN1 = L - 1 + LMD2 + LMD1
        LM1N1 = L + LD * (M - 2) + LMD1
        UB = U(LMN1)
        VB = V(LMN1)
        PB = P(LMN1)
        ROB = RO(LMN1)
        ASB = GAMMA * PB / ROB
        IF (M .NE. 1) CYCLE

        DUDX = (UB - U(L1MN1)) * DXR
        DPDX = (PB - P(L1MN1)) * DXR
        DRODX = (ROB - RO(L1MN1)) * DXR
        DVDY = (4.0 * V(L, 2, N1) - V(L, 3, N1)) * 0.5 * DYR
        V(LMN3) = 0.0

        URHS = -UB * DUDX - DPDX / ROB
        RORHS = -UB * DRODX - ROB * DUDX - (1 + NDIM) * ROB * BE * DVDY
        PRHS = -UB * DPDX + ASB * (RORHS + UB * DRODX)
        
        ! Updating solution
        U(LMN3) = U(LMN1) + URHS * DT
        P(LMN3) = P(LMN1) + PRHS * DT
        RO(LMN3) = RO(LMN1) + RORHS * DT
        IF (P(LMN3) .LE. 0.0) P(LMN3) = PLOW * PC
        IF (RO(LMN3) .LE. 0.0) RO(LMN3) = RLOW / G

      END DO
    END DO
    RETURN

  END SUBROUTINE INTER

END MODULE inter_module
