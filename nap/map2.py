
def MAP(IP, L, M, XWI, YW, NXNY, DY, LD1):
    BE = 1.0 / YW[L]
    if IP == 0:
        return BE, None, None, None, None, None, None
    Y = (M - 1) * DY
    AL = BE * Y * NXNY[L]
    DE = -BE * Y * XWI[L]
    if IP == 1:
        return BE, Y, AL, DE, None, None, None
    BE1 = 1.0 / YW[LD1]
    AL1 = BE1 * Y * NXNY[LD1]
    DE1 = -BE1 * Y * XWI[LD1]
    return BE, Y, AL, DE, BE1, AL1, DE1
