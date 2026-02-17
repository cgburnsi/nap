import numpy as np

def trapp(u, v, npart, x0a, cps, ci3, gamma, beta, k0):
    """
    Numerical integration using the trapezoidal method.

    Parameters:
        u (float): Lower limit of integration.
        v (float): Upper limit of integration.
        npart (int): Number of partitions for integration.
        x0a (float): Parameter for CPXF function.
        cps (float): Parameter for CPXF function.
        ci3 (float): Parameter for RHETF function.
        gamma (float): Parameter for RHETF function.
        beta (float): Parameter for RHETF function.
        k0 (float): Parameter for RHETF function.

    Returns:
        float: The computed Riemann sum (RIESUM).
    """

    # Define RHETF function
    def rhetf(a, b, c, d, e, n):
        return e * a**(1 - n) * b**n * np.exp(c * d * (1 - b / a) / (1 + d * (1 - b / a)))

    # Define FOXI1 function
    def foxi1(x, r):
        return x**2 * r

    # Define CPXF function
    def cpxf(x, y, z):
        return (x - y) / (1 - y) * z

    n = npart - 1
    part = npart
    h = (v - u) / part
    uph = u + h
    sum_ = 0.0

    cpx1 = cpxf(u, x0a, cps)
    cpx2 = cpxf(v, x0a, cps)

    rhet1 = rhetf(ci3, cpx1, gamma, beta, k0, 1)
    rhet2 = rhetf(ci3, cpx2, gamma, beta, k0, 1)

    # Calculate first and last terms of Riemann sum
    trm1 = foxi1(u, rhet1) / 2.0
    trm2 = foxi1(v, rhet2) / 2.0

    # Loop to calculate intermediate terms
    for _ in range(n):
        cpx = cpxf(uph, x0a, cps)
        rhet = rhetf(ci3, cpx, gamma, beta, k0, 1)
        sum_ += foxi1(uph, rhet)
        uph += h

    # Final Riemann sum
    riesum = h * (trm1 + sum_ + trm2)
    return riesum

