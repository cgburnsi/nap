"""
main.py  —  1D Finite Volume Method Solver
==========================================
This version keeps the existing 1D Euler numerics but stores solver state in
structured NumPy arrays using the shared dtypes in `src/dtypes/solver_dtypes.py`.
"""

import numpy as np
import matplotlib.pyplot as plt

from src.dtypes.solver_dtypes import CONSERVATIVE_DTYPE, PRIMITIVE_DTYPE


# =============================================================================
# 1. CONSTANTS AND PARAMETERS
# =============================================================================

GAMMA = 1.4
R_GAS = 287.058

X_START = 0.0
X_END = 1.0
NCELLS = 400

W_LEFT = (1.000, 0.0, 1.0)
W_RIGHT = (0.125, 0.0, 0.1)
DIAPHRAGM = 0.5

CFL = 0.5
T_END = 0.20
SCHEME = "muscl"
LIMITER = "van_leer"

BC_LEFT = "transmissive"
BC_RIGHT = "transmissive"

PLOT_ANALYTICAL = True
PRINT_INTERVAL = 100

LIMITER_MINMOD = 0
LIMITER_VAN_LEER = 1
LIMITER_SUPERBEE = 2


# =============================================================================
# 2. STATE HELPERS
# =============================================================================

def zeros_conservative(size):
    """Allocate conservative states with the shared structured dtype."""
    return np.zeros(size, dtype=CONSERVATIVE_DTYPE)


def zeros_primitive(size):
    """Allocate primitive states with the shared structured dtype."""
    return np.zeros(size, dtype=PRIMITIVE_DTYPE)


def copy_conservative(U):
    """Return a writable copy of a conserved scalar or array."""
    return np.array(U, dtype=CONSERVATIVE_DTYPE, copy=True)


def primitive_state(rho, u, p, v=0.0):
    """Construct one primitive state compatible with the shared dtype."""
    W = zeros_primitive(())
    W["rho"] = rho
    W["u"] = u
    W["v"] = v
    W["p"] = p
    return populate_primitive_derived(W)


def conservative_state(rho, rhou, rhoE, rhov=0.0):
    """Construct one conserved state compatible with the shared dtype."""
    U = zeros_conservative(())
    U["rho"] = rho
    U["rhou"] = rhou
    U["rhov"] = rhov
    U["rhoE"] = rhoE
    return U


def primitive_component_matrix(W):
    """Dense [rho, u, p] matrix for limiter algebra and plotting."""
    return np.column_stack([W["rho"], W["u"], W["p"]])


def populate_primitive_derived(W):
    """Fill derived primitive fields only where the state is initialized."""
    mask = W["rho"] > 0.0
    W["a"] = 0.0
    W["T"] = 0.0
    if np.any(mask):
        W["a"][mask] = sound_speed(W["rho"][mask], W["p"][mask])
        W["T"][mask] = temperature(W["rho"][mask], W["p"][mask])
    return W


def conservative_add_scaled(base, scale, delta):
    """Field-wise base + scale * delta for structured conservative arrays."""
    out = copy_conservative(base)
    for field in CONSERVATIVE_DTYPE.names:
        out[field] = base[field] + scale * delta[field]
    return out


def conservative_linear_combo(shape, terms):
    """Field-wise linear combination for structured conservative arrays."""
    out = zeros_conservative(shape)
    for weight, state in terms:
        for field in CONSERVATIVE_DTYPE.names:
            out[field] += weight * state[field]
    return out


# =============================================================================
# 3. EQUATION OF STATE
# =============================================================================

def pressure_from_internal_energy(rho, e):
    """p = (gamma - 1) * rho * e"""
    return (GAMMA - 1.0) * rho * e


def sound_speed(rho, p):
    """a = sqrt(gamma * p / rho)"""
    return np.sqrt(GAMMA * p / rho)


def temperature(rho, p):
    """T = p / (rho * R)"""
    return p / (rho * R_GAS)


def total_energy(rho, u, p, v=0.0):
    """
    Total energy per unit volume.
    E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2)
    """
    return p / (GAMMA - 1.0) + 0.5 * rho * (u**2 + v**2)


# =============================================================================
# 4. CONSERVED <-> PRIMITIVE VARIABLE CONVERSIONS
# =============================================================================

def prim_to_cons(W):
    """Convert primitive structured states into conserved structured states."""
    W_in = np.asarray(W, dtype=PRIMITIVE_DTYPE)
    rho = W_in["rho"]
    u = W_in["u"]
    v = W_in["v"]
    p = W_in["p"]

    U = zeros_conservative(W_in.shape)
    U["rho"] = rho
    U["rhou"] = rho * u
    U["rhov"] = rho * v
    U["rhoE"] = total_energy(rho, u, p, v)
    return U


def cons_to_prim(U):
    """Convert conserved structured states into primitive structured states."""
    U_in = np.asarray(U, dtype=CONSERVATIVE_DTYPE)
    rho = U_in["rho"]
    rhou = U_in["rhou"]
    rhov = U_in["rhov"]
    rhoE = U_in["rhoE"]

    u = rhou / rho
    v = rhov / rho
    kinetic = 0.5 * (rhou**2 + rhov**2) / rho
    p = (GAMMA - 1.0) * (rhoE - kinetic)

    W = zeros_primitive(U_in.shape)
    W["rho"] = rho
    W["u"] = u
    W["v"] = v
    W["p"] = p
    return populate_primitive_derived(W)


def check_physical(U, location=""):
    """Density and pressure must remain positive."""
    W = cons_to_prim(U)
    if np.any(W["rho"] <= 0.0) or np.any(W["p"] <= 0.0):
        raise ValueError(
            f"Non-physical state detected {location}: "
            f"min(rho)={W['rho'].min():.3e}, min(p)={W['p'].min():.3e}"
        )


# =============================================================================
# 5. PHYSICAL FLUX
# =============================================================================

def euler_flux(U):
    """1D Euler physical flux stored in the conservative structured dtype."""
    U_in = np.asarray(U, dtype=CONSERVATIVE_DTYPE)
    rho = U_in["rho"]
    rhou = U_in["rhou"]
    rhoE = U_in["rhoE"]
    u = rhou / rho
    p = (GAMMA - 1.0) * (rhoE - 0.5 * rhou**2 / rho)

    F = zeros_conservative(U_in.shape)
    F["rho"] = rhou
    F["rhou"] = rhou * u + p
    F["rhov"] = 0.0
    F["rhoE"] = u * (rhoE + p)
    return F


def max_wave_speed(U):
    """Maximum signal speed max(|u| + a) across the domain."""
    W = cons_to_prim(U)
    return float(np.max(np.abs(W["u"]) + W["a"]))


# =============================================================================
# 6. HLLC RIEMANN SOLVER
# =============================================================================

def hllc_flux(U_L, U_R):
    """HLLC numerical flux for a single interface."""
    W_L = cons_to_prim(U_L)
    W_R = cons_to_prim(U_R)

    rho_L = W_L["rho"]
    u_L = W_L["u"]
    p_L = W_L["p"]
    rho_R = W_R["rho"]
    u_R = W_R["u"]
    p_R = W_R["p"]
    a_L = W_L["a"]
    a_R = W_R["a"]
    E_L = U_L["rhoE"]
    E_R = U_R["rhoE"]

    rho_bar = 0.5 * (rho_L + rho_R)
    a_bar = 0.5 * (a_L + a_R)
    p_star = max(0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * rho_bar * a_bar, 0.0)

    q_L = 1.0 if p_star <= p_L else np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (p_star / p_L - 1))
    q_R = 1.0 if p_star <= p_R else np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (p_star / p_R - 1))

    S_L = u_L - a_L * q_L
    S_R = u_R + a_R * q_R

    if S_L >= 0.0:
        return euler_flux(U_L)
    if S_R <= 0.0:
        return euler_flux(U_R)

    num = p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)
    denom = rho_L * (S_L - u_L) - rho_R * (S_R - u_R)
    S_star = num / (denom + 1e-300)

    def star_state(rho, u, E, p, S):
        factor = rho * (S - u) / (S - S_star + 1e-300)
        e_star = E / rho + (S_star - u) * (S_star + p / (rho * (S - u + 1e-300)))
        return conservative_state(
            rho=factor,
            rhou=factor * S_star,
            rhoE=factor * e_star,
        )

    if S_star >= 0.0:
        U_star_L = star_state(rho_L, u_L, E_L, p_L, S_L)
        F_L = euler_flux(U_L)
        return conservative_add_scaled(F_L, S_L, conservative_add_scaled(U_star_L, -1.0, U_L))

    U_star_R = star_state(rho_R, u_R, E_R, p_R, S_R)
    F_R = euler_flux(U_R)
    return conservative_add_scaled(F_R, S_R, conservative_add_scaled(U_star_R, -1.0, U_R))


def hllc_flux_batch(U_L, U_R):
    """Vectorized HLLC across all faces using structured arrays."""
    rho_L = U_L["rho"]
    rhou_L = U_L["rhou"]
    E_L = U_L["rhoE"]
    rho_R = U_R["rho"]
    rhou_R = U_R["rhou"]
    E_R = U_R["rhoE"]

    u_L = rhou_L / rho_L
    u_R = rhou_R / rho_R
    p_L = (GAMMA - 1.0) * (E_L - 0.5 * rhou_L**2 / rho_L)
    p_R = (GAMMA - 1.0) * (E_R - 0.5 * rhou_R**2 / rho_R)
    a_L = np.sqrt(GAMMA * p_L / rho_L)
    a_R = np.sqrt(GAMMA * p_R / rho_R)

    rho_bar = 0.5 * (rho_L + rho_R)
    a_bar = 0.5 * (a_L + a_R)
    p_star = np.maximum(0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * rho_bar * a_bar, 0.0)

    q_L = np.where(
        p_star <= p_L,
        1.0,
        np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (p_star / p_L - 1)),
    )
    q_R = np.where(
        p_star <= p_R,
        1.0,
        np.sqrt(1 + (GAMMA + 1) / (2 * GAMMA) * (p_star / p_R - 1)),
    )

    S_L = u_L - a_L * q_L
    S_R = u_R + a_R * q_R

    num = p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)
    denom = rho_L * (S_L - u_L) - rho_R * (S_R - u_R)
    S_star = num / (denom + 1e-300)

    def star_state(rho, u, E, p, S):
        factor = rho * (S - u) / (S - S_star + 1e-300)
        e_star = E / rho + (S_star - u) * (S_star + p / (rho * (S - u + 1e-300)))
        U_star = zeros_conservative(rho.shape)
        U_star["rho"] = factor
        U_star["rhou"] = factor * S_star
        U_star["rhoE"] = factor * e_star
        return U_star

    U_starL = star_state(rho_L, u_L, E_L, p_L, S_L)
    U_starR = star_state(rho_R, u_R, E_R, p_R, S_R)

    F_L = euler_flux(U_L)
    F_R = euler_flux(U_R)
    F_starL = zeros_conservative(rho_L.shape)
    F_starR = zeros_conservative(rho_L.shape)

    for field in CONSERVATIVE_DTYPE.names:
        F_starL[field] = F_L[field] + S_L * (U_starL[field] - U_L[field])
        F_starR[field] = F_R[field] + S_R * (U_starR[field] - U_R[field])

    F = zeros_conservative(rho_L.shape)
    mask_left = S_L >= 0.0
    mask_star_left = (~mask_left) & (S_star >= 0.0)
    mask_star_right = (~mask_left) & (~mask_star_left) & (S_R >= 0.0)
    mask_right = ~(mask_left | mask_star_left | mask_star_right)

    for field in CONSERVATIVE_DTYPE.names:
        F[field][mask_left] = F_L[field][mask_left]
        F[field][mask_star_left] = F_starL[field][mask_star_left]
        F[field][mask_star_right] = F_starR[field][mask_star_right]
        F[field][mask_right] = F_R[field][mask_right]

    return F


# =============================================================================
# 7. SLOPE LIMITERS
# =============================================================================

def minmod(r):
    return np.maximum(0.0, np.minimum(1.0, r))


def van_leer(r):
    return (r + np.abs(r)) / (1.0 + np.abs(r) + 1e-300)


def superbee(r):
    return np.maximum(0.0, np.maximum(np.minimum(2.0 * r, 1.0), np.minimum(r, 2.0)))


def limiter_name_to_code(limiter_name):
    if limiter_name == "minmod":
        return LIMITER_MINMOD
    if limiter_name == "van_leer":
        return LIMITER_VAN_LEER
    if limiter_name == "superbee":
        return LIMITER_SUPERBEE
    raise ValueError(f"Unknown limiter '{limiter_name}'")


def apply_limiter(limiter_code, r):
    if limiter_code == LIMITER_MINMOD:
        return minmod(r)
    if limiter_code == LIMITER_VAN_LEER:
        return van_leer(r)
    if limiter_code == LIMITER_SUPERBEE:
        return superbee(r)
    raise ValueError(f"Unknown limiter code {limiter_code}")


# =============================================================================
# 8. MUSCL RECONSTRUCTION
# =============================================================================

def primitive_from_dense_triplets(W_dense):
    """Convert a dense [rho, u, p] matrix into primitive structured states."""
    W = zeros_primitive(W_dense.shape[0])
    W["rho"] = W_dense[:, 0]
    W["u"] = W_dense[:, 1]
    W["p"] = W_dense[:, 2]
    return populate_primitive_derived(W)


def muscl_reconstruct(U, limiter_code):
    """Reconstruct left and right face states at all interfaces."""
    ncells = U.shape[0]
    nfaces = ncells + 1

    W = cons_to_prim(U)
    W_dense = primitive_component_matrix(W)
    W_L_dense = np.zeros((nfaces, 3))
    W_R_dense = np.zeros((nfaces, 3))

    for i in range(1, nfaces - 1):
        if i >= 2:
            dW_left_im1 = W_dense[i - 1] - W_dense[i - 2]
        else:
            dW_left_im1 = np.zeros(3)
        dW_right_im1 = W_dense[i] - W_dense[i - 1]

        r_im1 = dW_left_im1 / (dW_right_im1 + 1e-300 * np.sign(dW_right_im1 + 1e-300))
        slope_im1 = apply_limiter(limiter_code, r_im1) * dW_right_im1
        W_L_dense[i] = W_dense[i - 1] + 0.5 * slope_im1

        dW_left_i = W_dense[i] - W_dense[i - 1]
        if i < nfaces - 2:
            dW_right_i = W_dense[i + 1] - W_dense[i]
        else:
            dW_right_i = np.zeros(3)

        r_i = dW_left_i / (dW_right_i + 1e-300 * np.sign(dW_right_i + 1e-300))
        slope_i = apply_limiter(limiter_code, r_i) * dW_right_i
        W_R_dense[i] = W_dense[i] - 0.5 * slope_i

    U_L = prim_to_cons(primitive_from_dense_triplets(W_L_dense))
    U_R = prim_to_cons(primitive_from_dense_triplets(W_R_dense))
    return U_L, U_R


def first_order_states(U):
    """First-order piecewise-constant face states."""
    ncells = U.shape[0]
    nfaces = ncells + 1

    U_L = zeros_conservative(nfaces)
    U_R = zeros_conservative(nfaces)
    U_L[1:-1] = U[:-1]
    U_R[1:-1] = U[1:]
    return U_L, U_R


# =============================================================================
# 9. BOUNDARY CONDITIONS
# =============================================================================

def apply_boundary_conditions(U_L, U_R, U, bc_left, bc_right):
    """Fill ghost states at the left and right boundary faces."""
    if bc_left == "transmissive":
        U_L[0] = U[0]
        U_R[0] = U[0]
    elif bc_left == "reflective":
        U_ghost = copy_conservative(U[0])
        U_ghost["rhou"] = -U[0]["rhou"]
        U_L[0] = U_ghost
        U_R[0] = U[0]
    else:
        raise ValueError(f"Unknown left boundary condition '{bc_left}'")

    if bc_right == "transmissive":
        U_L[-1] = U[-1]
        U_R[-1] = U[-1]
    elif bc_right == "reflective":
        U_ghost = copy_conservative(U[-1])
        U_ghost["rhou"] = -U[-1]["rhou"]
        U_L[-1] = U[-1]
        U_R[-1] = U_ghost
    else:
        raise ValueError(f"Unknown right boundary condition '{bc_right}'")

    return U_L, U_R


# =============================================================================
# 10. RESIDUAL ASSEMBLY
# =============================================================================

def compute_residual(U, dx, scheme, limiter_code):
    """Compute dU/dt = -R(U) in structured conservative form."""
    if scheme == "muscl":
        U_L, U_R = muscl_reconstruct(U, limiter_code)
    else:
        U_L, U_R = first_order_states(U)

    U_L, U_R = apply_boundary_conditions(U_L, U_R, U, BC_LEFT, BC_RIGHT)
    F_face = hllc_flux_batch(U_L, U_R)

    R = zeros_conservative(U.shape)
    for field in CONSERVATIVE_DTYPE.names:
        R[field] = (F_face[field][1:] - F_face[field][:-1]) / dx
    return R


# =============================================================================
# 11. SSP-RK3 TIME INTEGRATOR
# =============================================================================

def ssp_rk3_step(U, dt, dx, scheme, limiter_code):
    """Advance one time step with SSP-RK3 using structured conservative state."""
    R1 = compute_residual(U, dx, scheme, limiter_code)
    U1 = conservative_add_scaled(U, -dt, R1)

    R2 = compute_residual(U1, dx, scheme, limiter_code)
    U2_rhs = conservative_add_scaled(U1, -dt, R2)
    U2 = conservative_linear_combo(U.shape, [(0.75, U), (0.25, U2_rhs)])

    R3 = compute_residual(U2, dx, scheme, limiter_code)
    U3_rhs = conservative_add_scaled(U2, -dt, R3)
    return conservative_linear_combo(U.shape, [(1.0 / 3.0, U), (2.0 / 3.0, U3_rhs)])


# =============================================================================
# 12. CFL CONDITION
# =============================================================================

def compute_dt(U, dx, cfl, t, t_end):
    """Compute a stable explicit time step."""
    s_max = max_wave_speed(U)
    dt = cfl * dx / s_max
    return min(dt, t_end - t)


# =============================================================================
# 13. INITIAL CONDITIONS
# =============================================================================

def setup_mesh():
    """Build a uniform 1D mesh."""
    dx = (X_END - X_START) / NCELLS
    x_faces = np.linspace(X_START, X_END, NCELLS + 1)
    x_cells = 0.5 * (x_faces[:-1] + x_faces[1:])
    return x_cells, dx


def sod_initial_condition(x_cells):
    """Initialize the Sod shock tube using structured conservative state."""
    U = zeros_conservative(len(x_cells))
    x_diaphragm = X_START + DIAPHRAGM * (X_END - X_START)

    for i, x in enumerate(x_cells):
        rho, u, p = W_LEFT if x < x_diaphragm else W_RIGHT
        U[i] = prim_to_cons(primitive_state(rho, u, p))

    return U


# =============================================================================
# 14. ANALYTICAL SOLUTION
# =============================================================================

def sod_analytical(x_cells, t):
    """Exact Sod shock tube solution returned as a dense [rho, u, p] matrix."""
    from scipy.optimize import brentq

    rho_L, u_L, p_L = W_LEFT
    rho_R, u_R, p_R = W_RIGHT
    a_L = sound_speed(rho_L, p_L)
    a_R = sound_speed(rho_R, p_R)

    g = GAMMA
    g1 = (g - 1) / (2 * g)
    g2 = (g + 1) / (2 * g)
    g3 = 2 * g / (g - 1)
    g4 = 2 / (g - 1)
    g5 = 2 / (g + 1)
    g6 = (g - 1) / (g + 1)

    def f_L(p):
        if p <= p_L:
            return g4 * a_L * ((p / p_L) ** g1 - 1)
        A = g5 / rho_L
        B = g6 * p_L
        return (p - p_L) * np.sqrt(A / (p + B))

    def f_R(p):
        if p <= p_R:
            return g4 * a_R * ((p / p_R) ** g1 - 1)
        A = g5 / rho_R
        B = g6 * p_R
        return (p - p_R) * np.sqrt(A / (p + B))

    def f(p):
        return f_L(p) + f_R(p) + (u_R - u_L)

    p_min = min(p_L, p_R) * 1e-6
    p_max = max(p_L, p_R) * 100
    p_star = brentq(f, p_min, p_max, xtol=1e-12)
    u_star = 0.5 * (u_L + u_R) + 0.5 * (f_R(p_star) - f_L(p_star))

    if p_star <= p_L:
        rho_starL = rho_L * (p_star / p_L) ** (1 / g)
    else:
        rho_starL = rho_L * (p_star / p_L + g6) / (g6 * p_star / p_L + 1)

    if p_star <= p_R:
        rho_starR = rho_R * (p_star / p_R) ** (1 / g)
    else:
        rho_starR = rho_R * (p_star / p_R + g6) / (g6 * p_star / p_R + 1)

    if p_star <= p_L:
        S_HL = u_L - a_L
        a_starL = a_L * (p_star / p_L) ** g1
        S_TL = u_star - a_starL
    else:
        S_L = u_L - a_L * np.sqrt(g2 * p_star / p_L + g1)

    if p_star <= p_R:
        S_HR = u_R + a_R
        a_starR = a_R * (p_star / p_R) ** g1
        S_TR = u_star + a_starR
    else:
        S_R = u_R + a_R * np.sqrt(g2 * p_star / p_R + g1)

    S_contact = u_star
    x_diaphragm = X_START + DIAPHRAGM * (X_END - X_START)
    xi = (x_cells - x_diaphragm) / (t + 1e-300)

    W_exact = np.zeros((len(x_cells), 3))

    for i, s in enumerate(xi):
        if p_star <= p_L:
            if s <= S_HL:
                W_exact[i] = [rho_L, u_L, p_L]
            elif s <= S_TL:
                u_fan = (2 / (g + 1)) * (a_L + (g - 1) / 2 * u_L + s)
                a_fan = u_fan - s
                rho_fan = rho_L * (a_fan / a_L) ** g4
                p_fan = p_L * (a_fan / a_L) ** g3
                W_exact[i] = [rho_fan, u_fan, p_fan]
            elif s <= S_contact:
                W_exact[i] = [rho_starL, u_star, p_star]
            elif s <= (S_R if p_star > p_R else S_TR):
                W_exact[i] = [rho_starR, u_star, p_star]
            else:
                W_exact[i] = [rho_R, u_R, p_R]
        else:
            if s <= S_L:
                W_exact[i] = [rho_L, u_L, p_L]
            elif s <= S_contact:
                W_exact[i] = [rho_starL, u_star, p_star]
            elif s <= (S_R if p_star > p_R else S_TR):
                W_exact[i] = [rho_starR, u_star, p_star]
            else:
                W_exact[i] = [rho_R, u_R, p_R]

    return W_exact


# =============================================================================
# 15. MAIN SOLVER LOOP
# =============================================================================

def run_solver():
    """Run the simulation and return the final structured state."""
    print("=" * 60)
    print("  PyFVM — 1D Euler Solver")
    print("=" * 60)
    print(f"  Cells:   {NCELLS}")
    print(f"  Scheme:  {SCHEME}" + (f" + {LIMITER}" if SCHEME == "muscl" else ""))
    print(f"  CFL:     {CFL}")
    print(f"  t_end:   {T_END}")
    print(f"  BC:      left={BC_LEFT}, right={BC_RIGHT}")
    print("=" * 60)

    x_cells, dx = setup_mesh()
    U = sod_initial_condition(x_cells)
    limiter_code = limiter_name_to_code(LIMITER)

    t = 0.0
    iteration = 0
    time_hist = [t]
    mass_hist = [np.sum(U["rho"]) * dx]

    while t < T_END:
        dt = compute_dt(U, dx, CFL, t, T_END)
        U = ssp_rk3_step(U, dt, dx, SCHEME, limiter_code)

        t += dt
        iteration += 1

        mass = np.sum(U["rho"]) * dx
        mass_hist.append(mass)
        time_hist.append(t)

        if iteration % PRINT_INTERVAL == 0:
            s_max = max_wave_speed(U)
            print(
                f"  iter={iteration:5d}  t={t:.4f}  dt={dt:.2e}  "
                f"max_speed={s_max:.3f}  mass={mass:.6f}"
            )

        check_physical(U, f"at iteration {iteration}")

    print(f"\n  Done. {iteration} iterations, t={t:.4f}")
    print(f"  Mass conservation error: {abs(mass_hist[-1] - mass_hist[0]):.2e}")

    return x_cells, dx, U, np.array(time_hist), np.array(mass_hist)


# =============================================================================
# 16. POST-PROCESSING
# =============================================================================

def plot_results(x_cells, U, t_final):
    """Plot rho, u, p with optional analytical solution overlay."""
    W = cons_to_prim(U)
    W_dense = primitive_component_matrix(W)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(
        f"Sod Shock Tube  |  t={t_final:.3f}  |  "
        f"N={NCELLS}  |  {SCHEME}"
        + (f"+{LIMITER}" if SCHEME == "muscl" else ""),
        fontsize=12,
    )

    labels = ["Density  rho", "Velocity  u", "Pressure  p"]
    colors = ["steelblue", "darkorange", "forestgreen"]

    if PLOT_ANALYTICAL:
        W_exact = sod_analytical(x_cells, t_final)

    for ax, col, label, var_idx in zip(axes, colors, labels, range(3)):
        ax.plot(x_cells, W_dense[:, var_idx], color=col, linewidth=1.8, label="Numerical")

        if PLOT_ANALYTICAL:
            ax.plot(x_cells, W_exact[:, var_idx], "k--", linewidth=1.2, label="Exact", alpha=0.7)

        ax.set_xlabel("x")
        ax.set_ylabel(label)
        ax.set_xlim(X_START, X_END)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
    return fig


def plot_conservation(time_hist, mass_hist):
    """Plot mass conservation over time."""
    fig, ax = plt.subplots(figsize=(7, 3))
    error = np.abs(mass_hist - mass_hist[0])
    ax.semilogy(time_hist, error + 1e-16)
    ax.set_xlabel("Time")
    ax.set_ylabel("|Mass error|")
    ax.set_title("Mass Conservation")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    return fig


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__" or True:
    x, dx, U, time_hist, mass_hist = run_solver()

    W = cons_to_prim(U)
    rho = W["rho"]
    u = W["u"]
    p = W["p"]
    T = W["T"]
    a = W["a"]
    M = u / a

    fig_solution = plot_results(x, U, T_END)
    fig_conservation = plot_conservation(time_hist, mass_hist)
