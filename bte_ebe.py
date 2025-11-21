"""
Track the KEK BTe lattice
=============================================

"""

################################################################################
# Required Packages
################################################################################
import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

################################################################################
# User Parameters
################################################################################

########################################
# Lattice Files
########################################
xsuite_lattice_path     = "lattices/bte.py"
xsuite_optics_path      = "lattices/bte_import_optics.py"

########################################
# Initial Conditions
########################################
betx_init               = 24.097096226706224
bety_init               = 17.11998683241826
alfx_init               = .19232400480757933
alfy_init               = -2.131413440451613

########################################
# Tracking Parameters
########################################
ELE_STOP                = "injp"

################################################################################
# Load Initial distribution (from Andrea, Ocelot with Wakes)
################################################################################
initial_distribution    = np.load('999_amplitudes/andrea_dist_wake.npy', allow_pickle = True)
initial = xt.Table({
    "name":     np.arange(initial_distribution[0].shape[0]),
    "x":        +1 * initial_distribution[0],
    "px":       +1 * initial_distribution[1],
    "y":        +1 * initial_distribution[2],
    "py":       +1 * initial_distribution[3],
    "zeta":     -1 * initial_distribution[4],
    "delta":    +1 * initial_distribution[5] - np.mean(initial_distribution[5]),
    "state":    np.full(initial_distribution[0].shape, 1)})

################################################################################
# Load Lattice
################################################################################
env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]

################################################################################
# Twiss (4D)
################################################################################
tw = line.twiss4d(
    start   = xt.START,
    end     = xt.END,
    betx    = betx_init,
    bety    = bety_init,
    alfx    = alfx_init,
    alfy    = alfy_init)

################################################################################
# Build Particles
################################################################################
tracked_beam    = xt.Particles(
    p0c             = line.particle_ref.p0c,
    q0              = line.particle_ref.q0,
    mass0           = line.particle_ref.mass0,
    x               = initial.x,
    px              = initial.px,
    y               = initial.y,
    py              = initial.py,
    zeta            = initial.zeta,
    delta           = initial.delta)
initial_distribution = tracked_beam.copy()

################################################################################
# Track
################################################################################
line.track(
    particles               = tracked_beam,
    ele_stop                = ELE_STOP,
    turn_by_turn_monitor    = "ONE_TURN_EBE")
moni = line.record_last_track

################################################################################
# Calculate Statistical Emittances
################################################################################
def calculate_statistical_emittance(
    a, pa,
    delta           = None,
    plane           = "x",                 # 'x', 'y', or 'z' (longitudinal)
    particle_axis   = 0,           # your case: particles first
    element_axis    = 1,            # your case: elements second
    eps_tol         = 0.0,               # floor on determinant
):
    """
    Per-element statistical emittance and Twiss parameters from (a, pa),
    with correct handling of (x, px) vs (x, x') conventions and optional
    dispersion removal for transverse planes.

    Parameters
    ----------
    a, pa : ndarray
        Phase-space samples. For transverse: a=x or y [m], pa is either px/p0
        or angle x' [rad]. For longitudinal: a=zeta [m or s-equivalent], pa=delta.
        Shape can be (..., n_particles, n_elements, ...); set particle_axis/element_axis.
    delta : ndarray or None
        Relative momentum deviation per particle (same shape as a/pa).
        - Transverse: if provided, px is converted to angle via pa/(1+delta)
          and dispersion (D, D') is estimated and subtracted.
        - Longitudinal: ignored (pa is already delta).
        - If delta is None in transverse, pa is assumed to be angle already.
    plane : {'x','y','z'}
        Which plane to use. 'z' = longitudinal (zeta, delta).
    particle_axis, element_axis : int
        Axes in the input that correspond to particles and elements.
    eps_tol : float
        Non-negative floor applied to determinant before sqrt (useful to avoid
        tiny negative round-off).

    Returns
    -------
    out : dict of ndarrays with shape (n_elements,)
        epsilon, beta, alpha, gamma
        abar, pabar        # centroids (pabar is angle for transverse if delta provided)
        Dx, Dxp            # dispersion and angle-dispersion (zeros for 'z' plane or if not removed)
        dispersion_removed # bool per element
        mask_var_delta     # True where var(delta) > 0 was available (transverse only)
    """
    if plane not in ("x", "y", "z"):
        raise ValueError("plane must be 'x', 'y', or 'z'")

    a  = np.asarray(a,  dtype=float)
    pa = np.asarray(pa, dtype=float)
    if a.shape != pa.shape:
        raise ValueError("a and pa must have identical shapes")

    if particle_axis == element_axis:
        raise ValueError("particle_axis and element_axis must be different")

    # --- Reorder so: elements axis -> 0, particles axis -> 1
    axes = list(range(a.ndim))
    order = [element_axis, particle_axis] + [ax for ax in axes if ax not in (element_axis, particle_axis)]
    a  = np.transpose(a,  order)
    pa = np.transpose(pa, order)

    # Collapse any extra dims (require 2D)
    a  = a.reshape(a.shape[0], a.shape[1], -1)
    pa = pa.reshape(pa.shape[0], pa.shape[1], -1)
    if a.shape[2] != 1:
        raise ValueError("Inputs must be 2D (elements×particles) after reordering.")
    a  = a[..., 0]   # shape (n_ele, n_par)
    pa = pa[..., 0]

    n_ele, n_par = a.shape

    # Prepare delta with same reordering if provided
    if delta is not None:
        delta = np.asarray(delta, dtype=float)
        if delta.shape != np.transpose(np.asarray(delta), order).shape:
            # harmless; we just want to try the same transform
            pass
        delta = np.transpose(delta, order)
        delta = delta.reshape(n_ele, n_par, -1)
        if delta.shape[2] != 1:
            raise ValueError("delta must be 2D and match a/pa after reordering.")
        delta = delta[..., 0]
        if delta.shape != a.shape:
            raise ValueError("delta must match a/pa shape after reordering.")

    # ---------------- Longitudinal plane: (zeta, delta) ----------------
    if plane == "z":
        # No dispersion removal; pa is delta already.
        abar  = np.mean(a,  axis=1)
        pabar = np.mean(pa, axis=1)

        da  = a  - abar[:,  None]
        dpa = pa - pabar[:, None]

        aa = np.mean(da * da,   axis=1)
        pp = np.mean(dpa * dpa, axis=1)
        ap = np.mean(da * dpa,  axis=1)

        det = aa * pp - ap**2
        det = np.maximum(det, eps_tol if eps_tol > 0 else 0.0)
        eps = np.sqrt(det)

        with np.errstate(divide='ignore', invalid='ignore'):
            beta  = np.where(eps > 0, aa / eps, np.nan)
            alpha = np.where(eps > 0, -ap / eps, np.nan)
            gamma = np.where(eps > 0, pp / eps, np.nan)

        return dict(
            epsilon=eps, beta=beta, alpha=alpha, gamma=gamma,
            abar=abar, pabar=pabar,
            Dx=np.zeros(n_ele), Dxp=np.zeros(n_ele),
            dispersion_removed=np.zeros(n_ele, dtype=bool),
            mask_var_delta=np.zeros(n_ele, dtype=bool),
        )

    # --------------- Transverse planes: (x, px) or (y, py) ---------------
    # Decide which 'pa' to use downstream:
    # - If delta is provided: convert px/p0 -> angle via (1+delta)
    # - If delta is None: assume 'pa' is already angle
    if delta is not None:
        pa_eff = pa / (1.0 + delta)  # x' ≈ (px/p0)/(1+δ)
    else:
        pa_eff = pa  # assumed to be angle already

    # Centroids (per element, over particles)
    abar  = np.mean(a,      axis=1)
    pabar = np.mean(pa_eff, axis=1)

    ac   = a      - abar[:,  None]
    pac  = pa_eff - pabar[:, None]

    # Defaults
    Dx  = np.zeros(n_ele, dtype=float)   # position dispersion
    Dxp = np.zeros(n_ele, dtype=float)   # angle dispersion (D')
    dispersion_removed = np.zeros(n_ele, dtype=bool)
    mask_var_delta = np.zeros(n_ele, dtype=bool)

    # If delta provided, estimate (D, D') via correlations and subtract
    if delta is not None:
        dbar = np.mean(delta, axis=1)
        dc   = delta - dbar[:, None]

        var_d = np.mean(dc*dc, axis=1)
        mask_var_delta = var_d > 0

        if np.any(mask_var_delta):
            Dx[mask_var_delta]  = (np.mean(ac [mask_var_delta] * dc[mask_var_delta], axis=1)
                                   / var_d[mask_var_delta])
            Dxp[mask_var_delta] = (np.mean(pac[mask_var_delta] * dc[mask_var_delta], axis=1)
                                   / var_d[mask_var_delta])

        # Subtract dispersion (in the SAME variables used for the covariances)
        a_beta  = a      - Dx[:,  None]  * delta
        pa_beta = pa_eff - Dxp[:, None] * delta
        dispersion_removed = mask_var_delta.copy()
    else:
        a_beta, pa_beta = a, pa_eff

    # Moments on betatron coordinates
    ab   = np.mean(a_beta,  axis=1)
    pab  = np.mean(pa_beta, axis=1)
    da   = a_beta  - ab[:,  None]
    dpa  = pa_beta - pab[:, None]

    aa = np.mean(da  * da,  axis=1)   # <a^2>
    pp = np.mean(dpa * dpa, axis=1)   # <pa^2>
    ap = np.mean(da  * dpa, axis=1)   # <a pa>

    det = aa * pp - ap**2
    det = np.maximum(det, eps_tol if eps_tol > 0 else 0.0)
    eps = np.sqrt(det)

    with np.errstate(divide='ignore', invalid='ignore'):
        beta  = np.where(eps > 0, aa / eps, np.nan)
        alpha = np.where(eps > 0, -ap / eps, np.nan)
        gamma = np.where(eps > 0, pp / eps, np.nan)

    return dict(
        epsilon=eps, beta=beta, alpha=alpha, gamma=gamma,
        abar=abar, pabar=pabar,            # pabar is angle mean if delta was provided
        Dx=Dx, Dxp=Dxp,                    # zeros if no delta given
        dispersion_removed=dispersion_removed,
        mask_var_delta=mask_var_delta,
    )


import numpy as np

def transverse_eigen_emittances(
    x, px, y, py, delta=None, *,
    particle_axis=0, element_axis=1,
    remove_dispersion=True,  # subtract D, D' using correlations with delta
    use_angle=True,          # convert px,py -> x',y' via /(1+delta) if delta provided
    eps_tol=0.0              # nonnegative floor for numerical stability
):
    """
    Compute per-element transverse symplectic eigen-emittances (eps1>=eps2)
    from particle samples of (x, px, y, py) with optional dispersion removal.

    Parameters
    ----------
    x, px, y, py : ndarray, shape (..., n_particles, n_elements, ...)
        Particle coordinates. Units: x,y [m]; px,py either px/p0 or angles.
    delta : ndarray or None, same shape as x
        Relative momentum deviation. If provided and use_angle=True, px,py are
        converted to angles via /(1+delta). If remove_dispersion=True, D and D'
        are estimated per element and subtracted.
    particle_axis, element_axis : int
        Which axes correspond to particles and elements.
    eps_tol : float
        Floor applied to covariance-derived determinants/eigs.

    Returns
    -------
    dict with ndarrays of shape (n_elements,)
        eps1, eps2       : transverse eigen-emittances (linear invariants)
        Dx, Dxp, Dy, Dyp : dispersions used (zeros if not removed or no delta)
        n_used           : number of valid particles used per element
    """

    # -- helpers --
    def _reorder(A):
        A = np.asarray(A, float)
        axes = list(range(A.ndim))
        order = [element_axis, particle_axis] + [ax for ax in axes if ax not in (element_axis, particle_axis)]
        B = np.transpose(A, order).reshape(A.shape[element_axis], A.shape[particle_axis], -1)
        if B.shape[2] != 1:
            raise ValueError("Inputs must be 2D (elements×particles) after reordering.")
        return B[..., 0]  # (n_ele, n_part)

    def _symplectic_eigs(Sigma4):
        # J for 4D
        J = np.array([[0,1,0,0],
                      [-1,0,0,0],
                      [0,0,0,1],
                      [0,0,-1,0]], float)
        s = np.linalg.svd(J @ Sigma4, compute_uv=False)
        s = np.sort(s)[::-1]
        # for a covariance, singular values appear in (e1,e1,e2,e2)
        e1 = s[0]
        e2 = s[2] if s.size >= 3 else s[0]
        return max(e1, eps_tol), max(e2, eps_tol)

    # -- reorder inputs to (n_ele, n_part) --
    X  = _reorder(x);  PX = _reorder(px)
    Y  = _reorder(y);  PY = _reorder(py)
    if delta is not None:
        D = _reorder(delta)
    n_ele, n_part = X.shape

    # -- choose angle variables --
    if delta is not None and use_angle:
        Xp = PX / (1.0 + D)
        Yp = PY / (1.0 + D)
    else:
        Xp = PX
        Yp = PY

    # -- optional dispersion removal (per element) --
    Dx  = np.zeros(n_ele); Dxp = np.zeros(n_ele)
    Dy  = np.zeros(n_ele); Dyp = np.zeros(n_ele)

    if (delta is not None) and remove_dispersion:
        for i in range(n_ele):
            # finite mask
            m = np.isfinite(X[i]) & np.isfinite(Xp[i]) & np.isfinite(Y[i]) & np.isfinite(Yp[i]) & np.isfinite(D[i])
            if m.sum() < 3:
                continue
            d  = D[i, m]
            dc = d - d.mean()
            vard = np.mean(dc*dc)
            if vard <= 0:
                continue
            x  = X [i, m];  xp = Xp[i, m]
            y  = Y [i, m];  yp = Yp[i, m]
            Dx[i]  = np.mean((x  - x.mean())  * dc) / vard
            Dxp[i] = np.mean((xp - xp.mean()) * dc) / vard
            Dy[i]  = np.mean((y  - y.mean())  * dc) / vard
            Dyp[i] = np.mean((yp - yp.mean()) * dc) / vard
            # subtract dispersion (betatron coords)
            X [i, m] = x  - Dx[i]  * d
            Xp[i, m] = xp - Dxp[i] * d
            Y [i, m] = y  - Dy[i]  * d
            Yp[i, m] = yp - Dyp[i] * d

    # -- build 4D covariances and get eigen-emittances per element --
    eps1 = np.full(n_ele, np.nan); eps2 = np.full(n_ele, np.nan)
    n_used = np.zeros(n_ele, dtype=int)

    for i in range(n_ele):
        m = np.isfinite(X[i]) & np.isfinite(Xp[i]) & np.isfinite(Y[i]) & np.isfinite(Yp[i])
        n = int(m.sum()); n_used[i] = n
        if n < 4:
            continue
        # zero-mean columns [x, x', y, y']
        U = np.vstack([
            X [i, m] - np.mean(X [i, m]),
            Xp[i, m] - np.mean(Xp[i, m]),
            Y [i, m] - np.mean(Y [i, m]),
            Yp[i, m] - np.mean(Yp[i, m]),
        ]).T
        Sigma4 = (U.T @ U) / n
        e1, e2 = _symplectic_eigs(Sigma4)
        eps1[i], eps2[i] = e1, e2

    return dict(eps1=eps1, eps2=eps2, Dx=Dx, Dxp=Dxp, Dy=Dy, Dyp=Dyp, n_used=n_used)





line["volt_acw"] = 0
line.configure_bend_model(edge = "suppressed")
tt      = line.get_table(attr = True)
tt_sext = tt.rows[tt.element_type == "Sextupole"]
line.set(
    tt_sext,
    k2  = 0,
    k2s = 0)

emitt_x_dict    = calculate_statistical_emittance(
    a               = moni.x,
    pa              = moni.px,
    delta           = moni.delta,
    plane           = "x",
    particle_axis   = 0,
    element_axis    = 1)
emitt_y_dict    = calculate_statistical_emittance(
    a               = moni.y,
    pa              = moni.py,
    delta           = moni.delta,
    plane           = "y",
    particle_axis   = 0,
    element_axis    = 1)
emitt_z_dict    = calculate_statistical_emittance(
    a               = moni.x,
    pa              = moni.px,
    delta           = None,
    plane           = "z",
    particle_axis   = 0,
    element_axis    = 1)

# Arrays shaped (n_particles, n_elements)
out = transverse_eigen_emittances(moni.x, moni.px, moni.y, moni.py, moni.delta,
                                  particle_axis=0, element_axis=1,
                                  remove_dispersion=True, use_angle=True)

emitx, emity = out["eps1"], out["eps2"]  # should be flat for linear optics


################################################################################
# Plot
################################################################################
fig, axs = plt.subplots(3, 1, figsize = (10, 8), sharex = True)
axs[0].plot(tw.s, emitt_x_dict["epsilon"], label = "emitt_x")
axs[1].plot(tw.s, emitt_y_dict["epsilon"], label = "emitt_y")
axs[2].plot(tw.s, emitt_z_dict["epsilon"], label = "emitt_z")
plt.show()

fig, axs = plt.subplots(2, 1, figsize = (10, 8), sharex = True)
axs[0].plot(tw.s, emitx, label = "emitt_x")
axs[1].plot(tw.s, emity, label = "emitt_y")
plt.show()

fig, axs = plt.subplots(2, 1, figsize = (10, 8), sharex = True)
axs[0].plot(tw.s, emitt_x_dict["beta"], label = "Calc")
axs[0].plot(tw.s, tw.betx, label = "Twiss")
axs[1].plot(tw.s, emitt_y_dict["beta"], label = "Calc")
axs[1].plot(tw.s, tw.bety, label = "Twiss")

fig, axs = plt.subplots(2, 1, figsize = (10, 8), sharex = True)
axs[0].plot(tw.s, emitt_x_dict["alpha"], label = "Calc")
axs[0].plot(tw.s, tw.alfx, label = "Twiss")
axs[1].plot(tw.s, emitt_y_dict["alpha"], label = "Calc")
axs[1].plot(tw.s, tw.alfy, label = "Twiss")

fig, axs = plt.subplots(2, 1, figsize = (10, 8), sharex = True)
axs[0].plot(tw.s, emitt_x_dict["Dx"], label = "Calc")
axs[0].plot(tw.s, tw.dx, label = "Twiss")
axs[1].plot(tw.s, emitt_y_dict["Dx"], label = "Calc")
axs[1].plot(tw.s, tw.dy, label = "Twiss")

fig, axs = plt.subplots(2, 1, figsize = (10, 8), sharex = True)
axs[0].plot(tw.s, emitt_x_dict["Dxp"], label = "Calc")
axs[0].plot(tw.s, tw.dpx, label = "Twiss")
axs[1].plot(tw.s, emitt_y_dict["Dxp"], label = "Calc")
axs[1].plot(tw.s, tw.dpy, label = "Twiss")
plt.show()