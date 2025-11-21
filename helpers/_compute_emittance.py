
import numpy as np

######################################### Sigma Matrix & Emittances Function ######################
def compute_emittance_and_optics(
    x, px, y, py, zeta, delta, particle_axis=0, eps_tol=0.0):
    """
    Estimate optics functions and emittances from tracking coordinates.

    Parameters
    ----------
    x, px, y, py, zeta, delta : ndarray
        Arrays with identical shape. One axis enumerates particles.
        Example: shape (n_particles, n_elements) or (n_turns, n_particles, n_elements).
    particle_axis : int, optional
        Index of the axis corresponding to particles (default 0).
    eps_tol : float, optional
        Numerical floor for emittance / variances. Below this we return eps=0
        and beta/alpha/gamma = NaN at that point.

    Returns
    -------
    results : dict of ndarrays
        Each array has shape equal to the input shape with the particle axis
        removed (i.e. per element / observation point).

        Horizontal:
            'epsx' : geometric emittance in x
            'betx', 'alfx', 'gamx' : Twiss parameters in x
            'Dx', 'Dpx' : dispersion and its slope

        Vertical:
            'epsy' : geometric emittance in y
            'bety', 'alfy', 'gamy' : Twiss parameters in y
            'Dy', 'Dpy' : vertical dispersion and slope

        Longitudinal:
            'epsz' : “longitudinal emittance” from (zeta, delta)
            'betz', 'alfz', 'gamz' : Twiss-like parameters in longitudinal plane

    Notes
    -----
    - Transverse dispersion removal is done via:
        Dx   = cov(a, delta) / var(delta)
        Dpx  = cov(pa, delta) / var(delta)
        aβ   = (a - <a>) - Dx * (delta - <delta>)
        p_aβ = (pa - <pa>) - Dpx * (delta - <delta>)

    Then Twiss from the 2×2 covariance of (aβ, p_aβ):
        eps   = sqrt(σ_aa σ_pp - σ_ap^2)
        beta  =  σ_aa / eps
        alpha = -σ_ap / eps
        gamma =  σ_pp / eps

    - For longitudinal plane, no “dispersion removal” is applied; we simply
    treat (zeta, delta) as a canonical pair and compute Twiss from their
    covariance matrix.
    """
    # Ensure numpy arrays of float
    x     = np.asarray(x,     float)
    px    = np.asarray(px,    float)
    y     = np.asarray(y,     float)
    py    = np.asarray(py,    float)
    zeta  = np.asarray(zeta,  float)
    delta = np.asarray(delta, float)

    # Move particle axis to front so shape is (n_particles, n_points...)
    def move(a):
        return np.moveaxis(a, particle_axis, 0)

    X = move(x)
    PX = move(px)
    Y = move(y)
    PY = move(py)
    Z = move(zeta)
    D = move(delta)

    n_part = X.shape[0]
    other_shape = X.shape[1:]
    n_points = int(np.prod(other_shape))

    def alloc():
        return np.zeros(n_points, float)

    # Allocate flat arrays for results
    epsx = alloc(); betx = alloc(); alfx = alloc(); gamx = alloc(); Dx = alloc(); Dpx = alloc()
    epsy = alloc(); bety = alloc(); alfy = alloc(); gamy = alloc(); Dy = alloc(); Dpy = alloc()
    epsz = alloc(); betz = alloc(); alfz = alloc(); gamz = alloc()

    # Flatten “points” to 1D index, keep particles as first axis
    Xf  = X.reshape(n_part, n_points)
    PXf = PX.reshape(n_part, n_points)
    Yf  = Y.reshape(n_part, n_points)
    PYf = PY.reshape(n_part, n_points)
    Zf  = Z.reshape(n_part, n_points)
    Df  = D.reshape(n_part, n_points)

    for i in range(n_points):
        xx = Xf[:,  i]; pxx = PXf[:, i]
        yy = Yf[:,  i]; pyy = PYf[:, i]
        zz = Zf[:,  i]; dd  = Df[:,  i]

        # ---------- Horizontal plane ----------
        ac  = xx - xx.mean()
        pac = pxx - pxx.mean()
        dc  = dd - dd.mean()

        s_ad = np.mean(ac * dc)
        s_pd = np.mean(pac * dc)
        s_dd = np.mean(dc * dc)

        if s_dd > eps_tol:
            Dx_i  = s_ad / s_dd
            Dpx_i = s_pd / s_dd
        else:
            Dx_i  = 0.0
            Dpx_i = 0.0

        # Remove dispersion
        ab  = ac  - Dx_i  * dc
        pab = pac - Dpx_i * dc

        s_bb  = np.mean(ab * ab)
        s_bp  = np.mean(ab * pab)
        s_ppb = np.mean(pab * pab)

        det = s_bb * s_ppb - s_bp**2
        eps_i = np.sqrt(max(det, 0.0))

        if eps_i > eps_tol:
            betx_i = s_bb / eps_i
            alfx_i = -s_bp / eps_i
            gamx_i = s_ppb / eps_i
        else:
            eps_i  = 0.0
            betx_i = alfx_i = gamx_i = np.nan

        epsx[i] = eps_i; betx[i] = betx_i; alfx[i] = alfx_i; gamx[i] = gamx_i
        Dx[i]   = Dx_i;  Dpx[i]  = Dpx_i

        # ---------- Vertical plane ----------
        ac  = yy - yy.mean()
        pac = pyy - pyy.mean()
        # reuse dc, s_dd

        s_ad = np.mean(ac * dc)
        s_pd = np.mean(pac * dc)
        # s_dd as above

        if s_dd > eps_tol:
            Dy_i  = s_ad / s_dd
            Dpy_i = s_pd / s_dd
        else:
            Dy_i  = 0.0
            Dpy_i = 0.0

        # Remove vertical dispersion (usually tiny, but we treat it symmetrically)
        ab  = ac  - Dy_i  * dc
        pab = pac - Dpy_i * dc

        s_bb  = np.mean(ab * ab)
        s_bp  = np.mean(ab * pab)
        s_ppb = np.mean(pab * pab)

        det = s_bb * s_ppb - s_bp**2
        eps_i = np.sqrt(max(det, 0.0))

        if eps_i > eps_tol:
            bety_i = s_bb / eps_i
            alfy_i = -s_bp / eps_i
            gamy_i = s_ppb / eps_i
        else:
            eps_i  = 0.0
            bety_i = alfy_i = gamy_i = np.nan

        epsy[i] = eps_i; bety[i] = bety_i; alfy[i] = alfy_i; gamy[i] = gamy_i
        Dy[i]   = Dy_i;  Dpy[i]  = Dpy_i

        # ---------- Longitudinal plane (zeta, delta) ----------
        zc = zz - zz.mean()
        dc = dd - dd.mean()

        s_zz = np.mean(zc * zc)
        s_zd = np.mean(zc * dc)
        s_dd = np.mean(dc * dc)

        det = s_zz * s_dd - s_zd**2
        eps_i = np.sqrt(max(det, 0.0))

        if eps_i > eps_tol:
            betz_i = s_zz / eps_i
            alfz_i = -s_zd / eps_i
            gamz_i = s_dd / eps_i
        else:
            eps_i  = 0.0
            betz_i = alfz_i = gamz_i = np.nan

        epsz[i] = eps_i; betz[i] = betz_i; alfz[i] = alfz_i; gamz[i] = gamz_i

    # Reshape back to original “points” shape (everything except particle axis)
    def unflatten(a):
        return a.reshape(other_shape)

    results = dict(
        epsx = unflatten(epsx),
        betx = unflatten(betx),
        alfx = unflatten(alfx),
        gamx = unflatten(gamx),
        Dx   = unflatten(Dx),
        Dpx  = unflatten(Dpx),

        epsy = unflatten(epsy),
        bety = unflatten(bety),
        alfy = unflatten(alfy),
        gamy = unflatten(gamy),
        Dy   = unflatten(Dy),
        Dpy  = unflatten(Dpy),

        epsz = unflatten(epsz),
        betz = unflatten(betz),
        alfz = unflatten(alfz),
        gamz = unflatten(gamz),
    )

    return results