"""
Check if given points are within a DA frontier curve.
=============================================

"""

################################################################################
# Required Packages
################################################################################
import numpy as np

########################################
# Fit Ellipse to Survivor Boundary
########################################
def interpolate_edge(v0, v1, z0, z1, threshold = 0.0):
    if not (np.isfinite(z0) and np.isfinite(z1)) or z0 == z1:
        return 0.5 * (v0 + v1)
    t = np.clip((threshold - z0) / (z1 - z0), 0.0, 1.0)
    return v0 + t * (v1 - v0)

def calculate_boundary_points(X, Y, Z, turn_eps=0.0, x_tol='auto'):
    """
    Extract a clean frontier from positive-quadrant grids using *columns only*.

    Returns
    -------
    xb, yb : 1D arrays, sorted by xb, with yb non-increasing (no loops).
    """
    X   = np.asarray(X, float)
    Y   = np.asarray(Y, float)
    Z   = np.asarray(Z, float)

    assert X.shape == Y.shape == Z.shape
    
    rows, cols  = Z.shape
    thresh      = np.nanmax(Z) - float(turn_eps)

    # 1) Per column: take highest y that survives; interpolate to the crossing
    xs, ys = [], []
    for j in range(cols):
        zj      = Z[:, j]
        surv    = np.where(zj >= thresh)[0]
        if surv.size == 0:
            continue
        i = surv.max()
        if i < rows - 1 and np.isfinite(zj[i+1]):
            x_cross = interpolate_edge(X[i, j], X[i+1, j], Z[i, j], Z[i+1, j], thresh)
            y_cross = interpolate_edge(Y[i, j], Y[i+1, j], Z[i, j], Z[i+1, j], thresh)
        else:
            x_cross = X[i, j]
            y_cross = Y[i, j]
        xs.append(x_cross); ys.append(y_cross)

    if not xs:
        return np.array([]), np.array([])

    xs = np.asarray(xs); ys = np.asarray(ys)

    # 2) Sort by x (monotonic param for plotting)
    order   = np.argsort(xs)
    xs, ys  = xs[order], ys[order]

    # 3) Merge near-duplicate x's (keep the *max y* as the frontier)
    if x_tol == 'auto':
        # estimate step from median column spacing
        uniq = np.unique(np.round(xs, 12))
        if uniq.size > 1:
            step = np.median(np.diff(uniq))
        else:
            step = (xs.max() - xs.min()) / max(len(xs) - 1, 1) if len(xs) > 1 else 1.0
        x_tol = 0.25 * abs(step)

    merged_x = [xs[0]]
    merged_y = [ys[0]]
    for x, y in zip(xs[1:], ys[1:]):
        if abs(x - merged_x[-1]) <= x_tol:
            # take the upper envelope at this x-bin
            if y > merged_y[-1]:
                merged_y[-1] = y
        else:
            merged_x.append(x); merged_y.append(y)

    xs = np.asarray(merged_x); ys = np.asarray(merged_y)

    # 4) Enforce non-increasing y(x) to avoid tiny ups causing loops
    out_y       = np.empty_like(ys)
    cur         = ys[0]
    out_y[0]    = cur
    for k in range(1, len(ys)):
        cur         = min(cur, ys[k])
        out_y[k]    = cur

    return xs, out_y