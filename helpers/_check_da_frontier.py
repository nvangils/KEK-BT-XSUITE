"""
Check if given points are within a DA frontier curve.
=============================================

"""


#################################### Required Packages ###################################

import numpy as np


# Interpolation and Comparison Functions
################################################################################
def interpolate_frontier_y(xs, ys, xq, extrapolate = "linear"):
    """
    Vectorized interpolation of the frontier curve y_b(x) at arbitrary xq.

    Parameters
    ----------
    xs, ys : 1D arrays defining the boundary curve (must be sorted ascending in xs)
    xq     : array-like (any shape)
    extrapolate : {'linear', 'edge'}
        - 'linear'  → extend linearly beyond xs range
        - 'edge'    → clamp to edge values

    Returns
    -------
    yb : ndarray of same shape as xq
    """
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    xq = np.asarray(xq, float)

    # 1D interpolate inside domain
    yb = np.interp(xq.ravel(), xs, ys).reshape(xq.shape)

    if extrapolate == "edge":
        return yb

    # linear extrapolation on left/right
    if xs.size > 1:
        m_left  = (ys[1] - ys[0]) / (xs[1] - xs[0])
        m_right = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
    else:
        m_left = m_right = 0.0

    left_mask  = xq < xs[0]
    right_mask = xq > xs[-1]

    if np.any(left_mask):
        yb[left_mask]  = ys[0]  + m_left  * (xq[left_mask]  - xs[0])
    if np.any(right_mask):
        yb[right_mask] = ys[-1] + m_right * (xq[right_mask] - xs[-1])
    return yb


def compare_to_frontier(xs, ys, x_test, y_test, extrapolate = "linear"):
    """
    Compare arbitrary test points (x_test, y_test) to a frontier curve (xs, ys).

    Parameters
    ----------
    xs, ys : 1D arrays defining the frontier curve (xs ascending)
    x_test, y_test : array-like
        Coordinates of test points. Any broadcastable shape is allowed.
    extrapolate : {'linear', 'edge'}
    return_margin : bool
        If True, return the signed distance (y_test - y_b(x_test))
        If False, return boolean (True if y_test > y_b(x_test))

    Returns
    -------
    result : ndarray
        Boolean mask (if return_margin=False) or float array (if True).
        Same broadcasted shape as x_test and y_test.
    """
    x_test  = np.asarray(x_test, float)
    y_test  = np.asarray(y_test, float)

    yb      = interpolate_frontier_y(
        xs          = xs, 
        ys          = ys, 
        xq          = x_test,
        extrapolate = extrapolate)
    margin  = y_test - yb
    within  = (margin < 0)
    return within