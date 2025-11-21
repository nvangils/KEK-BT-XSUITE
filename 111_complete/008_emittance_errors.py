"""
Tracks input distribution to INJP and calculates emittance along BTe for random quadrupole misalignments
User choices: radiation ON/OFF tapering ON/OFF ECS settings and injection reference energy of line
Output: emittances along BTe with random quadrupole errors (alignment)
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       18-11-2025

"""
######################################### Required Packages ########################################

import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from helpers._fit_ellipse_da import calculate_boundary_points
# from helpers._local_momentum_taper import taper_magnets_to_local_momentum
from helpers._local_momentum_taper import taper_magnets_to_local_momentum
from tqdm import tqdm

######################################## User Parameters ##########################################

Radiation   = True      # TRUE ENABLES RADIATION 
Tapering    = False     # TRUE ENABLES TAPERING 
ECS_voltage = 0         # SET A VALUE IN MV TO ENABLE ELSE DEFAULT IS 18MV
ECS_lag     = None      # SET A LAG IN DEGREES TO ENABLE ELSE DEFAULT IS 180° ZERO CROSSING IN XSUITE

######################################### Lattice Files ########################################

xsuite_lattice_path     = "../lattices/bte.py"
xsuite_optics_path      = "../lattices/bte_import_optics.py"

######################################## Particles Filepath ####################################

particles_filepath      = "../data/yoshimotosan.dat"

######################################### Initial Conditions ####################################

betx_init               = 20.54043709
bety_init               = 18.01558997
alfx_init               = 1.00306534
alfy_init               = 0.370038838

#################################### Synchrotron Injection Energy Offset #######################

SYNCHROTRON_INJ_E_SCALE = 1.0076            # [1], relative energy increase at injection 1.006 
REFERENCE_ENERGY_HER    = 7.00E9           # [GeV], reference energy of HER

##################################### Final Element ############################################

ELE_STOP                = "injp"

##################################### Aperture Parameters ######################################

# T. Mori:
#   Round cross-section with 47 mm diameter at drifts and Q-magnets
#   Rectangle with about 30 mm height and 50-100 mm width at the dipoles
RADIUS_BP_APERTURE      = 0.0235           # [m],      Beam pipe radius at drifts and Q-magnets
HEIGHT_BEND_APERTURE    = 0.03             # [m],      Bend aperture height
WIDTH_BEND_APERTURE     = 0.05             # [m],      Bend aperture width

######################################### Load Lattice ###########################################

######################################### Load from Environment ##################################

env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]
if ECS_voltage is not None:
    line["volt_acw"] = ECS_voltage * 1E6
else :
    line["volt_acw"] = 18E6
if ECS_lag is not None:
    line["acw_phase"]  = ECS_lag
else :
    line["acw_phase"] = 180


######################################### Adjust Reference Energy #################################

line.particle_ref.p0c[0]   *= SYNCHROTRON_INJ_E_SCALE


######################################### Install Apertures #######################################

#################################### Bounding Aperture Helper Function ############################

def add_bounding_apertures(name, idx, aperture, aper_dict, idx_aper):
    idx_aper.append(idx-1)
    aper_dict[name + '_aper_entry'] = aperture.copy()
    idx_aper.append(idx)
    aper_dict[name + '_aper_exit'] = aperture.copy()
    
######################################### Insert Apertures Function ################################

def insert_apertures(line, idx_aper, aper_dict):
    shift = 0
    for idx, (aper_name, aper) in zip(idx_aper, aper_dict.items()):
        adjusted_idx = idx + shift  # Adjust index to account for prior insertions
        line.insert_element(
            at      = adjusted_idx + 1,
            name    = aper_name,
            element = aper)
        shift += 1  # Each insertion shifts future indices by 1

######################################### Generate the apertures ##################################

bp_aperture     = xt.LimitEllipse(
    a       = RADIUS_BP_APERTURE,
    b       = RADIUS_BP_APERTURE)
bend_aperture   = xt.LimitRect(
    min_x   = -1 * WIDTH_BEND_APERTURE / 2,
    max_x   = +1 * WIDTH_BEND_APERTURE / 2,
    min_y   = -1 * HEIGHT_BEND_APERTURE / 2,
    max_y   = +1 * HEIGHT_BEND_APERTURE / 2)

######################################### Get a clean table ######################################

tt  = line.get_table(attr = True)

######################################### Get bends and other active elements #####################

tt_bend     = tt.rows[
    (tt.element_type == "Bend") |
    (tt.element_type == "ThickSliceBend") |
    (tt.element_type == "RBend") |
    (tt.element_type == "ThickSliceRBend")]

tt_active   = tt.rows[
    (tt.element_type != 'Drift') & 
    (tt.element_type != 'Marker') &
    (tt.element_type != 'Bend') &
    (tt.element_type != "ThickSliceBend") &
    (tt.element_type != 'RBend') &
    (tt.element_type != "ThickSliceRBend") &
    (tt.element_type != '')]


###################### Get a list of which types of elements need aperture ########################

needs_bp_aperture   = np.unique(tt_active.element_type)
needs_bend_aperture = np.unique(tt_bend.element_type)


######################################## Calculate aperture positions  ############################

idx_aper    = []
aper_dict   = {}
for idx, (ele_name, ele_type) in enumerate(zip(tt.name, tt.element_type)):
    if ele_type in needs_bp_aperture:
        add_bounding_apertures(
            name        = ele_name,
            idx         = idx,
            aperture    = bp_aperture,
            aper_dict   = aper_dict,
            idx_aper    = idx_aper)
    if ele_type in needs_bend_aperture:
        add_bounding_apertures(
            name        = ele_name,
            idx         = idx,
            aperture    = bend_aperture,
            aper_dict   = aper_dict,
            idx_aper    = idx_aper)

####################### Install the apertures ######################################################

insert_apertures(line, idx_aper, aper_dict)


############################################# Check the apertures ##################################

aper_check = line.check_aperture()


####################### Setup #######################################################################

####################### Get clean table ######################

tt  = line.get_table(attr = True)

######################################### Build Particles ########################################

######################################### Load Particles Distribution ############################

initial_distribution    = np.loadtxt(particles_filepath)

######################################### Create Initial Table ###################################
initial_distribution    = xt.Table({
    "name":     np.arange(initial_distribution[:, 0].shape[0]),
    "x":        +1 * initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
    "px":       +1 * initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
    "y":        +1 * initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
    "py":       +1 * initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
    "zeta":     +1 * initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
    "delta":    +1 * initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
    "state":    np.full(initial_distribution[:, 0].shape, 1)})

######################################### Generate bunch ########################################

initial_beam    = xt.Particles(
    p0c             = line.particle_ref.p0c,
    q0              = line.particle_ref.q0,
    mass0           = line.particle_ref.mass0,
    x               = initial_distribution.x,
    px              = initial_distribution.px,
    y               = initial_distribution.y,
    py              = initial_distribution.py,
    zeta            = initial_distribution.zeta,
    delta           = initial_distribution.delta)
tracked_beam    = initial_beam.copy()

################################## Setup for tracking #######################################


################################## Enable RADIATION #########################################

########################### Taper to local momentum deviation ########################################
tw0 = line.twiss(
    start   = xt.START,
    end     = "injp",
    betx    = betx_init,
    bety    = bety_init,
    alfx    = alfx_init,
    alfy    = alfy_init)

if Radiation == True:
    line.configure_radiation(model = "mean")
else:
    pass

tw1 = line.twiss(
    start   = xt.START,
    end     = "injp",
    betx    = betx_init,
    bety    = bety_init,
    alfx    = alfx_init,
    alfy    = alfy_init)

if Tapering == True:
    taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
else:
    pass

tw2 = line.twiss(
    start   = xt.START,
    end     = "injp",
    betx    = betx_init,
    bety    = bety_init,
    alfx    = alfx_init,
    alfy    = alfy_init)

tw0.plot("x y")
plt.title("Before Radiation and Tapering")
tw1.plot("x y")
plt.title("After Radiation, Before Tapering")
tw2.plot("x y")
plt.title("After Radiation and Tapering")

######################################## Randomly offset quads ########################
quad_names      = tt.rows[tt.element_type == "Quadrupole"].name
quad_parents    = tt.rows[tt.element_type == "ThickSliceQuadrupole"].parent_name
quad_names      = list(set(list(set(quad_names)) + list(set(quad_parents))))
quad_names      = [name for name in quad_names if name is not None]

for quad_name in quad_names:
    line[quad_name].shift_y += np.random.normal(loc = 0.0, scale = 1E-4)

######################################## Quantum Synrad ########################
if Radiation == True:
    line.configure_radiation(model = "quantum")
else:
    pass

######################################## Track ##############################################

line.track(particles = tracked_beam, ele_stop  = ELE_STOP, turn_by_turn_monitor = "ONE_TURN_EBE")
moni = line.record_last_track

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

    for i in tqdm(range(n_points)):
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

results = compute_emittance_and_optics(
    x               = moni.x[tracked_beam.state == 1],
    px              = moni.px[tracked_beam.state == 1],
    y               = moni.y[tracked_beam.state == 1],
    py              = moni.py[tracked_beam.state == 1],
    zeta            = moni.zeta[tracked_beam.state == 1],
    delta           = moni.delta[tracked_beam.state == 1],
    particle_axis   = 0,)

gemitt_x    = results["epsx"]
gemitt_y    = results["epsy"]
gemitt_z    = results["epsz"]

betx        = results["betx"]
bety        = results["bety"]

dx          = results["Dx"]
dpx         = results["Dpx"]
dy          = results["Dy"]
dpy         = results["Dpy"]

######################################### Compute Emittances ######################################
print("Initial Emittances (mrad m):")
print(gemitt_x[0], gemitt_y[0], gemitt_z[0])
print("Final Emittances (mrad m):")
print(gemitt_x[-1], gemitt_y[-1], gemitt_z[-1])

# ######################################### Check optical function calculation ###################################
fig, axs = plt.subplots(1, 2, figsize=(18, 6))
axs[0].plot(tw0.s, tw2.betx, label = "Twiss")
axs[0].plot(tw0.s, betx[:tw0.s.shape[0]], linestyle = "--", label = "Tracking")
axs[1].plot(tw0.s, tw2.bety)
axs[1].plot(tw0.s, bety[:tw0.s.shape[0]], linestyle = "--")
axs[0].set_ylabel("Betx [m]")
axs[1].set_ylabel("Bety [m]")
axs[0].set_xlabel("s [m]")
for ax in axs:
    ax.legend()
    ax.grid()
fig.suptitle("Comparison of Twiss and Tracked Beta Functions")

fig, axs = plt.subplots(1, 2, figsize=(18, 6))
axs[0].plot(tw2.s, tw2.dx, label = "Twiss")
axs[0].plot(tw2.s, dx[:tw2.s.shape[0]], linestyle = "--", label = "Tracking")
axs[1].plot(tw2.s, tw2.dy)
axs[1].plot(tw2.s, dy[:tw2.s.shape[0]], linestyle = "--")
axs[0].set_ylabel("Dx [m]")
axs[1].set_ylabel("Dy [m]")
axs[0].set_xlabel("s [m]")
for ax in axs:
    ax.legend()
    ax.grid()
fig.suptitle("Comparison of Twiss and Tracked Dispersion Functions")

fig, axs = plt.subplots(1, 2, figsize=(18, 6))
axs[0].plot(tw2.s, tw2.dpx, label = "Twiss")
axs[0].plot(tw2.s, dpx[:tw2.s.shape[0]], linestyle = "--", label = "Tracking")
axs[1].plot(tw2.s, tw2.dpy)
axs[1].plot(tw2.s, dpy[:tw2.s.shape[0]], linestyle = "--")
axs[0].set_ylabel("Dpx []")
axs[1].set_ylabel("Dpy []")
axs[0].set_xlabel("s [m]")
for ax in axs:
    ax.legend()
    ax.grid()
fig.suptitle("Comparison of Twiss and Tracked Dispersion Slope Functions")

# ######################################### Plot Beam Phase Space ###################################
fig, axs = plt.subplots(1, 3, figsize=(18, 6))
axs[0].plot(tw0.s[:-1], gemitt_x[:tw0.s.shape[0]-1])
axs[1].plot(tw0.s[:-1], gemitt_y[:tw0.s.shape[0]-1])
axs[2].plot(tw0.s[:-1], gemitt_z[:tw0.s.shape[0]-1])
axs[0].set_ylabel("Horizontal Emittance [m rad]")
axs[1].set_ylabel("Vertical Emittance [m rad]")
axs[2].set_ylabel("Longitudinal Emittance [m rad]")
axs[0].set_xlabel("s [m]")

fig, axs = plt.subplots(1, 3, figsize=(18, 6))
axs[0].plot(tw0.s[:-1], gemitt_x[:tw0.s.shape[0]-1] / gemitt_x[0])
axs[1].plot(tw0.s[:-1], gemitt_y[:tw0.s.shape[0]-1] / gemitt_y[0])
axs[2].plot(tw0.s[:-1], gemitt_z[:tw0.s.shape[0]-1] / gemitt_z[0])
axs[0].set_ylabel("Horizontal Emittance Ratio [1]")
axs[1].set_ylabel("Vertical Emittance Ratio [1]")
axs[2].set_ylabel("Longitudinal Emittance Ratio [1]")
axs[0].set_xlabel("s [m]")
plt.show()