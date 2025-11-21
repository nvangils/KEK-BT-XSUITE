
"""
Check the tracked beam against the previously computed DA for BTE lattice.
This maps the BTe tracked beam at INJP and compares it to the DA at the septum (everything in DA space)
Reference frame transform moves beam at INJP into reference frame of bump.
User choices: emittance blowup at the end of the line, radiation ON/OFF tapering ON/OFF ECS settings (LAG SCAN) and injection type
betatron/synchrotron or hybrid 
Output: 
MAXIMAL survival fraction within DA of HER  scanning the ECS voltage and ECS lag
Optimiser steps show dependence of survivability on ECS voltage and ECS lag 
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       18-11-2025

"""
######################################### Required Packages ########################################
import xtrack as xt
import xobjects as xo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import sys
from scipy.optimize import differential_evolution
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from helpers._calculate_action_amplitude import calculate_action_amplitude
from helpers._fit_ellipse_da import calculate_boundary_points
from helpers._check_da_frontier import compare_to_frontier
from helpers._local_momentum_taper import taper_magnets_to_local_momentum
from helpers._compute_emittance import compute_emittance_and_optics

######################################## User Parameters ##########################################
#Emittance blowup:
EMITTANCE_X_MEASURED = 140e-6/14000	 # MEASURED EMITTANCE AT BTE2 IN NORMALISED
EMITTANCE_Y_MEASURED = 140e-6/14000	 # MEASURED EMITTANCE AT BTE2 IN NORMALISED
#Line configuration:
Radiation = True 	                 # TRUE ENABLES RADIATION 
Tapering = True 	                 # TRUE ENABLES TAPERING 
#ECS settings:
ECS_lag = None                       # SET A VALUE IN DEGREES TO ENABLE ELSE DEFAULT IS 108° ZERO CROSSING IN XSUITE
set_voltage_bounds = (18E6, 100E6)   # VOLTAGE BOUNDS FOR THE OPTIMISER in VOLTS
set_lag_bounds = (175, 185)          # LAG BOUNDS FOR THE OPTIMISER in DEGREES
aperture_check = True                # TRUE ENABLES APERTURE CHECKING

#Injection type
#Betatron injection: 
BETATRON_OFFSET= 7.6e-3              # SET TO 0 IF PURE SYNCHROTRON INJECTION
BETATRON_ANGLE_OFFSET= 0.6e-3        # SET TO 0 IF PURE SYNCHROTRON INJECTION
#Synchrotron injection 
SYNCHROTRON_OFFSET_BT= 0.00165       # IN % WITH RESPECT TO  HER ENERGY. SET TO 0 FOR BETATRON INJECTION

######################################### Lattice Files ########################################

xsuite_lattice_path     = "../lattices/bte.py"
xsuite_optics_path      = "../lattices/bte_import_optics.py"

######################################## Particles Filepath ####################################

#particles_filepath      = "../data/input_electrons_wake.npy"
particles_filepath      = "../data/yohsimotosan.dat"


######################################### DA Filepaths ##########################################

##For 200_8 mm lattices uncomment below 
da_xy_filepath         = "../data/da_xy_folded_her_oide_200_8_her.npz"
da_zx_filepath         = "../data/da_zx_folded_her_oide_200_8_her.npz"
da_zy_filepath         = "../data/da_zy_folded_her_oide_200_8_her.npz"

##For 100_3 mm lattices uncomment below 
# da_xy_filepath         = "../data/da_xy_folded_her_oide_100_3_her.npz"
# da_zx_filepath         = "../data/da_zx_folded_her_oide_100_3_her.npz"
# da_zy_filepath         = "../data/da_zy_folded_her_oide_100_3_her.npz"

######################################### Initial Conditions ####################################
betx_init               = 20.540437095810788 
bety_init               = 18.015589974012826 
alfx_init               = 1.0030653048224993 
alfy_init               = .37003884861294073 
#################################### Synchrotron Injection Energy Offset #######################

SYNCHROTRON_INJ_E_SCALE = 1.000 + SYNCHROTRON_OFFSET_BT          # Relative energy so that increase at injection is 1.006 
REFERENCE_ENERGY_HER    = 7.00E9                                 # [GeV], reference energy of HER

#################################### Lattice Injection Optics Parameters #######################

##For 100_3 mm lattices uncomment below 
# BETX_HER_INJECTIO       = 99.9986 
# BETY_HER_INJECTIO       = 2.55671
# ALFX_HER_INJECTIO       = -8.50447 
# ALFY_HER_INJECTIO       = -1.47663
# DX_HER_INJECTIO         = -1.6
# DY_HER_INJECTIO         = 0.0
# DDX_HER_INJECTIO        = 26.7511
# DDY_HER_INJECTIO        = 0.00048039
# DDPX_HER_INJECTIO       = 5.2278e-05
# DDPY_HER_INJECTIO       = -0.000376
# DPX_HER_INJECTIO        = -0.120117  
# DPY_HER_INJECTIO        = -4.12801e-07
# BETS0_HER               = 8.056774638483907

##For 200_8 mm lattices uncomment below 
BETX_HER_INJECTIO       = 99.9987
BETY_HER_INJECTIO       = 2.55671
ALFX_HER_INJECTIO       = -8.50448
ALFY_HER_INJECTIO       = -1.47663
DX_HER_INJECTIO         = -1.60 
DY_HER_INJECTIO         = 0.0
DDX_HER_INJECTIO        = 23.70  
DDY_HER_INJECTIO        = -0.0001
DDPX_HER_INJECTIO       = 1.96375
DDPY_HER_INJECTIO       = -0.0004
DPX_HER_INJECTIO        = -0.120
DPY_HER_INJECTIO        = 0.0
BETS0_HER               = 8.0346

##################################### Final Element ############################################

ELE_STOP                = "injp"

#################################### Injection Parameters ######################################

X_SEPTUM_CUT            = 35.8E-3 * -1           # [m],      Septum edge position, w.r.t. HER design orbit (INJECTIO)
X_INJECTION_POINT       = 36.9E-3 * -1           # [m],      INJP x offset, w.r.t. HER design orbit (INJECTIO)
PX_INJECTION_POINT      = 2.5E-3 * -1           # [rad],    INJP px offset, w.r.t. HER design orbit (INJECTIO)
X_ORBIT_BUMP            = X_INJECTION_POINT - DX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1) -1/2*(DDX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1)**2)-BETATRON_OFFSET
PX_ORBIT_BUMP           = PX_INJECTION_POINT - DPX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1)-1/2*(DDPX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1)**2) -BETATRON_ANGLE_OFFSET

##################################### Aperture Parameters ######################################

# T. Mori:
#   Round cross-section with 47 mm diameter at drifts and Q-magnets
#   Rectangle with about 30 mm height and 50-100 mm width at the dipoles
RADIUS_BP_APERTURE      = 0.0235           # [m],      Beam pipe radius at drifts and Q-magnets
HEIGHT_BEND_APERTURE    = 0.03             # [m],      Bend aperture height
WIDTH_BEND_APERTURE     = 0.05             # [m],      Bend aperture width


######################################### Tracking Context #####################################

N_THREADS             = 10                # Number of CPU threads for tracking

######################################### Load DA Data ##########################################

da_xy   = np.load(file = da_xy_filepath, allow_pickle = True)
da_zx   = np.load(file = da_zx_filepath, allow_pickle = True)
da_zy   = np.load(file = da_zy_filepath, allow_pickle = True)

######################################### Convert Format ########################################

NORM_AX_GRID_XY = da_xy["FOLDED_NORM_A_GRID_HER"]
NORM_AY_GRID_XY = da_xy["FOLDED_NORM_B_GRID_HER"]
DA_XY           = da_xy["FOLDED_SURVIVING_HER"]
GEMITT_X_XY     = da_xy["GEMITT_X"]
GEMITT_Y_XY     = da_xy["GEMITT_Y"]
GEMITT_Z_XY     = da_xy["GEMITT_Z"]

NORM_AZ_GRID_ZX = da_zx["FOLDED_NORM_A_GRID_HER"]
NORM_AX_GRID_ZX = da_zx["FOLDED_NORM_B_GRID_HER"]
DA_ZX           = da_zx["FOLDED_SURVIVING_HER"]
GEMITT_X_ZX     = da_zx["GEMITT_X"]
GEMITT_Y_ZX     = da_zx["GEMITT_Y"]
GEMITT_Z_ZX     = da_zx["GEMITT_Z"]

NORM_AZ_GRID_ZY = da_zy["FOLDED_NORM_A_GRID_HER"]
NORM_AY_GRID_ZY = da_zy["FOLDED_NORM_B_GRID_HER"]
DA_ZY           = da_zy["FOLDED_SURVIVING_HER"]
GEMITT_X_ZY     = da_zy["GEMITT_X"]
GEMITT_Y_ZY     = da_zy["GEMITT_Y"]
GEMITT_Z_ZY     = da_zy["GEMITT_Z"]

######################################### Ensure Consistency of Emittances ######################

assert np.isclose(GEMITT_X_XY, GEMITT_X_ZX) and np.isclose(GEMITT_X_XY, GEMITT_X_ZY)
assert np.isclose(GEMITT_Y_XY, GEMITT_Y_ZX) and np.isclose(GEMITT_Y_XY, GEMITT_Y_ZY)
assert np.isclose(GEMITT_Z_XY, GEMITT_Z_ZX) and np.isclose(GEMITT_Z_XY, GEMITT_Z_ZY)
GEMITT_X_HER   = GEMITT_X_XY
GEMITT_Y_HER   = GEMITT_Y_XY
GEMITT_Z_HER   = GEMITT_Z_XY

######################################### Calculate DA Boundary Points #################################

boundary_a_xy, boundary_b_xy  = calculate_boundary_points(
    NORM_AX_GRID_XY, NORM_AY_GRID_XY, DA_XY)
boundary_a_zx, boundary_b_zx  = calculate_boundary_points(
    NORM_AZ_GRID_ZX, NORM_AX_GRID_ZX, DA_ZX)
boundary_a_zy, boundary_b_zy  = calculate_boundary_points(
    NORM_AZ_GRID_ZY, NORM_AY_GRID_ZY, DA_ZY)

######################################### Load Lattice ###########################################
######################################### Load from Environment ##################################

env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]

######################################### Adjust Reference Energy #################################

line.particle_ref.p0c[0] *= SYNCHROTRON_INJ_E_SCALE


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


######################################### Load Particles Distribution ############################

#initial_distribution    = np.load(particles_filepath, allow_pickle = True)
initial_distribution    = np.loadtxt(particles_filepath)

####### Set up bt to co bump at injection point ###########
bt_to_co_bump_rotation      = xt.YRotation(
    angle   = -1 * np.rad2deg(PX_INJECTION_POINT - PX_ORBIT_BUMP)) #angular rotation
bt_to_co_bump_translation   = xt.XYShift(
    dx      = -1 * (X_INJECTION_POINT - X_ORBIT_BUMP)) #N.B. sign convention should be negative

######################################### ECS Lag and Voltage Optimization ###########################################


##################### Create storage for optimisation points ###################

# store history of evaluations
eval_history = {
    "lag":   [],
    "volt":  [],
    "surv":  []   # survival fraction in [0,1]
}

######################################### 2D Survival Function ###################################
# survival_results = []

def compute_survival(lag, voltage, record = True):
    """
    Set env fields (lag, voltage), track a copy of initial_beam,
    compute amplitudes and check against DA frontiers.
    Returns survival fraction in [0,1].

    NOTE: If your voltage/environment variable is named differently,
    change the env key "v_acw" to the correct one.
    """
    # set environment knobs
    env["lag_acw"] = float(lag)
    env["volt_acw "]   = float(voltage)   # <-- change this key if needed
    print(f"Evaluating lag={lag:.4f}, volt={voltage:.4e}")

    initial_distribution = np.loadtxt(particles_filepath)
    initial_distribution = xt.Table({
        "name":  np.arange(initial_distribution[:, 0].shape[0]),
        "x":     initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
        "px":    initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
        "y":     initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
        "py":    initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
        "zeta":  initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
        "delta": initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
        "state": np.full(initial_distribution[:, 0].shape, 1)
    })
    initial_beam = xt.Particles(
        p0c=line.particle_ref.p0c,
        q0=line.particle_ref.q0,
        mass0=line.particle_ref.mass0,
        x=initial_distribution.x, px=initial_distribution.px,
        y=initial_distribution.y, py=initial_distribution.py,
        zeta=initial_distribution.zeta, delta=initial_distribution.delta
    )

    tracked_beam = initial_beam.copy()
        
    ################################## Enable RADIATION #########################################
    line.configure_radiation(model = "mean") #need to reset for loop to work
    
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

    ######################################## Track ##############################################
    if Radiation == True:
        line.configure_radiation(model = "quantum")
    else:
        pass

    line.track(particles = tracked_beam, ele_stop  = ELE_STOP)

    ######################################## Perform Transforms to the HER ######################

    # Want these parameters in the frame of the HER stored beam CO bump at INJECTIO / INJP

    ########################## Optimise Orbit Bump #################
    processed_delta = (line.particle_ref.p0c[0] / REFERENCE_ENERGY_HER) * \
        (1 + tracked_beam.delta[tracked_beam.state == 1]) - 1

    mean_delta_at_inj       = np.mean(processed_delta)
    print(f"Mean delta at injection:    {mean_delta_at_inj:.6e}")

    X_ORBIT_BUMP            = X_INJECTION_POINT - DX_HER_INJECTIO * (mean_delta_at_inj) -1/2*(DDX_HER_INJECTIO * (mean_delta_at_inj)**2) - BETATRON_OFFSET
    PX_ORBIT_BUMP           = PX_INJECTION_POINT - DPX_HER_INJECTIO * (mean_delta_at_inj)-1/2*(DDPX_HER_INJECTIO * (mean_delta_at_inj)**2) -BETATRON_ANGLE_OFFSET
    print(f"Optimised X_ORBIT_BUMP:     {X_ORBIT_BUMP*1E3:.3f} mm")
    print(f"Optimised PX_ORBIT_BUMP:    {PX_ORBIT_BUMP*1E3:.3f} mrad")

    ########################## Create Transforms from INJP to CO Bump at INJECTIO #################
    print(f"Mean x pre ref transform:   ", np.mean(tracked_beam.x[tracked_beam.state == 1]))
    print(f"Mean px pre ref transform:  ", np.mean(tracked_beam.px[tracked_beam.state == 1]))

    bt_to_co_bump_rotation      = xt.YRotation(
        angle   = -1 * np.rad2deg(PX_INJECTION_POINT - PX_ORBIT_BUMP)) #angular rotation
    bt_to_co_bump_translation   = xt.XYShift(
        dx      = -1 * (X_INJECTION_POINT - X_ORBIT_BUMP)) #N.B. sign convention should be negative

    ######################################## Process Results ####################################
        
    bt_to_co_bump_rotation.track(tracked_beam)
    bt_to_co_bump_translation.track(tracked_beam)

    print(np.mean(tracked_beam.px[tracked_beam.state == 1]))

        ######################################### Process Coordinates #################################

    print(f"Proportion that survived up to injection:        {np.sum(tracked_beam.state == 1) / tracked_beam.state.size:.4f}")

    processed_x     = tracked_beam.x[tracked_beam.state == 1]
    processed_px    = tracked_beam.px[tracked_beam.state == 1]
    processed_y     = tracked_beam.y[tracked_beam.state == 1]
    processed_py    = tracked_beam.py[tracked_beam.state == 1]
    processed_zeta  = tracked_beam.zeta[tracked_beam.state == 1]# type: ignore
    processed_delta = (line.particle_ref.p0c[0] / REFERENCE_ENERGY_HER) * \
        (1 + tracked_beam.delta[tracked_beam.state == 1]) - 1

        
    tracked_amplitudes  = xt.Table({
        "name":     np.arange(processed_x.shape[0]),
        "Ax":       calculate_action_amplitude(
            a           = processed_x,
            pa          = processed_px,
            delta       = processed_delta,
            beta        = BETX_HER_INJECTIO,
            alfa        = ALFX_HER_INJECTIO, #ASSUMED OF RING
            da          = DX_HER_INJECTIO,
            dpa         = DPX_HER_INJECTIO,
            gemitt_a    = GEMITT_X_HER,
            dda= DDX_HER_INJECTIO,
            ddpa= DDPX_HER_INJECTIO),
        "Ay":       calculate_action_amplitude(
            a           = processed_y,
            pa          = processed_py,
            delta       = processed_delta,
            beta        = BETY_HER_INJECTIO,
            alfa        = ALFY_HER_INJECTIO,  #ASSUMED OF RING
            da          = DY_HER_INJECTIO,
            dpa         = DPY_HER_INJECTIO,
            gemitt_a    = GEMITT_Y_HER,
            dda= DDY_HER_INJECTIO,
            ddpa= DDPY_HER_INJECTIO),
        "Az":       calculate_action_amplitude(
            a           = processed_zeta,
            pa          = processed_delta,
            delta       = processed_delta,
            beta        = BETS0_HER,
            alfa        = 0.0,
            da          = 0.0,
            dpa         = 0.0,
            gemitt_a    = GEMITT_Z_HER,
            dda= 0.0,
            ddpa= 0.0)})


    ######### RESCALE THE EMITTANCE ACCORDING TO MEASURED BLOWUP ALONG THE LINE #################

    if Radiation == True:
        results = compute_emittance_and_optics(
            x               = tracked_beam.x[tracked_beam.state == 1],
            px              = tracked_beam.px[tracked_beam.state == 1],
            y               = tracked_beam.y[tracked_beam.state == 1],
            py              = tracked_beam.py[tracked_beam.state == 1],
            zeta            = tracked_beam.zeta[tracked_beam.state == 1],
            delta           = tracked_beam.delta[tracked_beam.state == 1],
            particle_axis   = 0,)

        gemitt_x    = results["epsx"]
        gemitt_y    = results["epsy"]
        gemitt_z    = results["epsz"]

        ######################################### Compute Emittances ######################################
        EMITTANCE_X_TRACKED = gemitt_x
        EMITTANCE_Y_TRACKED = gemitt_y
        

        print(f"Tracked Emittance X at INJP:    {EMITTANCE_X_TRACKED*1E6:.3f} um")
        print(f"Tracked Emittance Y at INJP:    {EMITTANCE_Y_TRACKED*1E6:.3f} um")
        tracked_amplitudes  = xt.Table({
            "name":     np.arange(processed_x.shape[0]),
            "Ax":       tracked_amplitudes["Ax"] * np.sqrt(EMITTANCE_X_MEASURED / EMITTANCE_X_TRACKED),
            "Ay":       tracked_amplitudes["Ay"] * np.sqrt(EMITTANCE_Y_MEASURED / EMITTANCE_Y_TRACKED),
            "Az":       tracked_amplitudes["Az"]*1})
        
        print("x rescale factor: ", np.sqrt(EMITTANCE_X_MEASURED / EMITTANCE_X_TRACKED))
        print("y rescale factor: ", np.sqrt(EMITTANCE_Y_MEASURED / EMITTANCE_Y_TRACKED))

    else:
        pass
    

    ######################################### XY Test ############################################
    within_xy = compare_to_frontier(
        xs          = boundary_a_xy,
        ys          = boundary_b_xy,
        x_test      = tracked_amplitudes.Ax,
        y_test      = tracked_amplitudes.Ay)
    print(f"Proportion that survived to injection within XY DA:        {np.sum(within_xy) / within_xy.size:.4f}")

    ######################################### ZX Test ############################################

    within_zx = compare_to_frontier(
        xs          = boundary_a_zx,
        ys          = boundary_b_zx,
        x_test      = tracked_amplitudes.Az,
        y_test      = tracked_amplitudes.Ax)
    print(f"Proportion that survived to injection within ZX DA:        {np.sum(within_zx) / within_zx.size:.4f}")

    ######################################### ZY Test ############################################

    within_zy = compare_to_frontier(
        xs          = boundary_a_zy,
        ys          = boundary_b_zy,
        x_test      = tracked_amplitudes.Az,
        y_test      = tracked_amplitudes.Ay)
    print(f"Proportion that survived to injection ZY DA:        {np.sum(within_zy) / within_zy.size:.4f}")

    ##################################### Septum Test ############################################

    assert np.abs(PX_ORBIT_BUMP) < 1E-1, \
        "We assume the CO bump is parallel to the stored beam orbit and therefore the septum"

    within_septum   = processed_x < (X_SEPTUM_CUT - X_ORBIT_BUMP)
    print(f"Proportion within Septum:       {np.sum(within_septum) / within_septum.size:.4f}")

    ############################### Within all test ##############################################

    within  = within_xy * within_zx * within_zy * within_septum
    print(f"Proportion within All:          {np.sum(within) / within.size:.4f}")

    #survival_fraction = np.sum(within) / within.size
    #survival_results.append((lag, voltage, survival_fraction))
    #print(f"Lag {lag:.2f}, Volt {voltage/1e6:.2f} MV → Survival: {survival_fraction*100:.3f}%")
    print(f"Proportion within All:          {np.sum(within) / within.size:.4f}")
    surv = np.sum(within)/within.size
    print(f"Lag {lag:.2f}, Volt {voltage/1e6:.2f} MV → Survival: {surv:.3f}%")

    
    if record:
        eval_history["lag"].append(float(lag))
        eval_history["volt"].append(float(voltage))
        eval_history["surv"].append(float(surv))

    return surv
    #return -survival_fraction  # negative because we will minimize

###################### Objective wrapper for optimisers #######################

def objective_vector(x):
    """
    x: array-like [lag, voltage]
    returns: negative survival (because optimisers minimise)
    """
    lag = float(x[0])
    volt = float(x[1])
    surv = compute_survival(lag, volt, record = True)
    return float(-surv)


######################################### Run 2D Optimizer ########################################
bounds = [set_lag_bounds,  # lag [deg]
          set_voltage_bounds]  # voltage [V], example range
# === Define global tracking ===
global eval_counter, total_evals, pbar
eval_counter = 0
pbar = None

# Estimate total number of evaluations (popsize * dim * (maxiter+1))
# This is approximate since DE may stop earlier.
dim = 2
maxiter = 10
popsize = 8
total_evals = popsize * dim * (maxiter + 1)

# === Wrap your objective ===
def objective_vector_logged(x):
    global eval_counter, pbar
    val = objective_vector(x)

    # Initialize progress bar on first call
    if pbar is None:
        from tqdm.auto import tqdm as tqdm_auto
        pbar = tqdm_auto(total=total_evals, desc="Evaluating", ncols=80)

    eval_counter += 1
    pbar.update(1)

    # Optional: show occasional debug prints
    if eval_counter % 20 == 0 or eval_counter < 10:
        print(f"Eval {eval_counter:4d}: lag={x[0]:.4f}, volt={x[1]:.4e}, obj={val:.5f}")

    return val

# === Callback per generation ===
def de_callback(xk, convergence):
    print(f"\n--- Generation complete ---")
    print(f"Best so far: lag={xk[0]:.4f}, volt={xk[1]:.4e} | Convergence={convergence:.4e}")
    print("-----------------------------\n")
    return False  # return True to stop early

print(" Starting differential_evolution (2D search over lag, voltage)...\n")

de_res = differential_evolution(
    func=objective_vector_logged,
    bounds=bounds,
    strategy="best1bin",
    maxiter=maxiter,
    popsize=popsize,
    tol=1e-3,
    polish=False,
    disp=True,
    callback=de_callback,
)

# Close progress bar cleanly
if pbar is not None:
    pbar.close()
# Local refinement starting from DE result
x0 = de_res.x
print(f"DE result: lag={x0[0]:.6f}, volt={x0[1]:.6f}, surv={-de_res.fun:.6f}")

from scipy.optimize import differential_evolution, minimize

# refine with bounded L-BFGS-B (note: objective is noisy/non-smooth but this often helps)
min_res = minimize(lambda xx: objective_vector(xx), x0=x0, method="L-BFGS-B", bounds=bounds, options={"maxiter":10})

optimal = min_res.x
max_survival = -min_res.fun

print(f"Optimal Lag:            {optimal[0]:.6f} degrees")
print(f"Optimal Voltage:        {optimal[1]:.6f} (env key 'v_acw' assumed)")
print(f"Max Survival Fraction:  {max_survival * 100:.6f} %")


######### PLOT OPTIMISATION ###########

# === Extract evaluation history arrays ===
path_lags = np.array(eval_history["lag"])
path_volt = np.array(eval_history["volt"])/1e6
path_surv = np.array(eval_history["surv"]) * 100.0  # convert to %

# === Clip values to [0, 100] range ===
path_surv = np.clip(path_surv, 0, 100)

# === Figure ===
fig, ax = plt.subplots(figsize=(9, 6))

# --- Scatter plot colored by survival ---
sc = ax.scatter(
    path_lags, path_volt,
    c=path_surv,
    cmap="plasma",
    s=80,
    edgecolor="k",
    linewidth=0.1,
    alpha=0.9,
    label="Evaluated Points"
)

# --- Highlight best point ---
ax.scatter(
    [optimal[0]], [optimal[1]/1e6],
    color="limegreen",
    edgecolor="k",
    s=250,
    marker="*",
    zorder=5,
    label=f"Optimum ({optimal[0]:.1f}, {optimal[1]/1e6:.1e})"
)

# --- Colorbar ---
cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Survival Fraction [%]", fontsize=12)

# --- Labels and layout ---
ax.set_xlabel("ECS Lag [deg]", fontsize=13)
ax.set_ylabel("ECS Voltage [MV]", fontsize=13)
ax.set_title(
    f"Optimisation Landscape (Evaluated Points Only)\n"
    f"Max survival = {max_survival*100:.1f}% at lag={optimal[0]:.2f} deg, volt={optimal[1]/1e6:.1e} MV \n \n Radiation is {Radiation} and Tapering is {Tapering}",
    fontsize=14, weight="bold"
)
ax.legend(fontsize=10, loc="best", frameon=True, facecolor="white", edgecolor="gray")
ax.grid(alpha=0.25)
ax.set_facecolor("white")

plt.tight_layout()
plt.show()

