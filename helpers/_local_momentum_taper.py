"""
Taper to local momentum deviation (compensate for acceleration)
=============================================

"""

################################################################################
# Required Packages
################################################################################
import xtrack as xt
import xobjects as xo

import numpy as np
from scipy.constants import c as clight
import matplotlib.pyplot as plt

################################################################################
# Main Function
################################################################################
def taper_magnets_to_local_momentum(
        line:           xt.Line,
        ele_start:      str         = None,
        ele_stop:       str         = None,
        delta0:         float       = 0.0,
        atol_delta:     float       = 1E-9,
        max_iter:       int         = 100,
        **kwargs):

    ########################################
    # Assertions
    ########################################
    # Check for a valid tracker
    assert line._has_valid_tracker(), \
        "Line does not have a valid tracker. Please build the tracker first."

    # Check that the tracker is CPU context
    assert isinstance(line._context, xo.ContextCpu), \
        "Only CPU context is supported"
    
    # Check that there is a particle reference with |q0| = 1
    assert line.particle_ref is not None, \
        "Particle reference is not set"
    assert np.abs(line.particle_ref.q0) == 1, \
        "Only |q0| = 1 is supported (for now)"

    # Enforce no repeated elements
    assert len(set(line.element_names)) == len(line.element_names), \
        "Line must not contain repeated elements to use `taper_to_local_momentum(...)`. "

    # Enforce no Cavity Slices
    assert "SliceCavity" not in list(set(line.tracker._tracker_data_base._line_table.element_type)), \
        "Element type 'SliceCavity' is not supported for local momentum tapering."

    ########################################
    # Get the reference particle
    ########################################
    reference_particle  = line.particle_ref.copy()

    ########################################
    # Test if compensation is needed
    ########################################
    test_particle       = reference_particle.copy()
    test_particle.delta = line.attr._cache['delta_taper'].multisetter.get_values()[0]
    line.track(
        particles               = test_particle,
        turn_by_turn_monitor    = 'ONE_TURN_EBE')
    moni                = line.record_last_track
    reference_delta     = moni.delta[0, :]
    if test_particle.state[0] > 0 and np.all(abs(reference_delta) < atol_delta):
        return

    ########################################
    # Tapering Objects
    ########################################
    delta_taper_setter  = line.attr._cache['delta_taper'].multisetter

    ########################################
    # Elements to taper
    ########################################
    delta_taper0        = delta_taper_setter.get_values()
    mask_delta_taper    = line.attr._cache['delta_taper'].mask

    ########################################
    # Converge
    ########################################
    with xt.line._preserve_track_flags(line):

        i_iter = 0
        delta_taper = delta_taper0.copy()

        while True:
            print("Iteration:", i_iter)

            test_particle       = reference_particle.copy()
            test_particle.delta = delta0
            line.track(
                particles               = test_particle,
                ele_start               = ele_start,
                ele_stop                = ele_stop,
                turn_by_turn_monitor    = 'ONE_TURN_EBE')
            moni                = line.record_last_track
            
            delta       = moni.delta[0, :-1]
            delta_taper = delta_taper_setter.get_values()
            
            # fig = plt.figure(figsize = (10, 6))
            # plt.plot(delta[mask_delta_taper], c = "r", label = "Tracked Delta")
            # plt.plot(delta_taper, c = "b", label = "Tapering")
            # plt.legend()
            # fig = plt.figure(figsize = (10, 6))
            # plt.plot(abs(delta[mask_delta_taper] - delta_taper))
            # plt.show()
            
            print("Test particle state: ", test_particle.state[0])
            print("Max Delta Difference: ", np.max(abs(delta[mask_delta_taper] - delta_taper)))
            
            if test_particle.state[0] != 1:
                raise RuntimeError("Particle lost during tapering process")

            if test_particle.state[0] > 0 and \
                    np.all(abs(delta[mask_delta_taper] - delta_taper) < atol_delta):
                break
            delta_taper_setter.set_values(delta[mask_delta_taper])

            i_iter += 1
            if i_iter > max_iter:
                raise RuntimeError("Maximum number of iterations reached")
