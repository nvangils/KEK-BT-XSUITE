"""
Convert the KEK BTp 6D lattice comparison 
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       27-10-2025
"""

#################################### Required Packages ###################################

import sad2xs as s2x
import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

from helpers._sad_helpers import twiss_sad
from helpers._misc_helpers import create_comparison_plots

#################################### User Parameters ####################################


SAD_LATTICE_PATH    = '../lattices/bte.sad'
OUTPUT_DIRECTORY    = '../lattices'
OUTPUT_FILENAME     = 'bte'
LINE_NAME           = 'BTE'

#################################### Load Reference Data #########################


tw6d_sad    = twiss_sad(
    lattice_filename        = SAD_LATTICE_PATH,
    line_name               = LINE_NAME,
    method                  = "6d",
    closed                  = False,
    transport_line          = True,
    reverse_element_order   = False,
    reverse_bend_direction  = False,
    additional_commands     = "")

##################################### Convert Lattice ####################################


line    = s2x.convert_sad_to_xsuite(
    sad_lattice_path            = SAD_LATTICE_PATH,
    line_name                   = LINE_NAME,
    excluded_elements           = [],
    user_multipole_replacements = None,
    reverse_element_order       = False,
    reverse_bend_direction      = False,
    reverse_charge              = True,
    output_directory            = OUTPUT_DIRECTORY,
    output_filename             = OUTPUT_FILENAME,
    output_header               = "SuperKEKB BTE")
line.replace_all_repeated_elements()

######################################### Install ref energy shifts ########################################

######################################### Get Cavity Table ########################################

env     = line.env
tt      = line.get_table(attr = True)
tt_cav  = tt.rows[tt.element_type == "Cavity"]

######################################### Create Installation Lists ########################################

ref_energy_increases    = []
for cav in tt_cav.name:
    env.new(
        name        = f"{cav}_dE",
        parent      = "ReferenceEnergyIncrease",
        Delta_p0c   = 0)
    ref_energy_increases.append(env.place(
        name        = f"{cav}_dE",
        at          = 0,
        from_       = cav,
        from_anchor = "end"))

#################################### Insert RES and ZetaShift Elements ########################################

line.insert(ref_energy_increases, s_tol = 1E-6)
######################################### Get updated tables ########################################

tt                  = line.get_table(attr = True)
tt_cav              = tt.rows[tt.element_type == "Cavity"]
tt_ref_energy_shift = tt.rows[tt.element_type == "ReferenceEnergyIncrease"]

######################################### Tune ref energy shifts ########################################

particle_energy = line.particle_ref.p0c
for ref_energy_shift in tt_ref_energy_shift.name:

    tw_test = line.twiss(
    method              = "6d",
    start               = xt.START,
    end                 = xt.END,
    betx                = tw6d_sad.betx[0],
    bety                = tw6d_sad.bety[0],
    alfx                = tw6d_sad.alfx[0],
    alfy                = tw6d_sad.alfy[0])

    line[ref_energy_shift].Delta_p0c    = particle_energy * tw_test["delta", ref_energy_shift]
    particle_energy = particle_energy * (1 + tw_test["delta", ref_energy_shift])

ref_increase    = 0
for ref_energy_shift in tt_ref_energy_shift.name:
    ref_increase += line[ref_energy_shift].Delta_p0c
print(f"Total Cavity Voltage:       {np.sum(tt_cav.voltage) / 1E9} GV")
print(f"Reference Energy Increase:  {ref_increase / 1E9} GeV")

######################################### Twiss for Comparison ########################################

######################################### Get table ########################################

tt = line.get_table(attr = True)

######################################### Twiss ########################################

tw6d    = line.twiss(
    start               = xt.START,
    end                 = xt.END,
    betx                = tw6d_sad.betx[0],
    bety                = tw6d_sad.bety[0],
    alfx                = tw6d_sad.alfx[0],
    alfy                = tw6d_sad.alfy[0])

######################################### Comparisons ########################################

for tw, tw_sad, mode in zip([ tw6d], [tw6d_sad], [ "6D"]):
    create_comparison_plots(tw, tw_sad, suptitle = mode, zero_tol = 1E-10)

######################################### Show plots ########################################

plt.show()
