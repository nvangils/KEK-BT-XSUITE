"""
Convert the KEK BTe 4D lattice comparison
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       15-10-2025
"""

#################################### Required Packages ###################################


import sad2xs as s2x
import xtrack as xt
import matplotlib.pyplot as plt

from helpers._sad_helpers import twiss_sad
from helpers._misc_helpers import create_comparison_plots

#################################### User Parameters ####################################


SAD_LATTICE_PATH    = '../lattices/bte.sad'
OUTPUT_DIRECTORY    = '../lattices'
OUTPUT_FILENAME     = 'bte'
LINE_NAME           = 'BTE'

#################################### Load Reference Data #########################


tw4d_sad    = twiss_sad(
    lattice_filename        = SAD_LATTICE_PATH,
    line_name               = LINE_NAME,
    method                  = "4d",
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

#################################### Get table ########################################

tt = line.get_table(attr = True)

##################################### Twiss ###########################################

tw4d    = line.twiss4d(
    start               = xt.START,
    end                 = xt.END,
    betx                = tw4d_sad.betx[0],
    bety                = tw4d_sad.bety[0],
    alfx                = tw4d_sad.alfx[0],
    alfy                = tw4d_sad.alfy[0])


##################################### Comparisons #####################################


for tw, tw_sad, mode in zip([tw4d], [tw4d_sad], ["4D"]):
    create_comparison_plots(tw, tw_sad, suptitle = mode, zero_tol = 1E-10)

##################################### Show plots ######################################

plt.show()
