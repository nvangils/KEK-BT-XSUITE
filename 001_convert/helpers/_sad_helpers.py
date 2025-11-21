"""
(Unofficial) SAD to XSuite Converter
"""


###################################### Required Packages #####################################

import os
import subprocess
import numpy as np
import tfs
import xtrack as xt


###################################### SAD Twiss Print Function #####################################

def generate_twiss_print_function():
    """
    Generate a symmetric log-spaced array of length n_points:
      - Positive values: logspace from 10**upper_power down to 10**lower_power
      - (If n_points is odd) a zero in the center
      - Negative values: the mirror image of the positive side
    
    Args:
      lower_power (float): exponent for the smallest magnitude (linthresh), e.g. -12 → 1e-12
      upper_power (float): exponent for the largest magnitude,       e.g. -9  → 1e-9
      n_points    (int):   total length of the array
    
    Returns:
      numpy.ndarray of shape (n_points,)
    """

    TWISS_COMMAND = f"""
! -------------------------------- PRINT TWISS OF THE RING ---------------------
SaveTwissFile[filename_]:=Module[
{{fn, pos}},
fn=OpenWrite[filename];    ! Use OpenAppend[] if you do not wish to overwrite file
$FORM="12.10";
WriteString[fn, "@ ",
                StringFill["TIME"," ", 20],
                "%s ",
                "\\"",
                StringFill[DateString[]," ",-20],
                "\\"",
                "\\n"]; 
WriteString[fn, "@ ",
                StringFill["LENGTH"," ", 20],
                "%le",
                StringFill[ToString[LINE["LENG","$$$"]]," ",-22],
                "\\n"]; 
WriteString[fn, "@ ",
                StringFill["Q1"," ", 20],
                "%le",
                StringFill[ToString[Twiss["NX","$$$"]/(2*Pi)]," ",-22],
                "\\n"]; 
WriteString[fn, "@ ",
                StringFill["Q2"," ", 20],
                "%le",
                StringFill[ToString[Twiss["NY","$$$"]/(2*Pi)]," ",-22],
                "\\n"]; 
WriteString[fn, "@ ",
                StringFill["BETXMAX"," ", 20],
                "%le",
                StringFill[ToString[Max[Twiss["BX","*"]]]," ",-22],
                "\\n"]; 
WriteString[fn, "@ ",
                StringFill["BETYMAX"," ", 20],
                "%le",
                StringFill[ToString[Max[Twiss["BY","*"]]]," ",-22],
                "\\n"]; 
WriteString[fn, "* ",
                StringFill["NAME"," ", 20]," ",
                StringFill["KEYWORD"," ", -12],"    ",
                StringFill["S"," ", -12],"    ",
                StringFill["L"," ", -12],"    ",
                StringFill["BETX"," ", -12],"    ",
                StringFill["BETY"," ", -12],"    ",
                StringFill["ALFX"," ", -12],"    ",
                StringFill["ALFY"," ", -12],"    ",
                StringFill["MUX"," ", -12],"    ",
                StringFill["MUY"," ", -12],"    ",
                StringFill["DX"," ", -12],"    ",
                StringFill["DX_DC"," ", -12],"    ",
                StringFill["DY"," ", -12],"    ",
                StringFill["DY_DC"," ", -12],"    ",
                StringFill["DPX"," ", -12],"    ",
                StringFill["DPX_DC"," ", -12],"    ",
                StringFill["DPY"," ", -12],"    ",
                StringFill["DPY_DC"," ", -12],"    ",
                StringFill["X"," ", -12],"    ",
                StringFill["PX"," ", -12],"    ",
                StringFill["Y"," ", -12],"    ",
                StringFill["PY"," ", -12],"    ",
                StringFill["DZ"," ", -12],"    ",
                StringFill["DELTA"," ", -12],"    ",
                StringFill["GEO_X"," ", -12],"    ",
                StringFill["GEO_Y"," ", -12],"    ",
                StringFill["GEO_Z"," ", -12],"    ",
                StringFill["K0L"," ", -12],"    ",
                StringFill["K1L"," ", -12],"    ",
                StringFill["K2L"," ", -12],"    ",
                StringFill["BZ"," ", -12],"    ",
                StringFill["R1"," ", -12],"    ",
                StringFill["R2"," ", -12],"    ",
                StringFill["R3"," ", -12],"    ",
                StringFill["R4"," ", -12],
                "\\n"];
WriteString[fn, "$ ",
                StringFill["%s"," ", 20]," ",
                StringFill["%s"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],"    ",
                StringFill["%le"," ", -12],
                StringFill["%le"," ", -12],
                "\\n"];
pos=LINE["POSITION","*{{^$$$}}"]; ! Getting positions of elements 
Do[
    WriteString[fn,     " ",
                        StringFill[StringJoin["\\"",LINE["NAME",pos[i]],"\\""]," ", 21]," ",
                        StringFill[StringJoin["\\"",LINE["TYPENAME",pos[i]],"\\""]," ", -12],"    ",
                        LINE["LENG",pos[i]],"    ",
                        LINE["L",pos[i]],"    ",
                        Twiss["BX",pos[i]],"    ",
                        Twiss["BY",pos[i]],"    ",
                        Twiss["AX",pos[i]],"    ",
                        Twiss["AY",pos[i]],"    ",
                        Twiss["NX",pos[i]]/(2*Pi),"    ",
                        Twiss["NY",pos[i]]/(2*Pi),"    ",
                        Twiss["PEX",pos[i]],"    ",
                        Twiss["EX",pos[i]],"    ",
                        Twiss["PEY",pos[i]],"    ",
                        Twiss["EY",pos[i]],"    ",
                        Twiss["PEPX",pos[i]],"    ",
                        Twiss["EPX",pos[i]],"    ",
                        Twiss["PEPY",pos[i]],"    ",
                        Twiss["EPY",pos[i]],"    ",
                        Twiss["DX",pos[i]],"    ",
                        Twiss["DPX",pos[i]],"    ",
                        Twiss["DY",pos[i]],"    ",
                        Twiss["DPY",pos[i]],"    ",
                        Twiss["DZ",pos[i]],"    ",
                        Twiss["DDP",pos[i]],"    ",
                        LINE["GX",pos[i]],"    ",
                        LINE["GY",pos[i]],"    ",
                        LINE["GZ",pos[i]],"    ",
                        LINE["K0",pos[i]],"    ",
                        LINE["K1",pos[i]],"    ",
                        LINE["K2",pos[i]],"    ",
                        LINE["BZ",pos[i]],"    ",
                        Twiss["R1",pos[i]],"    ",
                        Twiss["R2",pos[i]],"    ",
                        Twiss["R3",pos[i]],"    ",
                        Twiss["R4",pos[i]],
                        "\\n"
                ]
    ,{{i,Length[pos]}}
    ];
Close[fn];
];
    """

    return TWISS_COMMAND


###################################### Closed Ring 4D Twiss Function #####################################

def twiss_sad(
        lattice_filename:       str,
        line_name:              str,
        method:                 str,
        closed:                 bool,
        transport_line:         bool    = False,
        reverse_element_order:  bool    = False,
        reverse_bend_direction: bool    = False,
        rf_enabled:             bool    = True,
        radiation:              bool    = False,
        rad_compensation:       bool    = False,
        rad_taper:              bool    = False,
        delta0:                 float   = 0.0,
        additional_commands:    str     = ''):
    """
    Generate a SAD command to compute the twiss parameters of a lattice.
    """

    
    ###################################### Assertions #####################################
    
    assert method in ["4d", "4D", "6d", "6D"], "Method must be '4D' or '6D'"

    
    ###################################### Configure Settings #####################################
    
    if method == '4D' or method == '4d':
        METHOD_FLAG = 'CALC4D;'
    else:
        METHOD_FLAG = 'CALC6D;'

    if closed:
        CLOSED_FLAG = 'CELL;'
    else:
        CLOSED_FLAG = 'INS;'
    
    if transport_line:
        TRANSPORT_FLAG = 'TRPT;'
    else:
        TRANSPORT_FLAG = 'NOTRPT;'

    if rf_enabled:
        RF_FLAG     = 'RFSW;'
    else:
        RF_FLAG     = 'NORFSW;'

    if radiation:
        RAD_FLAG    = """RAD;
NOFLUC;"""
    else:
        RAD_FLAG    = 'NORAD;'

    if rad_compensation:
        RADCOMP_FLAG    = 'RADCOD;'
    else:
        RADCOMP_FLAG    = 'NORADCOD;'
    if rad_taper:
        TAPER_FLAG      = 'RADTAPER;'
    else:
        TAPER_FLAG      = ''

    
    ###################################### Generate the twiss command #####################################
    
    SAD_COMMAND = f"""
FFS;

GetMAIN["./{lattice_filename}"];

USE {line_name};

{additional_commands};
DP0 = {delta0};
{RF_FLAG}
{RAD_FLAG}
{RADCOMP_FLAG}
{TAPER_FLAG}
{CLOSED_FLAG}
{TRANSPORT_FLAG}
{METHOD_FLAG}
COD;
CALC;
SAVE ALL;

{generate_twiss_print_function()}

SaveTwissFile["./temporary_sad_twiss.tfs"];

abort;
"""

    ########################################
    # Write and execute the SAD command
    ########################################
    with open("temporary_sad_twiss.sad", "w") as f:
        f.write(SAD_COMMAND)

    try:
        subprocess.run(
            ["sad", "temporary_sad_twiss.sad"],
            capture_output  = True,
            text            = True,
            timeout         = 30)
    except subprocess.TimeoutExpired:
        print("SAD Twiss timed out at 30s")
        print("This typically means the process failed, often due to a filepath issue.")
        print("Ensure that there are no spaces or special characters in the full path to the lattice file.")

        os.remove("temporary_sad_twiss.sad")
        raise subprocess.TimeoutExpired(cmd = ["sad", "temporary_sad_twiss.sad"], timeout = 30)

    ########################################
    # Read the data
    ########################################
    sad_twiss   = tfs.read("temporary_sad_twiss.tfs")

    ########################################
    # Remove temporary files
    ########################################
    # os.remove("temporary_sad_twiss.sad")
    # os.remove("temporary_sad_twiss.tfs")

    ########################################
    # Convert to TwissTable
    ########################################
    s_idx       = np.argsort(np.array(sad_twiss['S']), kind = "stable")
    tw_sad      = xt.TwissTable({
        "name":     np.array(sad_twiss['NAME'])[s_idx],
        "s":        np.array(sad_twiss['S'])[s_idx],
        "x":        np.array(sad_twiss['X'])[s_idx],
        "px":       np.array(sad_twiss['PX'])[s_idx],
        "y":        np.array(sad_twiss['Y'])[s_idx],
        "py":       np.array(sad_twiss['PY'])[s_idx],
        "zeta":     np.array(sad_twiss['DZ'])[s_idx],
        "delta":    np.array(sad_twiss['DELTA'])[s_idx],
        "betx":     np.array(sad_twiss['BETX'])[s_idx],
        "bety":     np.array(sad_twiss['BETY'])[s_idx],
        "alfx":     np.array(sad_twiss['ALFX'])[s_idx],
        "alfy":     np.array(sad_twiss['ALFY'])[s_idx],
        "dx":       np.array(sad_twiss['DX'])[s_idx],
        "dpx":      np.array(sad_twiss['DPX'])[s_idx],
        "dy":       np.array(sad_twiss['DY'])[s_idx],
        "dpy":      np.array(sad_twiss['DPY'])[s_idx],
        "mux":      np.array(sad_twiss['MUX'])[s_idx],
        "muy":      np.array(sad_twiss['MUY'])[s_idx],
        "R1":       np.array(sad_twiss['R1'])[s_idx],
        "R2":       np.array(sad_twiss['R2'])[s_idx],
        "R3":       np.array(sad_twiss['R3'])[s_idx],
        "R4":       np.array(sad_twiss['R4'])[s_idx]})
    tw_sad['qx']            = sad_twiss['Q1']
    tw_sad['qy']            = sad_twiss['Q2']
    tw_sad['circumference'] = sad_twiss['LENGTH']
    
    ########################################
    # Element Order Reversal
    ########################################
    if reverse_element_order:
    
        tw_sad.s        = tw_sad.s[-1] - tw_sad.s
        tw_sad.x        *= +1
        tw_sad.px       *= +1 * -1
        tw_sad.y        *= +1
        tw_sad.py       *= +1 * -1
        tw_sad.zeta     = tw_sad.zeta[-1] - tw_sad.zeta
        tw_sad.delta    *= +1
        tw_sad.betx     *= +1
        tw_sad.bety     *= +1
        tw_sad.alfx     *= +1 * -1
        tw_sad.alfy     *= +1 * -1
        tw_sad.dx       *= +1
        tw_sad.dpx      *= +1 * -1
        tw_sad.dy       *= +1
        tw_sad.dpy      *= +1 * -1
        tw_sad.mux      = tw_sad.mux[-1] - tw_sad.mux
        tw_sad.muy      = tw_sad.muy[-1] - tw_sad.muy
        tw_sad.R1       *= +1
        tw_sad.R2       *= +1
        tw_sad.R3       *= +1
        tw_sad.R4       *= +1

        tw_sad.name     = np.flip(tw_sad.name)
        tw_sad.s        = np.flip(tw_sad.s)
        tw_sad.x        = np.flip(tw_sad.x)
        tw_sad.px       = np.flip(tw_sad.px)
        tw_sad.y        = np.flip(tw_sad.y)
        tw_sad.py       = np.flip(tw_sad.py)
        tw_sad.zeta     = np.flip(tw_sad.zeta)
        tw_sad.delta    = np.flip(tw_sad.delta)
        tw_sad.betx     = np.flip(tw_sad.betx)
        tw_sad.bety     = np.flip(tw_sad.bety)
        tw_sad.alfx     = np.flip(tw_sad.alfx)
        tw_sad.alfy     = np.flip(tw_sad.alfy)
        tw_sad.dx       = np.flip(tw_sad.dx)
        tw_sad.dpx      = np.flip(tw_sad.dpx)
        tw_sad.dy       = np.flip(tw_sad.dy)
        tw_sad.dpy      = np.flip(tw_sad.dpy)
        tw_sad.mux      = np.flip(tw_sad.mux)
        tw_sad.muy      = np.flip(tw_sad.muy)
        tw_sad.R1       = np.flip(tw_sad.R1)
        tw_sad.R2       = np.flip(tw_sad.R2)
        tw_sad.R3       = np.flip(tw_sad.R3)
        tw_sad.R4       = np.flip(tw_sad.R4)

    ########################################
    # Bend Direction Reversal
    ########################################
    if reverse_bend_direction:
    
        tw_sad.x        *= -1
        tw_sad.px       *= -1
        tw_sad.y        *= +1
        tw_sad.py       *= +1
        tw_sad.zeta     *= +1
        tw_sad.delta    *= +1
        tw_sad.betx     *= +1
        tw_sad.bety     *= +1
        tw_sad.alfx     *= +1
        tw_sad.alfy     *= +1
        tw_sad.dx       *= -1
        tw_sad.dpx      *= -1
        tw_sad.dy       *= +1
        tw_sad.dpy      *= +1
        tw_sad.mux      *= +1
        tw_sad.muy      *= +1
        tw_sad.R1       *= +1
        tw_sad.R2       *= +1
        tw_sad.R3       *= +1
        tw_sad.R4       *= +1

    ########################################
    # Return the TwissTable
    ########################################
    return tw_sad

################################################################################
# Rebuild SAD lattice (post GEO for solenoid)
################################################################################
def rebuild_sad_lattice(
        lattice_filename:       str = 'sad_lattice.sad',
        line_name:              str = 'TEST_LINE',
        additional_commands:    str = '',
        output_filename:        str | None = None):
    """
    Generate a SAD command to compute the twiss parameters of a lattice.
    """
    
    ########################################
    # Generate the twiss command
    ########################################
    if output_filename is None:
        output_filename = lattice_filename

    SAD_COMMAND = f"""
FFS;

GetMAIN["./{lattice_filename}"];  (* Input your lattice file before modification *)

USE {line_name};

{additional_commands};

INS;    ! Transfer Line
CALC;   ! Twiss
SAVE ALL;

of=OpenWrite["./{output_filename}"];
WriteString[of, "MOMENTUM = "//MOMENTUM//";\\n"];
WriteString[of, "FSHIFT = 0;\\n"];
FFS["output "//of//" type"];                     (* Write element definition *)
WriteBeamLine[of, ExtractBeamLine[], Format->"MAIN", Name->{{"{line_name}"}}];  (* Write lattice order *)
Close[of];

abort;
"""

    ########################################
    # Write and execute the SAD command
    ########################################
    with open("temporary_sad_twiss.sad", "w") as f:
        f.write(SAD_COMMAND)

    try:
        subprocess.run(
            ["sad", "temporary_sad_twiss.sad"],
            capture_output  = True,
            text            = True,
            timeout         = 30)
    except subprocess.TimeoutExpired:
        print("SAD Twiss timed out at 30s")
        print("This typically means the process failed, often due to a filepath issue.")
        print("Ensure that there are no spaces or special characters in the full path to the lattice file.")

        os.remove("temporary_sad_twiss.sad")
        raise subprocess.TimeoutExpired(cmd = ["sad", "temporary_sad_twiss.sad"], timeout = 30)

    ########################################
    # Remove temporary files
    ########################################
    os.remove("temporary_sad_twiss.sad")
