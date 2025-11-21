"""
SuperKEKB BTE
================================================================================
Converted using the SAD2XS Converter

"""

################################################################################
# Import Packages
################################################################################
import xtrack as xt
import numpy as np

################################################################################
# Create or Get Environment
################################################################################
env = xt.get_environment(verbose = True)
env.vars.default_to_zero = True

########################################
# Key Global Variables
########################################
env["mass0"]    = 510998.95
env["p0c"]      = 7000008020.508931
env["q0"]       = 1.0
env["fshift"]   = 0.0

########################################
# Reference Particle
########################################
env.particle_ref    = xt.Particles(
    mass0   = env["mass0"],
    p0c     = env["p0c"],
    q0      = env["q0"])

################################################################################
# Import lattice
################################################################################

############################################################
# Drifts
############################################################
env.new(name = 'ldmp1a', parent = xt.Drift, length = 1.782975)
env.new(name = 'lmblh', parent = xt.Drift, length = 0.125)
env.new(name = 'ldmp1b', parent = xt.Drift, length = 0.6300695382524946)
env.new(name = 'lpmh', parent = xt.Drift, length = 0.1)
env.new(name = 'lblcc', parent = xt.Drift, length = 0.22389)
env.new(name = 'lce0a', parent = xt.Drift, length = 2.13389)
env.new(name = 'lqmd', parent = xt.Drift, length = 0.10735)
env.new(name = 'lqz', parent = xt.Drift, length = 0.18735)
env.new(name = 'lce01a', parent = xt.Drift, length = 0.0743464892023)
env.new(name = 'ld300', parent = xt.Drift, length = 0.3)
env.new(name = 'lce1a', parent = xt.Drift, length = 1.0628798402556)
env.new(name = 'lbell', parent = xt.Drift, length = 0.09)
env.new(name = 'lblcd', parent = xt.Drift, length = 0.21625)
env.new(name = 'lce2a1', parent = xt.Drift, length = 0.36625)
env.new(name = 'ld260', parent = xt.Drift, length = 0.26)
env.new(name = 'lqmc', parent = xt.Drift, length = 0.1069)
env.new(name = 'lbvqxa', parent = xt.Drift, length = 0.2772)
env.new(name = 'lblc', parent = xt.Drift, length = 0.1872)
env.new(name = 'ld500', parent = xt.Drift, length = 0.5)
env.new(name = 'lbva', parent = xt.Drift, length = 1.3722)
env.new(name = 'lbqx2', parent = xt.Drift, length = 0.4341)
env.new(name = 'lce21', parent = xt.Drift, length = 1.4922883866788)
env.new(name = 'lpmlh', parent = xt.Drift, length = 0.125)
env.new(name = 'ld350', parent = xt.Drift, length = 0.35)
env.new(name = 'lmz', parent = xt.Drift, length = 0.1)
env.new(name = 'lbqx0', parent = xt.Drift, length = 0.545325)
env.new(name = 'lce3a', parent = xt.Drift, length = 0.5388768489311)
env.new(name = 'lbqx', parent = xt.Drift, length = 0.43455)
env.new(name = 'lce2b1', parent = xt.Drift, length = 0.2772)
env.new(name = 'ld210', parent = xt.Drift, length = 0.21)
env.new(name = 'ld80', parent = xt.Drift, length = 0.08)
env.new(name = 'lchbh', parent = xt.Drift, length = 0.125)
env.new(name = 'lce5a', parent = xt.Drift, length = 2.9996938123642)
env.new(name = 'lcbbh', parent = xt.Drift, length = 0.125)
env.new(name = 'lmblh2', parent = xt.Drift, length = 0.6659979999999219)
env.new(name = 'lce5b', parent = xt.Drift, length = 2.820695812364244)
env.new(name = 'ld150', parent = xt.Drift, length = 0.15)
env.new(name = 'lce6a', parent = xt.Drift, length = 3.7041498113021)
env.new(name = 'lblce', parent = xt.Drift, length = 0.223425)
env.new(name = 'lbqx1a', parent = xt.Drift, length = 0.383425)
env.new(name = 'lmbsh', parent = xt.Drift, length = 0.1)
env.new(name = 'lbqx1b', parent = xt.Drift, length = 0.293425)
env.new(name = 'l4x', parent = xt.Drift, length = 1.2673062892642)
env.new(name = 'l4y1', parent = xt.Drift, length = 0.4016546983106)
env.new(name = 'lgv', parent = xt.Drift, length = 0.075)
env.new(name = 'ld100', parent = xt.Drift, length = 0.1)
env.new(name = 'l41', parent = xt.Drift, length = 4.7113996875748)
env.new(name = 'l42', parent = xt.Drift, length = 4.0213996875748)
env.new(name = 'ld250', parent = xt.Drift, length = 0.25)
env.new(name = 'l43', parent = xt.Drift, length = 4.7113996875748)
env.new(name = 'l44', parent = xt.Drift, length = 2.1739609875748)
env.new(name = 'lbqx14', parent = xt.Drift, length = 0.578425)
env.new(name = 'lqzc', parent = xt.Drift, length = 0.1869)
env.new(name = 'l45', parent = xt.Drift, length = 0.518425)
env.new(name = 'l4b', parent = xt.Drift, length = 0.76685)
env.new(name = 'lbqxc', parent = xt.Drift, length = 0.470775)
env.new(name = 'l46', parent = xt.Drift, length = 5.6266)
env.new(name = 'l47bell', parent = xt.Drift, length = 0.0883)
env.new(name = 'l47d', parent = xt.Drift, length = 2.117)
env.new(name = 'l48d1', parent = xt.Drift, length = 1.3765)
env.new(name = 'l48d2', parent = xt.Drift, length = 1.3253)
env.new(name = 'lw4d1', parent = xt.Drift, length = 0.1686)
env.new(name = 'lw4d2', parent = xt.Drift, length = 0.1948)
env.new(name = 'lw4d3', parent = xt.Drift, length = 6.1435516691721)
env.new(name = 'lmwh', parent = xt.Drift, length = 0.125)
env.new(name = 'lw2', parent = xt.Drift, length = 12.5729516691721)
env.new(name = 'lw1', parent = xt.Drift, length = 12.5729516691721)
env.new(name = 'lmbth', parent = xt.Drift, length = 0.125)
env.new(name = 'lw3', parent = xt.Drift, length = 12.5729516691721)
env.new(name = 'lwa1', parent = xt.Drift, length = 0.2872159127557)
env.new(name = 'lbs1a', parent = xt.Drift, length = 5.6694027571257)
env.new(name = 'ld310', parent = xt.Drift, length = 0.271095)
env.new(name = 'lx1', parent = xt.Drift, length = 0.328445)
env.new(name = 'lj1', parent = xt.Drift, length = 0.2378687700911)
env.new(name = 'lbb1', parent = xt.Drift, length = 0.27219)
env.new(name = 'lbq13a', parent = xt.Drift, length = 0.2259501292978)
env.new(name = 'lchbsh', parent = xt.Drift, length = 0.025)
env.new(name = 'ld100c', parent = xt.Drift, length = 0.0869)
env.new(name = 'lqmac', parent = xt.Drift, length = 0.1369)
env.new(name = 'lbq12a', parent = xt.Drift, length = 0.2262670634871)
env.new(name = 'lbq11a', parent = xt.Drift, length = 0.2461116432598)
env.new(name = 'lqma', parent = xt.Drift, length = 0.13735)
env.new(name = 'lbq13', parent = xt.Drift, length = 0.3628501292978)
env.new(name = 'lk1', parent = xt.Drift, length = 4.5173215836445)
env.new(name = 'ls1', parent = xt.Drift, length = 0.448445)
env.new(name = 'ls21', parent = xt.Drift, length = 1.39)
env.new(name = 'ls2a', parent = xt.Drift, length = 3.7210631900934)
env.new(name = 'ls3', parent = xt.Drift, length = 7.0693656873955)
env.new(name = 'ls4a', parent = xt.Drift, length = 6.0246656873932)
env.new(name = 'ls5a', parent = xt.Drift, length = 0.9527058505509)
env.new(name = 'ld320', parent = xt.Drift, length = 0.278635)
env.new(name = 'lbb20', parent = xt.Drift, length = 0.26727)
env.new(name = 'lbq23', parent = xt.Drift, length = 0.385535)
env.new(name = 'lbq22a', parent = xt.Drift, length = 0.228635)
env.new(name = 'lbb2', parent = xt.Drift, length = 0.31727)
env.new(name = 'lbq21a', parent = xt.Drift, length = 0.5682251651292)
env.new(name = 'lbb2c', parent = xt.Drift, length = 0.313025)
env.new(name = 'lx7a', parent = xt.Drift, length = 3.6698602160724)
env.new(name = 'lbq33', parent = xt.Drift, length = 0.8070274261686)
env.new(name = 'lbq32a', parent = xt.Drift, length = 0.4501274261686)
env.new(name = 'lbb3', parent = xt.Drift, length = 0.3070788105584)
env.new(name = 'lbq31', parent = xt.Drift, length = 0.3958097869185)
env.new(name = 'lbq3a', parent = xt.Drift, length = 0.4258097869185)
env.new(name = 'lbq32b', parent = xt.Drift, length = 0.1101274261686)
env.new(name = 'lbq34a', parent = xt.Drift, length = 0.3758097869184)
env.new(name = 'lbq3', parent = xt.Drift, length = 0.7831597869185)
env.new(name = 'lm1a1', parent = xt.Drift, length = 0.2589079197818)
env.new(name = 'lmbbv2', parent = xt.Drift, length = 0.17500875)
env.new(name = 'lm1x1', parent = xt.Drift, length = 1.584114981163)
env.new(name = 'lmqqv1', parent = xt.Drift, length = 5.7042204866495)
env.new(name = 'lmqqv2', parent = xt.Drift, length = 5.7042204866495)
env.new(name = 'lmqqv3', parent = xt.Drift, length = 5.7042204866495)
env.new(name = 'lm1y1', parent = xt.Drift, length = 0.2768957301224)
env.new(name = 'lm2a', parent = xt.Drift, length = 0.84847836137)
env.new(name = 'lmqqh1', parent = xt.Drift, length = 2.8784799140679)
env.new(name = 'lcvbh', parent = xt.Drift, length = 0.125)
env.new(name = 'lmqqh2', parent = xt.Drift, length = 6.8623051133512)
env.new(name = 'lm2ba', parent = xt.Drift, length = 0.34889)
env.new(name = 'lmbbh', parent = xt.Drift, length = 0.386485)
env.new(name = 'lmbbh1', parent = xt.Drift, length = 0.66778)
env.new(name = 'lbq411', parent = xt.Drift, length = 0.38389)
env.new(name = 'lbq421', parent = xt.Drift, length = 0.3325948)
env.new(name = 'lbq43a1', parent = xt.Drift, length = 0.3202424)
env.new(name = 'lbq43a2', parent = xt.Drift, length = 0.048)
env.new(name = 'lm4aa', parent = xt.Drift, length = 0.34889)
env.new(name = 'lm4b1', parent = xt.Drift, length = 0.649185169638)
env.new(name = 'lmbt1h', parent = xt.Drift, length = 0.15)
env.new(name = 'lm4c1', parent = xt.Drift, length = 0.5288153937651)
env.new(name = 'lm4da1', parent = xt.Drift, length = 0.32089)
env.new(name = 'lm4da2', parent = xt.Drift, length = 0.048)
env.new(name = 'ld090', parent = xt.Drift, length = 0.09)
env.new(name = 'lm4ea', parent = xt.Drift, length = 0.4478524)
env.new(name = 'lmbbhc', parent = xt.Drift, length = 0.385705)
env.new(name = 'lm4f1', parent = xt.Drift, length = 0.3028524)
env.new(name = 'lmi1', parent = xt.Drift, length = 0.006994006715600007)
env.new(name = 'lmi2', parent = xt.Drift, length = 0.2204038939567)
env.new(name = 'lmi3', parent = xt.Drift, length = 0.2491443516818)
env.new(name = 'lse3', parent = xt.Drift, length = 0.42)
env.new(name = 'lse2', parent = xt.Drift, length = 0.42)
env.new(name = 'lse1', parent = xt.Drift, length = 0.42)
env.new(name = 'lfr049', parent = xt.Drift, length = 0.3069249999999084)
env.new(name = 'lfr050', parent = xt.Drift, length = 0.186925)
env.new(name = 'lfr051', parent = xt.Drift, length = 0.4828000000001691)
env.new(name = 'lfr052', parent = xt.Drift, length = 4.564894935148431)
env.new(name = 'lfr057', parent = xt.Drift, length = 0.192385)
env.new(name = 'lfr058', parent = xt.Drift, length = 0.120185)
env.new(name = 'lfr083', parent = xt.Drift, length = 4.493300000000057)
env.new(name = 'lfr059', parent = xt.Drift, length = 4.089707703966105)
env.new(name = 'lfr081', parent = xt.Drift, length = 0.2949999999999547)
env.new(name = 'lfr082', parent = xt.Drift, length = 0.29500000000000004)

############################################################
# Bends
############################################################

########################################
# Base Elements
########################################
env.new(name = 'hbend01000058596', parent = xt.Bend, length = 1.0000585961533)
env.new(name = 'hbend01000234413', parent = xt.Bend, length = 1.0002344134579)
env.new(name = 'hbend01536104245', parent = xt.Bend, length = 1.5361042457936593)
env.new(name = 'hbend01549187396', parent = xt.Bend, length = 1.5491873961201)
env.new(name = 'hbend01812220000', parent = xt.Bend, length = 1.81222)
env.new(name = 'hbend01812604425', parent = xt.Bend, length = 1.8126044257406)
env.new(name = 'hbend01813369682', parent = xt.Bend, length = 1.813369682202)
env.new(name = 'hbend02034013700', parent = xt.Bend, length = 2.0340137)
env.new(name = 'hbend02078583912', parent = xt.Bend, length = 2.078583912105)
env.new(name = 'hbend02083423192', parent = xt.Bend, length = 2.0834231922962)
env.new(name = 'hbend02085631902', parent = xt.Bend, length = 2.0856319020778)
env.new(name = 'hbend02292170390', parent = xt.Bend, length = 2.2921703908216)
env.new(name = 'vbend01906011158', parent = xt.Bend, length = 1.9060111581781, rot_s_rad = +np.pi/2)
env.new(name = 'vbend01906019066', parent = xt.Bend, length = 1.9060190663188, rot_s_rad = +np.pi/2)

########################################
# Cloned Elements
########################################
env.new(
    name                    = 'se4',
    parent                  = 'hbend01000058596',
    k0                      = 'k0_se4',
    h                       = -0.0374978027730003,
    edge_entry_angle        = -0.01875,
    edge_exit_angle         = -0.01875)
env.new(
    name                    = 'se3',
    parent                  = 'hbend01000058596',
    k0                      = 'k0_se3',
    h                       = -0.0374978027730003,
    edge_entry_angle        = -0.01875,
    edge_exit_angle         = -0.01875)
env.new(
    name                    = 'se2',
    parent                  = 'hbend01000058596',
    k0                      = 'k0_se2',
    h                       = -0.03534792874734828,
    edge_entry_angle        = -0.017675,
    edge_exit_angle         = -0.017675)
env.new(
    name                    = 'se1',
    parent                  = 'hbend01000234413',
    k0                      = 'k0_se1',
    h                       = -0.03948074518200142,
    edge_entry_angle        = -0.03949)
env.new(
    name                    = 'bc1p',
    parent                  = 'hbend01536104245',
    k0                      = 'k0_bc1p',
    h                       = -0.061120770523621186,
    edge_exit_angle         = -0.09388787510751445)
env.new(
    name                    = 'b0e',
    parent                  = 'hbend01549187396',
    k0                      = 'k0_b0e',
    h                       = 0.04816754779637654,
    edge_entry_angle        = 0.037310278974079514,
    edge_exit_angle         = 0.037310278974079514)
env.new(
    name                    = 'bhd0',
    parent                  = 'hbend01812220000',
    k0                      = 'k0_bhd0',
    h                       = 0.010631886393128528,
    edge_exit_angle         = 0.019267317159355382)
env.new(
    name                    = 'bh4e',
    parent                  = 'hbend01812604425',
    k0                      = 'k0_bh4e',
    h                       = -0.039636133211219254,
    edge_entry_angle        = -0.03592231523895,
    edge_exit_angle         = -0.03592231523895)
env.new(
    name                    = 'bh4ae',
    parent                  = 'hbend01813369682',
    k0                      = 'k0_bh4ae',
    h                       = -0.03950045406920005,
    edge_entry_angle        = -0.03581446292115,
    edge_exit_angle         = -0.03581446292115)
env.new(
    name                    = 'b1e',
    parent                  = 'hbend02034013700',
    k0                      = 'k0_b1e',
    h                       = -0.05004507680552988,
    edge_entry_angle        = -0.05089618592,
    edge_exit_angle         = -0.05089618592)
env.new(
    name                    = 'bh1e',
    parent                  = 'hbend02078583912',
    k0                      = 'k0_bh1e',
    h                       = 0.046356434033408694,
    edge_entry_angle        = 0.0481778690022,
    edge_exit_angle         = 0.0481778690022)
env.new(
    name                    = 'bh2e',
    parent                  = 'hbend02083423192',
    k0                      = 'k0_bh2e',
    h                       = -0.04377206742903278,
    edge_entry_angle        = -0.0455978702282,
    edge_exit_angle         = -0.0455978702282)
env.new(
    name                    = 'bh1ae',
    parent                  = 'hbend02085631902',
    k0                      = 'k0_bh1ae',
    h                       = 0.03451539039131686,
    edge_entry_angle        = 0.0359931996564,
    edge_exit_angle         = 0.0359931996564)
env.new(
    name                    = 'bh3e',
    parent                  = 'hbend02292170390',
    k0                      = 'k0_bh3e',
    h                       = 0.04440017127423098,
    edge_entry_angle        = 0.0508863789711,
    edge_exit_angle         = 0.0508863789711)
env.new(
    name                    = 'bv1ue',
    parent                  = 'vbend01906011158',
    k0                      = 'k0_bv1ue',
    h                       = -0.038839136587668785,
    edge_entry_angle        = -0.03701391385505,
    edge_exit_angle         = -0.03701391385505)
env.new(
    name                    = 'bv1de',
    parent                  = 'vbend01906011158',
    k0                      = 'k0_bv1de',
    h                       = 0.038839136587668785,
    edge_entry_angle        = 0.03701391385505,
    edge_exit_angle         = 0.03701391385505)
env.new(
    name                    = 'bv2ue',
    parent                  = 'vbend01906019066',
    k0                      = 'k0_bv2ue',
    h                       = -0.0392146299182919,
    edge_entry_angle        = -0.03737191615145,
    edge_exit_angle         = -0.03737191615145)
env.new(
    name                    = 'bv2de',
    parent                  = 'vbend01906019066',
    k0                      = 'k0_bv2de',
    h                       = 0.0392146299182919,
    edge_entry_angle        = 0.03737191615145,
    edge_exit_angle         = 0.03737191615145)

############################################################
# Correctors
############################################################

########################################
# Base Elements
########################################
env.new(name = 'hcorr00200000000', parent = xt.Bend, length = 0.2)
env.new(name = 'hcorr00225000000', parent = xt.Bend, length = 0.225)
env.new(name = 'hcorr00344400000', parent = xt.Bend, length = 0.3444)
env.new(name = 'hcorr00580000000', parent = xt.Bend, length = 0.58)
env.new(name = 'vcorr00075000000', parent = xt.Bend, length = 0.075, rot_s_rad = +np.pi/2)
env.new(name = 'vcorr00200000000', parent = xt.Bend, length = 0.2, rot_s_rad = +np.pi/2)
env.new(name = 'vcorr00288000000', parent = xt.Bend, length = 0.288, rot_s_rad = +np.pi/2)
env.new(name = 'vcorr00344400000', parent = xt.Bend, length = 0.3444, rot_s_rad = +np.pi/2)

########################################
# Cloned Elements
########################################
env.new(name = 'hx01e', parent = 'hcorr00200000000', k0 = 'k0_hx01e')
env.new(name = 'hx02e', parent = 'hcorr00200000000', k0 = 'k0_hx02e')
env.new(name = 'hx03e', parent = 'hcorr00200000000', k0 = 'k0_hx03e')
env.new(name = 'hx04e', parent = 'hcorr00200000000', k0 = 'k0_hx04e')
env.new(name = 'hx05e', parent = 'hcorr00200000000', k0 = 'k0_hx05e')
env.new(name = 'hx06e', parent = 'hcorr00200000000', k0 = 'k0_hx06e')
env.new(name = 'hx07e', parent = 'hcorr00200000000', k0 = 'k0_hx07e')
env.new(name = 'hx08e', parent = 'hcorr00200000000', k0 = 'k0_hx08e')
env.new(name = 'hx09e', parent = 'hcorr00200000000', k0 = 'k0_hx09e')
env.new(name = 'hx10e', parent = 'hcorr00200000000', k0 = 'k0_hx10e')
env.new(name = 'hw11e', parent = 'hcorr00200000000', k0 = 'k0_hw11e')
env.new(name = 'hw12e', parent = 'hcorr00200000000', k0 = 'k0_hw12e')
env.new(name = 'hw13e', parent = 'hcorr00200000000', k0 = 'k0_hw13e')
env.new(name = 'ht14e', parent = 'hcorr00200000000', k0 = 'k0_ht14e')
env.new(name = 'hm15e', parent = 'hcorr00200000000', k0 = 'k0_hm15e')
env.new(name = 'hm16e', parent = 'hcorr00200000000', k0 = 'k0_hm16e')
env.new(name = 'hm17e', parent = 'hcorr00200000000', k0 = 'k0_hm17e')
env.new(name = '-fzhinj1e1', parent = 'hcorr00225000000', k0 = 'k0_fzhinj1e1')
env.new(name = '-fzhinj1e2', parent = 'hcorr00225000000', k0 = 'k0_fzhinj1e2')
env.new(name = '-fzhinj1e3', parent = 'hcorr00225000000', k0 = 'k0_fzhinj1e3')
env.new(name = '-zhqi4e', parent = 'hcorr00344400000', k0 = 'k0_zhqi4e')
env.new(name = 'bh0e', parent = 'hcorr00580000000', k0 = 'k0_bh0e')
env.new(name = '-fzvkicke', parent = 'vcorr00075000000', k0 = 'k0_fzvkicke')
env.new(name = 'vx01e', parent = 'vcorr00200000000', k0 = 'k0_vx01e')
env.new(name = 'vx02e', parent = 'vcorr00200000000', k0 = 'k0_vx02e')
env.new(name = 'vx03e', parent = 'vcorr00200000000', k0 = 'k0_vx03e')
env.new(name = 'vx04e', parent = 'vcorr00200000000', k0 = 'k0_vx04e')
env.new(name = 'vx05e', parent = 'vcorr00200000000', k0 = 'k0_vx05e')
env.new(name = 'vx06e', parent = 'vcorr00200000000', k0 = 'k0_vx06e')
env.new(name = 'vx07e', parent = 'vcorr00200000000', k0 = 'k0_vx07e')
env.new(name = 'vw08e', parent = 'vcorr00200000000', k0 = 'k0_vw08e')
env.new(name = 'vw09e', parent = 'vcorr00200000000', k0 = 'k0_vw09e')
env.new(name = 'vw10e', parent = 'vcorr00200000000', k0 = 'k0_vw10e')
env.new(name = 'va11e', parent = 'vcorr00200000000', k0 = 'k0_va11e')
env.new(name = 'va12e', parent = 'vcorr00200000000', k0 = 'k0_va12e')
env.new(name = 'va13e', parent = 'vcorr00200000000', k0 = 'k0_va13e')
env.new(name = 'va14e', parent = 'vcorr00200000000', k0 = 'k0_va14e')
env.new(name = 'va15e', parent = 'vcorr00200000000', k0 = 'k0_va15e')
env.new(name = 'vt16e', parent = 'vcorr00200000000', k0 = 'k0_vt16e')
env.new(name = 'vt17e', parent = 'vcorr00200000000', k0 = 'k0_vt17e')
env.new(name = 'vb18e', parent = 'vcorr00200000000', k0 = 'k0_vb18e')
env.new(name = 'vb19e', parent = 'vcorr00200000000', k0 = 'k0_vb19e')
env.new(name = 'vc20e', parent = 'vcorr00200000000', k0 = 'k0_vc20e')
env.new(name = 'vc21e', parent = 'vcorr00200000000', k0 = 'k0_vc21e')
env.new(name = 'vm22e', parent = 'vcorr00200000000', k0 = 'k0_vm22e')
env.new(name = 'vm23e', parent = 'vcorr00200000000', k0 = 'k0_vm23e')
env.new(name = 'vm24e', parent = 'vcorr00200000000', k0 = 'k0_vm24e')
env.new(name = 'vm25e', parent = 'vcorr00200000000', k0 = 'k0_vm25e')
env.new(name = 'vm26e', parent = 'vcorr00200000000', k0 = 'k0_vm26e')
env.new(name = 'vm27e', parent = 'vcorr00200000000', k0 = 'k0_vm27e')
env.new(name = 'vm28e', parent = 'vcorr00200000000', k0 = 'k0_vm28e')
env.new(name = 'fy61e', parent = 'vcorr00288000000', k0 = 'k0_fy61e')
env.new(name = '-zvqi5e', parent = 'vcorr00344400000', k0 = 'k0_zvqi5e')

############################################################
# Quadrupoles
############################################################

########################################
# Base Elements
########################################
env.new(name = 'quad00376200000', parent = xt.Quadrupole, length = 0.3762)
env.new(name = 'quad00525300000', parent = xt.Quadrupole, length = 0.5253)
env.new(name = 'quad00826150000', parent = xt.Quadrupole, length = 0.82615)
env.new(name = 'quad01015230000', parent = xt.Quadrupole, length = 1.01523)

########################################
# Cloned Elements
########################################
env.new(name = 'qx3e', parent = 'quad00376200000', k1 = 'k1_qx3e')
env.new(name = 'qx4e', parent = 'quad00376200000', k1 = 'k1_qx4e')
env.new(name = 'qx5e', parent = 'quad00376200000', k1 = 'k1_qx5e')
env.new(name = 'qxf5e', parent = 'quad00376200000', k1 = 'k1_qxf5e')
env.new(name = 'qaf3e', parent = 'quad00376200000', k1 = 'k1_qaf3e')
env.new(name = 'qaf5e', parent = 'quad00376200000', k1 = 'k1_qaf5e')
env.new(name = 'qaf7e', parent = 'quad00376200000', k1 = 'k1_qaf7e')
env.new(name = 'qaf9e', parent = 'quad00376200000', k1 = 'k1_qaf9e')
env.new(name = 'qbf1e', parent = 'quad00376200000', k1 = 'k1_qbf1e')
env.new(name = 'qbf3e', parent = 'quad00376200000', k1 = 'k1_qbf3e')
env.new(name = 'qcf1e', parent = 'quad00376200000', k1 = 'k1_qcf1e')
env.new(name = 'qcf3e', parent = 'quad00376200000', k1 = 'k1_qcf3e')
env.new(name = 'qmf5e', parent = 'quad00376200000', k1 = 'k1_qmf5e')
env.new(name = 'qmf7e', parent = 'quad00376200000', k1 = 'k1_qmf7e')
env.new(name = 'qmf9e', parent = 'quad00376200000', k1 = 'k1_qmf9e')
env.new(name = 'qx1e', parent = 'quad00525300000', k1 = 'k1_qx1e')
env.new(name = 'qx2e', parent = 'quad00525300000', k1 = 'k1_qx2e')
env.new(name = 'qxd1e', parent = 'quad00525300000', k1 = 'k1_qxd1e')
env.new(name = 'qxf1e', parent = 'quad00525300000', k1 = 'k1_qxf1e')
env.new(name = 'qxd2e', parent = 'quad00525300000', k1 = 'k1_qxd2e')
env.new(name = 'qxf2e', parent = 'quad00525300000', k1 = 'k1_qxf2e')
env.new(name = 'qxd3e', parent = 'quad00525300000', k1 = 'k1_qxd3e')
env.new(name = 'qxf3e', parent = 'quad00525300000', k1 = 'k1_qxf3e')
env.new(name = 'qxd4e', parent = 'quad00525300000', k1 = 'k1_qxd4e')
env.new(name = 'qxf4e', parent = 'quad00525300000', k1 = 'k1_qxf4e')
env.new(name = 'qxd5e', parent = 'quad00525300000', k1 = 'k1_qxd5e')
env.new(name = 'qxd6e', parent = 'quad00525300000', k1 = 'k1_qxd6e')
env.new(name = 'qxf6e', parent = 'quad00525300000', k1 = 'k1_qxf6e')
env.new(name = 'qxd7e', parent = 'quad00525300000', k1 = 'k1_qxd7e')
env.new(name = 'qwfe', parent = 'quad00525300000', k1 = 'k1_qwfe')
env.new(name = 'qwde', parent = 'quad00525300000', k1 = 'k1_qwde')
env.new(name = 'qaf1e', parent = 'quad00525300000', k1 = 'k1_qaf1e')
env.new(name = 'qad2e', parent = 'quad00525300000', k1 = 'k1_qad2e')
env.new(name = 'qad4e', parent = 'quad00525300000', k1 = 'k1_qad4e')
env.new(name = 'qad6e', parent = 'quad00525300000', k1 = 'k1_qad6e')
env.new(name = 'qad8e', parent = 'quad00525300000', k1 = 'k1_qad8e')
env.new(name = 'qad10e', parent = 'quad00525300000', k1 = 'k1_qad10e')
env.new(name = 'qtf1e', parent = 'quad00525300000', k1 = 'k1_qtf1e')
env.new(name = 'qtd2e', parent = 'quad00525300000', k1 = 'k1_qtd2e')
env.new(name = 'qtf3e', parent = 'quad00525300000', k1 = 'k1_qtf3e')
env.new(name = 'qtd4e', parent = 'quad00525300000', k1 = 'k1_qtd4e')
env.new(name = 'qbd2e', parent = 'quad00525300000', k1 = 'k1_qbd2e')
env.new(name = 'qbd4e', parent = 'quad00525300000', k1 = 'k1_qbd4e')
env.new(name = 'qcd2e', parent = 'quad00525300000', k1 = 'k1_qcd2e')
env.new(name = 'qcd4e', parent = 'quad00525300000', k1 = 'k1_qcd4e')
env.new(name = 'qcf5e', parent = 'quad00525300000', k1 = 'k1_qcf5e')
env.new(name = 'qmd1e', parent = 'quad00525300000', k1 = 'k1_qmd1e')
env.new(name = 'qmf2e', parent = 'quad00525300000', k1 = 'k1_qmf2e')
env.new(name = 'qmf3e', parent = 'quad00525300000', k1 = 'k1_qmf3e')
env.new(name = 'qmd4e', parent = 'quad00525300000', k1 = 'k1_qmd4e')
env.new(name = 'qmd6e', parent = 'quad00525300000', k1 = 'k1_qmd6e')
env.new(name = 'qmd8e', parent = 'quad00525300000', k1 = 'k1_qmd8e')
env.new(name = 'qmd10e', parent = 'quad00525300000', k1 = 'k1_qmd10e')
env.new(name = 'qi4e', parent = 'quad00826150000', k1 = 'k1_qi4e')
env.new(name = 'qi5e', parent = 'quad01015230000', k1 = 'k1_qi5e')

############################################################
# Sextupoles
############################################################

########################################
# Base Elements
########################################
env.new(name = 'sext00100000000', parent = xt.Sextupole, length = 0.1)

########################################
# Cloned Elements
########################################
env.new(name = 'smf1e', parent = 'sext00100000000', k2 = 'k2_smf1e')
env.new(name = 'smf2e', parent = 'sext00100000000', k2 = 'k2_smf2e')

############################################################
# Cavities
############################################################
env.new(
    name        = 'acw',
    parent      = xt.Cavity,
    length      = 3.033,
    frequency   = 'freq_acw * (1 + fshift)',
    voltage     = 'volt_acw',
    lag         = 'lag_acw')

############################################################
# Reference Shifts
############################################################

########################################
# XYShifts
########################################
env.new(
    name        = 'cinjp_dxy',
    parent      = xt.XYShift,
    dx          = 'dx_cinjp_dxy',
    dy          = 'dy_cinjp_dxy')

########################################
# YRotations (CHI1)
########################################
env.new(
    name        = 'cinjp_chi1',
    parent      = xt.YRotation,
    angle       = 'chi1_cinjp_chi1')


############################################################
# Markers
############################################################
ALL_MARKERS = [
    'cbe', 'che', 'cve', 'endinjk1', 'f1', 'fd', 'fu', 'fvc', 'g1', 'g2', 'g5',
    'ga1i', 'ga2i', 'ga3i', 'ga3o', 'ga4i', 'gki', 'gko', 'gq', 'gsli', 'gti',
    'injp', 'jb0e', 'jb1e', 'jbc1p', 'jbh0e', 'jbh1ae', 'jbh1e', 'jbh2e',
    'jbh3e', 'jbh4ae', 'jbh4e', 'jbhd0', 'jbv1de', 'jbv1ue', 'jbv2de', 'jbv2ue',
    'jendinjk1', 'jfzhinj1e1', 'jfzhinj1e2', 'jfzhinj1e3', 'jinjp', 'jpqi4e',
    'jqad10e', 'jqad2e', 'jqad4e', 'jqad6e', 'jqad8e', 'jqaf1e', 'jqaf3e',
    'jqaf5e', 'jqaf7e', 'jqaf9e', 'jqbd2e', 'jqbd4e', 'jqbf1e', 'jqbf3e',
    'jqcd2e', 'jqcd4e', 'jqcf1e', 'jqcf3e', 'jqcf5e', 'jqi4e', 'jqi5e',
    'jqmd10e', 'jqmd1e', 'jqmd4e', 'jqmd6e', 'jqmd8e', 'jqmf2e', 'jqmf3e',
    'jqmf5e', 'jqmf7e', 'jqmf9e', 'jqtd2e', 'jqtd4e', 'jqtf1e', 'jqtf3e',
    'jqwde', 'jqwfe', 'jqx1e', 'jqx2e', 'jqx3e', 'jqx4e', 'jqx5e', 'jqxd1e',
    'jqxd2e', 'jqxd3e', 'jqxd4e', 'jqxd5e', 'jqxd6e', 'jqxd7e', 'jqxf1e',
    'jqxf2e', 'jqxf3e', 'jqxf4e', 'jqxf5e', 'jqxf6e', 'jse1', 'jse2', 'jse3',
    'jse4', 'mqi4e', 'mqi5e', 'mse1', 'mse10', 'mse11', 'mse12', 'mse13',
    'mse14', 'mse15', 'mse16', 'mse17', 'mse19', 'mse2', 'mse20', 'mse3',
    'mse4', 'mse5', 'mse6', 'mse7', 'mse8', 'mse9', 'mt', 'mw', 'pbte',
    'pkicker1', 'pqi4e', 's8bte1', 'sp61h1', 'spqad10e_m', 'spqad2e_m',
    'spqad6e_m', 'spqaf1e_m', 'spqaf3e_s', 'spqaf5e_s', 'spqaf7e_s',
    'spqaf9e_s', 'spqbd2e_m', 'spqbf1e_s', 'spqbf3e_s', 'spqcd2e_m',
    'spqcd4e_m', 'spqcf1e_s', 'spqcf3e_s', 'spqcf5e_m', 'spqmd10e_m',
    'spqmd1e_1m', 'spqmd1e_2m', 'spqmd1e_3m', 'spqmd4e_m', 'spqmd6e_m',
    'spqmd8e_m', 'spqmf2e_1m', 'spqmf2e_2m', 'spqmf3e_m', 'spqmf5e_s',
    'spqmf7e_s', 'spqmf9e_s', 'spqtd2e_m', 'spqtd4e_m', 'spqtf1e_m',
    'spqtf3e_m', 'spqwde_1m', 'spqwde_2m', 'spqwde_3m', 'spqwfe_1m',
    'spqwfe_2m', 'spqwfe_3m', 'spqx1e_m', 'spqx2e_m', 'spqx3e_s', 'spqx4e_s',
    'spqx5e_s', 'spqxd1e_m', 'spqxd2e_m', 'spqxd3e_m', 'spqxd4e_m', 'spqxd5e_m',
    'spqxd6e_m', 'spqxd7e_m', 'spqxf1e_m', 'spqxf2e_m', 'spqxf3e_m',
    'spqxf4e_m', 'spqxf5e_s', 'spqxf6e_m', 'starther']
for marker in ALL_MARKERS:
    env.new(name = marker, parent = xt.Marker)

############################################################
# Create Line
############################################################
env.new_line(
    name        = 'line',
    components  = [
        'pbte', 'bc1p', 'jbc1p', 'ldmp1a', 'bh0e', 'jbh0e', 'lmblh', 'lmblh',
        'ldmp1b', 'lpmh', 'sp61h1', 'lpmh', 'lblcc', 'bhd0', 'jbhd0', 'lce0a',
        'lmblh', 'spqx1e_m', 'lmblh', 'lqmd', 'qx1e', 'jqx1e', 'lqz', 'hx01e',
        'lce01a', 'ld300', 'lce01a', 'vx01e', 'lqz', 'qx2e', 'jqx2e', 'lqmd',
        'lmblh', 'spqx2e_m', 'lmblh', 'lce1a', 'lpmh', 'mse1', 'lpmh', 'lbell',
        'lblcd', 'b0e', 'jb0e', 'lce2a1', 'hx02e', 'ld260', 'lqmc', 'qx3e',
        'jqx3e', 'lqmc', 'lmblh', 'spqx3e_s', 'lmblh', 'lbvqxa', 'bv1ue',
        'jbv1ue', 'lblc', 'lbell', 'ld500', 'lbva', 'bv1de', 'jbv1de', 'lbqx2',
        'qx4e', 'jqx4e', 'lqmc', 'lmblh', 'spqx4e_s', 'lmblh', 'lce21', 'lpmlh',
        'mse2', 'lpmlh', 'ld350', 'hx03e', 'lmz', 'lmblh', 'spqx5e_s', 'lmblh',
        'lqmc', 'qx5e', 'jqx5e', 'lbqx0', 'b1e', 'jb1e', 'lce3a', 'lmblh',
        'spqxd1e_m', 'lmblh', 'lqmd', 'qxd1e', 'jqxd1e', 'lbqx', 'bv1ue',
        'jbv1ue', 'lblc', 'lbell', 'ld500', 'lbva', 'bv1de', 'jbv1de', 'lce2b1',
        'ld210', 'ld80', 'hx04e', 'lqz', 'qxf1e', 'jqxf1e', 'lqmd', 'lmblh',
        'spqxf1e_m', 'lmblh', 'lchbh', 'che', 'lchbh', 'lce5a', 'lcbbh', 'cbe',
        'lcbbh', 'lcbbh', 's8bte1', 'lcbbh', 'lbell', 'ld80', 'vx02e', 'lqz',
        'qxd2e', 'jqxd2e', 'lqmd', 'lmblh', 'spqxd2e_m', 'lmblh2', 'fy61e',
        'lce5b', 'hx05e', 'lmz', 'lmblh', 'spqxf2e_m', 'lmblh', 'lqmd', 'qxf2e',
        'jqxf2e', 'lqmd', 'ld150', 'lchbh', 'lchbh', 'lpmlh', 'mse3', 'lpmlh',
        'lce6a', 'lchbh', 'che', 'lchbh', 'lblce', 'b1e', 'jb1e', 'lbqx1a',
        'vx03e', 'lqz', 'g1', 'qxd3e', 'jqxd3e', 'lqmd', 'lmbsh', 'spqxd3e_m',
        'lmbsh', 'lbqx1b', 'b1e', 'jb1e', 'l4x', 'g2', 'l4y1', 'ld500', 'lbell',
        'lgv', 'ld100', 'hx06e', 'lqz', 'qxf3e', 'jqxf3e', 'lqmd', 'lmblh',
        'spqxf3e_m', 'lmblh', 'l41', 'vx04e', 'lqz', 'qxd4e', 'jqxd4e', 'lqmd',
        'lmblh', 'spqxd4e_m', 'lmblh', 'l42', 'ld250', 'lpmlh', 'mse4', 'lpmlh',
        'lbell', 'ld100', 'hx07e', 'lqz', 'qxf4e', 'jqxf4e', 'lqmd', 'lmblh',
        'spqxf4e_m', 'lmblh', 'l43', 'vx05e', 'lqz', 'qxd5e', 'jqxd5e', 'lqmd',
        'lmblh', 'spqxd5e_m', 'lmblh', 'l44', 'b1e', 'jb1e', 'lbqx14', 'hx08e',
        'lqzc', 'qxf5e', 'jqxf5e', 'lqmc', 'lmblh', 'spqxf5e_s', 'lmblh', 'l45',
        'b1e', 'jb1e', 'l4b', 'b1e', 'jb1e', 'lbqxc', 'qxd6e', 'jqxd6e', 'lqmd',
        'lmblh', 'spqxd6e_m', 'lmblh', 'lmz', 'vx06e', 'l46', 'hx09e', 'lqz',
        'qxf6e', 'jqxf6e', 'lqmd', 'lmblh', 'spqxf6e_m', 'lmblh', 'l47bell',
        'ld250', 'l47d', 'acw', 'ld250', 'l47bell', 'lmbsh', 'vx07e', 'lqz',
        'qxd7e', 'jqxd7e', 'lqmd', 'lmblh', 'spqxd7e_m', 'lmblh', 'l48d1',
        'acw', 'l48d2', 'hx10e', 'lqz', 'f1', 'qwfe', 'jqwfe', 'lqmd', 'lmblh',
        'spqwfe_1m', 'lmblh', 'lbell', 'lpmlh', 'mse5', 'lpmlh', 'ld250',
        'lw4d1', 'acw', 'lw4d2', 'acw', 'lw4d3', 'vw08e', 'lqz', 'qwde',
        'jqwde', 'lqmd', 'lmblh', 'spqwde_1m', 'lmblh', 'lbell', 'ld250',
        'lmwh', 'mw', 'lmwh', 'lw2', 'hw11e', 'lqz', 'qwfe', 'jqwfe', 'lqmd',
        'lmblh', 'spqwfe_2m', 'lmblh', 'lbell', 'ld250', 'lmwh', 'mw', 'lmwh',
        'lw1', 'vw09e', 'lqz', 'qwde', 'jqwde', 'lqmd', 'lmblh', 'spqwde_2m',
        'lmblh', 'lbell', 'ld250', 'lmwh', 'mw', 'lmwh', 'lw2', 'hw12e', 'lqz',
        'qwfe', 'jqwfe', 'lqmd', 'lmblh', 'spqwfe_3m', 'lmblh', 'lbell',
        'ld250', 'lmwh', 'mw', 'lmwh', 'lw1', 'vw10e', 'lqz', 'qwde', 'jqwde',
        'lqmd', 'lmblh', 'spqwde_3m', 'lmblh', 'ld250', 'lbell', 'lmbth', 'mt',
        'lmbth', 'lw3', 'hw13e', 'lqz', 'qaf1e', 'jqaf1e', 'lqmd', 'lmblh',
        'spqaf1e_m', 'lmblh', 'lpmh', 'mse6', 'lpmh', 'ld100', 'lwa1', 'ga1i',
        'bh1ae', 'jbh1ae', 'lbs1a', 'va11e', 'ld310', 'bh1e', 'jbh1e', 'lx1',
        'qad2e', 'jqad2e', 'lqmd', 'lmbsh', 'spqad2e_m', 'lmbsh', 'lj1', 'bh1e',
        'jbh1e', 'lbb1', 'bh1e', 'jbh1e', 'lbq13a', 'lchbsh', 'che', 'lchbsh',
        'ld100c', 'qaf3e', 'jqaf3e', 'lqmac', 'lmblh', 'spqaf3e_s', 'lmblh',
        'lbq12a', 'bh1e', 'jbh1e', 'lbb1', 'bh1e', 'jbh1e', 'lbq11a', 'va12e',
        'lqz', 'qad4e', 'jqad4e', 'lqma', 'lpmlh', 'mse7', 'lpmlh', 'lbq12a',
        'bh1e', 'jbh1e', 'lbb1', 'bh1e', 'jbh1e', 'lbq13', 'qaf5e', 'jqaf5e',
        'lqmac', 'lmblh', 'spqaf5e_s', 'lmblh', 'lbq12a', 'bh1e', 'jbh1e',
        'lbb1', 'bh1e', 'jbh1e', 'lbq11a', 'va13e', 'lqz', 'qad6e', 'jqad6e',
        'lqma', 'lmblh', 'spqad6e_m', 'lmblh', 'lbq12a', 'bh1e', 'jbh1e',
        'lbb1', 'bh1e', 'jbh1e', 'lbq13', 'qaf7e', 'jqaf7e', 'lqmac', 'lmblh',
        'spqaf7e_s', 'lmblh', 'lbq12a', 'bh1e', 'jbh1e', 'lbb1', 'bh1e',
        'jbh1e', 'lbq11a', 'va14e', 'lqz', 'qad8e', 'jqad8e', 'lqma', 'lpmlh',
        'mse8', 'lpmlh', 'lbq12a', 'bh1e', 'jbh1e', 'lbb1', 'bh1e', 'jbh1e',
        'lbq13', 'qaf9e', 'jqaf9e', 'lqmac', 'lmblh', 'spqaf9e_s', 'lmblh',
        'lbq12a', 'bh1e', 'jbh1e', 'lbb1', 'bh1e', 'jbh1e', 'lbq11a', 'va15e',
        'lqz', 'qad10e', 'jqad10e', 'lqma', 'lmblh', 'spqad10e_m', 'lmblh',
        'lbq12a', 'bh1e', 'jbh1e', 'gki', 'lk1', 'gko', 'bh1e', 'jbh1e', 'gti',
        'ls1', 'qtf1e', 'jqtf1e', 'lqmd', 'lmblh', 'spqtf1e_m', 'lmblh', 'ls21',
        'lpmlh', 'mse9', 'lpmlh', 'ls2a', 'vt16e', 'lmz', 'lmblh', 'spqtd2e_m',
        'lmblh', 'lqmd', 'qtd2e', 'jqtd2e', 'ls3', 'qtf3e', 'jqtf3e', 'lqmd',
        'lmblh', 'spqtf3e_m', 'lmblh', 'lmz', 'ht14e', 'ls4a', 'vt17e', 'lqz',
        'qtd4e', 'jqtd4e', 'lqmd', 'lmblh', 'spqtd4e_m', 'lmblh', 'ls5a',
        'lcbbh', 'cbe', 'lcbbh', 'lpmlh', 'mse10', 'lpmlh', 'ld320', 'ga2i',
        'bh2e', 'jbh2e', 'g5', 'lbb20', 'bh2e', 'jbh2e', 'lbq23', 'qbf1e',
        'jqbf1e', 'lqmc', 'lmblh', 'spqbf1e_s', 'lmblh', 'lbq22a', 'bh2e',
        'jbh2e', 'lbb2', 'bh2e', 'jbh2e', 'lbq21a', 'vb18e', 'lqz', 'qbd2e',
        'jqbd2e', 'lqmd', 'lmblh', 'spqbd2e_m', 'lmblh', 'lbq22a', 'bh2e',
        'jbh2e', 'lbb2', 'bh2e', 'jbh2e', 'lbq23', 'qbf3e', 'jqbf3e', 'lqmc',
        'lmblh', 'spqbf3e_s', 'lmblh', 'lbq22a', 'bh2e', 'jbh2e', 'lbb2',
        'bh2e', 'jbh2e', 'lbq21a', 'vb19e', 'lqz', 'qbd4e', 'jqbd4e', 'lqmd',
        'lpmlh', 'mse11', 'lpmlh', 'lbq22a', 'bh2e', 'jbh2e', 'lbb2c', 'bh3e',
        'jbh3e', 'lx7a', 'ga3i', 'bh3e', 'jbh3e', 'lbq33', 'qcf1e', 'jqcf1e',
        'lqmc', 'lmblh', 'spqcf1e_s', 'lmblh', 'lbq32a', 'bh3e', 'jbh3e',
        'lbb3', 'bh3e', 'jbh3e', 'lbq31', 'vc20e', 'lqz', 'qcd2e', 'jqcd2e',
        'lqmd', 'lmblh', 'spqcd2e_m', 'lmblh', 'lbq3a', 'bh3e', 'jbh3e', 'lbb3',
        'bh3e', 'jbh3e', 'lbq33', 'qcf3e', 'jqcf3e', 'lqmc', 'lmblh',
        'spqcf3e_s', 'lmblh', 'lbell', 'lpmlh', 'mse12', 'lpmlh', 'lbq32b',
        'bh3e', 'jbh3e', 'lbb3', 'bh3e', 'jbh3e', 'lbq31', 'vc21e', 'lqz',
        'qcd4e', 'jqcd4e', 'lqmd', 'lmblh', 'spqcd4e_m', 'lmblh', 'lbq34a',
        'bh3e', 'jbh3e', 'lbb3', 'bh3e', 'jbh3e', 'ga3o', 'lbq3', 'qcf5e',
        'jqcf5e', 'gsli', 'lqmd', 'lmbsh', 'spqcf5e_m', 'lmbsh', 'lm1a1',
        'bv2ue', 'jbv2ue', 'lmbbv2', 'fu', 'lmbbv2', 'bv2ue', 'jbv2ue', 'lm1x1',
        'ld500', 'ld100', 'vm22e', 'lqz', 'qmd1e', 'jqmd1e', 'lqmd', 'lmblh',
        'spqmd1e_1m', 'lmblh', 'lmqqv1', 'lchbh', 'che', 'lchbh', 'lpmlh',
        'mse13', 'lpmlh', 'ld100', 'hm15e', 'lqz', 'qmf2e', 'jqmf2e', 'lqmd',
        'lmblh', 'spqmf2e_1m', 'lmblh', 'lmqqv2', 'ld500', 'ld100', 'vm23e',
        'lqz', 'fvc', 'qmd1e', 'jqmd1e', 'lqmd', 'lmblh', 'spqmd1e_2m', 'lmblh',
        'lmqqv1', 'lchbh', 'che', 'lchbh', 'lpmlh', 'mse14', 'lpmlh', 'ld100',
        'hm16e', 'lqz', 'qmf2e', 'jqmf2e', 'lqmd', 'lmblh', 'spqmf2e_2m',
        'lmblh', 'lmqqv2', 'ld500', 'ld100', 'vm24e', 'lqz', 'qmd1e', 'jqmd1e',
        'lqmd', 'lmblh', 'spqmd1e_3m', 'lmblh', 'lmqqv3', 'ld500', 'ld100',
        'hm17e', 'lqz', 'qmf3e', 'jqmf3e', 'lqmd', 'lmbsh', 'spqmf3e_m',
        'lmbsh', 'lm1y1', 'bv2de', 'jbv2de', 'lmbbv2', 'fd', 'lmbbv2', 'bv2de',
        'jbv2de', 'lm2a', 'lmqqh1', 'lcvbh', 'cve', 'lcvbh', 'lpmlh', 'mse15',
        'lpmlh', 'ld100', 'vm25e', 'lqz', 'qmd4e', 'jqmd4e', 'lqmd', 'lmblh',
        'spqmd4e_m', 'lmblh', 'lmqqh2', 'lmwh', 'mw', 'lmwh', 'lpmlh', 'mse16',
        'lpmlh', 'ld100c', 'qmf5e', 'jqmf5e', 'lqmc', 'lmblh', 'spqmf5e_s',
        'lmblh', 'lm2ba', 'ga4i', 'bh4e', 'jbh4e', 'lmbbh', 'bh4e', 'jbh4e',
        'lmbbh1', 'bh4e', 'jbh4e', 'lbq411', 'vm26e', 'lqz', 'qmd6e', 'jqmd6e',
        'lqmd', 'lmblh', 'spqmd6e_m', 'lmblh', 'lbq421', 'bh4e', 'jbh4e',
        'lmbbh', 'bh4e', 'jbh4e', 'lbq43a1', 'smf1e', 'lbq43a2', 'lpmlh',
        'lpmlh', 'ld100c', 'qmf7e', 'jqmf7e', 'lqmc', 'lmblh', 'spqmf7e_s',
        'lmblh', 'lm4aa', 'bh4e', 'jbh4e', 'lmbbh', 'bh4e', 'jbh4e', 'lm4b1',
        'lcvbh', 'cve', 'lcvbh', 'lpmlh', 'mse17', 'lpmlh', 'ld250', 'lbell',
        'ld100', 'vm27e', 'lqz', 'qmd8e', 'jqmd8e', 'lqmd', 'lmblh',
        'spqmd8e_m', 'lmblh', 'ld100', 'lmbt1h', 'mt', 'lmbt1h', 'ld100',
        'lm4c1', 'bh4e', 'jbh4e', 'lmbbh', 'bh4e', 'jbh4e', 'lm4da1', 'smf2e',
        'lm4da2', 'ld090', 'ld100c', 'qmf9e', 'jqmf9e', 'lqmc', 'lmblh',
        'spqmf9e_s', 'lmblh', 'lm4ea', 'bh4ae', 'jbh4ae', 'lmbbhc', 'bh4ae',
        'jbh4ae', 'lm4f1', 'ld100', 'vm28e', 'lqz', 'qmd10e', 'jqmd10e', 'lqmd',
        'lmblh', 'spqmd10e_m', 'lmi1', 'gq', 'lmi2', 'lmi3', 'se4', 'jse4',
        'lse3', 'se3', 'jse3', 'lse2', 'se2', 'jse2', 'mse19', 'lse1', 'se1',
        'jse1', 'mse20', 'injp', 'jinjp', 'cinjp_dxy', 'cinjp_chi1', 'lfr049',
        'jpqi4e', 'qi4e', 'jqi4e', 'starther', 'lfr050', 'mqi4e', 'lfr051',
        '-zhqi4e', 'lfr052', 'mqi5e', 'lfr057', 'qi5e', 'jqi5e', 'lfr058',
        '-zvqi5e', 'lfr083', '-fzvkicke', 'lfr059', '-fzhinj1e1', 'jfzhinj1e1',
        'lfr081', '-fzhinj1e2', 'jfzhinj1e2', 'lfr082', '-fzhinj1e3',
        'jfzhinj1e3', 'endinjk1', 'jendinjk1'])
line = env.lines['line']
line.particle_ref = env.particle_ref.copy()

################################################################################
# Configure Modelling
################################################################################

########################################
# Set integrators
########################################
tt          = line.get_table()
tt_drift    = tt.rows[tt.element_type == "Drift"]
tt_bend     = tt.rows[tt.element_type == "Bend"]
tt_quad     = tt.rows[tt.element_type == "Quadrupole"]
tt_sext     = tt.rows[tt.element_type == "Sextupole"]
tt_oct      = tt.rows[tt.element_type == "Octupole"]
tt_mult     = tt.rows[tt.element_type == "Multipole"]
tt_sol      = tt.rows[tt.element_type == "Solenoid"]
tt_cavi     = tt.rows[tt.element_type == "Cavity"]

line.set(
    tt_drift,
    model               = "exact")
line.set(
    tt_bend,
    model               = "mat-kick-mat",
    integrator          = "uniform",
    num_multipole_kicks = 3)
line.set(
    tt_quad,
    model               = "mat-kick-mat",
    integrator          = "uniform",
    num_multipole_kicks = 3)
line.set(
    tt_sext,
    model               = "mat-kick-mat",
    integrator          = "yoshida4",
    num_multipole_kicks = 1)
line.set(
    tt_oct,
    model               = "mat-kick-mat",
    integrator          = "yoshida4",
    num_multipole_kicks = 1)
line.set(
    tt_mult,
    num_multipole_kicks = 3)
line.set(
    tt_sol,
    num_multipole_kicks = 3)
line.set(
    tt_cavi,
    model               = "drift-kick-drift-exact",
    integrator          = "yoshida4",
    absolute_time       = False)

########################################
# Set bend edges
########################################
line.configure_bend_model(edge = "full")

########################################
# Replace repeated elements
########################################
line.replace_all_repeated_elements()

################################################################################
# Offset markers
################################################################################

############################################################
# Get length of the line
############################################################
length   = line.get_length()

############################################################
# Offset marker locations
############################################################
MARKER_POSITIONS = {
    'pqi4e':                       [465.984121549171],
    'pkicker1':                    [483.189423463286]}

#########################################
# Install Markers
############################################################
marker_insertions   = []
for marker, insert_at_s_values in MARKER_POSITIONS.items():
    for insert_at_s in insert_at_s_values:
        if (length - insert_at_s) > 1.00E-09:
            marker_insertions.append(
                env.place(name = marker, at = insert_at_s))
        else:
            line.append_element(name = marker)
try:
    line.insert(marker_insertions, s_tol = 1.00E-09)
except AssertionError as err:
    print("Couldn't insert all the markers. Usually this is because of negative drifts")
    print(err)
    pass


########################################
# Replace repeated elements
########################################
line.replace_all_repeated_elements()
