"""
SuperKEKB BTP
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
env["p0c"]      = 3841286237.7791615
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
env.new(name = 'll3373', parent = xt.Drift, length = 0.1525)
env.new(name = 'll3374', parent = xt.Drift, length = 0.003)
env.new(name = 'll3376', parent = xt.Drift, length = 0.086)
env.new(name = 'lqf564', parent = xt.Drift, length = 0.204)
env.new(name = 'll3378', parent = xt.Drift, length = 0.086)
env.new(name = 'll3380', parent = xt.Drift, length = 0.0)
env.new(name = 'll3381', parent = xt.Drift, length = 0.0)
env.new(name = 'll3382', parent = xt.Drift, length = 0.086)
env.new(name = 'll3384', parent = xt.Drift, length = 0.0885086)
env.new(name = 'gv570', parent = xt.Drift, length = 0.078)
env.new(name = 'll3386', parent = xt.Drift, length = 0.10597)
env.new(name = 'll3393', parent = xt.Drift, length = 0.1319388)
env.new(name = 'll3397', parent = xt.Drift, length = 0.035)
env.new(name = 'll3398', parent = xt.Drift, length = 0.093377)
env.new(name = 'll3399', parent = xt.Drift, length = 0.0425)
env.new(name = 'ws572', parent = xt.Drift, length = 0.0)
env.new(name = 'll3401', parent = xt.Drift, length = 0.0425)
env.new(name = 'll3402', parent = xt.Drift, length = 0.0425)
env.new(name = 'sc572', parent = xt.Drift, length = 0.0)
env.new(name = 'll3404', parent = xt.Drift, length = 0.0425)
env.new(name = 'll3405', parent = xt.Drift, length = 0.0425)
env.new(name = 'll3406', parent = xt.Drift, length = 0.10597)
env.new(name = 'll3413', parent = xt.Drift, length = 0.1319388)
env.new(name = 'gv574', parent = xt.Drift, length = 0.1205)
env.new(name = 'ln001', parent = xt.Drift, length = 0.479025)
env.new(name = 'ln002', parent = xt.Drift, length = 0.381212)
env.new(name = 'ln003', parent = xt.Drift, length = 0.421212)
env.new(name = 'ln004', parent = xt.Drift, length = 2.0291)
env.new(name = 'ln005', parent = xt.Drift, length = 0.4)
env.new(name = 'ln006', parent = xt.Drift, length = 2.0)
env.new(name = 'ln007', parent = xt.Drift, length = 1.348151)
env.new(name = 'ln008', parent = xt.Drift, length = 0.215)
env.new(name = 'ln009', parent = xt.Drift, length = 0.35)
env.new(name = 'ln010', parent = xt.Drift, length = 0.225)
env.new(name = 'ln011', parent = xt.Drift, length = 0.0425)
env.new(name = 'sc584', parent = xt.Drift, length = 0.0)
env.new(name = 'ln012', parent = xt.Drift, length = 0.0425)
env.new(name = 'll3464', parent = xt.Drift, length = 0.2)
env.new(name = 'll3465', parent = xt.Drift, length = 0.003)
env.new(name = 'll3467', parent = xt.Drift, length = 0.086)
env.new(name = 'lqf584', parent = xt.Drift, length = 0.204)
env.new(name = 'll3469', parent = xt.Drift, length = 0.086)
env.new(name = 'll3471', parent = xt.Drift, length = 0.0)
env.new(name = 'll3472', parent = xt.Drift, length = 0.0)
env.new(name = 'll3473', parent = xt.Drift, length = 0.086)
env.new(name = 'll3475', parent = xt.Drift, length = 0.0885086)
env.new(name = 'gv584', parent = xt.Drift, length = 0.085)
env.new(name = 'll3477', parent = xt.Drift, length = 0.12357)
env.new(name = 'll3478e', parent = xt.Drift, length = 0.21660130193995)
env.new(name = 'l3481p', parent = xt.Drift, length = 0.2129794552984)
env.new(name = 'l3482', parent = xt.Drift, length = 0.09995700001888677)
env.new(name = 'l3483p', parent = xt.Drift, length = 0.19194344671905)
env.new(name = 'l3486p', parent = xt.Drift, length = 0.19194344671905)
env.new(name = 'l3487', parent = xt.Drift, length = 0.0425)
env.new(name = 'sc612', parent = xt.Drift, length = 0.0)
env.new(name = 'l3489', parent = xt.Drift, length = 0.0425)
env.new(name = 'l3490', parent = xt.Drift, length = 0.10708053029929795)
env.new(name = 'l3491p', parent = xt.Drift, length = 0.20995290201745)
env.new(name = 'l3492', parent = xt.Drift, length = 0.0)
env.new(name = 'sr613', parent = xt.Drift, length = 0.0)
env.new(name = 'l3497p', parent = xt.Drift, length = 0.18995290201745002)
env.new(name = 'l3498', parent = xt.Drift, length = 0.0425)
env.new(name = 'sp613', parent = xt.Drift, length = 0.0)
env.new(name = 'l3500', parent = xt.Drift, length = 0.3475441959651028)
env.new(name = 'l3501', parent = xt.Drift, length = 0.075)
env.new(name = 'sl613', parent = xt.Drift, length = 0.0)
env.new(name = 'l3503', parent = xt.Drift, length = 0.075)
env.new(name = 'l3504', parent = xt.Drift, length = 0.075)
env.new(name = 'sc613', parent = xt.Drift, length = 0.0)
env.new(name = 'l3506', parent = xt.Drift, length = 0.075)
env.new(name = 'l3507p', parent = xt.Drift, length = 0.18995290201745002)
env.new(name = 'l3510p', parent = xt.Drift, length = 0.20995290201745)
env.new(name = 'l3511', parent = xt.Drift, length = 0.10708053029929795)
env.new(name = 'l3512', parent = xt.Drift, length = 0.0425)
env.new(name = 'sc614', parent = xt.Drift, length = 0.0)
env.new(name = 'l3514', parent = xt.Drift, length = 0.0425)
env.new(name = 'l3515p', parent = xt.Drift, length = 0.19194344671905)
env.new(name = 'l3518p', parent = xt.Drift, length = 0.19194344671905)
env.new(name = 'l3519', parent = xt.Drift, length = 0.09995700001888677)
env.new(name = 'l3520p', parent = xt.Drift, length = 0.2129794552984)
env.new(name = 'l3523p', parent = xt.Drift, length = 0.2129794552984)
env.new(name = 'l3524', parent = xt.Drift, length = 0.12499554470159183)
env.new(name = 'l3526', parent = xt.Drift, length = 0.15)
env.new(name = 'l3528', parent = xt.Drift, length = 0.0805)
env.new(name = 'l3529', parent = xt.Drift, length = 0.0)
env.new(name = 'sp616', parent = xt.Drift, length = 0.0)
env.new(name = 'l3531', parent = xt.Drift, length = 0.1075)
env.new(name = 'l3533', parent = xt.Drift, length = 0.118)
env.new(name = 'l3534', parent = xt.Drift, length = 0.118)
env.new(name = 'l3536', parent = xt.Drift, length = 0.13)
env.new(name = 'l3537', parent = xt.Drift, length = 0.0575)
env.new(name = 'wm616', parent = xt.Drift, length = 0.0)
env.new(name = 'l3539', parent = xt.Drift, length = 0.07)
env.new(name = 'sc616', parent = xt.Drift, length = 0.0)
env.new(name = 'l3541', parent = xt.Drift, length = 0.0425)
env.new(name = 'l3542', parent = xt.Drift, length = 0.04)
env.new(name = 'gv616', parent = xt.Drift, length = 0.107)
env.new(name = 'l3544', parent = xt.Drift, length = 0.15451363373028146)
env.new(name = 'l3545', parent = xt.Drift, length = 0.035)
env.new(name = 'l3549', parent = xt.Drift, length = 0.035)
env.new(name = 'l3550', parent = xt.Drift, length = 0.140666)
env.new(name = 'l3551', parent = xt.Drift, length = 0.035)
env.new(name = 'l3555', parent = xt.Drift, length = 0.035)
env.new(name = 'l3556', parent = xt.Drift, length = 0.088)
env.new(name = 'gv618', parent = xt.Drift, length = 0.107)
env.new(name = 'l3558', parent = xt.Drift, length = 0.0830000000221193)
env.new(name = 'l3559', parent = xt.Drift, length = 0.0425)
env.new(name = 'sp618', parent = xt.Drift, length = 0.0)
env.new(name = 'l3561', parent = xt.Drift, length = 0.1075)
env.new(name = 'l3563', parent = xt.Drift, length = 0.11806872546977257)
env.new(name = 'l3564', parent = xt.Drift, length = 0.058)
env.new(name = 'l3566', parent = xt.Drift, length = 0.07)
env.new(name = 'l3567', parent = xt.Drift, length = 0.16268516054277254)
env.new(name = 'l3568', parent = xt.Drift, length = 0.07)
env.new(name = 'l3570', parent = xt.Drift, length = 0.058)
env.new(name = 'l3571', parent = xt.Drift, length = 0.11120092069177256)
env.new(name = 'lcq2bc', parent = xt.Drift, length = 0.1530000000221193)
env.new(name = 'lpmlih', parent = xt.Drift, length = 0.065)
env.new(name = 'l12', parent = xt.Drift, length = 0.252)
env.new(name = 'lcq31', parent = xt.Drift, length = 2.1317964757573)
env.new(name = 'ld250', parent = xt.Drift, length = 0.25)
env.new(name = 'lpmlh', parent = xt.Drift, length = 0.125)
env.new(name = 'lbell', parent = xt.Drift, length = 0.09)
env.new(name = 'lmblh', parent = xt.Drift, length = 0.125)
env.new(name = 'lqmc', parent = xt.Drift, length = 0.1065)
env.new(name = 'lqzc', parent = xt.Drift, length = 0.1865)
env.new(name = 'l2a1', parent = xt.Drift, length = 5.0631267939784)
env.new(name = 'l2b1', parent = xt.Drift, length = 4.4731267939784)
env.new(name = 'lcbbh', parent = xt.Drift, length = 0.125)
env.new(name = 'l2c1', parent = xt.Drift, length = 4.7631267939784)
env.new(name = 'lmz', parent = xt.Drift, length = 0.1)
env.new(name = 'lqcc', parent = xt.Drift, length = 0.1065)
env.new(name = 'lchbh', parent = xt.Drift, length = 0.125)
env.new(name = 'l2d1', parent = xt.Drift, length = 3.901672509765)
env.new(name = 'ld780', parent = xt.Drift, length = 0.753625)
env.new(name = 'lbqx1a', parent = xt.Drift, length = 0.9482825863412)
env.new(name = 'lbqx1b', parent = xt.Drift, length = 0.8046839302589)
env.new(name = 'l4x', parent = xt.Drift, length = 1.2561212769403)
env.new(name = 'l4y1', parent = xt.Drift, length = 1.0678775122092)
env.new(name = 'ld500', parent = xt.Drift, length = 0.5)
env.new(name = 'lqzd', parent = xt.Drift, length = 0.188)
env.new(name = 'lqmd', parent = xt.Drift, length = 0.108)
env.new(name = 'l41', parent = xt.Drift, length = 4.091148102808)
env.new(name = 'l42', parent = xt.Drift, length = 4.781148102808)
env.new(name = 'l4e1', parent = xt.Drift, length = 4.711148102808)
env.new(name = 'l44', parent = xt.Drift, length = 2.665416115171)
env.new(name = 'l45', parent = xt.Drift, length = 0.9947)
env.new(name = 'l46', parent = xt.Drift, length = 0.9359965768022)
env.new(name = 'l4b', parent = xt.Drift, length = 1.7506965769027)
env.new(name = 'lbqxb', parent = xt.Drift, length = 0.9633482729526)
env.new(name = 'l47', parent = xt.Drift, length = 5.6966)
env.new(name = 'l48', parent = xt.Drift, length = 6.0666)
env.new(name = 'l49', parent = xt.Drift, length = 5.8747963455943)
env.new(name = 'lw1', parent = xt.Drift, length = 10.651101424581)
env.new(name = 'lmwh', parent = xt.Drift, length = 0.125)
env.new(name = 'lw2', parent = xt.Drift, length = 10.651101424581)
env.new(name = 'lw3', parent = xt.Drift, length = 10.651101424581)
env.new(name = 'lmbth', parent = xt.Drift, length = 0.125)
env.new(name = 'lw4', parent = xt.Drift, length = 10.901101424581)
env.new(name = 'lwa1', parent = xt.Drift, length = 0.8335399845332)
env.new(name = 'lbs1', parent = xt.Drift, length = 7.1704076743589)
env.new(name = 'lx1', parent = xt.Drift, length = 0.8426808408165)
env.new(name = 'lj1', parent = xt.Drift, length = 0.7014546180649)
env.new(name = 'lbb1', parent = xt.Drift, length = 1.2993616690823)
env.new(name = 'lbq11a', parent = xt.Drift, length = 0.4845359633269)
env.new(name = 'lqmac', parent = xt.Drift, length = 0.1365)
env.new(name = 'lbq12a', parent = xt.Drift, length = 0.7348528970567)
env.new(name = 'lbq13a', parent = xt.Drift, length = 0.88469748946)
env.new(name = 'lbq13b', parent = xt.Drift, length = 0.263)
env.new(name = 'lqma', parent = xt.Drift, length = 0.138)
env.new(name = 'lbq12e', parent = xt.Drift, length = 0.489852902178)
env.new(name = 'lbq13', parent = xt.Drift, length = 1.14769748946)
env.new(name = 'lbq12b', parent = xt.Drift, length = 0.739852902178)
env.new(name = 'lbq12c', parent = xt.Drift, length = 0.4848528970567)
env.new(name = 'lbq12bc', parent = xt.Drift, length = 0.739052902178)
env.new(name = 'lka', parent = xt.Drift, length = 2.0728074498198)
env.new(name = 'lkb', parent = xt.Drift, length = 2.1295074498198)
env.new(name = 'lt01', parent = xt.Drift, length = 0.0098262379965)
env.new(name = 'lt1a', parent = xt.Drift, length = 0.3098262379965)
env.new(name = 'lt2', parent = xt.Drift, length = 0.529675)
env.new(name = 'lt21', parent = xt.Drift, length = 5.7224068168709)
env.new(name = 'lt3a', parent = xt.Drift, length = 4.4628158745934)
env.new(name = 'lt3b', parent = xt.Drift, length = 4.7628158745934)
env.new(name = 'lt4a', parent = xt.Drift, length = 2.9113731726384)
env.new(name = 'lt5', parent = xt.Drift, length = 0.561175)
env.new(name = 'lt61', parent = xt.Drift, length = 0.423175)
env.new(name = 'lt62', parent = xt.Drift, length = 0.5613055273338)
env.new(name = 'lbq21a', parent = xt.Drift, length = 0.3678922205778)
env.new(name = 'ld150', parent = xt.Drift, length = 0.138)
env.new(name = 'lbq22a', parent = xt.Drift, length = 0.2999266274684)
env.new(name = 'lbb2', parent = xt.Drift, length = 1.1510059824524)
env.new(name = 'lbq23a', parent = xt.Drift, length = 0.9106820722265)
env.new(name = 'lbq22b', parent = xt.Drift, length = 0.8649413571611)
env.new(name = 'lbq21e', parent = xt.Drift, length = 0.3878697203732)
env.new(name = 'lbq22c', parent = xt.Drift, length = 0.6378140055175)
env.new(name = 'l6a', parent = xt.Drift, length = 0.9637902506258)
env.new(name = 'ld415', parent = xt.Drift, length = 0.415)
env.new(name = 'l61a', parent = xt.Drift, length = 0.301895954985)
env.new(name = 'lbb3', parent = xt.Drift, length = 0.2660264397464)
env.new(name = 'lbq31a', parent = xt.Drift, length = 0.3170847864822)
env.new(name = 'lbq32a', parent = xt.Drift, length = 0.1943664346056)
env.new(name = 'lbq33', parent = xt.Drift, length = 0.3452506224009)
env.new(name = 'lm1a', parent = xt.Drift, length = 0.533114126045)
env.new(name = 'ls1a', parent = xt.Drift, length = 0.4695222767127)
env.new(name = 'ls11a', parent = xt.Drift, length = 4.800676496676401)
env.new(name = 'ls11b', parent = xt.Drift, length = 0.26)
env.new(name = 'ls13', parent = xt.Drift, length = 5.4606764966764)
env.new(name = 'ls14a', parent = xt.Drift, length = 5.1006764966764)
env.new(name = 'ls14b', parent = xt.Drift, length = 0.26)
env.new(name = 'ls12', parent = xt.Drift, length = 5.0006764966764)
env.new(name = 'ls1b', parent = xt.Drift, length = 0.4695222767127)
env.new(name = 'lm177a', parent = xt.Drift, length = 1.8736441555988)
env.new(name = 'lcvbh', parent = xt.Drift, length = 0.125)
env.new(name = 'lmqh1', parent = xt.Drift, length = 3.8545980534172)
env.new(name = 'lmqh2', parent = xt.Drift, length = 3.8545980534172)
env.new(name = 'lmqh3', parent = xt.Drift, length = 4.2545980534172)
env.new(name = 'ld100c', parent = xt.Drift, length = 0.088)
env.new(name = 'lm076b', parent = xt.Drift, length = 0.3612225)
env.new(name = 'lbb4', parent = xt.Drift, length = 0.4462705905835)
env.new(name = 'lbq41a', parent = xt.Drift, length = 0.460575)
env.new(name = 'lmbsh', parent = xt.Drift, length = 0.1)
env.new(name = 'lbq42a', parent = xt.Drift, length = 0.28187)
env.new(name = 'lm076a', parent = xt.Drift, length = 0.3312225)
env.new(name = 'lchbsh', parent = xt.Drift, length = 0.1)
env.new(name = 'lpmsh', parent = xt.Drift, length = 0.05)
env.new(name = 'lmzc', parent = xt.Drift, length = 0.088)
env.new(name = 'lm3a', parent = xt.Drift, length = 2.042319554274)
env.new(name = 'lm55a', parent = xt.Drift, length = 0.438416)
env.new(name = 'lm2', parent = xt.Drift, length = 0.55987)
env.new(name = 'lm1b', parent = xt.Drift, length = 1.9920950786619)
env.new(name = 'lm076c', parent = xt.Drift, length = 0.3812225)
env.new(name = 'lm046a', parent = xt.Drift, length = 0.30187)
env.new(name = 'lm146a', parent = xt.Drift, length = 0.5312225)
env.new(name = 'lmbt1h', parent = xt.Drift, length = 0.15)
env.new(name = 'ld100', parent = xt.Drift, length = 0.1)
env.new(name = 'lmi1a', parent = xt.Drift, length = 3.1084048484628)
env.new(name = 'ld1362', parent = xt.Drift, length = 1.35)
env.new(name = 'lmi3a', parent = xt.Drift, length = 0.8399201544954)
env.new(name = 'lsep1', parent = xt.Drift, length = 2.9533460728692)
env.new(name = 'lsqi6', parent = xt.Drift, length = 0.29714)
env.new(name = 'liqi6', parent = xt.Drift, length = 10.026178863617815)

############################################################
# Bends
############################################################

########################################
# Base Elements
########################################
env.new(name = 'hbend00800080265', parent = xt.Bend, length = 0.8000802659853)
env.new(name = 'hbend00800345567', parent = xt.Bend, length = 0.8003455671831)
env.new(name = 'hbend01051031867', parent = xt.Bend, length = 1.0510318674703)
env.new(name = 'hbend01052186955', parent = xt.Bend, length = 1.052186955893)
env.new(name = 'hbend01053066114', parent = xt.Bend, length = 1.0530661140669)
env.new(name = 'hbend01053415951', parent = xt.Bend, length = 1.0534159512182)
env.new(name = 'hbend01053454617', parent = xt.Bend, length = 1.0534546175758)
env.new(name = 'hbend01249749769', parent = xt.Bend, length = 1.2497497691714)
env.new(name = 'hbend01540814142', parent = xt.Bend, length = 1.5408141428297513)
env.new(name = 'hbend01547076259', parent = xt.Bend, length = 1.5470762593395)
env.new(name = 'hbend01580497288', parent = xt.Bend, length = 1.5804972885128)
env.new(name = 'hbend01836770243', parent = xt.Bend, length = 1.8367702431958)
env.new(name = 'hbend01881065986', parent = xt.Bend, length = 1.881065986441423)
env.new(name = 'hbend02095017481', parent = xt.Bend, length = 2.095017481049)
env.new(name = 'hbend03326168461', parent = xt.Bend, length = 3.326168461962433)
env.new(name = 'hbend03326168462', parent = xt.Bend, length = 3.326168462149516)
env.new(name = 'vbend01053873538', parent = xt.Bend, length = 1.0538735383142, rot_s_rad = +np.pi/2)
env.new(name = 'vbend01844887346', parent = xt.Bend, length = 1.8448873467866, rot_s_rad = +np.pi/2)

########################################
# Cloned Elements
########################################
env.new(
    name                    = 'sept2',
    parent                  = 'hbend00800080265',
    k0                      = 'k0_sept2',
    h                       = 0.05574440702530643,
    edge_entry_angle        = 0.0223,
    edge_exit_angle         = 0.0223)
env.new(
    name                    = 'sept1',
    parent                  = 'hbend00800345567',
    k0                      = 'k0_sept1',
    h                       = 0.06520308494200895,
    edge_entry_angle        = 0.052185)
env.new(
    name                    = 'b3p',
    parent                  = 'hbend01051031867',
    k0                      = 'k0_b3p',
    h                       = -0.09684993860842801,
    edge_entry_angle        = -0.05089618592,
    edge_exit_angle         = -0.05089618592)
env.new(
    name                    = 'bh1p',
    parent                  = 'hbend01052186955',
    k0                      = 'k0_bh1p',
    h                       = 0.09157662435875959,
    edge_entry_angle        = 0.0481778648075,
    edge_exit_angle         = 0.0481778648075)
env.new(
    name                    = 'b2p',
    parent                  = 'hbend01053066114',
    k0                      = 'k0_b2p',
    h                       = -0.08270349250490378,
    edge_entry_angle        = -0.04354612273595,
    edge_exit_angle         = -0.04354612273595)
env.new(
    name                    = 'bh1ap',
    parent                  = 'hbend01053415951',
    k0                      = 'k0_bh1ap',
    h                       = 0.06833618620541379,
    edge_entry_angle        = 0.0359932142971,
    edge_exit_angle         = 0.0359932142971)
env.new(
    name                    = 'bh1bp',
    parent                  = 'hbend01053454617',
    k0                      = 'k0_bh1bp',
    h                       = 0.034367576090286525,
    edge_entry_angle        = 0.0181023408636,
    edge_exit_angle         = 0.0181023408636)
env.new(
    name                    = 'bh2p',
    parent                  = 'hbend01249749769',
    k0                      = 'k0_bh2p',
    h                       = -0.09796271517999312,
    edge_entry_angle        = -0.0612144403418,
    edge_exit_angle         = -0.0612144403418)
env.new(
    name                    = 'bm611p',
    parent                  = 'hbend01540814142',
    k0                      = 'k0_bm611p',
    h                       = 0.10696134841633706,
    edge_exit_angle         = 0.16480755837603275)
env.new(
    name                    = 'bm616p',
    parent                  = 'hbend01540814142',
    k0                      = 'k0_bm616p',
    h                       = 0.10696134841633706,
    edge_entry_angle        = 0.16480755837603275)
env.new(
    name                    = 'b1p',
    parent                  = 'hbend01547076259',
    k0                      = 'k0_b1p',
    h                       = -0.08480036190388614,
    edge_entry_angle        = -0.06559631334245,
    edge_exit_angle         = -0.06559631334245)
env.new(
    name                    = 'bh1cp',
    parent                  = 'hbend01580497288',
    k0                      = 'k0_bh1cp',
    h                       = 0.10249719192136915,
    edge_entry_angle        = 0.08099826695595,
    edge_exit_angle         = 0.08099826695595)
env.new(
    name                    = 'bh4p',
    parent                  = 'hbend01836770243',
    k0                      = 'k0_bh4p',
    h                       = 0.045312528175647175,
    edge_entry_angle        = 0.0416143516985,
    edge_exit_angle         = 0.0416143516985)
env.new(
    name                    = 'bm612p',
    parent                  = 'hbend01881065986',
    k0                      = 'k0_bm612p',
    h                       = 0.11249205912579933,
    edge_entry_angle        = -0.16480755837603275,
    edge_exit_angle         = 0.37641254454233136)
env.new(
    name                    = 'bm615p',
    parent                  = 'hbend01881065986',
    k0                      = 'k0_bm615p',
    h                       = 0.11249205912579933,
    edge_entry_angle        = 0.37641254454233136,
    edge_exit_angle         = -0.16480755837603275)
env.new(
    name                    = 'bh3p',
    parent                  = 'hbend02095017481',
    k0                      = 'k0_bh3p',
    h                       = 0.10114441782991687,
    edge_entry_angle        = 0.1059496617321,
    edge_exit_angle         = 0.1059496617321)
env.new(
    name                    = 'bm614p',
    parent                  = 'hbend03326168461',
    k0                      = 'k0_bm614p',
    h                       = -0.1131670114869193,
    edge_exit_angle         = -0.37641254454233136)
env.new(
    name                    = 'bm613p',
    parent                  = 'hbend03326168462',
    k0                      = 'k0_bm613p',
    h                       = -0.11316701148055414,
    edge_entry_angle        = -0.37641254454233136)
env.new(
    name                    = 'bv1up',
    parent                  = 'vbend01053873538',
    k0                      = 'k0_bv1up',
    h                       = -0.05794454347405184,
    edge_entry_angle        = -0.0305331105285,
    edge_exit_angle         = -0.0305331105285)
env.new(
    name                    = 'bv1dp',
    parent                  = 'vbend01053873538',
    k0                      = 'k0_bv1dp',
    h                       = 0.05794454347405184,
    edge_entry_angle        = 0.0305331105285,
    edge_exit_angle         = 0.0305331105285)
env.new(
    name                    = 'bv2up',
    parent                  = 'vbend01844887346',
    k0                      = 'k0_bv2up',
    h                       = -0.09992257066529898,
    edge_entry_angle        = -0.0921729431394,
    edge_exit_angle         = -0.0921729431394)
env.new(
    name                    = 'bv2dp',
    parent                  = 'vbend01844887346',
    k0                      = 'k0_bv2dp',
    h                       = 0.09992257066529898,
    edge_entry_angle        = 0.0921729431394,
    edge_exit_angle         = 0.0921729431394)

############################################################
# Correctors
############################################################

########################################
# Base Elements
########################################
env.new(name = 'hcorr00200000000', parent = xt.Bend, length = 0.2)
env.new(name = 'hcorr00992261064', parent = xt.Bend, length = 0.992261064834161)
env.new(name = 'vcorr00200000000', parent = xt.Bend, length = 0.2, rot_s_rad = +np.pi/2)

########################################
# Cloned Elements
########################################
env.new(name = 'bx616', parent = 'hcorr00200000000', k0 = 'k0_bx616')
env.new(name = 'hx03p', parent = 'hcorr00200000000', k0 = 'k0_hx03p')
env.new(name = 'hx04p', parent = 'hcorr00200000000', k0 = 'k0_hx04p')
env.new(name = 'hx05p', parent = 'hcorr00200000000', k0 = 'k0_hx05p')
env.new(name = 'hx06p', parent = 'hcorr00200000000', k0 = 'k0_hx06p')
env.new(name = 'hx07p', parent = 'hcorr00200000000', k0 = 'k0_hx07p')
env.new(name = 'hx08p', parent = 'hcorr00200000000', k0 = 'k0_hx08p')
env.new(name = 'hx09p', parent = 'hcorr00200000000', k0 = 'k0_hx09p')
env.new(name = 'hw10p', parent = 'hcorr00200000000', k0 = 'k0_hw10p')
env.new(name = 'hw11p', parent = 'hcorr00200000000', k0 = 'k0_hw11p')
env.new(name = 'hw12p', parent = 'hcorr00200000000', k0 = 'k0_hw12p')
env.new(name = 'ht13p', parent = 'hcorr00200000000', k0 = 'k0_ht13p')
env.new(name = 'ht14p', parent = 'hcorr00200000000', k0 = 'k0_ht14p')
env.new(name = 'hb15p', parent = 'hcorr00200000000', k0 = 'k0_hb15p')
env.new(name = 'hb16p', parent = 'hcorr00200000000', k0 = 'k0_hb16p')
env.new(name = 'hm17p', parent = 'hcorr00200000000', k0 = 'k0_hm17p')
env.new(name = 'hm18p', parent = 'hcorr00200000000', k0 = 'k0_hm18p')
env.new(name = 'hm19p', parent = 'hcorr00200000000', k0 = 'k0_hm19p')
env.new(name = 'hm20p', parent = 'hcorr00200000000', k0 = 'k0_hm20p')
env.new(name = 'hm21p', parent = 'hcorr00200000000', k0 = 'k0_hm21p')
env.new(name = 'bp581', parent = 'hcorr00992261064', k0 = 'k0_bp581')
env.new(name = 'by616', parent = 'vcorr00200000000', k0 = 'k0_by616')
env.new(name = 'by618', parent = 'vcorr00200000000', k0 = 'k0_by618')
env.new(name = 'vx03p', parent = 'vcorr00200000000', k0 = 'k0_vx03p')
env.new(name = 'vx04p', parent = 'vcorr00200000000', k0 = 'k0_vx04p')
env.new(name = 'vx05p', parent = 'vcorr00200000000', k0 = 'k0_vx05p')
env.new(name = 'vx06p', parent = 'vcorr00200000000', k0 = 'k0_vx06p')
env.new(name = 'vx07p', parent = 'vcorr00200000000', k0 = 'k0_vx07p')
env.new(name = 'vx08p', parent = 'vcorr00200000000', k0 = 'k0_vx08p')
env.new(name = 'vx09p', parent = 'vcorr00200000000', k0 = 'k0_vx09p')
env.new(name = 'vw10p', parent = 'vcorr00200000000', k0 = 'k0_vw10p')
env.new(name = 'vw11p', parent = 'vcorr00200000000', k0 = 'k0_vw11p')
env.new(name = 'vw12p', parent = 'vcorr00200000000', k0 = 'k0_vw12p')
env.new(name = 'vw13p', parent = 'vcorr00200000000', k0 = 'k0_vw13p')
env.new(name = 'va14p', parent = 'vcorr00200000000', k0 = 'k0_va14p')
env.new(name = 'va15p', parent = 'vcorr00200000000', k0 = 'k0_va15p')
env.new(name = 'va16p', parent = 'vcorr00200000000', k0 = 'k0_va16p')
env.new(name = 'va17p', parent = 'vcorr00200000000', k0 = 'k0_va17p')
env.new(name = 'va18p', parent = 'vcorr00200000000', k0 = 'k0_va18p')
env.new(name = 'vt19p', parent = 'vcorr00200000000', k0 = 'k0_vt19p')
env.new(name = 'vt20p', parent = 'vcorr00200000000', k0 = 'k0_vt20p')
env.new(name = 'vb21p', parent = 'vcorr00200000000', k0 = 'k0_vb21p')
env.new(name = 'vb22p', parent = 'vcorr00200000000', k0 = 'k0_vb22p')
env.new(name = 'vc23p', parent = 'vcorr00200000000', k0 = 'k0_vc23p')
env.new(name = 'vc24p', parent = 'vcorr00200000000', k0 = 'k0_vc24p')
env.new(name = 'vc25p', parent = 'vcorr00200000000', k0 = 'k0_vc25p')
env.new(name = 'vm26p', parent = 'vcorr00200000000', k0 = 'k0_vm26p')
env.new(name = 'vm27p', parent = 'vcorr00200000000', k0 = 'k0_vm27p')
env.new(name = 'vm28p', parent = 'vcorr00200000000', k0 = 'k0_vm28p')
env.new(name = 'vm29p', parent = 'vcorr00200000000', k0 = 'k0_vm29p')
env.new(name = 'vm30p', parent = 'vcorr00200000000', k0 = 'k0_vm30p')
env.new(name = 'vm31p', parent = 'vcorr00200000000', k0 = 'k0_vm31p')
env.new(name = 'vm32p', parent = 'vcorr00200000000', k0 = 'k0_vm32p')
env.new(name = 'vm33p', parent = 'vcorr00200000000', k0 = 'k0_vm33p')

############################################################
# Quadrupoles
############################################################

########################################
# Base Elements
########################################
env.new(name = 'quad00204000000', parent = xt.Quadrupole, length = 0.204)
env.new(name = 'quad00374000000', parent = xt.Quadrupole, length = 0.374)
env.new(name = 'quad00387000000', parent = xt.Quadrupole, length = 0.387)
env.new(name = 'quad00460780000', parent = xt.Quadrupole, length = 0.46078)
env.new(name = 'quad00524000000', parent = xt.Quadrupole, length = 0.524)
env.new(name = 'quad00583720000', parent = xt.Quadrupole, length = 0.58372)

########################################
# Cloned Elements
########################################
env.new(name = 'qf564', parent = 'quad00204000000', k1 = 'k1_qf564')
env.new(name = 'qd564', parent = 'quad00204000000', k1 = 'k1_qd564')
env.new(name = 'qf584', parent = 'quad00204000000', k1 = 'k1_qf584')
env.new(name = 'qd584', parent = 'quad00204000000', k1 = 'k1_qd584')
env.new(name = 'qd616', parent = 'quad00374000000', k1 = 'k1_qd616')
env.new(name = 'qf616', parent = 'quad00374000000', k1 = 'k1_qf616')
env.new(name = 'qxd1p', parent = 'quad00387000000', k1 = 'k1_qxd1p')
env.new(name = 'qxf1p', parent = 'quad00387000000', k1 = 'k1_qxf1p')
env.new(name = 'qxd2p', parent = 'quad00387000000', k1 = 'k1_qxd2p')
env.new(name = 'qxf2p', parent = 'quad00387000000', k1 = 'k1_qxf2p')
env.new(name = 'qxd3p', parent = 'quad00387000000', k1 = 'k1_qxd3p')
env.new(name = 'qxd4p', parent = 'quad00387000000', k1 = 'k1_qxd4p')
env.new(name = 'qxf6p', parent = 'quad00387000000', k1 = 'k1_qxf6p')
env.new(name = 'qxd7p', parent = 'quad00387000000', k1 = 'k1_qxd7p')
env.new(name = 'qwfp', parent = 'quad00387000000', k1 = 'k1_qwfp')
env.new(name = 'qwdp', parent = 'quad00387000000', k1 = 'k1_qwdp')
env.new(name = 'qad1p', parent = 'quad00387000000', k1 = 'k1_qad1p')
env.new(name = 'qad3p', parent = 'quad00387000000', k1 = 'k1_qad3p')
env.new(name = 'qad5p', parent = 'quad00387000000', k1 = 'k1_qad5p')
env.new(name = 'qad7p', parent = 'quad00387000000', k1 = 'k1_qad7p')
env.new(name = 'qad9p', parent = 'quad00387000000', k1 = 'k1_qad9p')
env.new(name = 'qad11p', parent = 'quad00387000000', k1 = 'k1_qad11p')
env.new(name = 'qtf1p', parent = 'quad00387000000', k1 = 'k1_qtf1p')
env.new(name = 'qtf3p', parent = 'quad00387000000', k1 = 'k1_qtf3p')
env.new(name = 'qbd4p', parent = 'quad00387000000', k1 = 'k1_qbd4p')
env.new(name = 'qcd2p', parent = 'quad00387000000', k1 = 'k1_qcd2p')
env.new(name = 'qcd4p', parent = 'quad00387000000', k1 = 'k1_qcd4p')
env.new(name = 'qcd6p', parent = 'quad00387000000', k1 = 'k1_qcd6p')
env.new(name = 'qi7p', parent = 'quad00460780000', k1 = 'k1_qi7p')
env.new(name = 'qf618', parent = 'quad00524000000', k1 = 'k1_qf618')
env.new(name = 'qd618', parent = 'quad00524000000', k1 = 'k1_qd618')
env.new(name = 'qxf3p', parent = 'quad00524000000', k1 = 'k1_qxf3p')
env.new(name = 'qxf4p', parent = 'quad00524000000', k1 = 'k1_qxf4p')
env.new(name = 'qxd5p', parent = 'quad00524000000', k1 = 'k1_qxd5p')
env.new(name = 'qxf5p', parent = 'quad00524000000', k1 = 'k1_qxf5p')
env.new(name = 'qxd6p', parent = 'quad00524000000', k1 = 'k1_qxd6p')
env.new(name = 'qaf2p', parent = 'quad00524000000', k1 = 'k1_qaf2p')
env.new(name = 'qaf4p', parent = 'quad00524000000', k1 = 'k1_qaf4p')
env.new(name = 'qaf6p', parent = 'quad00524000000', k1 = 'k1_qaf6p')
env.new(name = 'qaf8p', parent = 'quad00524000000', k1 = 'k1_qaf8p')
env.new(name = 'qaf10p', parent = 'quad00524000000', k1 = 'k1_qaf10p')
env.new(name = 'qtd2p', parent = 'quad00524000000', k1 = 'k1_qtd2p')
env.new(name = 'qtd4p', parent = 'quad00524000000', k1 = 'k1_qtd4p')
env.new(name = 'qbf1p', parent = 'quad00524000000', k1 = 'k1_qbf1p')
env.new(name = 'qbd2p', parent = 'quad00524000000', k1 = 'k1_qbd2p')
env.new(name = 'qbf3p', parent = 'quad00524000000', k1 = 'k1_qbf3p')
env.new(name = 'qcf1p', parent = 'quad00524000000', k1 = 'k1_qcf1p')
env.new(name = 'qcf3p', parent = 'quad00524000000', k1 = 'k1_qcf3p')
env.new(name = 'qcf5p', parent = 'quad00524000000', k1 = 'k1_qcf5p')
env.new(name = 'qmf1p', parent = 'quad00524000000', k1 = 'k1_qmf1p')
env.new(name = 'qmd2p', parent = 'quad00524000000', k1 = 'k1_qmd2p')
env.new(name = 'qmd3p', parent = 'quad00524000000', k1 = 'k1_qmd3p')
env.new(name = 'qmf4p', parent = 'quad00524000000', k1 = 'k1_qmf4p')
env.new(name = 'qmd5p', parent = 'quad00524000000', k1 = 'k1_qmd5p')
env.new(name = 'qmf6p', parent = 'quad00524000000', k1 = 'k1_qmf6p')
env.new(name = 'qmd7p', parent = 'quad00524000000', k1 = 'k1_qmd7p')
env.new(name = 'qmf8p', parent = 'quad00524000000', k1 = 'k1_qmf8p')
env.new(name = 'qmd9p', parent = 'quad00524000000', k1 = 'k1_qmd9p')
env.new(name = 'qmf10p', parent = 'quad00524000000', k1 = 'k1_qmf10p')
env.new(name = 'qmd11p', parent = 'quad00524000000', k1 = 'k1_qmd11p')
env.new(name = 'qmf12p', parent = 'quad00524000000', k1 = 'k1_qmf12p')
env.new(name = 'qmd13p', parent = 'quad00524000000', k1 = 'k1_qmd13p')
env.new(name = 'qi6p', parent = 'quad00583720000', k1 = 'k1_qi6p')

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
env.new(
    name        = 'smd1p',
    parent      = 'sext00100000000',
    k2          = 'k2_smd1p',
    rot_s_rad   = '1.5707963267948966')
env.new(
    name        = 'smd2p',
    parent      = 'sext00100000000',
    k2          = 'k2_smd2p',
    rot_s_rad   = '1.5707963267948966')

############################################################
# Cavities
############################################################
env.new(
    name        = 'ak571',
    parent      = xt.Cavity,
    length      = 0.039,
    frequency   = 'freq_ak571 * (1 + fshift)',
    voltage     = 'volt_ak571',
    lag         = 'lag_ak571')
env.new(
    name        = 'ac571',
    parent      = xt.Cavity,
    length      = 0.5598368,
    frequency   = 'freq_ac571 * (1 + fshift)',
    voltage     = 'volt_ac571',
    lag         = 'lag_ac571')
env.new(
    name        = 'ac5712',
    parent      = xt.Cavity,
    length      = 1.3296124,
    frequency   = 'freq_ac5712 * (1 + fshift)',
    voltage     = 'volt_ac5712',
    lag         = 'lag_ac5712')
env.new(
    name        = 'ak572',
    parent      = xt.Cavity,
    length      = 0.039,
    frequency   = 'freq_ak572 * (1 + fshift)',
    voltage     = 'volt_ak572',
    lag         = 'lag_ak572')
env.new(
    name        = 'ac572',
    parent      = xt.Cavity,
    length      = 1.8894492,
    frequency   = 'freq_ac572 * (1 + fshift)',
    voltage     = 'volt_ac572',
    lag         = 'lag_ac572')
env.new(
    name        = 'ak573',
    parent      = xt.Cavity,
    length      = 0.039,
    frequency   = 'freq_ak573 * (1 + fshift)',
    voltage     = 'volt_ak573',
    lag         = 'lag_ak573')
env.new(
    name        = 'ac573',
    parent      = xt.Cavity,
    length      = 0.5598368,
    frequency   = 'freq_ac573 * (1 + fshift)',
    voltage     = 'volt_ac573',
    lag         = 'lag_ac573')
env.new(
    name        = 'ac5732',
    parent      = xt.Cavity,
    length      = 1.3296124,
    frequency   = 'freq_ac5732 * (1 + fshift)',
    voltage     = 'volt_ac5732',
    lag         = 'lag_ac5732')
env.new(
    name        = 'ak574',
    parent      = xt.Cavity,
    length      = 0.039,
    frequency   = 'freq_ak574 * (1 + fshift)',
    voltage     = 'volt_ak574',
    lag         = 'lag_ak574')
env.new(
    name        = 'ac574',
    parent      = xt.Cavity,
    length      = 1.8894492,
    frequency   = 'freq_ac574 * (1 + fshift)',
    voltage     = 'volt_ac574',
    lag         = 'lag_ac574')
env.new(
    name        = 'ak617',
    parent      = xt.Cavity,
    length      = 0.039,
    frequency   = 'freq_ak617 * (1 + fshift)',
    voltage     = 'volt_ak617',
    lag         = 'lag_ak617')
env.new(
    name        = 'ac617',
    parent      = xt.Cavity,
    length      = 1.8894492,
    frequency   = 'freq_ac617 * (1 + fshift)',
    voltage     = 'volt_ac617',
    lag         = 'lag_ac617')
env.new(
    name        = 'ak618',
    parent      = xt.Cavity,
    length      = 0.039,
    frequency   = 'freq_ak618 * (1 + fshift)',
    voltage     = 'volt_ak618',
    lag         = 'lag_ak618')
env.new(
    name        = 'ac618',
    parent      = xt.Cavity,
    length      = 1.8894492,
    frequency   = 'freq_ac618 * (1 + fshift)',
    voltage     = 'volt_ac618',
    lag         = 'lag_ac618')

############################################################
# Markers
############################################################
ALL_MARKERS = [
    'bs611', 'bs612', 'bs613', 'bs614', 'bs615', 'bs616', 'bx584', 'bx618',
    'by584', 'cbp', 'chp', 'cvp', 'fha121p', 'fha122p', 'fha123p', 'fha124p',
    'fha125p', 'fha126p', 'fhc161p', 'fhc162p', 'fhm201p', 'fhm202p', 'fhm203p',
    'fhx041p', 'fhx071p', 'fhx072p', 'ftcell', 'fvc', 'fvm271p', 'fvt201p',
    'g1', 'g2', 'g5', 'g6', 'ga1i', 'ga2i', 'ga3i', 'gki', 'gko', 'gsi', 'gso',
    'injp', 'mbm611p', 'mbm612p', 'mbm613p', 'mbm614p', 'mbm615p', 'mbm616p',
    'ms', 'msli', 'mt', 'mwp', 'p6sect', 'pac571', 'pac572', 'pac573', 'pac574',
    'pbt', 'pbtp', 'pinjax0', 'pm574', 'pstart', 'sp564', 'sp580', 'sp584',
    'spqad11p_a', 'spqad1p_a', 'spqad3p_a', 'spqad5p_a', 'spqad7p_a',
    'spqad9p_a', 'spqaf10p_k', 'spqaf2p_k', 'spqaf4p_k', 'spqaf6p_k',
    'spqaf8p_k', 'spqbd2p_k', 'spqbd4p_a', 'spqbf1p_k', 'spqbf3p_k',
    'spqcd2p_a', 'spqcd4p_a', 'spqcd6p_a', 'spqcf1p_k', 'spqcf3p_k',
    'spqcf5p_k', 'spqmd11p_k', 'spqmd13p_k', 'spqmd2p_1k', 'spqmd2p_2k',
    'spqmd3p_k', 'spqmd5p_k', 'spqmd7p_k', 'spqmd9p_k', 'spqmf10p_k',
    'spqmf12p_k', 'spqmf1p_1k', 'spqmf1p_2k', 'spqmf1p_3k', 'spqmf4p_k',
    'spqmf6p_k', 'spqmf8p_k', 'spqtd2p_k', 'spqtd4p_k', 'spqtf1p_a',
    'spqtf3p_a', 'spqwdp_1a', 'spqwdp_2a', 'spqwdp_3a', 'spqwfp_1a',
    'spqwfp_2a', 'spqwfp_3a', 'spqwfp_4a', 'spqxd1p_a', 'spqxd2p_a',
    'spqxd3p_a', 'spqxd4p_a', 'spqxd5p_k', 'spqxd6p_k', 'spqxd7p_a',
    'spqxf1p_a', 'spqxf2p_a', 'spqxf3p_k', 'spqxf4p_k', 'spqxf5p_k',
    'spqxf6p_a', 'sx571', 'sx573', 'sy571', 'sy573']
for marker in ALL_MARKERS:
    env.new(name = marker, parent = xt.Marker)

############################################################
# Create Line
############################################################
env.new_line(
    name        = 'line',
    components  = [
        'pstart', 'sp564', 'll3373', 'll3374', 'qf564', 'll3376', 'lqf564',
        'll3378', 'lqf564', 'll3380', 'll3381', 'll3382', 'qd564', 'll3384',
        'gv570', 'll3386', 'ak571', 'pac571', 'ac571', 'sx571', 'sy571',
        'ac5712', 'pac571', 'ak571', 'll3393', 'ak572', 'pac572', 'ac572',
        'pac572', 'ak572', 'll3397', 'll3398', 'll3399', 'ws572', 'll3401',
        'll3402', 'sc572', 'll3404', 'll3405', 'll3406', 'ak573', 'pac573',
        'ac573', 'sx573', 'sy573', 'ac5732', 'pac573', 'ak573', 'll3413',
        'ak574', 'pac574', 'ac574', 'pac574', 'ak574', 'pm574', 'gv574',
        'sp580', 'pbt', 'ln001', 'ln002', 'bp581', 'ln003', 'ln004', 'ln005',
        'ln006', 'ln007', 'ln008', 'by584', 'ln009', 'bx584', 'ln010', 'ln011',
        'sc584', 'ln012', 'sp584', 'll3464', 'll3465', 'qf584', 'll3467',
        'lqf584', 'll3469', 'lqf584', 'll3471', 'll3472', 'll3473', 'qd584',
        'll3475', 'gv584', 'll3477', 'p6sect', 'll3478e', 'pbtp', 'bm611p',
        'mbm611p', 'bs611', 'l3481p', 'l3482', 'l3483p', 'mbm612p', 'bm612p',
        'mbm612p', 'bs612', 'l3486p', 'l3487', 'sc612', 'l3489', 'l3490',
        'l3491p', 'l3492', 'sr613', 'mbm613p', 'bm613p', 'mbm613p', 'bs613',
        'l3497p', 'l3498', 'sp613', 'l3500', 'l3501', 'sl613', 'l3503', 'l3504',
        'sc613', 'l3506', 'l3507p', 'mbm614p', 'bm614p', 'mbm614p', 'bs614',
        'l3510p', 'l3511', 'l3512', 'sc614', 'l3514', 'l3515p', 'mbm615p',
        'bm615p', 'mbm615p', 'bs615', 'l3518p', 'l3519', 'l3520p', 'mbm616p',
        'bm616p', 'mbm616p', 'bs616', 'l3523p', 'l3524', 'bx616', 'l3526',
        'by616', 'l3528', 'l3529', 'sp616', 'l3531', 'qd616', 'l3533', 'l3534',
        'qf616', 'l3536', 'l3537', 'wm616', 'l3539', 'sc616', 'l3541', 'l3542',
        'gv616', 'l3544', 'l3545', 'ak617', 'ac617', 'ak617', 'l3549', 'l3550',
        'l3551', 'ak618', 'ac618', 'ak618', 'l3555', 'l3556', 'gv618', 'l3558',
        'l3559', 'sp618', 'l3561', 'qf618', 'l3563', 'l3564', 'by618', 'l3566',
        'l3567', 'l3568', 'bx618', 'l3570', 'l3571', 'qd618', 'ftcell',
        'lcq2bc', 'lpmlih', 'msli', 'lpmlih', 'l12', 'b1p', 'lcq31', 'ld250',
        'lpmlh', 'ms', 'lpmlh', 'lbell', 'lmblh', 'spqxd1p_a', 'lmblh', 'lqmc',
        'qxd1p', 'lqzc', 'vx03p', 'l2a1', 'lmblh', 'spqxf1p_a', 'lmblh', 'lqmc',
        'qxf1p', 'lqzc', 'hx03p', 'l2b1', 'ld250', 'lcbbh', 'cbp', 'lcbbh',
        'lbell', 'lmblh', 'spqxd2p_a', 'lmblh', 'lqmc', 'qxd2p', 'lqzc',
        'vx04p', 'l2c1', 'hx04p', 'lmz', 'lmblh', 'spqxf2p_a', 'lmblh', 'lqmc',
        'qxf2p', 'lqcc', 'lchbh', 'chp', 'lchbh', 'lpmlh', 'ms', 'lpmlh',
        'l2d1', 'lchbh', 'chp', 'lchbh', 'ld780', 'b2p', 'lbqx1a', 'vx05p',
        'lqzc', 'g1', 'qxd3p', 'lqmc', 'lmblh', 'spqxd3p_a', 'lmblh', 'lbqx1b',
        'b2p', 'l4x', 'g2', 'l4y1', 'ld500', 'lmz', 'hx05p', 'lqzd', 'qxf3p',
        'lqmd', 'lmblh', 'spqxf3p_k', 'lmblh', 'l41', 'ld250', 'lpmlh', 'ms',
        'lpmlh', 'lbell', 'lmz', 'vx06p', 'lqzc', 'qxd4p', 'lqmc', 'lmblh',
        'spqxd4p_a', 'lmblh', 'l42', 'hx06p', 'lqzd', 'qxf4p', 'lqmd', 'lmblh',
        'spqxf4p_k', 'lmblh', 'l4e1', 'vx07p', 'lqzd', 'qxd5p', 'lqmd', 'lmblh',
        'spqxd5p_k', 'lmblh', 'l44', 'b3p', 'l45', 'hx07p', 'lqzd', 'qxf5p',
        'lqmd', 'lmblh', 'spqxf5p_k', 'lmblh', 'l46', 'b3p', 'l4b', 'b3p',
        'lbqxb', 'qxd6p', 'lqmd', 'lmblh', 'spqxd6p_k', 'lmblh', 'lmz', 'vx08p',
        'l47', 'hx08p', 'lqzc', 'qxf6p', 'lqmc', 'lmblh', 'spqxf6p_a', 'lmblh',
        'l48', 'vx09p', 'lqzc', 'qxd7p', 'lqmc', 'lmblh', 'spqxd7p_a', 'lmblh',
        'l49', 'hx09p', 'lqzc', 'qwfp', 'lqmc', 'lmblh', 'spqwfp_1a', 'lmblh',
        'lbell', 'lpmlh', 'ms', 'lpmlh', 'ld250', 'lw1', 'vw10p', 'lqzc',
        'qwdp', 'lqmc', 'lmblh', 'spqwdp_1a', 'lmblh', 'lbell', 'ld250', 'lmwh',
        'mwp', 'lmwh', 'lw2', 'hw10p', 'lqzc', 'qwfp', 'lqmc', 'lmblh',
        'spqwfp_2a', 'lmblh', 'lbell', 'ld250', 'lmwh', 'mwp', 'lmwh', 'lw3',
        'vw11p', 'lqzc', 'qwdp', 'lqmc', 'lmblh', 'spqwdp_2a', 'lmblh', 'lbell',
        'ld250', 'lmwh', 'mwp', 'lmwh', 'lw2', 'hw11p', 'lqzc', 'qwfp', 'lqmc',
        'lmblh', 'spqwfp_3a', 'lmblh', 'lbell', 'ld250', 'lmwh', 'mwp', 'lmwh',
        'lw3', 'vw12p', 'lqzc', 'qwdp', 'lqmc', 'lmblh', 'spqwdp_3a', 'lmblh',
        'lbell', 'ld500', 'lw2', 'hw12p', 'lqzc', 'qwfp', 'lqmc', 'lmblh',
        'spqwfp_4a', 'lmblh', 'lbell', 'lmbth', 'mt', 'lmbth', 'lw4', 'vw13p',
        'lqzc', 'qad1p', 'lqmc', 'lmblh', 'spqad1p_a', 'lmblh', 'lbell',
        'lpmlh', 'ms', 'lpmlh', 'lwa1', 'ga1i', 'bh1ap', 'lbs1', 'bh1p', 'lx1',
        'qaf2p', 'lqmd', 'lmblh', 'spqaf2p_k', 'lmblh', 'lj1', 'bh1p', 'lbb1',
        'bh1p', 'lbq11a', 'va14p', 'lqzc', 'qad3p', 'lqmac', 'lmblh',
        'spqad3p_a', 'lmblh', 'lbq12a', 'bh1p', 'lbb1', 'bh1p', 'lbq13a', 'chp',
        'lbq13b', 'qaf4p', 'lqma', 'lpmlh', 'ms', 'lpmlh', 'lmblh', 'spqaf4p_k',
        'lmblh', 'lbq12e', 'bh1p', 'lbb1', 'bh1p', 'lbq11a', 'va15p', 'lqzc',
        'qad5p', 'lqmac', 'lmblh', 'spqad5p_a', 'lmblh', 'lbq12a', 'bh1p',
        'lbb1', 'bh1p', 'lbq13', 'qaf6p', 'lqma', 'lmblh', 'spqaf6p_k', 'lmblh',
        'lbq12b', 'bh1p', 'lbb1', 'bh1p', 'lbq11a', 'va16p', 'lqzc', 'qad7p',
        'lqmac', 'lmblh', 'spqad7p_a', 'lmblh', 'lpmlh', 'ms', 'lpmlh',
        'lbq12c', 'bh1p', 'lbb1', 'bh1p', 'lbq13', 'qaf8p', 'lqma', 'lmblh',
        'spqaf8p_k', 'lmblh', 'lbq12b', 'bh1p', 'lbb1', 'bh1p', 'lbq11a',
        'va17p', 'lqzc', 'qad9p', 'lqmac', 'lmblh', 'spqad9p_a', 'lmblh',
        'lbq12a', 'bh1p', 'lbb1', 'bh1p', 'lbq13', 'qaf10p', 'lqma', 'lmblh',
        'spqaf10p_k', 'lmblh', 'lbq12bc', 'bh1bp', 'gki', 'lka', 'va18p',
        'lqzc', 'qad11p', 'lqmc', 'lmblh', 'spqad11p_a', 'lmblh', 'lkb', 'gko',
        'lt01', 'bh1cp', 'lt1a', 'lmblh', 'spqtf1p_a', 'lmblh', 'lqmc', 'qtf1p',
        'lt2', 'bv1up', 'lt21', 'qtd2p', 'lqmd', 'lmblh', 'spqtd2p_k', 'lmblh',
        'lmz', 'vt19p', 'lt3a', 'ht13p', 'lqzc', 'qtf3p', 'lqmc', 'lmblh',
        'spqtf3p_a', 'lmblh', 'lt3b', 'vt20p', 'lqzd', 'qtd4p', 'lqmd', 'lmblh',
        'spqtd4p_k', 'lmblh', 'lbell', 'lpmlh', 'ms', 'lpmlh', 'lt4a', 'ht14p',
        'lmz', 'lmblh', 'spqbf1p_k', 'lmblh', 'lqmd', 'qbf1p', 'lt5', 'bv1dp',
        'lt61', 'g5', 'lt62', 'ga2i', 'bh2p', 'lbq21a', 'vb21p', 'ld150',
        'qbd2p', 'lqmd', 'lmblh', 'spqbd2p_k', 'lmblh', 'lbq22a', 'bh2p',
        'lbb2', 'bh2p', 'lbq23a', 'hb15p', 'lqzd', 'qbf3p', 'lqmd', 'lmblh',
        'spqbf3p_k', 'lmblh', 'lpmlh', 'ms', 'lpmlh', 'lbq22b', 'bh2p', 'lbb2',
        'bh2p', 'lbq21e', 'vb22p', 'lqzc', 'qbd4p', 'lqmc', 'lmblh',
        'spqbd4p_a', 'lmblh', 'lbq22c', 'bh2p', 'lbb2', 'bh2p', 'l6a', 'hb16p',
        'lqzd', 'qcf1p', 'lqmd', 'lmblh', 'spqcf1p_k', 'lmblh', 'ld415',
        'lpmlh', 'ms', 'lpmlh', 'l61a', 'ga3i', 'bh3p', 'lbb3', 'bh3p',
        'lbq31a', 'vc23p', 'lqzc', 'qcd2p', 'lqmc', 'lmblh', 'spqcd2p_a',
        'lmblh', 'lbq32a', 'bh3p', 'lbb3', 'bh3p', 'lbq33', 'qcf3p', 'lqmd',
        'lmblh', 'spqcf3p_k', 'lmblh', 'lbq32a', 'bh3p', 'lbb3', 'bh3p',
        'lbq31a', 'vc24p', 'lqzc', 'qcd4p', 'lqmc', 'lmblh', 'spqcd4p_a',
        'lmblh', 'lbq32a', 'bh3p', 'lbb3', 'bh3p', 'lbq33', 'qcf5p', 'lqmd',
        'lmblh', 'spqcf5p_k', 'lmblh', 'lbq32a', 'bh3p', 'lbb3', 'g6', 'bh3p',
        'lbq31a', 'vc25p', 'lqzc', 'qcd6p', 'lqmc', 'lmblh', 'spqcd6p_a',
        'lmblh', 'lbq32a', 'bh3p', 'gsi', 'lm1a', 'bv2up', 'ls1a', 'fvc',
        'qmf1p', 'lqmd', 'lmblh', 'spqmf1p_1k', 'lmblh', 'lmz', 'hm17p',
        'ld500', 'ls11a', 'smd1p', 'ls11b', 'vm26p', 'lqzd', 'qmd2p', 'lqmd',
        'lmblh', 'spqmd2p_1k', 'lmblh', 'lpmlh', 'ms', 'lpmlh', 'ld250', 'ls13',
        'hm18p', 'lqzd', 'qmf1p', 'lqmd', 'lmblh', 'spqmf1p_2k', 'lmblh',
        'lchbh', 'chp', 'lchbh', 'lpmlh', 'ms', 'lpmlh', 'ls14a', 'smd2p',
        'ls14b', 'vm27p', 'lqzd', 'qmd2p', 'lqmd', 'lmblh', 'spqmd2p_2k',
        'lmblh', 'ls12', 'lchbh', 'chp', 'lchbh', 'lpmlh', 'ms', 'lpmlh',
        'lbell', 'lmz', 'hm19p', 'lmz', 'lmblh', 'spqmf1p_3k', 'lmblh', 'lqmd',
        'qmf1p', 'ls1b', 'bv2dp', 'gso', 'lm177a', 'lcvbh', 'cvp', 'lcvbh',
        'lpmlh', 'ms', 'lpmlh', 'lmz', 'vm28p', 'lqzd', 'qmd3p', 'lqmd',
        'lmblh', 'spqmd3p_k', 'lmblh', 'lmqh1', 'ld500', 'lmz', 'hm20p', 'lqzd',
        'qmf4p', 'lqmd', 'lmblh', 'spqmf4p_k', 'lmblh', 'lmqh2', 'lcvbh', 'cvp',
        'lcvbh', 'ld250', 'lmz', 'vm29p', 'lqzd', 'qmd5p', 'lqmd', 'lmblh',
        'spqmd5p_k', 'lmblh', 'lmqh3', 'lmwh', 'mwp', 'lmwh', 'lpmlh', 'ms',
        'lpmlh', 'ld100c', 'qmf6p', 'lqmd', 'lmblh', 'spqmf6p_k', 'lmblh',
        'lm076b', 'bh4p', 'lbb4', 'bh4p', 'lbq41a', 'vm30p', 'lqzd', 'qmd7p',
        'lqmd', 'lmbsh', 'spqmd7p_k', 'lmbsh', 'lbq42a', 'bh4p', 'lbb4', 'bh4p',
        'lm076a', 'lchbsh', 'chp', 'lchbsh', 'lpmsh', 'ms', 'lpmsh', 'lmzc',
        'qmf8p', 'lqmd', 'lmblh', 'spqmf8p_k', 'lmblh', 'lm076b', 'bh4p',
        'lm3a', 'ld500', 'lmz', 'vm31p', 'lqzd', 'qmd9p', 'lqmd', 'lmblh',
        'spqmd9p_k', 'lmblh', 'lm55a', 'bh4p', 'lbb4', 'bh4p', 'lm2', 'qmf10p',
        'lqmd', 'lmblh', 'spqmf10p_k', 'lmblh', 'lchbh', 'chp', 'lchbh',
        'lpmlh', 'ms', 'lpmlh', 'lm1b', 'bh4p', 'lm076c', 'vm32p', 'ld150',
        'qmd11p', 'lqmd', 'lmbsh', 'spqmd11p_k', 'lmbsh', 'lm046a', 'bh4p',
        'lbb4', 'bh4p', 'lm146a', 'lmbt1h', 'mt', 'lmbt1h', 'ld100', 'lmz',
        'hm21p', 'lqzd', 'qmf12p', 'lqmd', 'lmblh', 'spqmf12p_k', 'lmblh',
        'lmi1a', 'lcvbh', 'cvp', 'lcvbh', 'lpmlh', 'ms', 'lpmlh', 'lmz',
        'vm33p', 'ld1362', 'qmd13p', 'lqmd', 'lmblh', 'spqmd13p_k', 'lmblh',
        'lmi3a', 'sept2', 'ms', 'lsep1', 'sept1', 'ms', 'injp', 'lsqi6', 'qi6p',
        'liqi6', 'qi7p'])
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
    'fhx041p':                     [77.049493235882],
    'fhx071p':                     [111.092162514335],
    'fhx072p':                     [113.893890958708],
    'fha121p':                     [234.113918571113],
    'fha122p':                     [249.578651399646],
    'fha123p':                     [261.565061813404],
    'fha124p':                     [273.551472227162],
    'fha125p':                     [285.537716471762],
    'fha126p':                     [292.196833562442],
    'fvt201p':                     [319.675480490862],
    'fhc161p':                     [353.434041612098],
    'fhc162p':                     [365.409232693881],
    'fvm271p':                     [409.757084281749],
    'fhm201p':                     [432.876573892591],
    'fhm202p':                     [441.241870803345],
    'fhm203p':                     [453.511632677790],
    'pinjax0':                     [480.752506560096]}

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
