"""The WaveBlocks Project

This file contains some ready made potentials in several variables
and with several separate energy levels. This is a pure data file
without any code. To load the potentials, use the methods of
:py:class:`BlockFactory`.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

####################################################################
# Potentials in one dimension and with one energy level (D=1, N=1) #
####################################################################

# Free particle
free_particle = {}
free_particle["variables"] = ["x"]
free_particle["potential"] = "c"
free_particle["defaults"] = {"c":"0"}
free_particle["number_levels"] = 1

# Simple harmonic potential 1D
quadratic = {}
quadratic["variables"] = ["x"]
quadratic["potential"] = "1/2 * sigma * x**2"
quadratic["defaults"] = {"sigma":"1/2"}
quadratic["number_levels"] = 1

# Perturbed harmonic potential
pert_quadratic = {}
pert_quadratic["variables"] = ["x"]
pert_quadratic["potential"] = "1/2 * sigma * x**2 + 1/2 * delta**2 * x**2"
pert_quadratic["defaults"] = {"sigma":0.05, "delta":0.2}
pert_quadratic["number_levels"] = 1

# A simple fourth order anharmonic potential
quartic = {}
quartic["variables"] = ["x"]
quartic["potential"] = "1/4 * sigma * x**4"
quartic["defaults"] = {"sigma":0.05}
quartic["number_levels"] = 1

quartic_reg = {}
quartic_reg["variables"] = ["x"]
quartic_reg["potential"] = "1/4 * sigma * x**4 + tau * sigma * x**2"
quartic_reg["defaults"] = {"sigma":0.05, "tau":1}
quartic_reg["number_levels"] = 1

# A potential consisting of a cosine wave
cos_osc = {}
cos_osc["variables"] = ["x"]
cos_osc["potential"] = "a * (1 - cos(b*x))"
cos_osc["defaults"] = {"a":0.07, "b":1.0}
cos_osc["number_levels"] = 1

# A potential consisting of a hyperbolic cosine
cosh_osc = {}
cosh_osc["variables"] = ["x"]
cosh_osc["potential"] = "a * cosh(b * x)"
cosh_osc["defaults"] = {"a":"1", "b":"1"}
cosh_osc["number_levels"] = 1

# A potential consisting of a hyperbolic cosine squared
# Its Taylor expansion is:
#   a  +  a b^2 x^4  +  1/3 a b^4 x^8  + ...
# which for the default parameters reduces to:
#   1  +  x^4  +  1/3 x^8  + ...
cosh_osc_sq = {}
cosh_osc_sq["variables"] = ["x"]
cosh_osc_sq["potential"] = "a * cosh(b * x**2)**2"
cosh_osc_sq["defaults"] = {"a":"1", "b":"1"}
cosh_osc_sq["number_levels"] = 1

# The Morse potential
morse = {}
morse["variables"] = ["x"]
morse["potential"] = "D * (1 - exp(-a*(x-x0)))**2"
morse["defaults"] = {"D":3.0, "a":0.5, "x0":0.0}
morse["number_levels"] = 1

morse_zero = {}
morse_zero["variables"] = ["x"]
morse_zero["potential"] = "D * (exp(-2*a*(x-x0)) - 2*exp(-a*(x-x0)))"
morse_zero["defaults"] = {"D":3.0, "a":0.5, "x0":0.0}
morse_zero["number_levels"] = 1

morse_zero_2 = {}
morse_zero_2["variables"] = ["x"]
morse_zero_2["potential"] = "l**2 * (exp(-2*(x-x0)) - 2*exp(-(x-x0)))"
morse_zero_2["defaults"] = {"l":1.0, "x0":0.0}
morse_zero_2["number_levels"] = 1

# A double well potential
double_well = {}
double_well["variables"] = ["x"]
double_well["potential"] = "sigma * (x**2 - 1)**2"
double_well["defaults"] = {"sigma":1.0}
double_well["number_levels"] = 1

double_well2 = {}
double_well2["variables"] = ["x"]
double_well2["potential"] = "a*x**4 - b*x**2"
double_well2["defaults"] = {"a":1.0, "b":1.0}
double_well2["number_levels"] = 1

# The Eckart potential
eckart = {}
eckart["variables"] = ["x"]
eckart["potential"] = "sigma * cosh(x/a)**(-2)"
eckart["defaults"] = {"sigma":100*3.8088*10**(-4), "a":1.0/(2.0*0.52918)}
eckart["number_levels"] = 1

# A smooth unit step like wall
wall = {}
wall["variables"] = ["x"]
wall["potential"] = "atan(sigma*x) + pi/2"
wall["defaults"] = {"sigma":10.0}
wall["number_levels"] = 1

# A narrow 'V'-like potential
v_shape = {}
v_shape["variables"] = ["x"]
v_shape["potential"] = "1/2 * sqrt(tanh(x)**2+4*delta**2)"
v_shape["defaults"] = {"delta":0.2}
v_shape["number_levels"] = 1

# Kratzer
kratzer = {}
kratzer["variables"] = ["x"]
kratzer["potential"] = "1/2 * (x**2 + b*(b-1) * x**(-2))"
kratzer["defaults"] = {"b":2.0}
kratzer["number_levels"] = 1

#####################################################################
# Potentials in two dimensions and with one energy level (D=2, N=1) #
#####################################################################

# Simple harmonic potential 2D
quadratic_2d = {}
quadratic_2d["variables"] = ["x", "y"]
quadratic_2d["potential"] = "1/2 * (sigmax * x**2 + sigmay * y**2)"
quadratic_2d["defaults"] = {"sigmax":"1/2", "sigmay":"1/2"}
quadratic_2d["number_levels"] = 1

# A simple fourth order anharmonic potential
quartic_2d = {}
quartic_2d["variables"] = ["x", "y"]
quartic_2d["potential"] = "sigmax * x**4 + sigmay * y**4"
quartic_2d["defaults"] = {"sigmax":"1", "sigmay":"1"}
quartic_2d["number_levels"] = 1

# A simple fourth order anharmonic potential
quartic_2d_rotsym = {}
quartic_2d_rotsym["variables"] = ["x", "y"]
quartic_2d_rotsym["potential"] = "sigmax**2 * x**4 + 2*sigmax*sigmay * x**2*y**2 + sigmay**2 * y**4"
quartic_2d_rotsym["defaults"] = {"sigmax":"1", "sigmay":"1"}
quartic_2d_rotsym["number_levels"] = 1

# A simple fourth order anharmonic potential
quartic_reg_2d = {}
quartic_reg_2d["variables"] = ["x", "y"]
quartic_reg_2d["potential"] = "sigmax * x**4 + sigmay * y**4 + taux * x**2 + tauy * y**2"
quartic_reg_2d["defaults"] = {"sigmax":"1", "sigmay":"1", "taux":"1", "tauy":"1"}
quartic_reg_2d["number_levels"] = 1

# A simple fourth order anharmonic potential
quartic_rotsym_reg_2d = {}
quartic_rotsym_reg_2d["variables"] = ["x", "y"]
quartic_rotsym_reg_2d["potential"] = "sigmax * x**4 + sigmay * y**4 + 2*sqrt(sigmax*sigmay)*x**2*y**2 + taux * x**2 + tauy * y**2"
quartic_rotsym_reg_2d["defaults"] = {"sigmax":"1", "sigmay":"1", "taux":"1", "tauy":"1"}
quartic_rotsym_reg_2d["number_levels"] = 1

# A simple sixth order anharmonic potential
sextic_2d = {}
sextic_2d["variables"] = ["x", "y"]
sextic_2d["potential"] = "sigmax * x**6 + sigmay * y**6"
sextic_2d["defaults"] = {"sigmax":"1", "sigmay":"1"}
sextic_2d["number_levels"] = 1

# A simple sixth order anharmonic potential
sextic_reg_2d = {}
sextic_reg_2d["variables"] = ["x", "y"]
sextic_reg_2d["potential"] = "sigmax * x**6 + sigmay * y**6 + taux * x**2 + tauy * y**2"
sextic_reg_2d["defaults"] = {"sigmax":"1", "sigmay":"1", "taux":"1", "tauy":"1"}
sextic_reg_2d["number_levels"] = 1

# A simple sixth order anharmonic potential
sextic_rotsym_reg_2d = {}
sextic_rotsym_reg_2d["variables"] = ["x", "y"]
sextic_rotsym_reg_2d["potential"] = "ax**3 * x**6 + ay**3 * y**6 + 3 * ax**2 * ay * x**4 * y**2 + 3 * ax * ay**2 * x**2 * y**4 + taux * x**2 + tauy * y**2"
sextic_rotsym_reg_2d["defaults"] = {"ax":"1", "ay":"1", "taux":"1", "tauy":"1"}
sextic_rotsym_reg_2d["number_levels"] = 1

# A potential consisting of a cosine wave part in 2D
cos_osc_2d = {}
cos_osc_2d["variables"] = ["x", "y"]
cos_osc_2d["potential"] = "ax * (1 - cos(bx*x)) + ay * (1 - cos(by*y))"
cos_osc_2d["defaults"] = {"ax":"1", "bx":"1", "ay":"1", "by":"1"}
cos_osc_2d["number_levels"] = 1

# A potential consisting of a cosine wave part in 2D
cos_osc_rotsym_2d = {}
cos_osc_rotsym_2d["variables"] = ["x", "y"]
cos_osc_rotsym_2d["potential"] = "-a * cos(b * sqrt(x**2+y**2)**2)"
cos_osc_rotsym_2d["defaults"] = {"a":"1", "b":"1"}
cos_osc_rotsym_2d["number_levels"] = 1

# A potential consisting of a cosine wave part in 2D
cos_osc_mul_2d = {}
cos_osc_mul_2d["variables"] = ["x", "y"]
cos_osc_mul_2d["potential"] = "-cos(a*x) * cos(b*y)"
cos_osc_mul_2d["defaults"] = {"a":"1", "b":"1"}
cos_osc_mul_2d["number_levels"] = 1

# A potential consisting of a cosine wave part in 2D
cos_osc_add_2d = {}
cos_osc_add_2d["variables"] = ["x", "y"]
cos_osc_add_2d["potential"] = "-(cos(a*x) + cos(b*y))"
cos_osc_add_2d["defaults"] = {"a":"1", "b":"1"}
cos_osc_add_2d["number_levels"] = 1

# A potential consisting of a hyperbolic cosine
cosh_osc_2d = {}
cosh_osc_2d["variables"] = ["x", "y"]
cosh_osc_2d["potential"] = "ax * (1 + cosh(bx*x)) + ay * (1 + cosh(by*y))"
cosh_osc_2d["defaults"] = {"ax":"1", "bx":"1", "ay":"1", "by":"1"}
cosh_osc_2d["number_levels"] = 1

# A potential consisting of a hyperbolic cosine
cosh_osc_rotsym_2d = {}
cosh_osc_rotsym_2d["variables"] = ["x", "y"]
cosh_osc_rotsym_2d["potential"] = "a * cosh(b * sqrt(x**2+y**2))"
cosh_osc_rotsym_2d["defaults"] = {"a":"1", "b":"1"}
cosh_osc_rotsym_2d["number_levels"] = 1

# A potential consisting of a hyperbolic cosine
cosh_osc_mul_2d = {}
cosh_osc_mul_2d["variables"] = ["x", "y"]
cosh_osc_mul_2d["potential"] = "cosh(a*x) * cosh(b*y)"
cosh_osc_mul_2d["defaults"] = {"a":"1", "b":"1"}
cosh_osc_mul_2d["number_levels"] = 1

# A potential consisting of a hyperbolic cosine
cosh_osc_add_2d = {}
cosh_osc_add_2d["variables"] = ["x", "y"]
cosh_osc_add_2d["potential"] = "cosh(a*x) + cosh(b*y)"
cosh_osc_add_2d["defaults"] = {"a":"1", "b":"1"}
cosh_osc_add_2d["number_levels"] = 1

# A quad well potential with 4 wells
quad_well = {}
quad_well["variables"] = ["x", "y"]
quad_well["potential"] = "ax*x**4 + ay*y**4 - bx*x**2 - by*y**2 -cx*x - cy*y"
quad_well["defaults"] = {"ax":1.0, "ay":1.0, "bx":3.0, "by":3.0, "cx":0.0, "cy":0.0}
quad_well["number_levels"] = 1

# A two dimensional double well potential
double_well_2d = {}
double_well_2d["variables"] = ["x", "y"]
double_well_2d["potential"] = "ax*x**4 + ay*y**4 - bx*x**2 - by*y**2 -cx*x - cy*y"
double_well_2d["defaults"] = {"ax":1.0, "ay":1.0, "bx":4.0, "by":0.0, "cx":0.0, "cy":0.0}
double_well_2d["number_levels"] = 1

# A two dimensional harmonic double well potential
double_well_harmonic_2d = {}
double_well_harmonic_2d["variables"] = ["x", "y"]
double_well_harmonic_2d["potential"] = "ax*x**4 + ay*y**4 - bx*x**2 - by*y**2 -cx*x - cy*y"
double_well_harmonic_2d["defaults"] = {"ax":1.0, "ay":0.0, "bx":4.0, "by":-1.0, "cx":0.0, "cy":0.0}
double_well_harmonic_2d["number_levels"] = 1

# A parabolic channel potential
channel_2d = {}
channel_2d["variables"] = ["x", "y"]
channel_2d["potential"] = "sigmax*x + 1/2*sigmay*y**2"
channel_2d["defaults"] = {"sigmax":"0.0", "sigmay":"0.45"}
channel_2d["number_levels"] = 1

# Another harmonic channel potential
harmonic_channel = {}
harmonic_channel["variables"] = ["x", "y"]
harmonic_channel["potential"] = "1/2*w**2*x**2 + sigma*y"
harmonic_channel["defaults"] = {"w":"1.0", "sigma":"-0.1"}
harmonic_channel["number_levels"] = 1

# A potential consisting of circular pit of radius R with steep walls
corral_rotsym_2d = {}
corral_rotsym_2d["variables"] = ["x", "y"]
corral_rotsym_2d["potential"] = "atan(sigma*(sqrt(x**2+y**2) - R)) + pi/2"
corral_rotsym_2d["defaults"] = {"sigma":"10", "R":"8"}
corral_rotsym_2d["number_levels"] = 1

# A potential consisting of circular pit of radius R with steep walls
circle_pit_2d = {}
circle_pit_2d["variables"] = ["x", "y"]
circle_pit_2d["potential"] = "atan(sigma*(sqrt(x**2+y**2) - R)) + pi/2"
circle_pit_2d["defaults"] = {"sigma":"10", "R":"8"}
circle_pit_2d["number_levels"] = 1

# A potential consisting of a ring like corral
corral_ring = {}
corral_ring["variables"] = ["x", "y"]
corral_ring["potential"] = "-sqrt(delta**2 + tanh(sqrt(x**2 + y**2) - R)**2 * tanh(sqrt(x**2 + y**2) + R)**2)/2"
corral_ring["defaults"] = {"delta":"1", "R":"3"}
corral_ring["number_levels"] = 1

# A potential consisting of a ring like valley
ring_valley = {}
ring_valley["variables"] = ["x", "y"]
ring_valley["potential"] = "sqrt(delta**2 + tanh(sqrt(x**2 + y**2) - R)**2 * tanh(sqrt(x**2 + y**2) + R)**2)/2"
ring_valley["defaults"] = {"delta":"1", "R":"3"}
ring_valley["number_levels"] = 1

# A potential consisting of plane with a Gaussian hill
gauss_hill_2d = {}
gauss_hill_2d["variables"] = ["x", "y"]
gauss_hill_2d["potential"] = "exp(-(sigmax*x**2+sigmay*y**2))"
gauss_hill_2d["defaults"] = {"sigmax":"1", "sigmay":"1"}
gauss_hill_2d["number_levels"] = 1

# A threefold morse potential
# This formula is for real inputs only due to atan2
morse_threefold = {}
morse_threefold["variables"] = ["x", "y"]
morse_threefold["potential"] = "(1 - exp((x**2+y**2) * (-sigma-(1-cos(3*atan2(y,x)))**2/16)))**2"
morse_threefold["defaults"] = {"sigma":"0.05"}
morse_threefold["number_levels"] = 1

# A threefold morse potential
# This formula has a singularity at (0,0)
morse_threefold_2 = {}
morse_threefold_2["variables"] = ["x", "y"]
morse_threefold_2["potential"] = """(1 - exp( (-16*sigma*(x**2+y**2)**3-(x*(x**2-3*y**2)-(x**2+y**2)**(3/2))**2)
                                              /(16*(x**2 + y**2)**2) ) )**2"""
morse_threefold_2["defaults"] = {"sigma":"0.05"}
morse_threefold_2["number_levels"] = 1

# A 2D tunneling example called the 'Eckart bottleneck potential'
eckart_bn = {}
eckart_bn["variables"] = ["x", "y"]
eckart_bn["potential"] = "v0*cosh(a*x)**(-2) + k/2*(1-sigma*exp(-l*x**2)) * y**2"
eckart_bn["defaults"] = {"v0":"0.425", "a":"1.3624", "k":"0.06784", "sigma":"0.5", "l":"0.25"}
eckart_bn["number_levels"] = 1

# Henon-Heiles 2D
henon_heiles = {}
henon_heiles["variables"] = ["x", "y"]
henon_heiles["potential"] = "1/2*a*(x**2 + y**2) + b*(x**2*y - y**3/3)"
henon_heiles["defaults"] = {"a":"1", "b":"1/2"}
henon_heiles["number_levels"] = 1

# Small well on top of high mountain
sine_maar = {}
sine_maar["variables"] = ["x", "y"]
sine_maar["potential"] = "sin(x**2+y**2) + alpha * exp(-sigma*(x**2+y**2))"
sine_maar["defaults"] = {"alpha":0.8, "sigma":1.0}
sine_maar["number_levels"] = 1

# Pullen-Edmonds
pullen_edmonds = {}
pullen_edmonds["variables"] = ["x", "y"]
pullen_edmonds["potential"] = "1/2*(b1*x**2 + b2*y**2) + a*x**2*y**2"
pullen_edmonds["defaults"] = {"a":"1", "b1":"1", "b2":"1"}
pullen_edmonds["number_levels"] = 1

#######################################################################
# Potentials in three dimensions and with one energy level (D=3, N=1) #
#######################################################################

# Simple harmonic potential 3D
quadratic_3d = {}
quadratic_3d["variables"] = ["x", "y", "z"]
quadratic_3d["potential"] = "1/2 * (sigmax * x**2 + sigmay * y**2 + sigmaz * z**2)"
quadratic_3d["defaults"] = {"sigmax":"1/2", "sigmay":"1/2", "sigmaz":"1/2"}
quadratic_3d["number_levels"] = 1

# A harmonic tube potential
harmonic_tube = {}
harmonic_tube["variables"] = ["x", "y", "z"]
harmonic_tube["potential"] = "1/2*wx**2*x**2 + 1/2*wy**2*y**2 + sigma*z"
harmonic_tube["defaults"] = {"wx":"1.0", "wy":"1.0", "sigma":"-0.1"}
harmonic_tube["number_levels"] = 1

######################################################################
# Potentials in four dimensions and with one energy level (D=4, N=1) #
######################################################################

# Simple harmonic potential 4D
quadratic_4d = {}
quadratic_4d["variables"] = ["x1", "x2", "x3", "x4"]
quadratic_4d["potential"] = "1/2 * (sigma1 * 1**2 + sigma2 * x2**2 + sigma3 * x3**2 + sigma4 * x4**2)"
quadratic_4d["defaults"] = {"sigma1":"1/2", "sigma2":"1/2", "sigma3":"1/2", "sigma4":"1/2"}
quadratic_4d["number_levels"] = 1

# A harmonic hyper-tube potential
harmonic_hypertube = {}
harmonic_hypertube["variables"] = ["x1", "x2", "x3", "x4"]
harmonic_hypertube["potential"] = "1/2*w1**2*x1**2 + 1/2*w2**2*x2**2 + 1/2*w3**2*x3**2 + sigma*x4"
harmonic_hypertube["defaults"] = {"w1":"1.0", "w2":"1.0", "w3":"1.0", "sigma":"-0.1"}
harmonic_hypertube["number_levels"] = 1



#####################################################################
# Potentials in one dimension and with two energy levels (D=1, N=2) #
#####################################################################

# Double harmonic potential for two components
two_quadratic = {}
two_quadratic["variables"] = ["x"]
two_quadratic["potential"] = [["1/2*sigma*x**2", "0"             ],
                              ["0",              "1/2*sigma*x**2"]]
two_quadratic["defaults"] = {"sigma":0.05}
two_quadratic["number_levels"] = 2

# Double quartic anharmonic potential for two components
two_quartic = {}
two_quartic["variables"] = ["x"]
two_quartic["potential"] = [["1/4*sigma*x**4", "0"             ],
                            ["0",              "1/8*sigma*x**4"]]
two_quartic["defaults"] = {"sigma":"1"}
two_quartic["number_levels"] = 2

# A potential with a single avoided crossing
delta_gap = {}
delta_gap["variables"] = ["x"]
delta_gap["potential"] = [["1/2 * tanh(x)", "delta"         ],
                          ["delta",         "-1/2 * tanh(x)"]]
delta_gap["defaults"] = {"delta":"0.2"}
delta_gap["number_levels"] = 2

# Diagonalized single avoided crossing
delta_gap_diag = {}
delta_gap_diag["variables"] = ["x"]
delta_gap_diag["potential"] = [["sqrt(delta**2 + tanh(x)**2/4)", "0"                             ],
                               ["0",                             "-sqrt(delta**2 + tanh(x)**2/4)"]]
delta_gap_diag["number_levels"] = 2

# A potential with two avoided crossings in series
two_crossings = {}
two_crossings["variables"] = ["x"]
two_crossings["potential"] = [["tanh(x-rho)*tanh(x+rho)/2", "delta/2"                   ],
                              ["delta/2",                   "-tanh(x-rho)*tanh(x+rho)/2"]]
two_crossings["defaults"] = {"rho":3.0}
two_crossings["number_levels"] = 2

######################################################################
# Potentials in two dimensions and with two energy levels (D=2, N=2) #
######################################################################

delta_gap_rotsym = {}
delta_gap_rotsym["variables"] = ["x", "y"]
delta_gap_rotsym["potential"] = [["tanh(sqrt(x**2 + y**2))/2",                      "delta"],
                                 ["delta",                     "-tanh(sqrt(x**2 + y**2))/2"]]
delta_gap_rotsym["number_levels"] = 2

conic = {}
conic["variables"] = ["x", "y"]
conic["potential"] = [["x",  "y"],
                      ["y", "-x"]]
conic["number_levels"] = 2

conic_avoided = {}
conic_avoided["variables"] = ["x", "y"]
conic_avoided["potential"] = [["x",                   "sqrt(y**2+delta**2)"],
                              ["sqrt(y**2+delta**2)", "-x"                 ]]
conic_avoided["defaults"] = {"delta":1.0}
conic_avoided["number_levels"] = 2

conic_avoided_c = {}
conic_avoided_c["variables"] = ["x", "y"]
conic_avoided_c["potential"] = [["x",           "y + I*delta"],
                                ["y - I*delta", "-x"         ]]
conic_avoided_c["number_levels"] = 2



########################################################################
# Potentials in three dimensions and with two energy levels (D=3, N=2) #
########################################################################



#######################################################################
# Potentials in one dimension and with three energy levels (D=1, N=3) #
#######################################################################

# Decoupled harmonic potentials for three components
three_quadratic = {}
three_quadratic["variables"] = ["x"]
three_quadratic["potential"] = [["1/2 * sigma * x**2", "0",                  "0"                 ],
                                ["0",                  "1/2 * sigma * x**2", "0"                 ],
                                ["0",                  "0",                  "1/2 * sigma * x**2"]]
three_quadratic["defaults"] = {"sigma":0.05}
three_quadratic["number_levels"] = 3

# A potential with three energy levels and multiple crossings
three_levels = {}
three_levels["variables"] = ["x"]
three_levels["potential"] = [["tanh(x+rho) + tanh(x-rho)", "delta1",       "delta2"         ],
                             ["delta1",                    "-tanh(x+rho)", "0"              ],
                             ["delta2",                    "0",            "1 - tanh(x-rho)"]]
three_levels["defaults"] = {"rho":3.0}
three_levels["number_levels"] = 3



########################################################################
# Potentials in two dimensions and with three energy levels (D=2, N=3) #
########################################################################



##########################################################################
# Potentials in three dimensions and with three energy levels (D=3, N=3) #
##########################################################################



######################################################################
# Potentials in one dimension and with four energy levels (D=1, N=4) #
######################################################################

# Decoupled harmonic potentials for four components
four_quadratic = {}
four_quadratic["variables"] = ["x"]
four_quadratic["potential"] = [["1/2 * sigma * x**2", "0",                  "0",                  "0"                 ],
                               ["0",                  "1/2 * sigma * x**2", "0",                  "0"                 ],
                               ["0",                  "0",                  "1/2 * sigma * x**2", "0"                 ],
                               ["0",                  "0",                  "0",                  "1/2 * sigma * x**2"]]
four_quadratic["defaults"] = {"sigma":0.05}
four_quadratic["number_levels"] = 4

# Harmonic and higher order anharmonic potentials for four components
four_powers = {}
four_powers["variables"] = ["x"]
four_powers["potential"] = [["1/2 * sigma * x**2", "0",                  "0",                  "0"                 ],
                            ["0",                  "1/4 * sigma * x**4", "0",                  "0"                 ],
                            ["0",                  "0",                  "1/6 * sigma * x**6", "0"                 ],
                            ["0",                  "0",                  "0",                  "1/8 * sigma * x**8"]]
four_powers["defaults"] = {"sigma":0.05}
four_powers["number_levels"] = 4



######################################################################
# Potentials in one dimension and with five energy levels (D=1, N=5) #
######################################################################

# Decoupled harmonic potential for five components
five_quadratic = {}
five_quadratic["variables"] = ["x"]
five_quadratic["potential"] = [["1/2 * sigma * x**2", "0",                  "0",                  "0",                  "0"                 ],
                               ["0",                  "1/2 * sigma * x**2", "0",                  "0",                  "0"                 ],
                               ["0",                  "0",                  "1/2 * sigma * x**2", "0",                  "0"                 ],
                               ["0",                  "0",                  "0",                  "1/2 * sigma * x**2", "0"                 ],
                               ["0",                  "0",                  "0",                  "0",                  "1/2 * sigma * x**2"]]
five_quadratic["defaults"] = {"sigma":0.05}
five_quadratic["number_levels"] = 5
