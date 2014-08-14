algorithm = "fourier"

T = 2 * 4.4
dt = 0.05

dimension = 1
ncomponents = 1

eps = 0.1

potential = "quartic"
sigma = 4.0

# The grid of our simulation domain
limits = [(-6.283185307179586, 6.283185307179586)]
number_nodes = [8192]

# The parameter set of the initial wavepacket
Q = [[1.0 ]]
P = [[1.0j]]
q = [[0.0 ]]
p = [[1.0 ]]
S = [[0.0 ]]

# What it takes to specify a wavepacket!
wp0 = {
    "type" : "HagedornWavepacket",
    "dimension" : 1,
    "ncomponents": 1,
    "eps" : eps,
    "Pi" : [q,p,Q,P,S],
    "basis_shapes" : [{
            "type" : "HyperbolicCutShape",
            "K" : 7,
            "dimension" : 1
            }],
    "coefficients" : [[ ((0,), 1.0) ]],
    "innerproduct" : {
        "type" : "HomogeneousInnerProduct",
        "delegate" : {
            "type" : "DirectHomogeneousQuadrature",
            'qr': {
                'type': 'TensorProductQR',
                'dimension': 1,
                'qr_rules': [{'dimension': 1, 'order': 14, 'type': 'GaussHermiteQR'}]
                }
            }
        }
    }

# Which wavepackets are initial values
initvals = [ wp0 ]

# How often do we write data to disk
write_nth = 1
