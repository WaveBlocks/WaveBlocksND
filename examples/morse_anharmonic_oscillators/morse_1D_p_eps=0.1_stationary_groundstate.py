algorithm = "hagedorn"
propagator = "semiclassical"
splitting_method = "Y4"

T = 12
dt = 0.01

dimension = 1
ncomponents = 1

eps = 0.1

potential ="morse_zero"
l = 1.0
x0 = 0.0

# The parameter set of the initial wavepacket
# Parameter values computed by 'ComputeGroundstate.py'
Q = [[0.84089642+0.j]]
P = [[0.+1.18920712j]]
q = [[0.0]]
p = [[0.0]]
S = [[0.0]]

# What it takes to specify a wavepacket!
wp0 = {
    "type" : "HagedornWavepacket",
    "dimension" : 1,
    "ncomponents": 1,
    "eps" : eps,
    "Pi" : [q,p,Q,P,S],
    "basis_shapes" : [{
        "type" : "HyperbolicCutShape",
        "K" : 10,
        "dimension" : 1
    }],
# Coefficient values computed by 'ComputeGroundstate.py'
    "coefficients" : [[ ((0,), 9.98929399e-01+0.j),
                        ((1,), 4.45551745e-02+0.j),
                        ((2,), 2.02921328e-03+0.j),
                        ((3,), 1.22086756e-02+0.j),
                        ((4,), 1.26458378e-03+0.j),
                        ((5,), 9.35094064e-05+0.j),
                        ((6,), 3.34407369e-04+0.j),
                        ((7,), 5.19805901e-05+0.j),
                        ((8,), 5.08918795e-06+0.j) ]],
    "quadrature" : {
        "type" : "HomogeneousQuadrature",
	'qr': {
            'type': 'TensorProductQR',
            'dimension': 1,
            'qr_rules': [{'dimension': 1, 'order': 14, 'type': 'GaussHermiteQR'}]
        }
    }
}

# Which wavepackets are initial values
initvals = [ wp0 ]

leading_component = 0

# How often do we write data to disk
write_nth = 5

matrix_exponential = "pade"
