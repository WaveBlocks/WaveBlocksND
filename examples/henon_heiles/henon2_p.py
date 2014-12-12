algorithm = "hagedorn"
propagator = "semiclassical"
splitting_method = "Y4"

T = 6
dt = 0.01

dimension = 2
ncomponents = 1

eps = 0.05

potential = "henon_heiles"
a = 1
b = 3

# The parameter set of the initial wavepacket
Q = [[1.0, 0.0],
     [0.0, 1.0]]

P = [[1.0j, 0.0 ],
     [0.0,  1.0j]]

q = [[ 0.06],
     [ 0.0]]

p = [[-0.01],
     [ 0.01]]

S = [[0.0]]

# What it takes to specify a wavepacket!
wp0 = {
    "type" : "HagedornWavepacket",
    "dimension" : 2,
    "ncomponents": 1,
    "eps" : eps,
    "Pi" : [q,p,Q,P,S],
    "basis_shapes" : [{
        "type" : "HyperCubicShape",
        "limits" : [10, 10],
        "dimension" : 2
    }],
    "coefficients" : [[ ((0,0), 1.0) ]],
    "innerproduct" : {
        "type" : "HomogeneousInnerProduct",
        "delegate" : {
            "type" : "DirectHomogeneousQuadrature",
            'qr': {
                'type': 'TensorProductQR',
                'dimension': 2,
                'qr_rules': [{'dimension': 1, 'order': 15, 'type': 'GaussHermiteQR'},
                             {'dimension': 1, 'order': 15, 'type': 'GaussHermiteQR'}],
            }
        }
    }
}

# Which wavepackets are initial values
initvals = [ wp0 ]

leading_component = 0

# How often do we write data to disk
write_nth = 5

matrix_exponential = "arnoldi"
arnoldi_steps = 175
