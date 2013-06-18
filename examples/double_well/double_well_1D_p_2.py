algorithm = "hagedorn"
propagator = "semiclassical"
splitting_method = "Y4"

T = 2 * 4.4
dt = 0.0005

dimension = 1
ncomponents = 1

eps = 0.1

potential = {}
potential["variables"] = ["x"]
potential["potential"] = "x**4 - x**2"

# The parameter set of the initial wavepacket
Q = [[1.0]]
P = [[1.0j]]
q = [[0.70710678118654746]]
p = [[-0.70710678118654746 + 2.0*0.1]]
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
        "K" : 32,
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
                'qr_rules': [{'dimension': 1, 'order': 36, 'type': 'GaussHermiteQR'}]
            }
        }
    }
}

# Which wavepackets are initial values
initvals = [ wp0 ]

leading_component = 0

# How often do we write data to disk
write_nth = 10

matrix_exponential = "pade"
