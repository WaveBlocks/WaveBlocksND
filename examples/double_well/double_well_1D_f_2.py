algorithm = "fourier"
propagator = "fourier"

T = 2 * 4.4
dt = 0.0005

dimension = 1
ncomponents = 1

eps = 0.1

potential = {}
potential["variables"] = ["x"]
potential["potential"] = "x**4 - x**2"

# The grid of our simulation domain
limits = [(-6.283185307179586, 6.283185307179586)]
number_nodes = [4096]

# The parameter set of the initial wavepacket
Q = [[1.0]]
P = [[1.0j]]
q = [[0.70710678118654746]]
p = [[-0.70710678118654746 + 0.2]]
S = [[0.0]]

# What it takes to specify a wavepacket!
wp0 = {
    "type": "HagedornWavepacket",
    "dimension": 1,
    "ncomponents": 1,
    "eps": eps,
    "Pi": [q, p, Q, P, S],
    "basis_shapes": [{
        "type": "HyperbolicCutShape",
        "K": 32,
        "dimension": 1
    }],
    "coefficients": [[((0,), 1.0)]],
    "innerproduct": {
        "type": "HomogeneousInnerProduct",
        "delegate": {
            "type": "DirectHomogeneousQuadrature",
            'qr': {
                'type': 'TensorProductQR',
                'dimension': 1,
                'qr_rules': [{'dimension': 1, 'order': 36, 'type': 'GaussHermiteQR'}]
            }
        }
    }
}

# Which wavepackets are initial values
initvals = [wp0]

leading_component = 0

# How often do we write data to disk
write_nth = 10
