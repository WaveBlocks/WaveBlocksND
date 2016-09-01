algorithm = "hagedorn"
propagator = "semiclassical"
splitting_method = "Y61"

T = 15
dt = 0.01

dimension = 1
ncomponents = 1

eps = 0.1104536

potential = 'morse_zero'
D = 0.0572
a = 0.983
x0 = 5.03855

# The parameter set of the initial wavepacket
# Sigma = PQ^{-1} = 1.3836
Q = [[1.0]]
P = [[1.0j]]
q = [[4.53]]
p = [[0.0]]
S = [[0.0]]

wp0 = {
    'type': 'HagedornWavepacket',
    'dimension': 1,
    'ncomponents': 1,
    'eps': eps,
    'Pi': [q, p, Q, P, S],
    'basis_shapes': [{
        'type': 'HyperbolicCutShape',
        'K': 32,
        'dimension': 1
    }],
    'coefficients': [[((0,), 1.0)]],
    'innerproduct': {
        'type': 'HomogeneousInnerProduct',
        'delegate': {
            'type': 'DirectHomogeneousQuadrature',
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

matrix_exponential = 'pade'
