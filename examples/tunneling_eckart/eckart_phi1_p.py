algorithm = "hagedorn"
propagator = "semiclassical"
splitting_method = "Y4"

T = 70
dt = 0.005

dimension = 1
ncomponents = 1

# Note: the eps in the paper is our eps**2
eps = 0.1530417681822

potential = "eckart"
sigma = 100 * 3.8008 * 10**(-4.0)
a =  1.0 / (2.0 * 0.52918)

# The parameter set of the initial wavepacket
Q = [[ 3.5355339059327 ]]
P = [[ 0.2828427124746j]]
q = [[-7.5589045088306 ]]
p = [[ 0.2478854736792 ]]
S = [[0.0]]

# What it takes to specify a wavepacket!
wp0 = {
    "type" : "HagedornWavepacket",
    "dimension" : dimension,
    "ncomponents": ncomponents,
    "eps" : eps,
    "Pi" : [q,p,Q,P,S],
    "basis_shapes" : [{
        "type" : "HyperbolicCutShape",
        "K" : 512,
        "dimension" : 1
    }],
    "coefficients" : [[ ((1,), 1.0) ]],
    "innerproduct" : {
        "type" : "HomogeneousInnerProduct",
        "delegate" : {
            "type" : "DirectHomogeneousQuadrature",
            'qr': {
                'type': 'TensorProductQR',
                'dimension': 1,
                'qr_rules': [{'dimension': dimension, 'order': 516, 'type': 'GaussHermiteQR'}]
            }
        }
    }
}

# Which wavepackets are initial values
initvals = [ wp0 ]

leading_component = 0

# How often do we write data to disk
write_nth = 20

matrix_exponential = "arnoldi"
arnoldi_steps = 15
