algorithm = "hagedorn"
propagator = "semiclassical"
splitting_method = "Y4"

T = 12
dt = 0.01

dimension = 1
ncomponents = 1

eps = 0.5

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
        "K" : 32,
        "dimension" : 1
    }],
# Coefficient values computed by 'ComputeGroundstate.py'
    "coefficients" : [[ ((0,), 9.71624580e-01+0.j),
                        ((1,), 2.17565922e-01+0.j),
                        ((2,), 4.94239114e-02+0.j),
                        ((3,), 6.91773925e-02+0.j),
                        ((4,), 3.24717929e-02+0.j),
                        ((5,), 1.15303880e-02+0.j),
                        ((6,), 1.12231124e-02+0.j),
                        ((7,), 6.93226427e-03+0.j),
                        ((8,), 3.12172192e-03+0.j),
                        ((9,), 2.52232926e-03+0.j),
                        ((10,), 1.78234459e-03+0.j),
                        ((11,), 9.46130107e-04+0.j),
                        ((12,), 6.92830580e-04+0.j),
                        ((13,), 5.20991520e-04+0.j),
                        ((14,), 3.12383872e-04+0.j),
                        ((15,), 2.18769288e-04+0.j),
                        ((16,), 1.68187826e-04+0.j),
                        ((17,), 1.10395188e-04+0.j),
                        ((18,), 7.64700398e-05+0.j),
                        ((19,), 5.88906431e-05+0.j),
                        ((20,), 4.12571591e-05+0.j),
                        ((21,), 2.88198122e-05+0.j),
                        ((22,), 2.20380938e-05+0.j),
                        ((23,), 1.61259169e-05+0.j),
                        ((24,), 1.14353642e-05+0.j),
                        ((25,), 8.61368186e-06+0.j),
                        ((26,), 6.43633879e-06+0.j),
                        ((27,), 4.57387421e-06+0.j),
                        ((28,), 3.23790802e-06+0.j),
                        ((29,), 2.33991946e-06+0.j),
                        ((30,), 1.55696339e-06+0.j) ]],
    "quadrature" : {
        "type" : "HomogeneousQuadrature",
	'qr': {
            'type': 'TensorProductQR',
            'dimension': 1,
            'qr_rules': [{'dimension': 1, 'order': 36, 'type': 'GaussHermiteQR'}]
        }
    }
}

# Which wavepackets are initial values
initvals = [ wp0 ]

leading_component = 0

# How often do we write data to disk
write_nth = 5

matrix_exponential = "pade"
