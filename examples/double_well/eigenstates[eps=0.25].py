# Compute with --nopq --noPQ

dimension = 1
ncomponents = 1

potential = "double_well"
b = 5.0

eigenstate_of_level = 0

eigenstates_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]

eps = 0.25

hawp_template = {
    "type" : "HagedornWavepacket",
    "dimension" : dimension,
    "ncomponents": 1,
    "eps" : eps,
    "basis_shapes" : [{
            "type" : "HyperCubicShape",
            "limits" : [128],
            "dimension" : 1
            }]
    }

innerproduct = {
    "type" : "HomogeneousInnerProduct",
    "delegate" : {
        "type" : "DirectHomogeneousQuadrature",
        'qr': {
            'type': 'TensorProductQR',
            'dimension': 1,
            'qr_rules': [{'dimension': 1, 'order': 132, 'type': 'GaussHermiteQR'}]
            }
        }
    }
