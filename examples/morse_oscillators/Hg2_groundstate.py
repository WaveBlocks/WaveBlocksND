dimension = 1
ncomponents = 1

potential = 'morse_zero'
D = 0.004164356975043448
a = 0.896696467191408
x0 = 5.542566723349327

eps = 0.048360430020609635

groundstate_of_level = 0

hawp_template = {
    'type': 'HagedornWavepacket',
    'dimension': dimension,
    'ncomponents': 1,
    'eps': eps,
    'basis_shapes': [{
        'type': 'HyperbolicCutShape',
        'K': 32,
        'dimension': 1
    }]
}

innerproduct = {
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
