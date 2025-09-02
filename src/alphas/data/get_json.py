import json

"""
This is just an example showing you how to easily grab the json data
for one of the ensembles as an example
"""

if __name__ == '__main__':
    # grab json data
    fn = 'f4l20t40b200m000_HISQ_pppa-info.json'
    with open(fn, 'r') as in_file: data = json.load(in_file) 

    # example: get dH
    dH = data['dH']

    # example: query available keys
    print(data.keys())

    """
    info:
    dH
    spatial plaquette
    temporal plaquette
    Re[spatial Polyakov loop]
    Re[temporal Polyakov loop]
    Im[spatial Polyakov loop]
    Im[temporal Polyakov loop]
    chiral condensate -- only available for strong coupling ensembles (beta_b = 10/g_0^2 >= 7.5)
    average CG iterations (fermion)
    average CG iterations (hasenbusch)
    cut -- an array of 1's and 0's: 1 means that the config is below the thermalization cut and 0 means that it is above

    ... the rest is just information about the ensemble that it grabbed from the json file in the ensemble directory...
    """
