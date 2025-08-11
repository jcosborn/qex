"""
Extracts HMC information from all ensembles and stores in nice little JSON files
"""

import os
import sys
import subprocess
import typing
import re
import json
import statistics

LOCATION = 'JLab'
FLAVOR = 'f4'
SMEAR_BC = 'HISQ_pppa'
CPATH = os.getcwd()

READCG = False
READPBP = False
HASENBUSCH = True
START_NEW_TRAJECTORY = True

COUPLING_CONVERT = {
    '200': '20.0',
    '180': '18.0',
    '160': '16.0',
    '140': '14.0',
    '120': '12.0',
    '100': '10.0',
    '900': '9.00',
    '850': '8.50',
    '800': '8.00',
    '750': '7.50',
    '725': '7.25',
    '700': '7.00'
}
MASS_CONVERT = {
    '000'  : '0.000',
    '0005' : '0.005',
    '00025': '0.0025',
    '0001' : '0.001'
}
DATA = {
    '24.24.24.48': {
        '200': ['000'],
        '180': ['000'],
        '160': ['000'],
        '140': ['000'],
        '120': ['000'],
        '100': ['000'],
        '900': ['000'],
        '850': ['000'],
        '800': ['000'],
        '750': ['0005', '00025', '0001'],
        '725': ['0005', '00025', '0001'],
        '700': ['0005', '00025', '0001']
    },
    '32.32.32.64': {
        '200': ['000'],
	'180': ['000'],
	'160': ['000'],
	'140': ['000'],
	'120': ['000'],
	'100': ['000'],
	'900': ['000'],
	'850': ['000'],
	'800': ['000'],
	'750': ['0005', '00025', '0001'],
        '725': ['0005',	'00025', '0001'],
        '700': ['0005', '00025', '0001']
    },
    '40.40.40.80': {
        '200': ['000'],
	'160': ['000'],
	'180': ['000'],
	'160': ['000'],
	'140': ['000'],
	'120': ['000'],
	'100': ['000'],
	'900': ['000'],
	'850': ['000'],
	'800': ['000'],
	'750': ['0005', '00025', '0001'],
        '725': ['0005',	'00025', '0001'],
        '700': ['0005', '00025', '0001']
    },
    '48.48.48.96': {
        '200': ['000'],
	'180': ['000'],
	'160': ['000'],
	'140': ['000'],
	'120': ['000'],
	'100': ['000'],
	'900': ['000'],
	'850': ['000'],
	'800': ['000'],
	'750': ['0005', '00025', '0001'],
        '725': ['0005', '00025', '0001'],
        '700': ['0005', '00025', '0001']
    }
}

LOCAL_DATA = {}
HCGKEY = 'average CG iterations (hasenbusch)'
FCGKEY = 'average CG iterations (fermion)'

def vol(volume: str) -> str: return ''.join(['l', volume.split('.')[0], 't', volume.split('.')[-1]])

def configuration(f: str) -> int: return int(f.split('_')[-1].replace('.log', ''))

def files(path: str) -> (int, list[str]):
    # get thermalization cut (if it exists)
    cut = 0
    cf = '/'.join([path, 'cut'])
    if os.path.isfile(cf):
        with open(cf, 'r') as in_file: cut = int(in_file.readlines()[0])

    # add file to list of files (if log file and above cut)
    result = []
    for	f in os.listdir(path):
        if not f.endswith('.log'): continue
        try:
            configuration(f)
            result.append(f)
        except ValueError: continue
    result.sort(key = configuration)
    return (cut, result)

def finished(content: list[str]) -> bool:
    return ' s] Total time (Init - Finalize): ' in ''.join(content)
    
def simple_measurements(data: dict[str, any], line: list[str]) -> bool:
    if not line: return False
    tag = line[0].replace(':', '')
    match tag:
        case 'MEASplaq':
            data['spatial plaquette'].append(float(line[2]))
            data['temporal plaquette'].append(float(line[4]))
            return True
        case 'MEASploop':
            data['Re[spatial Polyakov loop]'].append(float(line[2]))
            data['Im[spatial Polyakov loop]'].append(float(line[3]))
            data['Re[temporal Polyakov loop]'].append(float(line[5]))
            data['Im[temporal Polyakov loop]'].append(float(line[6]))
            return True
        case 'ACC' | 'REJ':
            dH = float(line[1].replace(',', ''))
            data['dH'].append(dH)
            data['acceptance'].append(1 if tag == 'ACC' else 0)
            return True
        case _: pass
    return False

def complicated_measurements(data: dict[str, any], line: list[str]) -> None:
    global HCGKEY, FCGKEY, LOCAL_DATA
    global READCG, READPBP, HASENBUSCH
    
    if not line: return None
    tag = line[0].replace(':', '')
    match tag:
        case 'kinetic':
            READCG = not READCG
            if READCG and (HCGKEY in LOCAL_DATA):
                data[HCGKEY].append(statistics.mean(LOCAL_DATA[HCGKEY]))
                data[FCGKEY].append(statistics.mean(LOCAL_DATA[FCGKEY]))
            if READCG: (LOCAL_DATA[HCGKEY], LOCAL_DATA[FCGKEY]) = ([], [])
        case 'stagSolve':
            if not READCG: return None
            if HASENBUSCH:
                LOCAL_DATA[HCGKEY].append(float(line[1]))
                HASENBUSCH = False
            else:
                LOCAL_DATA[FCGKEY].append(float(line[1]))
                HASENBUSCH = True
        case 'MEASpbp': pass
    

def catalogue(volume: str, coupling: str, mass: str) -> None:
    global START_NEW_TRAJECTORY
    
    # path information
    ensemble = ''.join([FLAVOR, vol(volume), 'b', coupling, 'm', mass, '_', SMEAR_BC])
    path = '/'.join([CPATH, volume, ensemble, ''])

    # data to be collected
    data = {
        'dH':                         [],
        'spatial plaquette':          [],
        'temporal plaquette':         [],
        'Re[spatial Polyakov loop]':  [],
        'Im[spatial Polyakov loop]':  [],
        'Re[temporal Polyakov loop]': [],
        'Im[temporal Polyakov loop]': [],
        'chiral condensate':          [],
        FCGKEY:                       [],
        HCGKEY:                       [],
        'acceptance':                 [],
        'cut':                        []
    }

    # data collection
    if not os.path.isdir(path): return
    (cut, logs) = files(path)
    for log in logs:
        with open(path + log, 'r') as in_file:
            lines = in_file.readlines()
            if not finished(lines): continue
            for line in lines:
                spln = line.split()
                simple = simple_measurements(data, spln)
                if not simple: complicated_measurements(data, spln)
                if 'kinetic:' in line:
                    if not START_NEW_TRAJECTORY: data['cut'].append(0 if configuration(log) > cut else 1)
                    START_NEW_TRAJECTORY = not START_NEW_TRAJECTORY
    data['running'] = ''.join([FLAVOR, vol(volume), 'b', coupling]) in subprocess.run(
        ['squeue', '--format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R"', '--me'],
        capture_output = True,
        text = True,
        check = True
    ).stdout

    # read in information about ensemble
    with open(path + ensemble + '.json', 'r') as in_file:
        info = json.loads(re.sub(r'(?<=:)\s*0(\d+)', r'"\g<0>"', in_file.read()))
        for key, value in info.items(): data[key] = value
        
    # save data in json format to disk
    dpath = '/'.join([CPATH, 'hmc', ''])
    if not os.path.isdir(dpath): os.mkdir(CPATH + 'hmc')
    with open(dpath + ensemble + '-info.json', 'w+') as out_file:
        json.dump(data, out_file, indent = 4)

    # tell user that you've done your job
    print('saved hmc info:', dpath + ensemble + '-info.json')
 
if __name__ == '__main__':
    for volume, vDATA in DATA.items():
        for coupling, cvDATA in vDATA.items():
            [*map(lambda mass: catalogue(volume, coupling, mass), cvDATA)]
    
