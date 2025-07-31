import sys,os
import warnings
import typing
import math
import statistics as stat
import pickle as pkl

from collections.abc import Iterator

# Warning handling
class ZeroDivisionAutocorrelationWarning(UserWarning):
    def __init__(self): print('zero division in autocorrelation time')

class NameCheckWarning(UserWarning):
    def __init__(self, file: str) -> None: print(file + ' failed name check')

class TemporalExtentWarning(UserWarning):
    def __init__(self) -> None: print('-t=<t> not specified, so setting t = 2*l')

class NoCutFileWarning(UserWarning):
    def __init__(self) -> None: print('no cut file found; setting cut to zero')
    
# Error handling
class UnspecifiedOptionError(Exception):
    def __init__(self, tag: str) -> None: 
        self.tag = tag
        self._message = "must specify -" + self.tag + "=<" + self.tag + ">"
        self._message += " in command line"
    def __str__(self) -> str: return repr(self._message)
    def mend(self, spatial: str) -> str: 
        TemporalExtentWarning
        return str(2*int(spatial))

class NoLogFilesError(Exception):
    def __init__(self, epath: str) -> None: 
        self._message = "No log files found in " + epath
    def __str__(self) -> str: return repr(self._message)

# Progress bar
def update_progress(progress: list[int] | list[float]) -> None:
    """Progress bar                                                                
    Creates progress bar                                                        
    Adapted from: https://stackoverflow.com/questions/3160699/python-progress-bar   
    """
    (length,status) = (10, "")
    if isinstance(progress,int): progress = float(progress)
    if progress < 0: (progress,status) = (0,"Halt...\r\n")
    if progress >= 1: (progress,status) = (1,"Done...\r\n")
    block = int(round(length*progress))
    (hashes,progress) = ("#"*block + "-"*(length-block),round(100*progress))
    text = "\rPercent: [{0}] {1}% {2}".format(hashes,progress,status)
    print(text,end='')

# Autocorrelation time estimation
def batch_data(data: list[float], n: int) -> Iterator[float]: 
    for i in range(0, len(data), n): yield data[i: i + n]

def batch_means(data: list[float], size: int) -> list[float]:
    batched_data = [*batch_data(data,size)]
    return [stat.mean(bdata) for bdata in batched_data]

def effective_sample_size(data: list[float], size: int) -> list[float]:
    return stat.variance(data)/stat.variance(batch_means(data,size))

def autocorrelation_time(data: list[str]) -> float:
    fdata = [*map((lambda x: float(x)),data)]
    batch_size = round(len(data)**(2./3.))
    try: return batch_size/effective_sample_size(fdata,batch_size)
    except ZeroDivisionError: pass
    ZeroDivisionAutocorrelationWarning
    return 1.0

# Reads command line and saves option states
def cmd_line(options: dict[str]) -> tuple[str]:
    option_tags = [*options.keys()]
    for arg in sys.argv:
        cmd_flag = arg.split('=')
        cmd_tag = cmd_flag[0].replace('-', '')
        if cmd_tag in option_tags: options[cmd_tag] = cmd_flag[-1]
    for option_tag, option_state in options.items():
        if option_state is None: 
            try: raise UnspecifiedOptionError(option_tag)
            except UnspecifiedOptionError as err: 
                if err.tag == 't': options['t'] = err.mend(options['l']); pass
    return (state for state in options.values())
    
# Find files
def configuration_number(filename: str) -> int:
    return int(filename.split('_')[-1].replace('.log', ''))

def find_files() -> list[str]:
    files = []
    for file in os.listdir(epath):
        if file.endswith('.log'):
            try: 
                fpath = '/'.join([epath,file])
                configuration_number(fpath)
                files.append(fpath)
            except ValueError: NameCheckWarning(file)
    if not files: raise NoLogFilesError
    files.sort(key = configuration_number)
    return files

# Get cut
def get_cut() -> int:
    cpath = '/'.join([epath,'cut'])
    if os.path.isfile(cpath): 
        with open(cpath,'r') as cut_file:
            cut = int(cut_file.readlines()[0])
            print('cut at:',cut)
            return cut
    NoCutFileWarning
    return 0

# Extract information from files
def extract_header(inline: str) -> None: 
    line = inline.split(':')
    tag = line[0]
    if tag in header.keys(): header[tag] = line[-1].replace('\n','')

def appendMeas(tag: str, member: str) -> None:
    measurements[tag].append(float(member))

def appendMC(tag: str, member: str) -> None:
    try: monte_carlo[tag].append(float(member))
    except stat.StatisticsError:
        print('not enough data in appendMC'); pass
    
def appendCG(state: bool, line: list[str]) -> None:
    (iters,gflops) = (int(line[1]),float(line[3].replace('Gf/s','')))
    if state: fCGIters.append(iters); fGFLOP.append(gflops)
    else: hCGIters.append(iters); hGFLOP.append(gflops)

def appendSmear(state: bool, smear: str, line: list[str]) -> None:
    gflops = float(line[2].replace('Gf/s',''))
    if state:
        match smear:
            case 'link': flGFLOP.append(gflops)
            case 'force': ffGFLOP.append(gflops)
            case _: pass

def appendAccept(acc: bool, line: list[str]) -> None:
    acceptance.append(acc)
    if acc: dH = float(line[1].replace(',',''))
    else: dH = 0.0
    measurements['dH'].append(dH)
    measurements['dH^2'].append(dH*dH)
    measurements['exp(-dH)'].append(math.exp(-dH))

def extract_measurement(line: list[str]) -> None:
    match line[0].replace(':',''):
        case 'MEASplaq':
            appendMeas('spatial plaquette',line[2])
            appendMeas('temporal plaquette',line[4])
            appendMeas('total plaquette',line[6])
        case 'MEASploop':
            appendMeas('Re[spatial Polyakov loop]',line[2])
            appendMeas('Im[spatial Polyakov loop]',line[3])
            appendMeas('Re[temporal Polyakov loop]',line[5])
            appendMeas('Im[temporal Polyakov loop]',line[6])
        #case 'MEASpbp': appendMeas('chiral condensate',line[-1])
        case _: pass

def extract_monte_carlo(line: list[str]) -> None:
    global measureForceState,measureFermionCGState
    global fCGIters,hCGIters,fGFLOP,hGFLOP,flGFLOP,ffGFLOP
    match line[0].replace(':',''):
        case 'kinetic': 
            if measureForceState:
                try:
                    appendMC('FF CG iterations [fermion]',stat.mean(fCGIters))
                    appendMC('FF CG iterations [Hasenbusch]',stat.mean(hCGIters))
                    appendMC('FF CG GFLOP rate [fermion]',stat.mean(fGFLOP))
                    appendMC('FF CG GFLOP rate [Hasenbusch]',stat.mean(hGFLOP))
                    appendMC('FF link smear GFLOP rate',stat.mean(flGFLOP))
                    appendMC('FF force smear GFLOP rate',stat.mean(ffGFLOP))
                except stat.StatisticsError: print('not enough data in appendMC'); pass
                (fCGIters,hCGIters,fGFLOP,hGFLOP) = ([],[],[],[])
                (flGFLOP,ffGFLOP) = ([],[])
            measureForceState = not measureForceState
        case 'stagSolve':
            if measureForceState:
                appendCG(measureFermionCGState,line)
                measureFermionCGState = not measureFermionCGState
        case 'linkSmear': appendSmear(measureForceState,'link',line)
        case 'forceSmear': appendSmear(measureForceState,'force',line)
        case 'ACC': appendAccept(1,line)
        case 'REJ': appendAccept(0,line)
        case _: pass
    try:
        if 'Total' in line[2]:
            timing['time'].append(float(line[-2]))
            timing['configs'] += 10
    except IndexError: pass
        
def extract_content(filename: str) -> None:
    with open(filename,'r') as in_file:
        for line in in_file.readlines(): 
            split_line = line.replace('\n','').split()
            if split_line:
                extract_header(line)
                extract_measurement(split_line)
                extract_monte_carlo(split_line)

# Summarize Monte Carlo information
def expectation_value(tag: str, data: list[float], prob = None) -> str:
    if ((tag != 'dH') and (tag != 'dH^2') and (tag != 'exp(-dH)')):
        (mean, svar, tau) = (stat.mean(data), stat.variance(data), autocorrelation_time(data))
        return '<'+tag+'> = '+str(mean)+'+/-'+str(math.sqrt(tau*svar/len(data)))
    elif ((tag == 'dH') or (tag == 'dH^2') or (tag == 'exp(-dH)')):
        Z = sum(prob)
        mean = sum(d*p for d, p in zip(data, prob)) / Z
        svar = sum((d - mean)*(d - mean)*p for d, p
                   in zip(data, prob)) / Z
        tau = autocorrelation_time(data)
        return '<'+tag+'> = '+str(mean)+'+/-'+str(math.sqrt(tau*svar))
        
def summarize() -> None:
    print('\n' + 50 * '-.')
    print('configurations:',timing['configs'])
    for tag,state in header.items(): print(tag+':',state)
    print(50 * '-.')
    for tag, data in measurements.items(): print(expectation_value(tag, data, measurements['exp(-dH)']))
    print(50 * '-.')
    for tag, data in monte_carlo.items(): print(tag + ' =',stat.mean(data))
    print('acceptance rate =',sum(acceptance)/len(acceptance))
    try: print('average time [s] =',stat.mean(timing['time']))
    except stat.StatisticsError: print('average time [s] = <not enough files>')
    print(50 * '-.')

def write_data() -> None:
    data = {'measurements': measurements, 'monte carlo': monte_carlo, 'information': timing}
    with open(rpath+'/hmc/'+ensemble+'.bin','wb+') as out_file: pkl.dump(data,out_file)
    
# Main execution
if __name__ == "__main__":
    print(50 * '-.')
    header = {
        'trajectory length': None,
        'beta': None,
        'mass': None,
        'hasenbusch mass': None
    }
    measurements = {
        'spatial plaquette': [],
        'temporal plaquette': [],
        'total plaquette': [],
        'Re[spatial Polyakov loop]': [],
        'Im[spatial Polyakov loop]': [],
        'Re[temporal Polyakov loop]': [],
        'Im[temporal Polyakov loop]': [],
        'dH': [],
        'dH^2': [],
        'exp(-dH)': []
    }
    monte_carlo = {
        'FF CG iterations [fermion]': [],
        'FF CG GFLOP rate [fermion]': [],
        'FF CG iterations [Hasenbusch]': [],
        'FF CG GFLOP rate [Hasenbusch]': [],
        'FF link smear GFLOP rate': [],
        'FF force smear GFLOP rate': [],
    }
    timing = {
        'core-hours': [],
        'time': [],
        'configs': 0
    }

    measureForceState = False
    measureFermionCGState = False
    (fCGIters,hCGIters,fGFLOP,hGFLOP) = ([],[],[],[])

    (flGFLOP,ffGFLOP) = ([],[])

    acceptance = []

    (l,t,b,m,postfix,save) = cmd_line(
        {'l': None, 't': None, 'b': None, 'm': None, 'postfix': None, 'save': 'True'}
    )
    if postfix is None: postfix = ''
    else: postfix = '-' + postfix
    
    volume = '.'.join([l,l,l,t])
    ensemble = ''.join(['f4l',l,'t',t,'b',b,'m',m,'_HISQ_pppa',postfix])

    rpath = os.getcwd()
    vpath = '/'.join([rpath,volume])
    epath = '/'.join([vpath,ensemble])

    files = find_files()

    cut = get_cut()
    for idx,file in enumerate(files):
        cfg = configuration_number(file)
        if cfg > cut: 
            extract_content(file)
            update_progress(idx/len(files))
    timing['cut'] = cut
    summarize()
    if save == 'True':
        if not os.path.isdir(rpath+'/hmc/'): os.mkdir(rpath+'/hmc/')
        write_data()
