"""
Data processing module for gauge flows.

Author: Curtis Peterson

~ Notes ~

Using cubic spline interpolation. 

Code for interpolation based off of:
- https://blog.scottlogic.com/2020/05/18/cubic-spline-in-python-and-alteryx.html
- https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation
"""

""" Imports """
# For specifying types
from typing import Tuple, List

# For figuring out where in spline interpolation is
import bisect

# For operating system-specific tasks
import os as os

# For getting user arguments
import sys as sys

""" Cubic spline interpolation code """

# Compute changes
def dx(x: List[float]) -> List[float]:
    # Return changes
    return [x[i + 1] - x[i] for i in range(len(x) - 1)]

# Create tridiagonal matrix
def tridiag_mat(n: List[float], h: List[float]) -> Tuple[List[float], List[float], List[float]]:
    """
    a_i = h_i / (h_i + h_{i+1}) for i in [0, n - 3]
    b_i = 2 for i in [0, n-1]
    c_i = h_i / (h_{i-1} + h_i) for i in [1, n - 2]
    """

    # Create A matrix
    A = [h[i] / (h[i] + h[i + 1]) for i in range(n - 2)] + [0.]

    # Create B matrix
    B = [2] * n

    # Create C matrix
    C = [0.] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(n - 2)]

    # Return A, B and C rows
    return A, B, C

# Create column vector of d_n coefficients
def d_vec(n: int, h: List[float], y: List[float]):
    """
    d_0 = 0
    d_i/6 = ([y_{i+1} - y_i]/h_{i+1} - [y_i - y_{i-1}]/h_i) / (h_i + h_{i+1})
    d_{n-1} = 0
    """
    # Function to compute d_i * (h_i + h_{i-1})/ 6
    d_i = lambda i: (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]

    # Return d vector
    return [0.] + [6. * d_i(ind) / (h[ind] + h[ind - 1]) for ind in range(1, n - 1)] + [0.]

# Solve tridiagonal system
def solve(A: List[float], B: List[float], C: List[float], D: List[float]):
    """
    Standard solve for tridiagonal system of equations. Taken directly from
    https://blog.scottlogic.com/2020/05/18/cubic-spline-in-python-and-alteryx.html
    """
    c_p = C + [0]
    d_p = [0] * len(B)
    X = [0] * len(B)

    c_p[0] = C[0] / B[0]
    d_p[0] = D[0] / B[0]
    for i in range(1, len(B)):
        c_p[i] = c_p[i] / (B[i] - c_p[i - 1] * A[i - 1])
        d_p[i] = (D[i] - d_p[i - 1] * A[i - 1]) / (B[i] - c_p[i - 1] * A[i - 1])

    X[-1] = d_p[-1]
    for i in range(len(B) - 2, -1, -1):
        X[i] = d_p[i] - c_p[i] * X[i + 1]

    return X

# Create spline interpolation function
def spline_func(x: List[float], y: List[float]):
    """
    Does spline interpolation and returns spline function
    """
    
    """ Do checks """
    # Define length of x
    n = len(x)

    # Check length
    if n < 3:
        # Raise error
        raise ValueError('x is too short')
    if n != len(y):
        # Raise error
        raise ValueError('x and y are not the same length')
    
    """ Construct tridiagonal s.o.e. and solve it """
    # Get changes
    h = dx(x)

    # Get rows of tridiagonal system
    A, B, C = tridiag_mat(n, h)

    # Get d vector
    D = d_vec(n, h, y)

    # Solve system of equations
    M = solve(A, B, C, D)

    """ Construct coefficients """
    # Define function for C_0
    C_0 = lambda i: (M[i+1] - M[i]) * h[i] * h[i] / 6.

    # Define function for C_1
    C_1 = lambda i: M[i] * h[i] * h[i] / 2.

    # Define function for C_2
    C_2 = lambda i: y[i+1] - y[i] - (M[i+1] + 2. * M[i]) * h[i] * h[i] / 6.

    # Define function for C_3
    C_3 = lambda i: y[i]

    """ Define spline function and return """
    # Define spline function
    def spline(x_val):
        # Get index of spline
        idx = min(bisect.bisect(x, x_val) - 1, n - 2)

        # Define Z
        Z = (x_val - x[idx]) / h[idx]

        # Return interpolation
        return Z * (C_2(idx) + Z * (C_1(idx) + C_0(idx) * Z)) + C_3(idx) 

    # Return spline function
    return spline

""" Convenience functions """

# For updating progress
def update_progress(progress):
    """Progress bar
    This little helper function creates a progress bar
    
    Taken from: https://stackoverflow.com/questions/3160699/python-progress-bar       
    """

    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block),
                                               round(progress*100), status)
    sys.stdout.write(text)
    sys.stdout.flush()

    # Return nothing 
    return None

""" For finding and processing data """

# For finding data
def find_files():
    # Internal function to get configuration
    def conf_keys(fn):
        return int(fn.split('_')[-1].replace('.log', ''))

    # Define out file name
    out_fn = flow_trans[options['flow']] + '_' + options['l'] + options['t']
    out_fn += '_' + options['b'][0] + '.' + options['b'][1:]

    # Create empty list of files
    files = []

    # Go through files
    for f in os.listdir(data_path):
        # Check if correct
        if (f.startswith(options['flow'] + '_')) and (out_fn not in f):
            # Try to grab data
            try:
                # Try to grab key
                test = conf_keys(data_path + f)

                # Add file
                files.append(data_path + f)
            except ValueError:
                # Print info
                print('Hi, I wasn\'t able to process', data_path + f)

                # Continue
                pass

    # Check if files found
    if files:
        # Sort files
        files.sort(key = conf_keys)

        # Return files
        return files
    else:
        # Otherwise, raise FileNotFoundError
        raise FileNotFoundError('No files with ' + options['flow'] + ' flow')

    # Return nothing
    return None

# Read file
def read_data(fn):
    # Open file                                                                             
    with open(fn, 'r') as in_file:
        # Define data arrays                                                                
        flts,t2EC,t2EP,Q,ReP,ImP = [],[],[],[],[],[]
        
        # Cycle through lines                                                               
        for line in in_file.readlines():
            # Inconsistency option
            inconsist = False

            # Flow verification
            if (('flow_act: ' in line) or ('c1: 'in line)) and (options['verif'] == 'true'):
                # Check if flow action is displayed
                if 'flow_act: ' in line:
                    # Check if C0p0 and action chose to be Wilson
                    if (options['flow'] != 'C0p0') and ('Wilson' in line):
                        # Print warning about inconsistent flows
                        print('Warning: inconsistent flows. Line:', line)

                        # Set inconsistency to true
                        inconsist = True
                    elif (options['flow'] == 'C0p0') and ('Wilson' not in line):
                        # Print warning about inconsistent flows
                        print('Warning: inconsistent flows. Line:', line)

                        # Set inconsistency to true
                        inconsist = True
                elif ('c1: ' in line) and (options['flow'] != 'C0p0'): 
                    # Split line
                    sl = line.split()

                    # Define flow
                    opt_flow = str(round(float(options['flow'][1:].replace('p', '.')), 2))
                    opt_flow = opt_flow.replace('.', 'p')

                    # Define file flow
                    fl_flow = str(round(float(sl[-1]), 2)).replace('.', 'p')

                    # Check that rounded flow name matches flow in file
                    if opt_flow != fl_flow:
                        # Print warning about inconsistent flows
                        print('Warning: inconsistent flows. Line:', line)

                        print(opt_flow, fl_flow)

                        # Set inconsistency to true
                        inconsist = True

            # Check if flow line                                                            
            if ('FLOW' in line) and ('SUMMARY' not in line) and (not inconsist):
                # Get split line                                                            
                sl = line.split()

                # Append flow time                                                          
                flts.append(float(sl[1]))

                # Append t^2 E(t) from clover                                               
                t2EC.append(float(sl[3]))

                # Append topological charge
                Q.append(float(sl[7]))

                # Check if zero                                                             
                if float(sl[1]) == 0.:
                    # Append zero                                                           
                    t2EP.append(float(sl[2]))
                else:
                    # Append t^2 E(t) from plaquette                                        
                    t2EP.append(float(sl[6]) / float(sl[1])**2.)

                # Real/imaginary parts of spatial/temporal Polyakov loop
                ResP = float(sl[-2])
                ImsP = float(sl[-1])
                RetP = float(sl[-4])
                ImtP = float(sl[-3])

                # Space/time-averaged real/imaginary Polyakov loop
                ReP.append(3.*0.25*ResP+0.25*RetP)
                ImP.append(3.*0.25*ImsP+0.25*ImtP)
    
    # Return data
    return flts,t2EC,t2EP,Q,ReP,ImP

# Process data
def process_data():
    """ Set things up """
    # Define file name
    #out_fn = flow_trans[options['flow']] + '_' + options['l'] + options['t'] 
    #out_fn += '_' + options['cpling'][0] + '.' + options['cpling'][1:]
    out_fn = 'f' + options['nf'] + 'l' + options['l'] + 't' + options['t']
    out_fn += 'stggb' + options['b'] + 'm' + options['m']
    out_fn += '-' + options['flow'] + 'flow_cut.dat'
    
    # Set cut to zero
    cut = 0

    """ Getting cut file """
    # Try to open cut file                                                       
    try:
        # Open cut file                                                    
        with open(file_path() + 'cut', 'r') as cut_file:
            # Save cut                                                           
            cut = int(cut_file.readlines()[0])

        # Tell user what cut will be set to
        print('Cut file found and cut set to', cut)
    except FileNotFoundError:
        # Tell user that cut will be kept to zero                                
        print('No cut file. Keeping at config. zero.')

    # Open out file
    with open(data_path + out_fn, 'w+') as out_file:
        # Cycle through files
        for fn_ind, fn in enumerate(files):
            """ Get data """
            # Get configuration
            config = fn.split('_')[-1].replace('.log', '')

            # Check that configuration greater than cut
            if int(config) >= cut:
                # Try to open file
                try:
                    # Read data
                    flts,t2EC,t2EP,Q,ReP,ImP = read_data(fn)
                except FileNotFoundError:
                    # Tell user that data does not exit
                    print(fn, 'does not exist.')

                    # Continue
                    pass

                """ Do interpolation """
                # Start try/catch statement
                try:
                    # Do interpolation of t2EC and t2EP
                    t2ECsp, t2EPsp = spline_func(flts, t2EC), spline_func(flts, t2EP)

                    # Do interpolation of topological charge
                    Qsp = spline_func(flts, Q)

                    # Do interpolation of Polyakov loop components
                    RePsp,ImPsp = spline_func(flts,ReP),spline_func(flts,ImP)
                    
                    """ Write data and finish up """
                    # Define initial flt
                    flt = 0.

                    # Cycle through dt flow times
                    while True:
                        # Create text to write to output file
                        text = [
                            config, str(round(flt, 2)), str(t2EPsp(flt)), 
                            str(0.0),
                            str(t2ECsp(flt)), str(Qsp(flt)),
                            str(RePsp(flt)), str(ImPsp(flt))
                        ]

                        # Write to out file
                        out_file.write(' '.join(text) + '\n')
                
                        # Check flow time
                        if flt + options['dt'] >= max(flts): break
                        else: flt += options['dt']
                        
                    # Update progress
                    update_progress((fn_ind + 1) / len(files))
                except ValueError:
                    # Tell user that something went wrong
                    print('Warning! Something wrong with', fn + '.', 
                          'Flow could be broken.')

                    # Continue
                    pass

    # Create symbolic link
    os.system('ln -sf ' + data_path + out_fn + ' ./flow/' + out_fn)

""" Main code """
# Main code
if __name__ == "__main__":
    """ Specify options and grab user input """
    # Specify options
    options = {
        'dt' : 0.05, 'flow' : 'C0p0', 'm' : None,
        'l' : None, 't' : None, 'b' : None,
        'nf' : '4', 'bc' : 'pppa', 'ext' : '',
        'verif': 'true'
    }

    for arg in sys.argv:
        option = arg.split('=')[0].replace('-', '')
        if option in options.keys(): options[option] = arg.split('=')[-1]
    
    """ Define a few convenient things """
    # Get current directory
    path = os.getcwd()

    # Check if time specified
    options['t'] = options['l'] if options['t'] is None else options['t']

    # Define volume
    vol = '.'.join([options['l'], options['l'], options['l'], options['t']])

    # Define ensemble
    ens = 'f' + options['nf'] + 'l' + options['l'] + 't' + options['t']
    ens += 'b' + options['b'] + 'm' + options['m'].replace('p','')
    ens += '_HISQ_' + options['bc'] + options['ext']
    
    # Define file path
    file_path = lambda : path + '/' + vol + '/'  + ens + '/'

    # Define path for saving data
    data_path = file_path() + 'flow/'

    # Define flow translation dictionary
    flow_trans = {'C0p0' : 'W'}

    # Check if flow not in dicationary
    if options['flow'] not in flow_trans.keys():
        # Put in dictionary
        flow_trans[options['flow']] = options['flow']

    """ Find and process data """
    # Find files
    files = find_files()

    # Process data
    process_data()

