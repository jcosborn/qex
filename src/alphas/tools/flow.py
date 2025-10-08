"""
Cubic spline:
- https://blog.scottlogic.com/2020/05/18/cubic-spline-in-python-and-alteryx.html
- https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation 
"""

import argparse
import pathlib
import bisect
import pickle
from typing import Callable as callable

# constants

CWD: pathlib.Path = pathlib.Path.cwd()
EPS: float = 0.01 # flow time epsilon that data will be saved in
LAST: dict[str, str] = {
    '20': '12.500',
    '24': '18.000',
    '32': '32.000',
    '40': '50.000',
    '48': '72.000'
}
COUPS: dict[str, str] = {
    '700': '7p00',
    '725': '7p25',
    '750': '7p50',
    '800': '8p00',
    '850': '8p00',
    '900': '9p00',
    '100': '10p0',
    '120': '12p0',
    '140': '14p0',
    '160': '16p0',
    '180': '18p0',
    '200': '20p0'
}
MASSES: dict[str, str] = {
    '000':   '0p000',
    '0001':  '0p001',
    '00025': '0p0025',
    '0005':  '0p005'
}

# procedures implementing cubic spline

def dx(x: list[float]) -> list[float]:
    return [x[i + 1] - x[i] for i in range(len(x) - 1)]

def dvec(n: int, h: list[float], y: list[float]) -> list[float]:
    # d_0 = 0
    # d_i/6 = ([y_{i+1} - y_i]/h_{i+1} - [y_i - y_{i-1}]/h_i) / (h_i + h_{i+1})
    # d_{n-1} = 0
    d: callable[int] = lambda i: (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]
    return [0.] + [6. * d(ind) / (h[ind] + h[ind - 1]) for ind in range(1, n - 1)] + [0.]

def solve(a: list[float], b: list[float], c: list[float], d: list[float]) -> list[float]:
    # from: https://blog.scottlogic.com/2020/05/18/cubic-spline-in-python-and-alteryx.html
    cp: list[float] = c + [0]
    dp: list[float] = [0] * len(b)
    x: list[float] = [0] * len(b)
    
    # calculate cp, dp
    cp[0] = c[0] / b[0]
    dp[0] = c[0] / b[0]
    for i in range(1, len(b)):
        cp[i] = cp[i] / (b[i] - cp[i - 1] * a[i - 1])
        dp[i] = (d[i] - dp[i - 1] * a[i - 1]) / (b[i] - cp[i - 1] * a[i - 1])

    # calculate solution & return
    x[-1] = dp[-1]
    for i in range(len(b) - 2, -1, -1): x[i] = dp[i] - cp[i] * x[i + 1]
    return x

def spline(x: list[float], y: list[float]) -> callable[[float], float]:
    n: int = len(x)
    if n < 3: raise ValueError('x is too short')
    if n != len(y): raise ValueError('x and y must be same length')
    
    # setup & solve tridiagonal system:
    # a_i = h_i / (h_i + h_{i+1}) for i in [0, n - 3]
    # b_i = 2 for i in [0, n-1]
    # c_i = h_i / (h_{i-1} + h_i) for i in [1, n - 2]
    h: list[float] = dx(x)
    a: list[float] = [h[i] / (h[i] + h[i + 1]) for i in range(n - 2)] + [0.]
    b: list[float] = [2] * n
    c: list[float] = [0.] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(n - 2)]
    d: list[float] = dvec(n, h, y)
    m: list[float] = solve(a, b, c, d)

    # construct spline coefficients
    c0: callable[int] = lambda i: (m[i+1] - m[i]) * h[i] * h[i] / 6.
    c1: callable[int] = lambda i: m[i] * h[i] * h[i] / 2.
    c2: callable[int] = lambda i: y[i+1] - y[i] - (m[i+1] + 2. * m[i]) * h[i] * h[i] / 6.
    c3: callable[int] = lambda i: y[i]
    
    # define spline as callable & return
    def spline_fcn(xv):
        idx = min(bisect.bisect(x, xv) - 1, n - 2)
        z = (xv - x[idx]) / h[idx]
        return z * (c2(idx) + z * (c1(idx) + c0(idx) * z)) + c3(idx)
    return spline_fcn

# data processing

def get_cut(ens_path: pathlib.Path) -> int:
    cut_file = ens_path / 'cut'
    if cut_file.is_file():
        with open(cut_file, 'r') as in_file:
            return int(in_file.readlines()[0])
    return 0
    
def cfg(fn: pathlib.Path) -> int:
    return int(str(fn).split('_')[-1].replace('.log', ''))

def files(pth: pathlib.Path, flow: str) -> list[str]:
    fns = []
    for fn in pth.glob(flow + '*.log'): fns.append(fn)
    fns.sort(key = cfg)
    return fns
        
def process_data(vol: str, ens: str, flw: str):
    ens_path: pathlib.Path = CWD / vol / ens
    data_path: pathlib.Path = ens_path / 'flow'
    flow_path: pathlib.Path = CWD / 'flow'
    flow_file: pathlib.Path = flow_path / ens
    cut: int = get_cut(ens_path)
    data: dict[str, any] = {
        'Ep': [],
        'Es': [],
        'Ec': [],
        'Q':  []
    }
    itr: int = 0
    with open(str(flow_file) + '-' + flw + 'flow.dat', 'w+') as out_file:
        for fn in files(data_path, flw):
            # check if configuration is below thermalization cut
            cfg: int = int(str(fn).split('_')[-1].replace('.log', ''))
            if (cfg <= cut): continue
            
            # get data
            last_time: str = ''
            ts: list[float] = []
            ps: list[float] = []
            rs: list[float] = []
            cs: list[float] = []
            qs: list[float] = []
            ipss: list[float] = []
            rpss: list[float] = []
            ipts: list[float] = []
            rpts: list[float] = []
            with open(fn, 'r') as in_file:
                lines: list[str] = in_file.readlines()
                last_time = lines[-1].split()[1]
                for line in lines:
                    # format:
                    # time plaq rect clov topo re(polys) im(polys) re(polyt) im(polyt)
                    # where polys is spatial Polyakov loop and polyt is temporal " "
                    # see: src/alphas/alphasflow.nim in milc-qcd fork of QEX
                    _, t, p, r, c, q, rps, ips, rpt, ipt = line.split()
                    ts.append(float(t))
                    ps.append(1. - float(p))
                    rs.append(1. - float(r))
                    cs.append(1. - float(c))
                    qs.append(float(q))
                    rpss.append(float(rps))
                    ipss.append(float(ips))
                    rpts.append(float(rpt))
                    ipts.append(float(ipt))

            # check to make sure unfinished flows don't contaminate outputs
            if last_time != LAST[vol.split('.')[0]]:
                print('WARNING!', fn, 'did not finish!')
                continue

            # prepare data dictionary for new member
            for key in data.keys():
                if key != 'flow_times': data[key].append([])

            # interpolate over data w/ spline
            pspline = spline(ts, ps)
            rspline = spline(ts, rs)
            cspline = spline(ts, cs)
            qspline = spline(ts, qs)
            rpsspline = spline(ts, rpss)
            ipsspline =	spline(ts, ipss)
            rptspline =	spline(ts, rpts)
            iptspline = spline(ts, ipts)

            # save interpolated data in prescribed step size
            ots: list[float] = []
            dt = max(ts) - min(ts)
            nt = int(round(dt / EPS))
            for ti in range(nt + 1):
                t = EPS * ti
                p: float = pspline(t)
                r: float = rspline(t)
                c: float = cspline(t)
                q: float = qspline(t)
                rps: float = rpsspline(t)
                ips: float = ipsspline(t)
                rpt: float = rptspline(t)
                ipt: float = iptspline(t)
                
                # write data to .dat file
                out_data = [
                    str(cfg),
                    f"{t:.2f}",
                    f"{p:.16f}",
                    f"{r:.16f}",
                    f"{c:.16f}",
                    f"{q:.16f}",
                    f"{rps:.16f}",
                    f"{ips:.16f}",
                    f"{rpt:.16f}",
                    f"{ipt:.16f}"
                ]
                out_file.write(' '.join(out_data) + '\n')

                # save data to data dictionary
                if 'flow_times' not in data.keys(): ots.append(t)
                data['Ep'][itr].append(p)
                data['Es'][itr].append(r)
                data['Ec'][itr].append(c)
                data['Q'][itr].append(q)
            if 'flow_times' not in data.keys(): data['flow_times'] = ots
            itr += 1

    # save data dictionary to disk
    b: str = ens.split('b')[-1].split('m')[0]
    m: str = ens.split('m')[-1].split('_')[0]
    fn: str = '_'.join([COUPS[b], vol.replace('.', 'l'), MASSES[m], flw]) + '.bin'
    with open(flow_path / fn, 'wb+') as out_file: pickle.dump(data, out_file)
        
# main execution

if __name__ == '__main__':
    # parse arguments
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description = 'Tool for setting up gradient flow workflows',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-ns',
        '--spatial-extent',
        help = 'Spatial extent of ensemble volume, <ns>',
        required = True
    )
    parser.add_argument(
        '-nt',
        '--temporal-extent',
        help = 'Temporal extent of ensemble volume, <nt>',
        required = True
    )
    parser.add_argument(
        '-e',
        '--ensemble',
        help = 'Ensemble name; format: "f4l<ns>t<nt>b<beta_b>m<mass>_HISQ_pppa" or "all"',
        required = True
    )
    parser.add_argument(
        '-f',
        '--flow',
        help = 'Flow name; e.g., wilson',
        required = True
    )
    args: argparse.Namespace = parser.parse_args()

    # get result of argument parse
    ns: str = args.spatial_extent
    nt: str = args.temporal_extent
    ensemble: str = args.ensemble
    flow: str = args.flow
    
    # get data
    volume = '.'.join([ns, ns, ns, nt])
    process_data(volume, ensemble, flow)
