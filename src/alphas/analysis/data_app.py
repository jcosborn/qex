import os
import json
import datetime
import statistics as stat

import numpy as np
import streamlit as st
import pyerrors as pe
import gvar as gv

import plot

### information ###

DATASETS = {
    '20': {
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
    '24': {
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
    '32': {
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
    '40': {
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
    '48': {
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

INTEGRATOR = {'2MN': '2nd-order Omelyan'}
VOLUMES = {
    "20": r"$20^3 \times 40$",
    "24": r"$24^3 \times 48$",
    "32": r"$32^3 \times 64$",
    "40": r"$40^3 \times 80$",
    "48": r"$48^3 \times 96$"
}
BETAS = {
    "700": "7.00",
    "725": "7.25",
    "750": "7.50",
    "800": "8.00",
    "850": "8.50",
    "900": "9.00",
    "100": "10.0",
    "120": "12.0",
    "140": "14.0",
    "160": "16.0",
    "180": "18.0",
    "200": "20.0"
}
MASSES = {
    "000"  : "0.0000",
    "0005" : "0.0050",
    "00025": "0.0025",
    "0001" : "0.0010"
}

ENSVOL = {
    "20": "l20t40",
    "24": "l24t48",
    "32": "l32t64",
    "40": "l40t40",
    "48": "l48t96"
}

### helper procedures ###

def weighted(data, expdh):
    try:
        Z = sum(expdh)
        m = sum(data*expdh)/Z
        e = np.sqrt(sum(expdh*(data - m)**2.)/Z)
        return gv.gvar(m, e)
    except ValueError:
        tdata = list(data)
        texpdh = list(expdh)
        if len(tdata) > len(texpdh):
            while len(tdata) != len(texpdh): del tdata[-1]
        elif len(data) < len(expdh):
            while len(tdata) != len(texpdh): del texpdh[-1]
        tdata = np.array(tdata)
        texpdh = np.array(texpdh)
        Z = sum(texpdh)
        m = sum(tdata*texpdh)/Z
        e = np.sqrt(sum(texpdh*(tdata - m)**2.)/Z)
        return gv.gvar(m, e)

def proper(data):
    obs = pe.Obs([data], ['dummy'])
    obs.gamma_method()
    return (
        gv.gvar(obs.e_tauint['dummy'], obs.e_dtauint['dummy']), 
        gv.gvar(obs.value, obs.dvalue)
    )

### page layout ###

# page configuration
st.set_page_config(page_title = 'ensembles', page_icon = ':atom_symbol:')

# page title
st.markdown(r"# Fermilab Lattice and MILC Collaboration $\boldsymbol{\alpha_{\mathrm{s}}(M_{Z})}$ ensembles")

# navigation information
st.markdown("""
""")

# Construct side bar
st.sidebar.markdown("""
## Ensemble Parameters
""")

### side bar ###

# volume selector
volume = st.sidebar.selectbox(
    r"$N_{\mathrm{s}} \equiv L/a$", 
    options = [*VOLUMES.keys()]
)

# bare gauge coupling selector
coupling = st.sidebar.selectbox(
    r"$\beta_{b} \equiv 6/g_{0}^2$", 
    options = [*DATASETS[volume].keys()][::-1],
    format_func = lambda beta: BETAS[beta]
)

# bare gauge coupling selector
mass = st.sidebar.selectbox(
    r"$am_{\mathrm{f}}$", 
    options = DATASETS[volume][coupling][::-1],
    format_func = lambda m: MASSES[m]
)

### display ensemble information ###

ensemble = ''.join(['f4', ENSVOL[volume], 'b', coupling, 'm', mass, '_HISQ_pppa'])
fn = ensemble + '-info.json'

if os.path.exists('../data/' + fn):
    ### side bar display ###
    with open('../data/' + fn, 'r') as in_file: data = json.load(in_file)
    st.sidebar.markdown(
    f"""
    ## Summary \n
    - configurations: {len(data['dH']) - sum(data['cut'])}
    - running: {data['running']}
    - last update: {datetime.datetime.fromtimestamp(os.path.getctime('../data/' + fn))}
    
    ## Molecular Dynamics
    - trajectory length: {data['hmc']['trajectory-length']}
    - outer integrator 
      - fields: fermion and Hasenbusch
      - algorithm: {INTEGRATOR[data['fermion']['integrator']]}
      - steps: {data['fermion']['steps']}
    - inner integrator
      - fields: gauge
      - algorithm: {INTEGRATOR[data['gauge']['integrator']]}
      - steps: {data['gauge']['steps']}*

    ## Action
    - Hasenbusch mass: {data['action']['hasenbusch-mass']}
    - plaquette coefficient: 1
    - rectangle coefficent: -1/20
    - tadpole: 1
    - fat7 (level 1)
      - one link: 1/8
      - three link: 1/16
      - five link: 1/64
      - seven link: 1/384
      - lepage: 0.0
    - fat7 (level 2)
      - one link: 1
      - three link: 1/16
      - five link: 1/64
      - seven link: 1/384
      - lepage: -1/8
    - naik: -1/24
    - unitary projection: Cayley-Hamilton

    --------

    *per outer integrator gauge field update 
    """)

    ### main page plots ###

    cfgs = [*range(len(data['dH']))]

    # fcn(dH)
    dH = np.array(data['dH'])
    dH2 = dH*dH
    expdH = np.exp(-dH)

    # cg iterations
    hcg = data['average CG iterations (hasenbusch)']
    fcg = data['average CG iterations (fermion)']

    def plot_hmc_observable(h, meas, lbl):
        cfg_cut = [cfg for idx, cfg in enumerate(cfgs) if not data['cut'][idx]] 
        meas_cut = [m for idx, m in enumerate(meas) if not data['cut'][idx]]
        h.scatter(
            [cfg for idx, cfg in enumerate(cfgs) if data['cut'][idx]], 
            [m for idx, m in enumerate(meas) if data['cut'][idx]], 
            **{'facecolor': 'none', 'edgecolor': 'grey', 'alpha': 0.1}
        )
        h.axis.axvline(min(cfg_cut), color = 'k')
        try:
            h.scatter(
                cfg_cut, meas_cut, 
                **{'facecolor': 'none', 'edgecolor': 'grey', 'alpha': 0.5}
            )
        except ValueError:
            if len(cfg_cut) > len(meas_cut):
                while len(cfg_cut) != len(meas_cut): del cfg_cut[-1]
            elif len(cfg_cut) < len(meas_cut):
                while len(cfg_cut) != len(meas_cut): 
                    del dH[-1]
                    del dH2[-1]
                    del expdH[-1]
            h.scatter(
                cfg_cut, meas_cut, 
                **{'facecolor': 'none', 'edgecolor': 'grey', 'alpha': 0.5}
            )
        expdh_cut = [m for idx, m in enumerate(expdH) if not data['cut'][idx]]
        exp = weighted(np.array(meas_cut), np.array(expdh_cut))
        x = [min(cfg_cut), max(cfg_cut) - 1]
        y = [exp, exp]
        h.fill_between(x, y, color = 'magenta', alpha = 0.1)
        h.axis.axhline(exp.mean, color = 'magenta', alpha = 0.5)
        h.set_title('$\\overline{' + lbl + '} =' + str(exp) + '$') 
    
    # information
    st.markdown("""
    ## Hamiltonian Monte Carlo metrics
    """)

    # cumulative acceptance rate
    dhp = plot.Plot()
    allacc = [1 for idx, _ in enumerate(data['acceptance']) if not data['cut'][idx]]
    accs = [m for idx, m in enumerate(data['acceptance']) if not data['cut'][idx]]
    cfg_cut = [cfg for idx, cfg in enumerate(cfgs) if not data['cut'][idx]] 
    nrm = np.cumsum(allacc)
    accr = np.cumsum(accs)
    cum_acc_rate = accr/nrm
    dhp.line(cfg_cut, cum_acc_rate , color = 'grey')
    dhp.axis.axvline(min(cfg_cut), color = 'k')
    dhp.set_title(
        '$\\mathrm{acceptance \\ rate \\ = \\ ' + str(cum_acc_rate[-1]) + '}$'
    )
    dhp.decorate(
        xlim = [0, len(cfgs) - 1],
        ylim = [0., 1.],
        xlabel = '$\mathrm{configuration}$',
        ylabel = '$\mathrm{cumulative \\ acc \\ rate}$'
    )
    dhp.axis.axhline(cum_acc_rate[-1], color = 'magenta', alpha = 0.5)
    st.pyplot(fig = dhp.handle)

    # dH
    dhp = plot.Plot()
    plot_hmc_observable(dhp, dH, '\mathrm{d}\mathcal{H}')
    dhp.decorate(
        xlim = [0, len(cfgs) - 1],
        ylim = [-1., 1.],
        xlabel = '$\mathrm{configuration}$',
        ylabel = '$\mathrm{d}\mathcal{H}$'
    )
    st.pyplot(fig = dhp.handle)

    # dH^2
    dhp = plot.Plot()
    plot_hmc_observable(dhp, dH2, '\mathrm{d}\mathcal{H}^2')
    dhp.decorate(
        xlim = [0, len(cfgs) - 1],
        ylim = [-0.1, 2.0],
        xlabel = '$\mathrm{configuration}$',
        ylabel = '$\mathrm{d}\mathcal{H}^2$'
    )
    st.pyplot(fig = dhp.handle)

    # exp(-dH)
    dhp = plot.Plot()
    plot_hmc_observable(dhp, expdH, '\exp(-\mathrm{d}\mathcal{H})')
    dhp.decorate(
        xlim = [0, len(cfgs) - 1],
        ylim = [0.0, 2.0],
        xlabel = '$\mathrm{configuration}$',
        ylabel = '$\exp(-\mathrm{d}\mathcal{H})$'
    )
    st.pyplot(fig = dhp.handle)

    # average CG per trajectory (fermion)
    hcg_lbl = 'avg \\ CG \\ itns \\ per \\ traj'
    dhp = plot.Plot()
    plot_hmc_observable(dhp, fcg, '\\mathrm{' + hcg_lbl + '}')
    dhp.decorate(
        xlabel = '$\mathrm{configuration}$',
        ylabel = '$\\mathrm{' + hcg_lbl + ' \\ (fermion)}$'
    )
    st.pyplot(fig = dhp.handle)

    # information
    st.markdown("""
    ## Observables
    Errors and autocorrelation times calculated from 
    implementation of the $\Gamma$-method in [pyerrors](https://github.com/fjosw/pyerrors). 
    - pyerrors: https://doi.org/10.1016/j.cpc.2023.108750
    - $\Gamma$-method: https://doi.org/10.1016/S0010-4655(03)00467-3
    """)

    def plot_observable(h, meas, lbl):
        cfg_cut = [cfg for idx, cfg in enumerate(cfgs) if not data['cut'][idx]] 
        meas_cut = [m for idx, m in enumerate(meas) if not data['cut'][idx]]
        h.scatter(
            [cfg for idx, cfg in enumerate(cfgs) if data['cut'][idx]], 
            [m for idx, m in enumerate(meas) if data['cut'][idx]], 
            **{'facecolor': 'none', 'edgecolor': 'grey', 'alpha': 0.1}
        )
        h.axis.axvline(min(cfg_cut), color = 'k')
        h.scatter(
            cfg_cut, meas_cut, 
            **{'facecolor': 'none', 'edgecolor': 'grey', 'alpha': 0.5}
        )
        (tau, exp) = proper(meas)
        x = [min(cfg_cut), max(cfg_cut) - 1]
        y = [exp, exp]
        h.fill_between(x, y, color = 'magenta', alpha = 0.1)
        h.axis.axhline(exp.mean, color = 'magenta', alpha = 0.5)
        ttl = '$\\overline{' + lbl + '} =' + str(exp) + '$, '
        ttl += '$\\tau = ' + str(tau) + '$'
        h.set_title(ttl)

    splaq = np.array(data['spatial plaquette'])
    tplaq = np.array(data['temporal plaquette'])
    plaq = 0.5*(splaq + tplaq)

    srpoly = np.array(data['Re[spatial Polyakov loop]'])
    trpoly = np.array(data['Re[temporal Polyakov loop]'])

    sipoly = np.array(data['Im[spatial Polyakov loop]'])
    tipoly = np.array(data['Im[temporal Polyakov loop]'])

    rpoly = 0.5*(srpoly + trpoly)
    ipoly = 0.5*(sipoly + tipoly)

    poly = np.sqrt(rpoly*rpoly + ipoly*ipoly)

    # plaquette
    dhp = plot.Plot()
    plot_observable(dhp, plaq, '\mathrm{plaquette}')
    dhp.decorate(
        xlim = [0, len(cfgs) - 1],
        ylim = [0., 1.],
        xlabel = '$\mathrm{configuration}$',
        ylabel = '$\mathrm{plaquette}$'
    )
    st.pyplot(fig = dhp.handle)

    # Polyakov loop
    dhp = plot.Plot()
    dhp.scatter(
        rpoly, ipoly, **{'facecolor': 'none', 'edgecolor': 'grey', 'alpha': 0.5}
    )
    (tau, exp) = proper(poly)
    ttl = '$\\overline{|\\mathrm{Polyakov}|} =' + str(exp) + '$, '
    ttl += '$\\tau = ' + str(tau) + '$'
    dhp.set_title(ttl)
    dhp.decorate(
        xlim = [-0.25, 0.25],
        ylim = [-0.25, 0.25],
        xlabel = '$\Re\mathrm{Polyakov}$',
        ylabel = '$\Im\mathrm{Polyakov}$'
    )
    st.pyplot(fig = dhp.handle)

    # scatter plot of Polyakov loop in argand plane


else: st.image("ensemble_not_found.png")   
