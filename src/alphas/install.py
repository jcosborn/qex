#!/usr/bin/python3

import argparse
import typing
import sys
import os

MDEVOLVE = 'https://github.com/jxy/MDevolve'

NIMV = 'nim-2.2.4'
NIMP = '-linux_x64'
NIM = 'https://nim-lang.org/download/' + NIMV + NIMP + '.tar.xz'

def dest(path, dir) -> str: return '/'.join([path, dir])

def isdir(path, dir) -> bool: return os.path.isdir(dest(path, dir))

def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog = 'install',
        description = 'Installs MILC collaboration alphas software'
    )
    p.add_argument('-p', '--prefix', required = True, help = 'desired install path')
    p.add_argument('-q', '--qex', required = True, help = 'path to QEX')
    p.add_argument(
        '--compile_gradient_flow', 
        help = 'compile gradient flow', 
        type = bool,
        default = False
    )
    p.add_argument(
        '--run_hmc_regression_test', 
        help = 'runs alphashmc in regression mode after compilation', 
        type = bool,
        default = False
    )
    return p.parse_args()

mkdir = lambda path, dir: os.mkdir(dest(path, dir))
cd = lambda path, dir: os.chdir(dest(path, dir))
wget = lambda url: os.system('wget ' + url)
clone = lambda url: os.system('git clone ' + url)
tarx = lambda file: os.system('tar -xvf ' + file + '.tar.xz')
targ = lambda file: os.system('tar -xvf ' + file + '.tar.gz')
export = lambda name, path: os.system('export ' + name + '=' + path)

install_qop_qio = lambda qex: os.system(dest(qex, 'bootstrap-travis'))
nimble_install = lambda: os.system('nimble install')

system = lambda cmd: os.system(cmd)

symlink = lambda source, dest: os.system('ln -sf ' + source + ' ' + dest)

def install_nim(nim: str):
    wget(NIM)
    tarx(nim + NIMP)

def install_mdevolve(mdevolve: str, nimexec: str, nimbleexec: str):
    if not isdir('./', 'MDevolve'): clone(MDEVOLVE)
    cd('./', 'MDevolve')
    nimble_install()

def configure(qex: str, qmp: str, qio: str):
    system(' '.join([
        dest(qex,'configure'),
        'qmpdir:' + qmp,
        'qiodir:' + qio
    ]))

def install_hmc(build: str, bin: str):
    system('make alphashmc')
    dustbin = dest(dest(build, 'bin'), 'alphashmc')
    symlink(dustbin, bin)

def install_gradient_flow(build: str, bin: str): 
    system('make alphasflow')
    dustbin = dest(dest(build, 'bin'), 'alphasflow')
    symlink(dustbin, bin)

def hmc_regress(): system('./bin/alphashmc')

if __name__ == '__main__':
    args = args()
    path = args.prefix
    qex = args.qex
    compile_gradient_flow = args.compile_gradient_flow
    run_hmc_regression_test = args.run_hmc_regression_test

    if not isdir(path, 'build'): mkdir(path, 'build')
    if not isdir(path, 'deps'): mkdir(path, 'deps')
    if not isdir(path, 'bin'): mkdir(path, 'bin')

    deps = dest(path, 'deps')
    build = dest(path, 'build')
    bin = dest(path, 'bin')

    mdevolve = dest(deps, 'mdevolve')
    nim = dest(deps, NIMV)
    nimpath = dest(nim, 'bin')
    nimexec = dest(nimpath, 'nim')
    nimbleexec = dest(nimpath, 'nimble')

    export('nim', nimexec)
    export('nimble', nimbleexec)

    cd(path, 'deps')
    install_qop_qio(qex)
    qmp = dest(deps, 'qmp')
    qio = dest(deps, 'qio')

    install_nim(nim)
    install_mdevolve(mdevolve, nimexec, nimbleexec)

    cd(path, 'build')
    configure(qex, qmp, qio)

    install_hmc(build, bin)
    if compile_gradient_flow: install_gradient_flow(build, bin)

    if run_hmc_regression_test: hmc_regress()

    # Final printout reminding folks how to use both binary files
    


    

    