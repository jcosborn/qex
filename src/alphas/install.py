#!/usr/bin/python3

# ... I wrote this before I knew pathlib existed. Oh well...
import argparse
import subprocess
import sys
import os

MDEVOLVE = 'https://github.com/jxy/MDevolve'

NIMV = 'nim-2.2.2'
NIMP = '-linux_x64'
NIM = 'https://nim-lang.org/download/' + NIMV + NIMP + '.tar.xz'
GRID = 'https://github.com/ctpeterson/Grid-HISQ.git'
SPACK = 'https://github.com/spack/spack.git'

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
        '--grid-backend',
        help = 'compile with Grid [https://github.com/paboyle/Grid] backend',
        type = bool,
        default = False
    )
    p.add_argument(
        '--system',
        help = 'target system for Grid compilation',
        type = str,
        choices = ['local', 'lq1', 'perlmutter'],
        default = 'local'
    )
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
    p.add_argument(
        '--build-cpus',
        help = 'number of CPUs to use for building (make -j<n>)',
        type = str,
        default = '4'
    )
    return p.parse_args()

run = lambda cmd: os.system(cmd)
mkdir = lambda path, dir: os.mkdir(dest(path, dir))
cd = lambda path, dir: os.chdir(dest(path, dir))
wget = lambda url: os.system('wget ' + url)
clone = lambda url: os.system('git clone ' + url)
tarx = lambda file: os.system('tar -xvf ' + file + '.tar.xz')
targ = lambda file: os.system('tar -xvf ' + file + '.tar.gz')
mv = lambda source, dest: os.system('mv ' + source + ' ' + dest)

def subsystem(cmd: str) -> str:
    return subprocess.check_output(cmd, shell = True).decode('utf-8').strip()

install_qop_qio = lambda qex: os.system(dest(qex, 'bootstrap-travis'))
nimble_install = lambda nimble, nim: os.system(nimble + ' install --nim=' + nim)

system = lambda cmd: os.system(cmd)

symlink = lambda source, dest: os.system('ln -sf ' + source + ' ' + dest)

def install_nim(nim: str):
    wget(NIM)
    tarx(nim + NIMP)

def install_mdevolve(mdevolve: str, nimexec: str, nimbleexec: str):
    if not isdir(deps, 'MDevolve'): clone(MDEVOLVE)
    cd(deps, 'MDevolve')
    nimble_install(nimbleexec, nimexec)

def install_grid(deps: str, machine: str, build_cpus: str) -> str:
    cd(deps, '')
    if not isdir(deps, 'Grid'): 
        clone(GRID)
        mv(deps + '/Grid-HISQ', deps + '/Grid')
    cd(deps, 'Grid')

    # bootstrap (ensures that Eigen is also installed)
    system('./bootstrap.sh')
    if not isdir(deps, 'Grid/build'): mkdir(deps, 'Grid/build')

    # install spack for getting Grid dependencies
    clone(SPACK)
    spack = dest(deps, dest('Grid', 'spack'))
    setup_env = dest(spack, 'share/spack/setup-env.sh')
    run('chmod u+x ' + setup_env)
    run(setup_env)

    # install FFTW w/ Spack... takes a while
    run(spack + '/bin/spack install -v -j ' + build_cpus + ' fftw')
    fftw = subsystem(
        "echo `" + spack + "/bin/spack find --paths fftw | grep ^fftw | awk '{print $2}'`"
    )
    
    # configure
    cd(deps, 'Grid/build')
    grid = dest(dest(deps, 'Grid'), 'build')
    config = '--prefix=' + grid + ' '
    config += '--with-fftw=' + fftw + ' '
    if machine == 'local':
        config += '--enable-simd=GEN '
        config += '--enable-comms=mpi-auto '
    elif machine == 'lq1':
        #config += '--enable-simd=AVX512 '
        config += '--enable-simd=AVX '
        config += '--enable-comms=mpi-auto '
        #config += '--enable-shm=shmget '
        #config += '--enable-shmpath=/dev/hugepages '
    config += '--disable-fermion-reps '
    config += '--disable-gparity '
    system('../configure ' + config)

    # make & make install (make distributed over "build_cpus" CPUs)
    system('make -j' + build_cpus)
    system('make install -j' + build_cpus)

    # return path
    return grid

def configure(qex: str, qmp: str, qio: str, grid: str, nim: str, grid_backend: bool):
    os.environ['NIM'] = nim
    if grid_backend:
        grid_config = grid + '/grid-config'
        include = dest(grid, 'include')
        flags = '-I' + include + ' ' 
        flags += subsystem(grid_config + ' --cxxflags') + ' '
        ld = subsystem(grid_config + ' --ldflags') + ' '
        libs = subsystem(grid_config + ' --libs') + ' '
        libs += ' -L' + dest(grid, 'lib') + ' -lGrid' + ' '
        system(' '.join([
            dest(qex, 'configure'),
            'qmpdir:' + qmp,
            'qiodir:' + qio,
            'griddir:' + grid,
            'cppflagsAlways:' + '"' + flags + ld + libs + '"'
        ]))
    else:
        system(' '.join([
            dest(qex, 'configure'),
            'qmpdir:' + qmp,
            'qiodir:' + qio
        ]))

def install_hmc(build: str, bin: str):
    system('make cpp alphashmc')
    dustbin = dest(dest(build, 'bin'), 'alphashmc')
    symlink(dustbin, bin)

def install_gradient_flow(build: str, bin: str): 
    system('make cpp alphasflow')
    dustbin = dest(dest(build, 'bin'), 'alphasflow')
    symlink(dustbin, bin)

def hmc_regress(): system('./bin/alphashmc')

if __name__ == '__main__':
    args = args()
    path = args.prefix
    qex = args.qex
    compile_gradient_flow = args.compile_gradient_flow
    run_hmc_regression_test = args.run_hmc_regression_test
    grid_backend = args.grid_backend

    if not isdir(path, 'build'): mkdir(path, 'build')
    if not isdir(path, 'deps'): mkdir(path, 'deps')
    if not isdir(path, 'bin'): mkdir(path, 'bin')

    deps = dest(path, 'deps')
    build = dest(path, 'build')
    bin = dest(path, 'bin')

    cd(path, 'deps')
    install_qop_qio(qex)
    qmp = dest(deps, 'qmp')
    qio = dest(deps, 'qio')

    nim = dest(deps, NIMV)
    install_nim(nim)
    nimpath = dest(nim, 'bin')
    nimexec = dest(nimpath, 'nim')
    nimbleexec = dest(nimpath, 'nimble')

    mdevolve = dest(deps, 'mdevolve')
    install_mdevolve(mdevolve, nimexec, nimbleexec)

    if grid_backend: grid = install_grid(deps, args.system, args.build_cpus)
    else: grid = ''

    cd(path, 'build')
    configure(qex, qmp, qio, grid, nimexec, grid_backend)

    install_hmc(build, bin)
    if compile_gradient_flow: install_gradient_flow(build, bin)

    if run_hmc_regression_test: hmc_regress()

    # Final printout reminding folks how to use both binary files

