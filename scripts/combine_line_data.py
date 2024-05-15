#!/usr/bin/env python3

import numpy as np
import pandas as pd
import h5py
from glob import glob
from tqdm import tqdm

import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='TODO')
parser.add_argument('-resolution', type=int, default=128,
                    help='Resolution to use')
parser.add_argument('-basepath', type=str, default='output/line_test_run',
                    help='Base path for input files')
parser.add_argument('-max_run_index', type=int, default=10,
                    help='Load data up to this run index')
parser.add_argument('-h5fname', type=str, default='dataset',
                    help='HDF5 dataset base name (without extension)')
args = parser.parse_args()

file_pattern = f'{args.basepath}*_sigmaz_{args.resolution}.npz'
files = sorted(glob(file_pattern))

tmp = np.load(files[0])
z = tmp['arr_0']
prev_file = ''
run_indices = []
h5fname = f'{args.h5fname}_{args.resolution}.hdf5'

with h5py.File(h5fname, 'w') as h5f:
    dset = h5f.create_dataset('z_line', data=z)

    for f in tqdm(files):
        run_index = int(f.replace(args.basepath, '').split('_')[0])

        if run_index > args.max_run_index:
            continue

        tmp = np.load(f)
        run_indices.append(run_index)
        cycle = tmp['cycle']
        h5path = f'run_{run_index}'

        # Load log file if it is a new run
        fname = f'{args.basepath}{run_index}_log.txt'
        new_run = fname != prev_file

        if new_run:
            log = pd.read_csv(fname, delim_whitespace=True)
            prev_file = fname
            dset = h5f.create_dataset(f'{h5path}/radius', data=log['x.2'])
            dset = h5f.create_dataset(f'{h5path}/voltage', data=log['voltage'])
            dset = h5f.create_dataset(f'{h5path}/time', data=log['time'])

        dset = h5f.create_dataset(f'{h5path}/sigmaz_{cycle}',
                                  data=tmp['uniform_data'])

        # Load line potential and interpolate to right grid
        fname = f'{args.basepath}{run_index}_line_{cycle:06d}.txt'
        tmp = np.loadtxt(fname)

        phi = tmp[:, 2]
        z_phi = tmp[:, 1]
        phi_z = np.interp(z, z_phi, phi)
        dset = h5f.create_dataset(f'{h5path}/phiz_{cycle}', data=phi_z)

        # Load rhs
        rhs = np.load(f.replace('sigmaz', 'rhs'))['uniform_data']
        dset = h5f.create_dataset(f'{h5path}/rhs_{cycle}', data=rhs)

        # Load sigma
        sigma = np.load(f.replace('sigmaz', 'sigma'))['uniform_data']
        dset = h5f.create_dataset(f'{h5path}/sigma_{cycle}', data=sigma)

        if new_run:
            tmp = np.load(f.replace('sigmaz', 'rhs'))
            h5f.create_dataset(f'{h5path}/r_grid', data=tmp['arr_0'])
            h5f.create_dataset(f'{h5path}/z_grid', data=tmp['arr_1'])

    dset = h5f.create_dataset('run_indices', data=run_indices)

print(f'Wrote {h5fname}')
