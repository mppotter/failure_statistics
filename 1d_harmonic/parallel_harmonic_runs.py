import numpy as np
from mpi4py import MPI
import os
from compile_single_sim_data import end_strain_from_log
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--temperature', type=float)
parser.add_argument('--strain_rate', type=float)
args = parser.parse_args()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

N = 50
temp = args.temperature
sr = args.strain_rate

if rank == 0:
    os.system(f"mkdir N{N}_sr{sr:g}_T{temp:.4f}")
comm.Barrier()

os.system(f"mkdir N{N}_sr{sr:g}_T{temp:.4f}/{rank}")
os.chdir(f"N{N}_sr{sr:g}_T{temp:.4f}/{rank}")

# Init run
os.system(f"python3 ../../generate_harmonic.py {N} --strain 0.000")
num_sims = 1024 // 32
f_strains = []

for i in range(num_sims):
    os.system("python3 ../../single_harmonic_chain_file.py " +
              f"--temperature {temp} --strain_rate {sr:g} " +
              f"--dimensionality 1 --total_integrated_atoms {N} " +
              "--harmonic_filename harmonic_potential.txt " +
              "--data_file final.lmp")

    os.system("flux run -n 1 -c 1 /usr/workspace/potter21/lammps/build/lmp -in in.harmonic " +
              "-screen out")
    os.system(f"mv log.tensile log.{i}")

    f_strains.append(end_strain_from_log(f"log.{i}"))

np.save("f_strains.npy", np.asarray(f_strains))
