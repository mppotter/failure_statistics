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
    os.system(f"mkdir N{N}_sr{sr:g}_T{temp:.1f}")
comm.Barrier()

os.system(f"mkdir N{N}_sr{sr:g}_T{temp:.1f}/{rank}")
os.chdir(f"N{N}_sr{sr:g}_T{temp:.1f}/{rank}")

# Init run
os.system(f"python3 ../../generate_tensile_Cu.py {N} --strain 0.000")
num_sims = 1024 // 32
f_strains = []

for i in range(num_sims):
    if not os.path.isfile(f"log.{i}"):
        os.system("python3 ../../single_copper_chain_file.py " +
                  f"--temperature {temp} --strain_rate {sr:g} " +
                  f"--dimensionality 3 --total_integrated_atoms {N} " +
                  "--meam_file ../../Cu.meam --library_file ../../library.meam " +
                  "--data_file final.lmp")

        os.system("flux run -n 1 -c 1 /usr/workspace/potter21/lammps/build/lmp -in in.copper " +
                  "-screen out")
        os.system(f"mv log.tensile log.{i}")

    f_strains.append(end_strain_from_log(f"log.{i}"))

np.save("f_strains.npy", np.asarray(f_strains))
