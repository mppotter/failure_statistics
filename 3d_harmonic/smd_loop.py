import os
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

temps = np.array([0.00001, 0.0003, 0.0020, 0.0060, 0.0180, 0.0300])
strains = np.arange(0.000, 0.070, 0.001)

os.system("mkdir smd")

T = temps[3]

for strain in strains:
    print(strain)
    logname = f"smd/log.T{T:.5f}_s{strain:.3f}_0"
    if os.path.isfile(logname):
        print(f"{logname} exists, skipping")
        continue
    os.system(f"python3 generate_smd_harmonic.py 50 --strain {strain}")
    os.system(f"python3 smd_harmonic.py --temperature {T:.5f} --ensemble nvt " +
              f"--dimensionality 3 --total_integrated_atoms 50 " +
              f"--log_name smd/log.T{T:.5f}_s{strain:.3f}_0 " +
              f"--strain {strain:.3f} --data_file final.lmp " +
              f"--rank {rank}")
    os.system(f"flux run -n 1 -c 1 /usr/workspace/potter21/lammps/build/lmp -in in.harmonic_smd_{rank} -screen out")
    print(f"T: {T:.5f}, Strain: {strain:.3f}")
