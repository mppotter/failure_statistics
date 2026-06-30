import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--temperature', type=float)
parser.add_argument('--strain_rate', type=float)
args = parser.parse_args()

T = args.temperature
sr = args.strain_rate

f_strains = []

for i in range(32):
    print(i)
    single_f_strains = np.load(f"N50_sr{sr:.4f}_T{T:.1f}/{i}/f_strains.npy")
    f_strains.append(single_f_strains)

np.save(f"N50_sr{sr:.4f}_T{T:.1f}/f_strains.npy", np.asarray(f_strains).reshape(-1))
