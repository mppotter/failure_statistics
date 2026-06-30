import os
import numpy as np

strains = np.arange(0.000, 0.056, 0.001)
T = 200.0

for strain in strains:
    for i in range(1):
        os.system(f"bash run_copper_smd.sh {i} {strain:.3f} {T:.1f}")
        print(f"T: {T:.1f}, Strain: {strain:.3f}, i: {i}")
