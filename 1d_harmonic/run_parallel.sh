#!/bin/bash
#flux: -N 1
#flux: -t 60
#flux: -q pbatch

T=0.0191
sr=0.0005

source /p/vast1/potter21/vast_torch/bin/activate
flux run -n 32  python3 parallel_harmonic_runs.py --temperature $T --strain_rate $sr

python3 collect_results.py --temperature $T --strain_rate $sr
