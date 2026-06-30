# Failure Statistics of Single-Particle Chains

This repository contains the code, scripts, and reference files needed to reproduce the data reported in *Failure Statistics of Single-Particle Chains*.

The paper examines three single-particle chain systems: a 1D truncated harmonic chain, a 3D truncated harmonic chain, and a 3D copper chain with interactions modeled using an MEAM potential. The repository is organized into three corresponding directories. Below, we briefly describe how the scripts and code can be used to reproduce the reported results.

## `1d_harmonic`

This directory contains the files for producing the 1D truncated harmonic system results. The basic entry point is `run_parallel.sh`, which invokes the Python script `parallel_harmonic_runs.py` across, by default, 32 ranks at a specified temperature and strain rate.

In the configuration provided here, this script runs 1,024 independent simulations and records the failure strain for each simulation. These 1,024 failure strains are automatically collated into a file named `f_strains.npy`.

Minor changes may be required to reflect the location of your LAMMPS executable and the job scheduler used on your system.

## `3d_harmonic`

As with `1d_harmonic`, the basic entry point is `run_parallel.sh`. By default, 1,024 simulations are performed across 32 ranks, and the resulting failure strains are stored in `f_strains.npy`.

We also provide the code and scripts needed to generate steered molecular dynamics (SMD) data for calculating the energy barrier. Running `smd_loop.py` in its current form will run approximately 70 simulations at a specific temperature across a range of strains.

For each strain, the energy barrier $\Delta E(\epsilon)$ can be obtained by integrating the force on the SMD spring, `f_4[4]`, with respect to the atomic separation, `f_4[6]`. The lower integration limit is the start of the simulation, and the upper integration limit is the point at which the force on the SMD spring crosses zero.

## `3d_copper`

As with the truncated harmonic systems, the basic entry point is `run_parallel.sh`. By default, 1,024 simulations are performed across 32 ranks, and the resulting failure strains are stored in `f_strains.npy`.

Similar to the 3D truncated harmonic system, `smd_loop.py` can be used to run the SMD simulations needed to calculate the energy barrier as a function of strain. Since this system requires a harmonic restoring force to be turned on only when a bond would otherwise break, we provide `Cu_breaking_rs.txt` as a lookup table for these critical bond lengths.

## Citation

A DOI will be added to this repository at the time of publication.

