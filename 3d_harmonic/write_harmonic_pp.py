import numpy as np
import matplotlib.pyplot as plt


def harmonic(r, k=20.0, r0=2.20, cutoff=2.75,
             grad=False):

    pot = (k/2)*(r - r0)**2
    d_pot = -k*(r - r0)

    pot = np.where(r > cutoff, 0.0, pot)
    # pot = np.where(r < inner_cut, 0.0, pot)

    d_pot = np.where(r > cutoff, 0.0, d_pot)
    # d_pot = np.where(r < inner_cut, 0.0, d_pot)

    if grad:
        return pot, d_pot
    else:
        return pot


def write_file(k, r0, cutoff, filename="harmonic_potential.txt"):
    header = \
    f"""# DATE: 2025-03-10  UNITS: lj CONTRIBUTOR: Mark Potter
# Harmonic potential w/ cutoff

COS_PP
N 5000 R 1.0 {1.40:.2f}
    """

    inds = np.arange(1, 5001, 1)
    r = np.linspace(1.0, 1.40, 5000)

    pot, force = harmonic(r, k=k, r0=r0, cutoff=cutoff,
                          grad=True)

    table_arr = np.stack([inds, r, pot, force], axis=1)
    fmt = ["%d", "%.7f", "%.7f", "%.7f"]

    np.savetxt(filename, table_arr, fmt=fmt,
               header=header, comments='')


if __name__ == "__main__":
    write_file(4.0, 3.92, 4.10)
