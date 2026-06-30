import numpy as np
import sim_cell
import argparse


def write_lammps_file(filename, sim_cell):
    header = ["# LAMMPS data file",
              "\n{} atoms".format(sim_cell.num_atoms),
              "{} atom types".format(sim_cell.num_types),
              "\n{} {} xlo xhi".format(sim_cell.xlo, sim_cell.xhi),
              "{} {} ylo yhi".format(sim_cell.ylo, sim_cell.yhi),
              "{} {} zlo zhi".format(sim_cell.zlo, sim_cell.zhi),
              # "\nMasses",
              # "\n1 1.0080",
              # "2 12.0107",
              # "3 15.9994",
              "\nAtoms # atomic\n",
              ]

    header_string = "\n".join(header)

    np.savetxt(filename, sim_cell.atoms, delimiter=" ", header=header_string, comments="",
               fmt=['%d', '%d', '%.5f', '%.5f', '%.5f',])


def add_coords(init_c, add_c):
    c = init_c
    c[1:] = init_c[1:] + add_c[1:]
    return c


def add_atom(atom, delta_y, atom_type):
    c = np.array([atom_type, 0.0, delta_y, 0.0])

    return add_coords(c, atom)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("num_atoms", type=int)
    parser.add_argument("--strain", type=float, default=0.0)
    parser.add_argument("--perturb_dist", type=float, default=0.0)
    args = parser.parse_args()
    num_atoms = args.num_atoms
    rng = np.random.default_rng()

    delta_y = 2.31

    # Generate the chain, one atom
    atoms_array = np.empty(shape=[num_atoms, 4])
    for i in range(0, num_atoms, 1):
        if i == 0:
            prev_atom = np.array([1, 5.0, 0.0, 5.0])
        else:
            prev_atom = atoms_array[i-1]

        # noise = rng.uniform(low=0.95, high=1.05)
        atom_type = 1
        atom = add_atom(prev_atom, delta_y, atom_type)
        atoms_array[i] = atom

    # Write the datafile
    pre_array = np.arange(1, num_atoms+1, 1).reshape(-1, 1)
    atoms = np.hstack([pre_array, atoms_array])

    yhi = (num_atoms) * delta_y

    # Apply a strain
    atoms[:, 3] = atoms[:, 3] * (1.0+args.strain)
    yhi = yhi * (1.0+args.strain)

    # Perturb atom 1 some distance
    atoms[0, 3] = atoms[0, 3] + args.perturb_dist

    num_atoms = atoms.shape[0]
    sc = sim_cell.SimCell(num_atoms, 2, 0.0, 10.0, 0.0, yhi, 0.0, 10.0, atoms)
    write_lammps_file("final.lmp", sc)
