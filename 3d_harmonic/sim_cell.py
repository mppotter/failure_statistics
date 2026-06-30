import numpy as np
from dataclasses import dataclass


@dataclass
class SimCell():
    num_atoms: int
    num_types: int
    xlo: float
    xhi: float
    ylo: float
    yhi: float
    zlo: float
    zhi: float
    atoms: np.ndarray


def assign_type(sim_cell, func):
    sim_cell.num_types = 4
    is_weakphase = sim_cell.atoms[:, 1] <= 2
    assignments = func(sim_cell.atoms[:, 2], sim_cell.atoms[:, 3],
                       sim_cell.atoms[:, 4]).astype(int)
    factors = 2*np.logical_and(is_weakphase, assignments) +\
        -2*np.logical_and(np.logical_not(is_weakphase),
                          np.logical_not(assignments))
    sim_cell.atoms[:, 1] += factors


def load_lammps_file(filename):
    file = open(filename, 'r')
    num_atoms, num_types = None, None
    xlo, xhi = None, None
    ylo, yhi = None, None
    zlo, zhi = None, None
    atoms = None

    line_num = 0

    for line in file:
        line = line.split()
        line_num += 1

        if len(line) == 2 and line[1] == 'atoms':
            num_atoms = int(line[0])

        if len(line) == 3 and ' '.join(line[1:]) == 'atom types':
            num_types = int(line[0])

        if len(line) == 4:
            if ' '.join(line[2:]) == 'xlo xhi':
                xlo, xhi = float(line[0]), float(line[1])
            if ' '.join(line[2:]) == 'ylo yhi':
                ylo, yhi = float(line[0]), float(line[1])
            if ' '.join(line[2:]) == 'zlo zhi':
                zlo, zhi = float(line[0]), float(line[1])

        if line and line[0] == 'Atoms':
            if line[2] != 'atomic':
                raise ValueError("LAMMPS file is not \'atomic\' style")
            atoms = np.genfromtxt(filename, skip_header=line_num,
                                  max_rows=num_atoms, usecols=(0, 1, 2, 3, 4,))
            break

    file.close()

    return SimCell(num_atoms, num_types, xlo, xhi, ylo, yhi, zlo, zhi, atoms)


def get_ts_line_number(filename, ts_lo, ts_hi):
    start_lines = []

    with open(filename, 'r') as fp:
        offset = 0
        for i, line in enumerate(fp):
            if line == "ITEM: TIMESTEP\n":
                line = next(fp)
                offset += 1
                # print(f"i={i}, line={line.strip()}")
                timestep = int(line.strip())
                if timestep >= ts_lo and timestep <= ts_hi:
                    start_lines.append(i + offset)

    return start_lines


def load_dump_steps(filename, start_lines):
    fp = open(filename, 'r')

    step_nums = []
    timesteps = []
    offset = 1

    for i, line in enumerate(fp):
        if i + offset not in start_lines:
            continue

        # Keep adding lines to header_info until atom info section is reached
        header = []
        for j in range(8):
            line = line.split()
            header.append(line)
            offset += 1
            line = next(fp)

        step_num = int(header[1][0])

        num_atoms, num_types = int(header[3][0]), None
        xlo, xhi = float(header[5][0]), float(header[5][1])
        ylo, yhi = float(header[6][0]), float(header[6][1])
        zlo, zhi = float(header[7][0]), float(header[7][1])

        atoms = np.loadtxt(filename, skiprows=i+offset,
                           max_rows=num_atoms)

        timestep = SimCell(num_atoms, num_types, xlo,
                           xhi, ylo, yhi, zlo, zhi, atoms)

        step_nums.append(step_num)
        timesteps.append(timestep)

    fp.close()

    return step_nums, timesteps


def load_dump_file(filename, first_ts=-1, last_ts=np.inf, step_nums=False):
    line_nums = get_ts_line_number(filename, first_ts, last_ts)
    steps, timesteps = load_dump_steps(filename, line_nums)

    if step_nums:
        return steps, timesteps

    return timesteps


def write_lammps_file(filename, sim_cell):
    header = ["# LAMMPS data file",
              "{} atoms".format(sim_cell.num_atoms),
              "{} atom types".format(sim_cell.num_types),
              "{} {} xlo xhi".format(sim_cell.xlo, sim_cell.xhi),
              "{} {} ylo yhi".format(sim_cell.ylo, sim_cell.yhi),
              "{} {} zlo zhi".format(sim_cell.zlo, sim_cell.zhi),
              "",
              "Atoms # atomic",
              ""]

    header_string = "\n".join(header)

    np.savetxt(filename, sim_cell.atoms, delimiter=" ", header=header_string, comments="",
               fmt=['%d', '%d', '%.15f', '%.15f', '%.15f',])


def load_npy_files(atoms_filename, header_filename):
    meta = np.load(header_filename)
    atoms = np.load(atoms_filename)

    return SimCell(int(meta[0, 0]), int(meta[0, 1]), meta[1, 0], meta[1, 1], meta[2, 0],
                   meta[2, 1], meta[3, 0], meta[3, 1], atoms)


def write_npy_files(atoms_filename, header_filename, sim_cell):
    header_array = np.array([[sim_cell.num_atoms, sim_cell.num_types,],
                             [sim_cell.xlo, sim_cell.xhi,],
                             [sim_cell.ylo, sim_cell.yhi,],
                             [sim_cell.zlo, sim_cell.zhi,],
                             ], dtype=float)

    np.save(header_filename, header_array)
    np.save(atoms_filename, sim_cell.atoms)
