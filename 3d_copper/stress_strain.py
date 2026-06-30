import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


def thermo_word_dict(header):
    d = {}
    for i, thermo_word in enumerate(header):
        d[thermo_word] = i
    return d


def last_N_lines(filename, N):
    bufsize = 8192
    fsize = os.stat(filename).st_size
    iter = 0

    with open(filename) as f:
        if bufsize > fsize:
            bufsize = fsize - 1
        fetched_lines = []

        while True:
            iter += 1
            if fsize-bufsize * iter < 0:
                return fetched_lines[-N:]

            f.seek(fsize-bufsize * iter)
            fetched_lines.extend(f.readlines())

            if len(fetched_lines) >= N or f.tell() == 0:
                return fetched_lines[-N:]


def find_run_end(filename):
    last_lines = [line.split() for line in last_N_lines(filename, 50)]

    for i, line in enumerate(last_lines):
        if line and line[0] == 'Loop':
            last_ts = last_lines[i-1][0]
            if last_ts == 'Fix':
                # only seems to happen when last ts == fix frequency
                last_ts = last_lines[i-2][0]
            return int(last_ts)


def filtered_log_lines(filename):
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("Fix halt condition"):
                continue
            yield line


def read_logfile(filename, num_steps=None, ts=None, header=None, header_line=None):
    i = 1

    with open(filename, 'r') as file:
        for line in file:
            line = line.split()
            if len(line) > 0 and line[0] == 'Step' and header_line is None:
                header_line = i
                header = line
                start_step = int(next(file).split()[0])
                next_step = int(next(file).split()[0])
                if ts is None:
                    ts = next_step - start_step
                break
            i += 1

    if num_steps is None:
        end_step = find_run_end(filename)
        num_steps = int((end_step - start_step) / ts)

    d = thermo_word_dict(header)
    if len(header) == 14:
        print("Warning: fix halt message may be included in data")
    data = np.genfromtxt(filtered_log_lines(filename), max_rows=num_steps+1,
                         skip_header=header_line, invalid_raise=False)
    return d, data


def read_sparse_logfile(filename, find_break_steps=False):
    data = []

    with open(filename, 'r') as file:
        # Find the header line
        for line in file:
            line = line.split()
            if len(line) > 0 and line[0] == 'Step':
                header = line
                break

        file.seek(0)

        # Read data following any header line
        break_steps = []
        in_data = False
        for line in file:
            line = line.split()
            if len(line) > 0:
                if line[0] == 'Step':
                    in_data = True
                    continue
                if line[0] == 'Loop':
                    in_data = False
                    break_steps.append(data[-1][0])
                    continue
            if in_data:
                a = np.array(line).astype(float)
                data.append(a)

        data = np.array(data)

    d = thermo_word_dict(header)

    if find_break_steps:
        break_steps = np.array(break_steps).astype(int)
        return break_steps, d, data

    return d, data


def eng_stress_strain(Pyy, Ly, Lx0, Lz0):
    eng_strain = (Ly / Ly[0]) - 1
    eng_stress = -Pyy  # / (Lx0*Lz0)

    return eng_stress, eng_strain


def stress_indices(stress, factor):
    max_stress = np.max(stress)

    offset_greater = np.roll(stress, 1) >= (max_stress*factor)
    lesser = stress <= (max_stress*factor)

    return np.where(np.logical_and(offset_greater, lesser))[0]


def failure_strain(stress, strain, factor=0.5):
    ind = stress_indices(stress, factor)
    return strain[ind[-1]]


def toughness(stress, strain):
    return np.trapz(stress, strain)


def yield_strength(stress, strain):
    ind = stress_indices(stress, 1.0)
    return stress[ind[-1]]


def norm_post_yield_toughness(stress, strain):
    yield_ind = stress_indices(stress, 1.0)[-1]
    max_stress = np.max(stress)

    return toughness(stress[yield_ind:], strain[yield_ind:]) / max_stress


def last_strain(strain):
    # Useful for sims with a CN halt condition

    return strain[-1]


def stress_strain_from_logfile(logfile, ts=None):
    td, data = read_logfile(logfile, ts=ts)

    Pyy, Ly = data[:, td['Pyy']], data[:, td["Ly"]]
    Lx0, Lz0 = data[0, td['Lx']], data[0, td['Lz']]
    stress, strain = eng_stress_strain(Pyy, Ly, Lx0, Lz0)

    return stress, strain


def stress_strain_td_data(td, data):
    Pyy, Ly = data[:, td['Pyy']], data[:, td["Ly"]]
    Lx0, Lz0 = data[0, td['Lx']], data[0, td['Lz']]
    stress, strain = eng_stress_strain(Pyy, Ly, Lx0, Lz0)

    return stress, strain


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile")
    parser.add_argument("plot_name")
    args = parser.parse_args()

    stress, strain = stress_strain_from_logfile(args.logfile)

    fs = last_strain(strain)
    print(f"Failure strain: {fs}")
    print(f"Failure r: {(2**(1/6))*(1+fs)}")

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(strain, stress)

    ax.set_xlabel("Strain")
    ax.set_ylabel("Stress (RU)")

    # ax.set_xlim(0.0, 0.11)

    plt.tight_layout()
    plt.savefig(args.plot_name, dpi=400)
