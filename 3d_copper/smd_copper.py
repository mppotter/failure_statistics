import argparse
import numpy as np
import write_smd_harmonic_file as harmonic


def system_init(w, data_file):
    cmds = [
        "units          metal",
        "atom_style     atomic",
        "boundary       p p p",
        f"read_data      {data_file}",
    ]

    return cmds


def dim_conditions(total_integrated_atoms, dimensionality, ensemble):
    if total_integrated_atoms == 0 and dimensionality != 3:
        raise ValueError("Need number of atoms to adjust DOF for temp.")

    if dimensionality != 3:
        if dimensionality == 1:
            factor = 2*total_integrated_atoms
            dim_cond = [
                "fix            one_dim all setforce 0.0 NULL 0.0",
                "compute        one_temp all temp",
                f"compute_modify one_temp extra/dof {factor}",
            ]
        elif dimensionality == 2:
            factor = 1*total_integrated_atoms
            dim_cond = [
                "fix            one_dim all setforce 0.0 NULL NULL",
                "compute        one_temp all temp",
                f"compute_modify one_temp extra/dof {factor}",
            ]
        else:
            raise ValueError(
                f"Dimensionality must be 1, 2, or 3; not {dimensionality}")
        if ensemble == 'nvt':
            fix_modify = "fix_modify    1 temp one_temp"
    else:
        dim_cond = []
        fix_modify = ''

    return dim_cond, fix_modify


def break_r_lookup(strain, break_r_file):
    arr = np.loadtxt(break_r_file)
    strains = arr[:, 0]
    break_rs = arr[:, 1]

    offset_greater = np.roll(strains, -1) >= strain
    lesser = strains <= strain

    ind = np.where(np.logical_and(offset_greater, lesser))[0][0]

    return break_rs[ind]


def set_coeff(meam_file, library_file, table_file, break_r, strain):
    if strain < -0.020:
        # i.e. never true
        F_0 = -0.200
    else:
        F_0 = 0.000

    # break_r = break_r * 1.00
    print(break_r)

    harmonic.write_file(2.0, break_r, 5.00, F_0=F_0)
    print(f"Break r: {break_r:.3f}")

    cmds = [
        "pair_style     hybrid/overlay meam table linear 5000",
        f"pair_coeff     * * meam {library_file} Cu Cu2 Cu3 Cu4 {meam_file} Cu Cu2 Cu3 Cu4",
        f"pair_coeff     1 2 table {table_file} COS_PP",
        f"pair_coeff     2 3 table {table_file} COS_PP",
        f"pair_coeff     1 4 table {table_file} COS_PP",
    ]

    return cmds


# def set_coeff(meam_file, library_file, table_file, break_r, strain):
#     cmds = [
#         "pair_style     meam",
#         f"pair_coeff    * * {library_file} Cu Cu2 {meam_file} Cu Cu2",
#     ]
# 
#     return cmds


def set_neighbors():
    cmds = [
        "neighbor       0.3 bin",
        "neigh_modify   every 10 delay 0 check no",
    ]

    return cmds


def set_thermo_style():
    style = "thermo_style   custom step temp pe etotal vol press pxx pyy pzz lx ly lz"

    cmds = [
        "thermo         20",
        style,
    ]

    return cmds


def init_velocities(temp, not_random, dimensionality):
    if not_random:
        seed = 488404201
    else:
        rng = np.random.default_rng()
        seed, = rng.integers(size=1, low=1, high=(2**31 - 1))

    if dimensionality == 1:
        return [f"velocity       all create {3*temp} {seed} temp one_temp dist gaussian",
                "velocity       all set 0.0 NULL 0.0",
                ]
    elif dimensionality == 2:
        return [f"velocity       all create {2*temp} {seed} temp one_temp dist gaussian",
                "velocity       all set 0.0 NULL 0.0",
                ]
    else:
        return [f"velocity       all create {temp} {seed}",
                ]


def eng_minimize():
    cmds = [
        "fix            1 all nve",
        "min_style      fire",
        "minimize       1.0e-4 1.0e-6 1000 10000",
        "unfix          1",
    ]

    return cmds


def nvt_relax(start_temp, end_temp, num_steps):
    cmds = [
        f"fix            1 all nvt temp {start_temp} {end_temp} 1.0",
        f"run            {num_steps}",
        "unfix          1",
    ]

    return cmds


def npt_relax(start_temp, end_temp, num_steps):
    cmds = [
        "variable       x_press equal pxx",
        "variable       y_press equal pyy",
        "variable       z_press equal pzz",
        f"fix            1 all npt temp {start_temp} {start_temp} 1.0"
        " x ${x_press} 0.0 1.0 y ${y_press} 0.0 1.0 z ${z_press} 0.0 1.0",
        f"run            {num_steps}",
        "unfix          1",
    ]

    return cmds


def smd_pull(w, couple_vel, ensemble, log_name, strain, break_r,
             fix_modify=None):
    # cutoff = break_r * 1.10
    # halt_condition = [
    #     "group          type1 type 1",
    #     f"compute        coord type1 coord/atom cutoff {cutoff} group type1",
    #     "compute        min_coord type1 reduce min c_coord",
    # ]

    # halt_condition.append("variable       min_coord equal c_min_coord")
    # halt_condition.append(
    #     "fix            10 all halt 20 v_min_coord < 1 error continue")

    if ensemble == 'nve':
        fix_cmd = ["fix            1 all nve",
                   ]
    elif ensemble == 'nvt':
        fix_cmd = [f"fix            1 all nvt temp {w} {w} 0.1",
                   fix_modify,
                   ]
    else:
        raise ValueError(f"Unrecognized ensemble: {ensemble}")

    # equil_sep = 2**(1/6) * (1+strain)
    equil_sep = 0.0

    # Carbon ids
    lo_id = 25
    hi_id = 26
    num_steps = int(200_000)
    cmds = [
        f"log            {log_name}",
        f"group          center_lo id {lo_id}",
        f"group          center_hi id {hi_id}",
        "fix            4 center_lo smd cvel 2.0 0.001 couple " +
        f"center_hi auto auto auto {equil_sep:.4f}",
        "thermo_style   custom step temp pe vol lx ly lz"
        " pxx pyy f_4[1] f_4[2] f_4[3] f_4[4] f_4[5] f_4[6] f_4[7]",
        # f"fix            4 all deform 1 y erate {strain_rate}",
        f"run            {int(num_steps)}",
    ]

    unfix = [
        # "unfix          10",
        "unfix          1",
        "unfix          4",
    ]

    if ensemble == 'nve':
        unfix.append("unfix     3")

    return [*fix_cmd, *cmds, *unfix]


def write_input(cmds, filename):
    with open(filename, 'w') as f:
        for cmd in cmds:
            f.write(cmd)
            f.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--temperature', default=100.0, type=float)
    parser.add_argument('--strain_rate', default=0.0001, type=float)
    parser.add_argument('--log_name', default='log.tensile')
    parser.add_argument('--not_random', action='store_true', default=False)
    parser.add_argument('--ensemble', default='nvt',
                        type=str, help="either nve or nvt")
    parser.add_argument('--couple_vel', type=float, default=0.0001)
    parser.add_argument('--data_file', type=str, default='final.lmp')
    parser.add_argument('--dimensionality', type=int, default=3)
    parser.add_argument('--total_integrated_atoms', type=int, default=-1)
    parser.add_argument('--init_strain', type=float, default=0.000)
    parser.add_argument('--meam_file', type=str)
    parser.add_argument('--library_file', type=str)
    parser.add_argument('--table_file', type=str)
    parser.add_argument('--strain', type=float)
    parser.add_argument('--break_r_file', type=str)
    args = parser.parse_args()

    start_temp, end_temp = args.temperature, args.temperature

    final_strain = 0.005

    break_r = break_r_lookup(args.strain, args.break_r_file)

    # Simulation boilerplate
    init = [
        *system_init(start_temp, args.data_file),
        *set_coeff(args.meam_file, args.library_file, args.table_file,
                   break_r, args.strain),
        *set_neighbors(),
        "timestep       0.010",
    ]

    # Change dimensionality of sim
    dim_cond, fix_modify = dim_conditions(args.total_integrated_atoms,
                                          args.dimensionality,
                                          args.ensemble, )

    # Minimize
    minimize = eng_minimize()

    # Initialize velocities
    init_vel = [
        *set_thermo_style(),
        *init_velocities(start_temp, args.not_random, args.dimensionality)
    ]

    # Uniaxial tension with cannonical ensemble
    smd = smd_pull(start_temp, args.couple_vel, args.ensemble,
                   args.log_name, args.init_strain, break_r,
                   fix_modify=fix_modify,
                   )

    dump = [
        "# dump           1 all custom 50 dump.copper_smd id type x y z fx fy fz vx vy vz",
    ]

    cmds = [*init, *dim_cond, *init_vel,
            'reset_timestep       0', *dump, *smd,]

    write_input(cmds, "in.copper_smd")
