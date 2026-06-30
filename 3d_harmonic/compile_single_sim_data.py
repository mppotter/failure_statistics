import numpy as np
import argparse
import stress_strain as ss


def recover_short_data(curve, length, continuation=False):
    # Note: assumes that -2 axis is time/strain, -1 axis is column
    reduced_shape = np.array(curve.shape)
    reduced_shape[-2] = length
    reduced_curve = np.zeros(shape=reduced_shape)

    if continuation:
        # "Continue" the data by extrapolating
        step_size = curve[..., 1, 0] - curve[..., 0, 0]
        # Case for a "constant" continuation
        if step_size == 0:
            return recover_short_data(curve, length, continuation=False)
        start = curve[..., 0, 0]
        end = length*step_size + start
        reduced_curve[..., 0] = np.arange(start, end, step_size)
        reduced_curve[..., :curve.shape[-2], 1:] = curve[..., 1:]
    else:
        # Data is frozen at last value
        reduced_curve[..., :curve.shape[-2], :] = curve[..., :]
        reduced_curve[..., curve.shape[-2]:, :] = curve[..., -1, :]

    return reduced_curve


def reduce_data(data_curve, num_samples, continuation=False):
    # Note: assumes that -2 axis is time/strain, -1 axis is column
    if len(data_curve.shape) == 1:
        data_curve = np.expand_dims(data_curve, 1)

    step = data_curve.shape[-2] / num_samples

    if step < 1:
        print("Warning: recovering short data")
        reduced_curve = recover_short_data(data_curve, num_samples, continuation)
    else:
        selection = np.arange(0, data_curve.shape[-2], step)
        selection = np.round(selection).astype(int)

        # Edge cases
        selection = selection[:num_samples]
        selection[-1] = data_curve.shape[-2] - 1

        reduced_curve = data_curve[..., selection, :]

    return reduced_curve


def first_coord_drop(arr):
    # Takes a single coord array, returns first index where arr < arr[0]
    def first(arr):
        ind = np.asarray(arr < arr[0]).nonzero()[0]
        if len(ind) == 0:
            return -1
        return ind[0]

    return first(arr)


def end_strain_from_log(logname):
    td, data = ss.read_logfile(logname)
    stress, strain = ss.stress_strain_td_data(td, data)

    coord = data[:, td['c_min_coord']]
    
    f_inds = first_coord_drop(coord)
    f_strain = strain[f_inds]
    
    return f_strain
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("batch_index", type=int)
    parser.add_argument("specimen_index", type=int)
    parser.add_argument("reduced_size", type=int)
    parser.add_argument("datatype", type=str, default='strain')
    args = parser.parse_args()

    # Load datas
    td, data = ss.read_logfile("log.tensile")
    stress, strain = ss.stress_strain_td_data(td, data)

    coord = data[:, td['c_min_coord']]

    # Create reduced r_avg data
    if args.datatype == 'all':
        # all returns [timestep, stress, strain]
        r_avg = data[:, [td['Step'], td['Pyy']]]
        r_avg = np.hstack([r_avg, strain.reshape(-1, 1)])

        reduced_timestep = reduce_data(data[:, td['Step']], args.reduced_size,
                                       continuation=True)
        reduced_stress = reduce_data(data[:, td['Pyy']], args.reduced_size,
                                     continuation=False)
        reduced_strain = reduce_data(strain, args.reduced_size,
                                     continuation=True)
        reduced_r_avg = np.hstack([reduced_timestep, reduced_stress,
                                   reduced_strain.reshape(-1, 1)])
    else:
        if args.datatype == 'strain':
            # Note: previously, this option recorded the bond length
            r_avg = strain
            continuation=True
        elif args.datatype == 'timestep':
            r_avg = data[:, td['Step']]
            continuation=True
        elif args.datatype == 'stress':
            r_avg = data[:, td['Pyy']]
            continuation=False
        else:
            raise ValueError(f"Datatype must be 'strain', 'timestep', " + \
                             f"'stress', or 'all'; not '{args.datatype}'.")
        reduced_r_avg = reduce_data(r_avg, args.reduced_size,
                                    continuation=continuation)

    # Create reduced binned coordination data
    reduced_coord = reduce_data(coord, args.reduced_size, continuation=False)

    # Calculate failure r_avgs
    f_inds = first_coord_drop(coord)
    f_r_avgs = r_avg[f_inds]
    print(f_r_avgs)

    np.save("spec_inds.npy", np.array([args.batch_index, args.specimen_index]))
    np.save("r_avg.npy", reduced_r_avg)
    np.save("coord_data.npy", reduced_coord)
    np.save("f_r_avg.npy", f_r_avgs)
