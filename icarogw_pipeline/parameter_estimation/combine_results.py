import os, json, pickle, re
from argparse import ArgumentParser
from tqdm import tqdm
import numpy as np

# ----------------------------------------------------------------------------------------- #
# Read a population folder, and combine all the individual events PE into a single dict     #
# ----------------------------------------------------------------------------------------- #

pars_pipeline_names = {
    'luminosity_distance': 'dL',
    'mass_1': 'm1d',
    'mass_2': 'm2d',
    'redshift': 'z',
    'mass_1_source': 'm1s',
    'mass_2_source': 'm2s',
}

def union_multiple_sets(*sets):
    res = set()
    for s in sets:
        res.update(s)
    return res


def pars_to_keep_fun(*options):

    categories = {
        'mass':     {'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'chirp_mass', 'total_mass', 'symmetric_mass_ratio', 'mass_ratio', 'chirp_mass_source', 'total_mass_source'},
        'distance': {'luminosity_distance', 'redshift', 'comoving_distance'},
        'spin':     {'chi_1', 'chi_2', 'a_1', 'cos_tilt_1', 'a_2', 'cos_tilt_2', 'phi_jl', 'phi_12', 'tilt_1', 'tilt_2', 'spin_1x', 'spin_1y', 'spin_1z', 'spin_2x', 'spin_2y', 'spin_2z', 'phi_1', 'phi_2', 'chi_eff', 'chi_1_in_plane', 'chi_2_in_plane', 'chi_p'},
        'skyloc':   {'ra', 'dec'},
        'sampling': {'log_likelihood', 'log_prior', 'L1_log_likelihood', 'H1_log_likelihood', 'V1_log_likelihood', 'K1_log_likelihood'},
        'detector': {'reference_frequency', 'waveform_approximant', 'minimum_frequency', 'L1_matched_filter_snr', 'L1_optimal_snr', 'H1_matched_filter_snr', 'H1_optimal_snr', 'V1_matched_filter_snr', 'V1_optimal_snr', 'K1_matched_filter_snr', 'K1_optimal_snr'},
        'other':    {'theta_jn', 'psi', 'phase', 'geocent_time', 'iota'},
    }

    for option in options:
        if not (option in categories or option in {'all', 'minimal'}) : raise KeyError(f"Unknown option '{option}' to specify parameters to keep. Please give one or more options from"+"\n"+"\n\t".join(list(categories.keys())))
        else: continue

    if 'all' in options:
        parameters_to_keep = union_multiple_sets(*[categories[option] for option in categories])
    elif 'minimal' in options:
        parameters_to_keep = {'mass_1', 'mass_2', 'luminosity_distance'}
    else:
        parameters_to_keep = union_multiple_sets(*[categories[option] for option in options])
    
    return parameters_to_keep


def rename(par):
    if par in pars_pipeline_names: return pars_pipeline_names[par]
    else:                          return par


def check_format(event_dirname, n):
    search = re.search("event_[0-9]{"+str(n)+"}(?P<extension>_?\w*)", event_dirname)
    return search.group('extension') == ''



def main():

    parser = ArgumentParser()
    parser.add_argument('-d', '--pop_dir', type=str, default = None)
    parser.add_argument('-n', '--ev_filename_format', type=int, default = 4)
    parser.add_argument('-s', '--snr_restriction', action="store_true", dest='snr_restriction', default = False)
    parser.add_argument('-t', '--snr_threshold', type=float, default = 12.)
    parser.add_argument('-p', '--parameters', nargs='+', default = 'all', help="Options: mass, distance, spin, skyloc, sampling, detector, other (options are stackable), all (equivalent to all of the previous ones stacked), minimal (saves only m1d, m2d, dL)")
    args = parser.parse_args()

    pop_dirname = os.path.basename(args.pop_dir)

    pe_dirpath = os.path.join(args.pop_dir, "parameter_estimation")
    event_dirname_list = [event_dirname for event_dirname in os.listdir(pe_dirpath) if 'event' in event_dirname and check_format(event_dirname, args.ev_filename_format)]
    N_events = len(event_dirname_list)

    print(f"\n * Combining PE results of the {N_events} events from population directory {pop_dirname}")

    if args.snr_restriction:
        injected_filepath = os.path.join(args.pop_dir, "population_observed.pickle")
        with open(injected_filepath, 'rb') as f:
            snr_arr = pickle.load(f)['snr']
    else:
        snr_arr = np.full(N_events, np.inf)
    
    pars_to_keep = pars_to_keep_fun(*args.parameters)
    print("\n * Only the value of the following parameters for each PE sample will be saved in the resulting file:\n\n\t"+"\n\t".join(sorted(list(pars_to_keep)))+"\n")

    combined_results = {}
    for i, event_dirname in tqdm(enumerate(sorted(event_dirname_list)), desc="Looping over events", unit=" event", total=N_events):

        ev_idx = int(re.search("event_(?P<id>[0-9]+)", event_dirname).group('id'))
        if snr_arr[ev_idx] < args.snr_threshold:
            continue

        try:
            event_result_filepath = os.path.join(pe_dirpath, event_dirname, "sampler/label_result.json")
            with open(event_result_filepath, 'r') as f:
                event_posterior_samples_dict = json.load(f)['posterior']['content']
            

            for par in (set(event_posterior_samples_dict.keys()) & pars_to_keep):

                if rename(par) not in combined_results:
                    combined_results[rename(par)] = {i: np.array(event_posterior_samples_dict[par])}
                else:
                    combined_results[rename(par)][i] = np.array(event_posterior_samples_dict[par])
        except FileNotFoundError:
            print(f"{event_dirname} has no available result.")
            continue
    
    snr_part = f"_snr{args.snr_threshold:.0f}" if args.snr_restriction else ""
    if 'all' in args.parameters:       pars_part = "_all"
    elif 'minimal' in args.parameters: pars_part = "_minimal" 
    else:                              pars_part = "_"+"-".join(list(args.parameters))
    combined_filename = f"combined_PE_samples_{pop_dirname}{snr_part}{pars_part}.pickle"
    combined_filepath = os.path.join(args.pop_dir, combined_filename)
    with open(combined_filepath, 'wb') as f:
        pickle.dump(combined_results, f, protocol = pickle.HIGHEST_PROTOCOL)

    print(f"\n * Finished ! \n")


# Execute the main function.
# Execute the main function.
if __name__=='__main__': 
    main()
