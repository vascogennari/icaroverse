import os, json, pickle, re
from optparse import OptionParser
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

def rename(par):
    if par in pars_pipeline_names: return pars_pipeline_names[par]
    else:                          return par

def check_format(event_dirname, n):
    search = re.search("event_[0-9]{"+str(n)+"}(?P<extension>_?\w*)", event_dirname)
    return search.group('extension') == ''



def main():

    parser = OptionParser()
    parser.add_option('-d', '--pop-dir', type='string', metavar='pop_dir', default = None)
    parser.add_option('-n', '--ev-filename-format', type='int', metavar='ev_filename_format', default = 4)
    parser.add_option('-s', '--snr-restriction', action="store_true", dest='snr_restriction', default = False)
    parser.add_option('-t', '--snr-threshold', type='float', metavar='snr_threshold', default = 12.)
    (opts, _) = parser.parse_args()

    pop_dirname = os.path.basename(opts.pop_dir)

    pe_dirpath = os.path.join(opts.pop_dir, "parameter_estimation")
    event_dirname_list = [event_dirname for event_dirname in os.listdir(pe_dirpath) if 'event' in event_dirname and check_format(event_dirname, opts.ev_filename_format)]
    N_events = len(event_dirname_list)

    print(f"\n * Combining PE results of the {N_events} events from population directory {pop_dirname}\n")

    if opts.snr_restriction:
        injected_filepath = os.path.join(opts.pop_dir, "population_observed.pickle")
        with open(injected_filepath, 'rb') as f:
            snr_arr = pickle.load(f)['snr']
    else:
        snr_arr = np.full(N_events, np.inf)

    combined_results = {}
    for i, event_dirname in tqdm(enumerate(sorted(event_dirname_list)), desc="Looping over events", unit=" event", total=N_events):

        ev_idx = int(re.search("event_(?P<id>[0-9]+)", event_dirname).group('id'))
        if snr_arr[ev_idx] < opts.snr_threshold:
            continue

        try:
            event_result_filepath = os.path.join(pe_dirpath, event_dirname, "sampler/label_result.json")
            with open(event_result_filepath, 'r') as f:
                event_posterior_samples_dict = json.load(f)['posterior']['content']
            

            for par in event_posterior_samples_dict:

                if rename(par) not in combined_results:
                    combined_results[rename(par)] = {i: np.array(event_posterior_samples_dict[par])}
                else:
                    combined_results[rename(par)][i] = np.array(event_posterior_samples_dict[par])
        except FileNotFoundError:
            print(f"{event_dirname} has no available result.")
            continue
    
    snr_part = f"_snr{opts.snr_threshold:.0f}" if opts.snr_restriction else ""
    combined_filename = f"combined_PE_samples_{pop_dirname}{snr_part}.pickle"
    combined_filepath = os.path.join(opts.pop_dir, combined_filename)
    with open(combined_filepath, 'wb') as f:
        pickle.dump(combined_results, f, protocol = pickle.HIGHEST_PROTOCOL)

    print(f"\n * Finished ! \n")


# Exerestrictione the main function.
if __name__=='__main__': 
    main()
