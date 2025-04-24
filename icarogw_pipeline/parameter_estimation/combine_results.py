import os, json, pickle
from optparse import OptionParser
from tqdm import tqdm
import numpy as np

# ----------------------------------------------------------------------------------------- #
# Read a population folder, and combine all the individual events PE into a single dict     #
# ----------------------------------------------------------------------------------------- #

rename = {
    'luminosity_distance': 'dL',
    'mass_1': 'm1d',
    'mass_2': 'm2d',
    'redshift': 'z',
    'mass_1_source': 'm1s',
    'mass_2_source': 'm2s',
}

def main():

    parser = OptionParser()
    parser.add_option('-d', '--pop-dir', type='string', metavar = 'pop_dir', default = None)
    (opts, _) = parser.parse_args()

    pop_dirname = os.path.basename(opts.pop_dir)

    pe_dirpath = os.path.join(opts.pop_dir, "parameter_estimation")
    event_dirname_list = os.listdir(pe_dirpath)
    N_events = len(event_dirname_list)

    print(f"\n * Combining PE results of the {N_events} events from population directory {pop_dirname}")

    combined_results = {}
    for i, event_dirname in tqdm(enumerate(sorted(event_dirname_list))):

        event_result_filepath = os.path.join(pe_dirpath, event_dirname, "sampler/label_result.json")
        with open(event_result_filepath, 'r') as f:
            event_posterior_samples_dict = json.load(f)['posterior']['content']

        for par in event_posterior_samples_dict:

            if rename[par] not in combined_results:
                combined_results[rename[par]] = {i: np.array(event_posterior_samples_dict[par])}
            else:
                combined_results[rename[par]][i] = np.array(event_posterior_samples_dict[par])
    
    combined_filename = f"combined_posterior_{pop_dirname}"
    with open(combined_filename, 'wb') as f:
        pickle.dump(combined_results, f, protocol = pickle.HIGHEST_PROTOCOL)


# Execute the main function.
if __name__=='__main__': 
    main()
