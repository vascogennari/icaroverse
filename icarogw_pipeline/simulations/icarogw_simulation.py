import os, sys, shutil, configparser
import numpy as np, pickle, pandas as pd, bilby
import matplotlib.pyplot as plt, seaborn as sns
from optparse import OptionParser
from scipy.integrate import simpson
from tqdm import tqdm
from multiprocessing import Process, Value, Lock, Manager
import json, re, time, datetime
import astropy
from _ctypes import PyObj_FromPtr  # see https://stackoverflow.com/a/15012814/355230

# Internal imports
sys.path.append('../')
import initialise, icarogw_runner as icarorun, snr_computation
import icarogw.simulation as icarosim
import icarogw.wrappers as icarowrap



######################
# Main function call #
######################

def main():

    # ------------------------------------------------ #
    # Read config file and initialise input parameters #
    # ------------------------------------------------ #

    parser = OptionParser(initialise.usage)
    parser.add_option(      '--config-file', type='string', metavar = 'config_file', default = None)
    parser.add_option('-n', '--n-processes', type='int',    metavar = 'n_processes', default = -1, help="Set the number of processes for parallelized injections generation from command line, if this should match some external structure (e.g. number of CPUs allocated to the simulation on a computing cluster job.)")
    (opts, _) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    # Initialise input parameters dictionary.
    input_pars = initialise.InitialiseOptions(Config)
    if opts.n_processes >= 0: input_pars['n-processes'] = opts.n_processes

    # Set output directory.
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])
    if not os.path.exists(os.path.join(input_pars['output'], 'plots')): os.makedirs(os.path.join(input_pars['output'], 'plots'))

    # Copy config file to output.
    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))

    # Deviate stdout and stderr to file.
    if not input_pars['screen-output']:
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_icarogw.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_icarogw.txt'), 'w')
    else: pass
    print('\n\n Starting  i c a r o g w  simulation\n\n')

    # Print run parameters.
    print(' * I will be running with the following parameters.\n')
    icarorun.print_dictionary(input_pars)

    # --------------------------------- #
    # Generate population or injections #
    # --------------------------------- #

    # Set additional input parameters for icarogw_runner.
    input_pars['zmax'] = input_pars['icarogw-sim-z-max']
    input_pars['ref-cosmology'] = {'H0': input_pars['truths']['H0'], 'Om0': input_pars['truths']['Om0']}

    # Initialise the model wrappers.
    tmp = icarorun.Wrappers(input_pars)
    m1w, m2w, rw, cw, ref_cosmo = tmp.return_Wrappers()
    input_pars['wrappers'] = {'m1w': m1w, 'm2w': m2w, 'rw': rw, 'cw': cw, 'ref-cosmo': ref_cosmo}

    # Initilialise the cosmology.
    cw.cosmology.build_cosmology(astropy.cosmology.FlatLambdaCDM(H0 = input_pars['truths']['H0'], Om0 = input_pars['truths']['Om0']))

    # Save and print true parameters.
    population_parameters = input_pars['wrappers']['m1w'].population_parameters + input_pars['wrappers']['rw'].population_parameters + input_pars['wrappers']['cw'].population_parameters
    if not input_pars['single-mass']: population_parameters += input_pars['wrappers']['m2w'].population_parameters
    print('\n * Using the following population parameters.\n')
    for key in population_parameters:
        if not key in input_pars['truths']: input_pars['truths'][key] = input_pars['all-truths'][key]
    print('\t{}'.format({key: input_pars['truths'][key] for key in population_parameters}))
    save_truths(input_pars['output'], {key: input_pars['truths'][key] for key in population_parameters})

    # Postprocess existing data.
    if input_pars['postprocessing']:
        try:
            samps_dict_astrophysical = pd.read_pickle(os.path.join(input_pars['output'], '{}_astrophysical.pickle').format(input_pars['run-type']))
            samps_dict_observed      = pd.read_pickle(os.path.join(input_pars['output'], '{}_observed.pickle'     ).format(input_pars['run-type']))
            print(    '\n * Reading existing {}.\n'.format(             input_pars['run-type']))
        except: print('\n * Existing {} not found. Exiting...\n'.format(input_pars['run-type']))

    # Generate new data.
    else:
        print('\n * Generating new {}.'.format(input_pars['run-type']))
        save_settings_pretty_json(input_pars['output'], input_pars)
        
        # Generate either a synthetic population or a set of injections for selection effects.
        if   input_pars['run-type'] == 'population':                            samps_dict_astrophysical, samps_dict_observed, strain_records = generate_population(input_pars)
        elif input_pars['run-type'] == 'injections' and input_pars['parallel']: samps_dict_astrophysical, samps_dict_observed = generate_injections_parallel(input_pars)
        elif input_pars['run-type'] == 'injections':                            samps_dict_astrophysical, samps_dict_observed = generate_injections(input_pars)
        else: raise ValueError('Invalid type of run. Please either select "population" or "injections".')

        with open(os.path.join(input_pars['output'], '{}_observed.pickle'     ).format(input_pars['run-type']), 'wb') as handle:
            pickle.dump(samps_dict_observed,      handle, protocol = pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(input_pars['output'], '{}_astrophysical.pickle').format(input_pars['run-type']), 'wb') as handle:
            pickle.dump(samps_dict_astrophysical, handle, protocol = pickle.HIGHEST_PROTOCOL)
        if input_pars['run-type'] == 'population' and input_pars['save-strain']:
            with open(os.path.join(input_pars['output'], 'strain_records.pickle').format(input_pars['run-type']), 'wb') as handle:
                pickle.dump(strain_records, handle, protocol = pickle.HIGHEST_PROTOCOL)

    # Plots section.
    plot_population(input_pars, samps_dict_astrophysical, samps_dict_observed)
    print(' * Finished.\n')




##########################
# Draw and detect events #
##########################

def generate_population(pars):
    '''
        Extract a set of samples from the specified probability distribution.

        The event values are sampled from the source frame distribution and
        then converted in the detector frame, assuming a reference cosmology.

        The information on the prior is not present here, and Jacobians for the
        source-detector conversion need to be handled in the running script.
    '''
    # Get distribution samples.
    m1s, m2s, zs, m1d, m2d, dL, _ = get_distribution_samples(pars)

    # Compute the SNR to select the detected events.
    SNR, idx_detected, _, additional_parameters, strain_records = compute_SNR(pars, m1s, m2s, zs, m1d, m2d, dL, save_strain=pars['save-strain'])
    clean_dict(additional_parameters, ['mass_1', 'mass_2', 'luminosity_distance', 'matched_filter_SNR'])
    
    # Save the number of detected events.
    print('\n * Number of detections: {}\n'.format(len(idx_detected)), flush = True)
    with open(os.path.join(pars['output'], 'number_detected_events.txt'), 'w') as f: f.write('{}'.format(len(idx_detected)))
    
    # Select the detected events.
    m1s_det = m1s[idx_detected]
    m2s_det = m2s[idx_detected]
    zs_det  = zs[ idx_detected]
    m1d_det = m1d[idx_detected]
    m2d_det = m2d[idx_detected]
    dL_det  = dL[ idx_detected]
    SNR_det = SNR[idx_detected]
    additional_parameters_det = {key: val[idx_detected] for key, val in additional_parameters.items()}
    
    # Save the population.
    samps_dict_observed      = {'m1s': m1s_det, 'm2s': m2s_det, 'z' : zs_det, 'm1d': m1d_det, 'm2d': m2d_det, 'dL': dL_det, 'snr': SNR_det, **additional_parameters_det}
    samps_dict_astrophysical = {'m1s': m1s,     'm2s': m2s,     'z' : zs,     'm1d': m1d,     'm2d': m2d,     'dL': dL,     'snr': SNR    , **additional_parameters    }
    
    return samps_dict_astrophysical, samps_dict_observed, strain_records


def generate_injections(pars):
    '''
        Generate injections for selection effects.

        The injection values are extracted in the source frame, but the hierarchical likelihood is expressed in the detector frame.
        We convert the injections from source to detector frame, and include this transformation in the prior weight.
    '''

    # Read checkpoint files if they exist.
    try:
        with open(os.path.join(pars['output'], 'injections_checkpoint_observed.pickle'     ), 'rb') as handle:
            samps_dict_observed      = pickle.load(handle)
        with open(os.path.join(pars['output'], 'injections_checkpoint_astrophysical.pickle'), 'rb') as handle:
            samps_dict_astrophysical = pickle.load(handle)
        inj_number_tmp = len(samps_dict_observed['m1d'])
        c = samps_dict_observed['n_batches_so_far']
        duration_cp = samps_dict_observed['run_time_so_far']
        print('\n * Reading existing injections from checkpoint files. Re-starting with {} generated and {} detected injections.\n'.format(len(samps_dict_astrophysical['m1d']), inj_number_tmp))
    except:
        # Initialize the dictionaries.
        samps_dict_observed      = {'m1s': [], 'm2s': [], 'z' : [], 'm1d': [], 'm2d': [], 'dL': [], 'snr': [], 'prior': []}
        samps_dict_astrophysical = {'m1s': [], 'm2s': [], 'z' : [], 'm1d': [], 'm2d': [], 'dL': [], 'snr': []}
        inj_number_tmp = 0
        c = 0
        duration_cp = 0.
    
    start_time = time.time() - duration_cp

    # Generate the injections.
    with tqdm(total = pars['injections-number'], initial = inj_number_tmp) as pbar:
        # Loop on banks of injections until the set number of detected injections is reached.
        while inj_number_tmp < pars['injections-number']:

            # Get distribution samples.
            if not pars['use-icarogw-sim-inj']:
                m1s, m2s, zs, m1d, m2d, dL, prior = get_distribution_samples(pars)

            else: # Use the icarogw.simulation implementation.
                # Draw the masses. Available options: PowerLaw, PowerLawPeak, MultiPeak.
                m1s, m2s, pdf_m = icarosim.generate_mass_inj(Nsamp = int(pars['injections-number-bank']), mass_model = pars['icarogw-sim-mass-model'], dic_param = pars['truths'])

                # Draw the luminosity distance.
                if   pars['icarogw-sim-draw-dL'] == 'uniform-dL':
                    dL, pdf_dL = icarosim.generate_dL_inj_uniform(  Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-z-max'])
                elif pars['icarogw-sim-draw-dL'] == 'uniform-z':
                    dL, pdf_dL = icarosim.generate_dL_inj_z_uniform(Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-z-max'])
                elif pars['icarogw-sim-draw-dL'] == 'uniform-volume':
                    dL, pdf_dL = icarosim.generate_dL_inj(          Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-z-max'])
                else:
                    raise ValueError('Unknown option for drawing the luminosity distance using icarogw.simulation. Exiting...')

                # Get the redshift.
                zs = icarosim.dl_to_z(dL)
                # Convert the masses to detector frame
                m1d = m1s * (1 + zs)
                m2d = m2s * (1 + zs)
                # Transform the prior from source to detector frame.
                # Only need to account for the masses (m1s, m2s) --> (m1d, m2d), as we draw directly from dL (i.e. no ddL/dz contribution).
                prior = (pdf_m * pdf_dL) / ((1 + zs)**2)

            # Compute the SNR to select the detected events.
            SNR, idx_detected, idx_softcut, additional_parameters, _ = compute_SNR(pars, m1s, m2s, zs, m1d, m2d, dL, save_strain=False)
            clean_dict(additional_parameters, ['mass_1', 'mass_2', 'luminosity_distance', 'matched_filter_SNR'])
            
            samps_dict_astrophysical['m1s'] = np.concatenate((samps_dict_astrophysical['m1s'], m1s[idx_softcut]))
            samps_dict_astrophysical['m2s'] = np.concatenate((samps_dict_astrophysical['m2s'], m2s[idx_softcut]))
            samps_dict_astrophysical['z'  ] = np.concatenate((samps_dict_astrophysical['z'  ], zs[ idx_softcut]))
            samps_dict_astrophysical['m1d'] = np.concatenate((samps_dict_astrophysical['m1d'], m1d[idx_softcut]))
            samps_dict_astrophysical['m2d'] = np.concatenate((samps_dict_astrophysical['m2d'], m2d[idx_softcut]))
            samps_dict_astrophysical['dL' ] = np.concatenate((samps_dict_astrophysical['dL' ], dL[ idx_softcut]))
            samps_dict_astrophysical['snr'] = np.concatenate((samps_dict_astrophysical['snr'], SNR[idx_softcut]))
            
            samps_dict_observed['m1s'  ] = np.concatenate((samps_dict_observed['m1s'  ], m1s[  idx_detected]))
            samps_dict_observed['m2s'  ] = np.concatenate((samps_dict_observed['m2s'  ], m2s[  idx_detected]))
            samps_dict_observed['z'    ] = np.concatenate((samps_dict_observed['z'    ], zs[   idx_detected]))
            samps_dict_observed['m1d'  ] = np.concatenate((samps_dict_observed['m1d'  ], m1d[  idx_detected]))
            samps_dict_observed['m2d'  ] = np.concatenate((samps_dict_observed['m2d'  ], m2d[  idx_detected]))
            samps_dict_observed['dL'   ] = np.concatenate((samps_dict_observed['dL'   ], dL[   idx_detected]))
            samps_dict_observed['snr'  ] = np.concatenate((samps_dict_observed['snr'  ], SNR[  idx_detected]))
            samps_dict_observed['prior'] = np.concatenate((samps_dict_observed['prior'], prior[idx_detected]))

            # Collect additional event parameters
            for key in additional_parameters:
                if key in samps_dict_astrophysical: 
                    samps_dict_astrophysical[key].append(additional_parameters[key][idx_softcut])
                    samps_dict_observed[     key].append(additional_parameters[key][idx_detected])
                else: 
                    samps_dict_astrophysical[key] = [additional_parameters[key][idx_softcut], ]
                    samps_dict_observed[     key] = [additional_parameters[key][idx_detected], ]


            number_detected_chunk = len(idx_detected)
            inj_number_tmp += number_detected_chunk
            c = c+1
            if inj_number_tmp > pars['injections-number']: pbar.update(pars['injections-number'])
            else:               pbar.update(number_detected_chunk)

            # Save injections to checkpoint files.
            if c % pars['inverse-checkpoint-rate'] == 0:
                samps_dict_observed['run_time_so_far'] = time.time() - start_time
                samps_dict_observed['n_batches_so_far'] = c
                with open(os.path.join(pars['output'], 'injections_checkpoint_observed.pickle'     ), 'wb') as handle:
                    pickle.dump(samps_dict_observed,      handle, protocol = pickle.HIGHEST_PROTOCOL)
                with open(os.path.join(pars['output'], 'injections_checkpoint_astrophysical.pickle'), 'wb') as handle:
                    pickle.dump(samps_dict_astrophysical, handle, protocol = pickle.HIGHEST_PROTOCOL)

    # Save the number of detected events.
    total_duration = time.time() - start_time
    print('\n * Generated {} injections (simulation duration : {}).'.format(pars['injections-number-bank'] * c, datetime.timedelta(seconds=total_duration)))
    with open(os.path.join(pars['output'], 'injections_number.txt'), 'w') as f:
        f.write('Generated: {}\n'.format(int(pars['injections-number-bank'] * c)))
        f.write('Detected:  {}'.format(len(samps_dict_observed['m1d'])))
    
    # Remove the checkpoint files.
    checkpoint_observed_path      = os.path.join(pars['output'], 'injections_parallel_checkpoint_observed.pickle'     )
    if os.path.exists(checkpoint_observed_path):      os.remove(checkpoint_observed_path)
    checkpoint_astrophysical_path = os.path.join(pars['output'], 'injections_parallel_checkpoint_astrophysical.pickle')
    if os.path.exists(checkpoint_astrophysical_path): os.remove(checkpoint_astrophysical_path)
    # Clean the result dict
    clean_dict(samps_dict_observed, ['run_time_so_far', 'n_batches_so_far'])

    return samps_dict_astrophysical, samps_dict_observed


def generate_injections_parallel(pars):
    '''
        Generate injections for selection effects. Launching multiple batches at once, within parallel processes.

        The injection values are extracted in the source frame, but the hierarchical likelihood is expressed in the detector frame.
        We convert the injections from source to detector frame, and include this transformation in the prior weight.
    '''
    # Read checkpoint files if they exist.
    samps_dict_observed, samps_dict_astrophysical, inj_number_cp, n_batches_cp, run_time_cp = load_checkpoint_parallel(pars)
    start_time = time.time() - run_time_cp

    with Manager() as manager:
        lock = Lock()  # Synchronization lock to avoid race conditions

        samps_dict_observed      = manager.dict(samps_dict_observed)
        samps_dict_astrophysical = manager.dict(samps_dict_astrophysical)
        inj_number = Value('i', inj_number_cp)
        n_batches  = Value('i', n_batches_cp)
        process_id = n_batches_cp

        processes = []  # List to manage active processes
        initial_counter_value = inj_number.value

        # Initialize the tqdm progress bar
        with tqdm(total=pars['injections-number'], desc="Injections generation progress", unit="event", initial=initial_counter_value) as pbar:

            while True:
                # Clean up completed processes
                processes = [p for p in processes if p.is_alive()]

                # Check the stopping condition
                if inj_number.value >= pars['injections-number']: 
                    break

                # Launch new processes if we haven't reached the max_process limit
                while len(processes) < pars['n-processes']:
                    process = Process(
                        target = worker_generate_injection_parallel, 
                        args   = (
                            lock, 
                            process_id, 
                            inj_number, 
                            n_batches,
                            samps_dict_observed, 
                            samps_dict_astrophysical, 
                            pars,
                        )
                    )
                    process.start()
                    processes.append(process)
                    process_id += 1  # Increment process ID
                    # Update the progress bar with the completed tasks
                    pbar.n = inj_number.value
                    pbar.refresh()

                # Throttle process creation to avoid overwhelming the system
                time.sleep(0.5)

                # Save checkpoint every k iterations
                if process_id % pars['inverse-checkpoint-rate'] == 0:
                    with lock: 
                        samps_dict_observed['run_time_so_far'] = time.time() - start_time
                        with open(os.path.join(pars['output'], 'injections_parallel_checkpoint_observed.pickle'), 'wb') as handle:
                            pickle.dump(dict(samps_dict_observed), handle, protocol = pickle.HIGHEST_PROTOCOL)
                        with open(os.path.join(pars['output'], 'injections_parallel_checkpoint_astrophysical.pickle'), 'wb') as handle:
                            pickle.dump(dict(samps_dict_astrophysical), handle, protocol = pickle.HIGHEST_PROTOCOL)

            # Wait for remaining processes to finish
            for process in processes:
                process.join()

            # Final update of the progress bar
            pbar.n = pars['injections-number']
            pbar.refresh()

        # Combine the results of all processes in a single dict
        samps_dict_observed_final, samps_dict_astrophysical_final = {}, {}
        available_pid = min([pid_ for pid_ in samps_dict_observed.keys() if isinstance(pid_, int)])
        for key in samps_dict_observed[available_pid]:

            samps_dict_observed_final[key] = np.concatenate(tuple(
                samps_dict_observed[_pid][key]
                for _pid in samps_dict_observed
                if isinstance(_pid, int)
            ))

            if key != 'prior':
                samps_dict_astrophysical_final[key] = np.concatenate(tuple(
                    samps_dict_astrophysical[_pid][key]
                    for _pid in samps_dict_astrophysical
                    if isinstance(_pid, int)
                ))

    # Save the number of detected events.
    total_duration = time.time() - start_time
    print(f"\n * Generated {pars['injections-number-bank'] * n_batches.value:.0f} injections (simulation duration : {datetime.timedelta(seconds=total_duration)}).")
    with open(os.path.join(pars['output'], 'injections_number.txt'), 'w') as f:
        f.write(f"Generated: {pars['injections-number-bank'] * n_batches.value:.0f}\n")
        f.write(f"Detected:  {inj_number.value:.0f}")

    # Remove the checkpoint files.
    checkpoint_observed_path      = os.path.join(pars['output'], 'injections_parallel_checkpoint_observed.pickle'     )
    if os.path.exists(checkpoint_observed_path):      os.remove(checkpoint_observed_path)
    checkpoint_astrophysical_path = os.path.join(pars['output'], 'injections_parallel_checkpoint_astrophysical.pickle')
    if os.path.exists(checkpoint_astrophysical_path): os.remove(checkpoint_astrophysical_path)
    # checkpoint_duration_path      = os.path.join(pars['output'], 'duration_checkpoint.pickle')
    # if os.path.exists(checkpoint_duration_path): os.remove(checkpoint_duration_path)

    return samps_dict_astrophysical_final, samps_dict_observed_final


def worker_generate_injection_parallel(lock, pid, inj_number, n_batches, samps_dict_observed, samps_dict_astrophysical, pars):
    """
    Individual injections batch worker for parallelized injections generation
    """
    # Get distribution samples.
    if not pars['use-icarogw-sim-inj']:
        np.random.seed(pid) # Otherwise all workers seem to start with the same seed hence generating the same (m1s, m2s, z) batches.
        m1s, m2s, zs, m1d, m2d, dL, prior = get_distribution_samples(pars)

    else: # Use the icarogw.simulation implementation.
        # Draw the masses. Available options: PowerLaw, PowerLawPeak, MultiPeak.
        m1s, m2s, pdf_m = icarosim.generate_mass_inj(Nsamp = int(pars['injections-number-bank']), mass_model = pars['icarogw-sim-mass-model'], dic_param = pars['truths'])

        # Draw the luminosity distance.
        if   pars['icarogw-sim-draw-dL'] == 'uniform-dL':
            dL, pdf_dL = icarosim.generate_dL_inj_uniform(  Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-z-max'])
        elif pars['icarogw-sim-draw-dL'] == 'uniform-z':
            dL, pdf_dL = icarosim.generate_dL_inj_z_uniform(Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-z-max'])
        elif pars['icarogw-sim-draw-dL'] == 'uniform-volume':
            dL, pdf_dL = icarosim.generate_dL_inj(          Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-z-max'])
        else:
            raise ValueError('Unknown option for drawing the luminosity distance using icarogw.simulation. Exiting...')

        # Get the redshift.
        zs = icarosim.dl_to_z(dL)
        # Convert the masses to detector frame
        m1d = m1s * (1 + zs)
        m2d = m2s * (1 + zs)
        # Transform the prior from source to detector frame.
        # Only need to account for the masses (m1s, m2s) --> (m1d, m2d), as we draw directly from dL (i.e. no ddL/dz contribution).
        prior = (pdf_m * pdf_dL) / ((1 + zs)**2)

    # Compute the SNR to select the detected events.
    SNR, idx_detected, idx_softcut, additional_parameters, _ = compute_SNR(pars, m1s, m2s, zs, m1d, m2d, dL, save_strain=False)
    clean_dict(additional_parameters, ['mass_1', 'mass_2', 'luminosity_distance', 'matched_filter_SNR'])
    n_det = len(idx_detected)

    # Build dictionaries for the current batch
    samps_batch_dict_astrophysical = {
        'm1s': m1s[idx_softcut], 
        'm2s': m2s[idx_softcut], 
        'z':   zs[ idx_softcut], 
        'm1d': m1d[idx_softcut], 
        'm2d': m2d[idx_softcut], 
        'dL':  dL[ idx_softcut], 
        'snr': SNR[idx_softcut], 
    }
    samps_batch_dict_observed = {
        'm1s':   m1s[  idx_detected], 
        'm2s':   m2s[  idx_detected], 
        'z':     zs[   idx_detected], 
        'm1d':   m1d[  idx_detected], 
        'm2d':   m2d[  idx_detected], 
        'dL':    dL[   idx_detected], 
        'snr':   SNR[  idx_detected], 
        'prior': prior[idx_detected], 
    }

    # Collect additional event parameters
    for key in additional_parameters:
        samps_batch_dict_astrophysical[key] = additional_parameters[key][idx_softcut]
        samps_batch_dict_observed[key]      = additional_parameters[key][idx_detected]

    with lock:

        if inj_number.value < pars['injections-number']: # Update the counter & the result dicts safely

            samps_dict_astrophysical[pid] = samps_batch_dict_astrophysical
            samps_dict_observed[pid] = samps_batch_dict_observed

            inj_number.value += n_det
            n_batches.value  += 1

            samps_dict_observed['n_det_so_far']     = inj_number.value
            samps_dict_observed['n_batches_so_far'] = n_batches.value

            print(f"Task {pid: 2d}: N_det this batch = {n_det:.0f}. Total : {inj_number.value: 9.0f}/{pars['injections-number']:.0f}")

        else:
            print(f"Task {pid}: already reached desired number of injections {pars['injections-number']}. No update.")


def load_checkpoint_parallel(pars):
    
    try:
        with open(os.path.join(pars['output'], 'injections_parallel_checkpoint_observed.pickle'     ), 'rb') as handle:
            samps_dict_observed          = pickle.load(handle)
        with open(os.path.join(pars['output'], 'injections_parallel_checkpoint_astrophysical.pickle'), 'rb') as handle: 
            samps_dict_astrophysical = pickle.load(handle)
        inj_number_cp = samps_dict_observed['n_det_so_far']
        n_batches_cp = samps_dict_observed['n_batches_so_far']
        run_time_cp = samps_dict_observed['run_time_so_far']
        print('\n * Reading existing injections from checkpoint files. Re-starting with {} generated and {} detected injections.\n'.format(n_batches_cp * pars['injections-number-bank'], inj_number_cp))
    except:
        # Initialize the dictionaries.
        samps_dict_observed      = {}
        samps_dict_astrophysical = {}
        inj_number_cp = 0
        n_batches_cp = 0
        run_time_cp = 0
        print('\n * No existing checkpoint file. Starting from zero simulated events.\n')

    return samps_dict_observed, samps_dict_astrophysical, inj_number_cp, n_batches_cp, run_time_cp



#######################
# Auxiliary functions #
#######################

def update_weights(w, val_dict):

    filt = {key: val_dict[key] for key in w.population_parameters}
    w.update(**filt)


def save_truths(path, dictionary):
    
    with open(os.path.join(path, 'true_parameters.txt'), 'w') as f:
        for key in dictionary.keys():
            max_len = len(max(dictionary.keys(), key = len))
            f.write('{}  {}\n'.format(key.ljust(max_len), dictionary[key]))


def read_truths(path):

    with open(os.path.join(path, 'true_parameters.txt'), 'r') as f:
        res = dict([line.strip().split('  ', 1) for line in f])
    res = {key: float(res[key]) for key in res.keys()}

    return res


def save_settings_pretty_json(path, dictionary):
    """Pretty JSON settings saving"""

    class NoIndent(object):
        """ Value wrapper."""
        def __init__(self, value):
            if not isinstance(value, (list, tuple, dict, np.ndarray)):
                raise TypeError('Only lists, tuples, dict, numpy.ndarray can be wrapped')
            self.value = value

    class MyEncoder(json.JSONEncoder):
        """
        Custom JSON encoder, only 1st level indented
        See https://stackoverflow.com/questions/42710879/write-two-dimensional-list-to-json-file
        """

        FORMAT_SPEC = '@@{}@@'  # Unique string pattern of NoIndent object ids.
        regex = re.compile(FORMAT_SPEC.format(r'(\d+)'))  # compile(r'@@(\d+)@@')

        def __init__(self, **kwargs):
            # Keyword arguments to ignore when encoding NoIndent wrapped values.
            ignore = {'cls', 'indent'}
            # Save copy of any keyword argument values needed for use here.
            self._kwargs = {k: v for k, v in kwargs.items() if k not in ignore}
            super(MyEncoder, self).__init__(**kwargs)

        def default(self, obj):
            return (self.FORMAT_SPEC.format(id(obj)) if isinstance(obj, NoIndent)
                        else super(MyEncoder, self).default(obj))
        
        def iterencode(self, obj, **kwargs):
            format_spec = self.FORMAT_SPEC  # Local var to expedite access.
            # Replace any marked-up NoIndent wrapped values in the JSON repr
            # with the json.dumps() of the corresponding wrapped Python object.
            for encoded in super(MyEncoder, self).iterencode(obj, **kwargs):
                match = self.regex.search(encoded)
                if match:
                    id = int(match.group(1))
                    no_indent = PyObj_FromPtr(id)
                    json_repr = json.dumps(no_indent.value, **self._kwargs)
                    # Replace the matched id string with json formatted representation
                    # of the corresponding Python object.
                    encoded = encoded.replace(
                                '"{}"'.format(format_spec.format(id)), json_repr)
                yield encoded

    dictionary_tosave = {}
    for key, value in dictionary.items():
        if   key != 'wrappers' and isinstance(value, (list, tuple, dict, np.ndarray)): dictionary_tosave[key] = NoIndent(value)
        elif key != 'wrappers':                                                        dictionary_tosave[key] = value
        else: pass # the wrapper entry in the input_pars 
    filename = os.path.join(path, 'analysis_settings.json')
    with open(filename, 'w') as f: json.dump(dictionary_tosave, f, indent=4, cls=MyEncoder)


def build_filter_subsample(N_evs, N_subset):

    filt = np.full(N_evs, False, dtype = bool)
    idxs = np.random.choice(filt.shape[0], N_subset, replace = False)
    for idx in idxs: filt[idx] = True
    return filt

def estimate_events_number(pars):

    if pars['observation-time'] < 0.:
        t_start, t_end = snr_computation.run_to_time_window(pars['observing-run'])
        pars['observation-time'] = (t_end - t_start) / (3600*24*365)

    print('\n * Estimating the number of astrophysical events from the rate, using R0 = {} [Gpc^-3 yr^-1] and Tobs = {} [yr].'.format(pars['R0'], pars['observation-time']))

    z_array = np.linspace(pars['bounds-z'][ 0], pars['bounds-z'][ 1], pars['N-points'])
    # Set the rate evolution.
    update_weights(pars['wrappers']['rw'], pars['truths'])
    # Project the rate on the light cone.
    tmp = pars['R0'] * pars['wrappers']['rw'].rate.evaluate(z_array) * pars['wrappers']['ref-cosmo'].dVc_by_dzdOmega_at_z(z_array) * 4*np.pi / (1+z_array)
    # Integrate in redshift and multiply by the observation time.
    events_number = round(simpson(tmp, z_array) * pars['observation-time'])

    print('\n * Drawing {} events from the population.'.format(events_number))
    return events_number

def clean_nans_in_pdf(pdf):

    if not np.all(np.isfinite(pdf)):
        pdf = np.nan_to_num(pdf, nan=0.0, posinf=0.0, neginf=0.0)
    return pdf

def get_distribution_samples(pars):
    '''
        Draw samples from the selected sitribution.
        Returns the source and detector frame samples.

        Whatever distribution is used, the output prior is expressed in the detector frame using the variables (m1d,m2d,dL).
        Please make sure to properly account for any variable transformation involved, by including the corresponding Jacobian.
    '''
    # Compute the number of astrophysical events, if required.
    if pars['estimate-events-number']: pars['events-number'] = estimate_events_number(pars)

    # Extract the number of events to generate, distinguishing between population and injections simulations.
    if   pars['run-type'] == 'population': N_events = pars['events-number'         ]
    elif pars['run-type'] == 'injections': N_events = pars['injections-number-bank']
    else: raise ValueError("Unknown run-type. Please choose between 'population' and 'injections'.")

    # Initialise the arrays.
    m1_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])
    m2_array = np.linspace(pars['bounds-m2'][0], pars['bounds-m2'][1], pars['N-points'])
    q_array  = np.linspace(pars['bounds-q' ][0], pars['bounds-q' ][1], pars['N-points'])
    z_array  = np.linspace(pars['bounds-z' ][0], pars['bounds-z' ][1], pars['N-points'])
    dL_array = np.linspace(pars['bounds-dL'][0], pars['bounds-dL'][1], pars['N-points'])

    # Set the sampler to draw events.
    if   pars['drawing-method'] == 'rejection-sampling'             : _sampler = rejection_sampling_1D
    elif pars['drawing-method'] == 'inverse-transform'              : _sampler = draw_samples_CDF_1D
    elif pars['drawing-method'] == 'deterministic-inverse-transform': _sampler = draw_stratified_samples_CDF_1D
    else: raise ValueError("Unknown drawing-method. Please choose between 'rejection-sampling', 'inverse-transform', 'deterministic-inverse-transform'.")

    # Set the sampler to draw events.
    if   pars['drawing-method'] == 'rejection-sampling'             : _sampler = rejection_sampling_1D
    elif pars['drawing-method'] == 'inverse-transform'              : _sampler = draw_samples_CDF_1D
    elif pars['drawing-method'] == 'deterministic-inverse-transform': _sampler = draw_stratified_samples_CDF_1D
    else: raise ValueError("Unknown drawing-method. Please choose between 'rejection-sampling', 'inverse-transform', 'deterministic-inverse-transform'.")

    # Rate evolution.
    update_weights(pars['wrappers']['rw'], pars['truths'])
    if not 'LuminosityProbability' in pars['model-rate']:
        tmp = pars['wrappers']['rw'].rate.evaluate(z_array)
        tmp = clean_nans_in_pdf(tmp)
        if not 'RedshiftProbability' in pars['model-rate']:
            tmp *= pars['wrappers']['ref-cosmo'].dVc_by_dzdOmega_at_z(z_array) * 4*np.pi / (1+z_array) # Convert from rate to probability distribution.
            zs, pdf_z = _sampler(z_array, tmp, N_events, 1)
            if pars['plot-astrophysical']: plot_injected_distribution(pars, z_array, pars['wrappers']['rw'], 'rate_evolution', rate_evolution = 1)
        else:
            zs, pdf_z = _sampler(z_array, tmp, N_events, 1)
            if pars['plot-astrophysical']: plot_injected_distribution(pars, z_array, pars['wrappers']['rw'], 'redshift_distribution', rate_evolution = 1, z_samps = zs)
        pdf_z_array = tmp
    else:
        tmp = pars['wrappers']['rw'].rate.evaluate(dL_array)
        tmp = clean_nans_in_pdf(tmp)
        dL, pdf_dL = _sampler(dL_array, tmp, N_events, 1)
        zs = pars['wrappers']['cw'].cosmology.dl2z(dL)
        if pars['plot-astrophysical']: plot_injected_distribution(pars, z_array, pars['wrappers']['rw'], 'redshift_distribution', rate_evolution = 1, z_samps = zs)
        pdf_z_array = tmp

    # Primary mass.
    update_weights(pars['wrappers']['m1w'], pars['truths'])
    if 'Redshift' in pars['model-primary']:
        m1s    = np.zeros(N_events)
        pdf_m1 = np.zeros(N_events)
        for i,z in tqdm(enumerate(zs),  total = len(zs)):
            tmp = pars['wrappers']['m1w'].pdf(m1_array, z)
            tmp = clean_nans_in_pdf(tmp)
            if not pars['drawing-method'] == 'deterministic-inverse-transform':
                # For redshift evolving distributions, we use the redshift samples to draw the masses.
                m1s[i], pdf_m1[i] = _sampler(m1_array, tmp,        1, seed = pars['seed'])
            else:
                m1s[i], pdf_m1[i] = _sampler(m1_array, tmp, N_events, seed = pars['seed'], quantile_index = i)
        if pars['plot-astrophysical']: plot_injected_distribution(pars, m1_array, pars['wrappers']['m1w'], 'm1z_redshift', redshift = True)
    else:
        tmp = pars['wrappers']['m1w'].pdf(m1_array)
        tmp = clean_nans_in_pdf(tmp)
        m1s, pdf_m1 = _sampler(m1_array, tmp, N_events, 2)
    pdf_m1_array = tmp

    # If required, remove the log10 contribution.
    if pars['log10-PDF']:
        m1s = np.power(10., m1s)
        pdf_m1 *= np.log10(np.e) / m1s   # Compute the Jacobian: |J_(log10(m1s))->(m1s)| = log10(e) / m1s.

    pdf_q_array = np.zeros(pars['N-points'])
    # Secondary mass.
    if not pars['single-mass']:
        if 'MassRatio' in pars['model-secondary']:
            update_weights(pars['wrappers']['m2w'], pars['truths'])
            tmp = pars['wrappers']['m2w'].pdf(q_array)
            tmp = clean_nans_in_pdf(tmp)
            qs, pdf_q = _sampler(q_array, tmp, N_events, 3)
            pdf_q_array = tmp

            # If required, remove the log10 contribution.
            if pars['log10-PDF']:
                qs = np.power(10., qs)
                pdf_q *= np.log10(np.e) / qs # Compute the Jacobian: |J_(log10(q))->(q)| = log10(e) / q.

            # If required, get m2 from q.
            if not pars['inverse-mass-ratio']:
                m2s = qs * m1s
                pdf_m2 = pdf_q / m1s         # Compute the Jacobian: |J_(m1s,q)->(m1s,m2s)| = 1/m1s, with q = m2/m1.
            else:
                m2s = m1s / qs
                pdf_m2 = pdf_q / m1s * qs**2 # Compute the Jacobian: |J_(m1s,q)->(m1s,m2s)| = q**2/m1s, with q=m1/m2.

            pdf_m1m2 = pdf_m1 * pdf_m2
            if pars['plot-astrophysical']: plot_injected_distribution(pars, q_array, pars['wrappers']['m2w'], 'q_source_frame', q_samps = qs)
        else:
            if 'Redshift' in pars['model-primary']:
                raise ValueError('The conditional secondary with redshift evolution in the primary is not implemented. Exiting...')
            pars['wrappers']['m2w'] = icarowrap.m1m2_conditioned_lowpass_m2(pars['wrappers']['m2w']) # Condition the secondary on the primary.
            update_weights(pars['wrappers']['m2w'], pars['truths'])
            m1s, m2s = pars['wrappers']['m2w'].prior.sample(N_events)
            pdf_m1m2 = pars['wrappers']['m2w'].prior.pdf(m1s, m2s)
    else:
        m2s = np.zeros(N_events)

    # Get detector frame quantities.
    m1d = m1s * (1 + zs)
    m2d = m2s * (1 + zs)
    if not 'LuminosityProbability' in pars['model-rate']:
        dL = pars['wrappers']['cw'].cosmology.z2dl(zs)

    if not pars['single-mass']:
        # Transform the prior from source to detector frame: |J_(m1s,m2s,z)->(m1d,m2d,dL)| = 1/ [(1+z)**2 * ddL/dz].
        prior = (pdf_m1m2 * pdf_z) / ((1 + zs)**2 * pars['wrappers']['ref-cosmo'].ddl_by_dz_at_z(zs))
    else:
        if not 'LuminosityProbability' in pars['model-rate']:
            # Transform the prior from source to detector frame: |J_(m1s,z)->(m1d,dL)| = 1/ [(1+z) * ddL/dz].
            prior = (pdf_m1 * pdf_z) / ((1 + zs)  * pars['wrappers']['ref-cosmo'].ddl_by_dz_at_z(zs))
        else:
            # Transform the prior from source to detector frame: |J_(m1s,dL)->(m1d,dL)| = 1/ (1+z).
            prior = (pdf_m1 * pdf_dL) / ((1 + zs))

    # Save injected population.
    data = np.column_stack((m1_array, q_array, z_array, pdf_m1_array, pdf_q_array, pdf_z_array))
    np.savetxt(os.path.join(pars['output'], 'true_population.txt'), data, delimiter="\t", header="m1s\tq\tz\tpdf_m1s\tpdf_q\tpdf_z")

    return m1s, m2s, zs, m1d, m2d, dL, prior


def rejection_sampling_1D(x, PDF, N, seed = None, quantile_index = None):
    '''
        Draw N samples from a distribution that follows the array PDF,
        using a rejection sampling algorithm.
    '''
    # Normalize PDF.
    PDF = PDF / np.sum(PDF)
    # Define proposal distribution (Uniform over range).
    x_min, x_max = np.min(x), np.max(x)
    M = np.max(PDF)  # Upper bound for rejection sampling.

    samples = []
    while len(samples) < N:
        # Sample uniformly from domain.
        x_s = np.random.uniform(x_min, x_max)
        # Compute probability density at x_s (interpolated).
        pdf_s = np.interp(x_s, x, PDF)
        # Accept or reject
        if np.random.uniform(0, M) < pdf_s: samples.append(x_s)

    # Convert to numpy array.
    samples = np.array(samples)
    # Compute corresponding probability densities.
    pdf_samples = np.interp(samples, x, PDF)

    return samples, pdf_samples


def draw_samples_CDF_1D(x, PDF, N, seed = None, quantile_index = None):
    '''
        Draw N samples from a distribution that follows the array PDF.

        Compute the cumulative distribution (CDF) and draw uniformly from [0,1].
        Interpolate the resulting samples to get a continuous extraction on x.
    '''
    # Normalize PDF to ensure it sums to 1.
    PDF = PDF / np.sum(PDF)
    # Compute CDF and ensure it starts at 0.
    CDF = np.cumsum(PDF)
    CDF = np.insert(CDF, 0, 0)  # Insert 0 at the beginning
    x_extended = np.insert(x, 0, x[0])  # Extend x for interpolation
    # Draw uniform samples from [0,1].
    tmp = np.random.uniform(0, 1, N)
    # Inverse transform sampling: interpolate using the CDF.
    samps = np.interp(tmp, CDF, x_extended)
    # Compute the derivative of the CDF (which is the PDF) using finite differences
    dCDF_dx = np.gradient(CDF, x_extended)  # Approximate derivative
    pdf_s = np.interp(samps, x_extended, dCDF_dx)  # Interpolate to get PDF values

    return samps, pdf_s


def draw_stratified_samples_CDF_1D(x, PDF, N, seed = 1, quantile_index = None):
    '''
    Generate N deterministic samples from a given 1D discrete PDF using stratified quantiles.
    
    Inputs:
        x    : array of x-values corresponding to the PDF
        PDF  : unnormalized PDF values over x
        N    : number of samples to generate
        seed : seed for random shuffling
        
    Returns:
        samps : array of smooth, noise-free samples
        pdf_s : interpolated PDF values at sample locations
    '''
    # Normalize the PDF to make it a proper probability distribution
    PDF = PDF / np.sum(PDF)
    # Compute the CDF and ensure it starts at 0
    CDF = np.cumsum(PDF)
    CDF = np.insert(CDF, 0, 0)
    x_extended = np.insert(x, 0, x[0])

    q = (np.arange(N) + 0.5) / N
    rng = np.random.default_rng(seed)
    q = q[rng.permutation(N)]
    if quantile_index is not None: q = q[quantile_index]

    # Inverse transform sampling via interpolation
    samps = np.interp(q, CDF, x_extended)
    # Compute the derivative (approximate PDF) using finite differences
    dCDF_dx = np.gradient(CDF, x_extended)
    pdf_s = np.interp(samps, x_extended, dCDF_dx)

    return samps, pdf_s



#######################
#      SNR utils      #
#######################

def compute_SNR(pars, m1s, m2s, zs, m1d, m2d, dL, save_strain=False):
    '''
    Compute the SNR with various methods for the input events

    Methods (in pars['SNR-method'])
    -------
    bilby: 
        Full matched filter SNR with bilby's Interferometer objects.
    pycbc:
        Optimal SNR with pycbc PSD and waveform objects.
        Optional: add unit centered gaussian fluctuation to mimic 
        matched filter SNR if pars['snr-pycbc-method']=='mf-fast'.
    proxy:
        Approximate SNR from inspiral leading order scaling with Mc and dL.
    flat-PSD:
        Approximation from redshift alone, assuming a flat PSD.
    lisabeta:
        For LISA analyses.

    Parameters
    ----------
    pars: dict
        At least 'SNR-method', and optional arguments.
    m1s: (m,) shape array-like
        source frame primary mass (used for proxy, flat-psd methods).
    m2s: (m,) shape array-like
        source frame secondary mass (used for proxy, flat-psd methods).
    zs:  (m,) shape array-like
        redshift (used for proxy, flat-psd methods).
    m1d: (m,) shape array-like
        detector frame primary mass (used for bilby, pycbc, lisabeta methods)
    m2d: (m,) shape array-like
        detector frame secondary mass (used for bilby, pycbc, lisabeta methods)
    dL:  (m,) shape array-like
        luminosity distance (for bilby, pycbc, lisabeta methods)
    save_strain:  bool
        Flag to save the strain (noise) data for each of the events.

    Returns
    -------
    SNR:                   (m,) shape array-like
        Signal-to-Noise ratio of input events
    idx_detected:          (m,) shape array-like
        indices of detected events
    additional_parameters: dict of (m,) shape array-like
        dictionary with arrays of additional parameters
    '''
    # Extract the number of events to generate, distinguishing between population and injections simulations
    if   pars['run-type'] == 'population': N_events = pars['events-number'         ]
    elif pars['run-type'] == 'injections': N_events = pars['injections-number-bank']
    else: raise ValueError("Unknown run-type. Please choose between 'population' and 'injections'.")
    strain_records = None

    # Use the full waveform to compute the SNR.
    if   pars['SNR-method'] == 'bilby':

        if pars['run-type'] == 'population': print('\n * Computing the MF SNR with the full waveform using bilby')
        bilby.core.utils.log.setup_logger(log_level=0)

        SNR = []
        additional_parameters = []
        if save_strain: strain_records = []

        bdp = snr_computation.BilbyDetectionPipeline(
            psd_dir             = pars['PSD-path'                     ],
            observing_run       = pars['observing-run'                ],
            reference_frequency = pars['snr-bilby-reference-frequency'],
            sampling_frequency  = pars['snr-bilby-sampling-frequency' ],
            approximant         = pars['snr-bilby-waveform'           ],
            precessing_apx      = pars['snr-bilby-precessing-wf'      ],
        )
 
        if   pars['run-type'] == 'population': iterator = tqdm(zip(m1d, m2d, dL), total=len(m1d), desc="MF SNR bilby", miniters=int(pars['events-number'])//200)
        elif pars['run-type'] == 'injections': iterator = zip(m1d, m2d, dL)

        for _m1, _m2, _dL in iterator:
            bdp.set_event_dict(
                {
                    'mass_1':              _m1,
                    'mass_2':              _m2,
                    'luminosity_distance': _dL,
                },
                set_duration_and_start_time=True
            )
            try: 
                bdp.set_ifos_list()
                bdp.set_waveform_generator()
                bdp.set_strain_data_from_psd()
                if save_strain:
                    strain_records.append({
                        'duration': bdp.duration,
                        'start_time': bdp.start_time,
                        'strain': [ifo.frequency_domain_strain for ifo in bdp.ifos_list],
                    })
                bdp.inject_signal()

            except RuntimeError as err: 
                print(f"\n{bdp.event_dict}\n")
                raise RuntimeError(err)

            event_dict = bdp.compute_matched_filter_SNR()

            SNR.append(event_dict['matched_filter_SNR'])
            additional_parameters.append(event_dict)
        
        SNR = np.array(SNR)
        idx_detected = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-cut'])
        idx_softcut  = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-soft-cut'])
        # additional_parameters is a list of single event dict, we convert it to a dict of arrays
        additional_parameters = pd.DataFrame(additional_parameters).to_dict(orient="list")
        for key in additional_parameters: additional_parameters[key] = np.array(additional_parameters[key])
        strain_records = np.array(strain_records)

    # Use the full waveform to compute the optimal SNR with bilby
    elif pars['SNR-method'] == 'pycbc':

        tmp_str = " approximate MF"*(pars['snr-pycbc-method']=='mf-fast') + " optimal"*(pars['snr-pycbc-method']=='opt')
        print(f'\n * Computing the{tmp_str} SNR with the full waveform using pycbc')

        SNR = np.zeros(N_events)

        try:
            detector_network = snr_computation.DetectorNetwork(
                observing_run = pars['observing-run'          ], 
                flow          = pars['snr-pycbc-f-low'        ], 
                delta_f       = pars['snr-pycbc-delta-f'      ],
                sample_rate   = pars['snr-pycbc-sampling-rate'],
                network       = pars['snr-pycbc-detectors'    ],
                psd_directory = pars['PSD-path'               ],
            )
        except AttributeError:
            raise ImportError("Please make sure PyCBC is properly installed if you wish to use it for SNR computation. (See https://pycbc.org/)")
        
        detector_network.load_psds()

        if   pars['run-type'] == 'population': iterator = tqdm(enumerate(zip(m1d, m2d, dL)), total=len(m1d), desc="Opt SNR pycbc")
        elif pars['run-type'] == 'injections': iterator = enumerate(zip(m1d, m2d, dL))

        for i, (_m1, _m2, _dL) in iterator:
            SNR[i] = detector_network.hit_network(
                m1=_m1, m2=_m2, dL=_dL,
                t_gps       = np.random.uniform(1240215503.0, 1240215503.0+3e7), # GW190425 trigtime. FIXME: Improve this.
                approximant = pars['snr-pycbc-waveform'  ],
                precessing  = pars['snr-pycbc-precession'],
                snr_method  = pars['snr-pycbc-method'    ],
            )

        idx_detected = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-cut'])
        idx_softcut  = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-soft-cut'])
        additional_parameters = {}


    # Use the quadrupole inspiral approximation to compute the SNR.
    elif pars['SNR-method'] == 'proxy':

        print(f'\n * Computing the SNR with the inspiral leading order approximation')

        theta        = icarosim.rvs_theta(N_events, 0., 1.4, pars['snr-proxy-theta-path']) # Average on extrinsic parameters.
        SNR, _, _    = icarosim.snr_samples(
            m1s, m2s, zs, theta = theta,
            numdet = pars['snr-proxy-N-detectors'  ],
            rho_s  = pars['snr-proxy-SNR-reference'],
            dL_s   = pars['snr-proxy-dL-reference' ],
            Md_s   = pars['snr-proxy-Mc-reference' ],
        )
        idx_detected = icarosim.snr_and_freq_cut(m1s, m2s, zs, SNR, snrthr = pars['SNR-cut'],      fgw_cut = pars['snr-proxy-fgw-cut'])
        idx_softcut  = icarosim.snr_and_freq_cut(m1s, m2s, zs, SNR, snrthr = pars['SNR-soft-cut'], fgw_cut = pars['snr-proxy-fgw-cut'])
        additional_parameters = {}


    # Use the flat PSD to compute the SNR.
    elif pars['SNR-method'] == 'flat-PSD':

        print('\n * Computing the SNR with the flat PSD.')
        SNR          = icarosim.snr_samples_flat(zs)
        idx_detected = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-cut'])
        idx_softcut  = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-soft-cut'])
        additional_parameters = {}

    # Use lisabeta to compute the SNR for LISA interferometer
    elif pars['SNR-method'] == 'lisabeta':

        print('\n * Computing the SNR with lisabeta')

        try:
            SNR = snr_computation.SNR_lisabeta(m1d, m1d/m2d, dL, Tobs = pars['observation-time'])
            idx_detected = snr_computation.cut_SNR(SNR)
            idx_softcut  = snr_computation.cut_SNR(SNR, snr_thr=pars['SNR-soft-cut'])
            additional_parameters = {}
        except AttributeError:
            raise ImportError("Please make sure lisabeta is properly installed if you wish to use it for SNR computation. (See https://pypi.org/project/lisabeta/)")

    else:
        raise ValueError('Unknown method to compute the SNR. Exiting...')

    return SNR, idx_detected, idx_softcut, additional_parameters, strain_records


def clean_dict(d, keys):
    """Remove dictionary entries"""
    for key in keys:
        if key in d: d.pop(key)



######################
# Plotting functions #
######################

def plot_population(pars, samps_dict_astrophysical, samps_dict_observed):

    colors = ['#0771AB', '#8F3A49']
    nbins  = 40
    alpha  = 0.5

    title = 'm1_source_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    m1_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])
    if not pars['log10-PDF']:
        if pars['plot-astrophysical']: plt.hist(samps_dict_astrophysical['m1s'], density = 1, bins = nbins, color = colors[0], alpha = alpha, label = 'Astrophysical')
        plt.hist(samps_dict_observed[     'm1s'], density = 1, bins = nbins, color = colors[1], alpha = alpha, label = 'Observed'     )
    else:
        if pars['plot-astrophysical']: plt.hist(np.log10(samps_dict_astrophysical['m1s']), density = 1, bins = nbins, color = colors[0], alpha = alpha, label = 'Astrophysical')
        plt.hist(np.log10(samps_dict_observed[     'm1s']), density = 1, bins = nbins, color = colors[1], alpha = alpha, label = 'Observed'     )
    update_weights(pars['wrappers']['m1w'], pars['truths'])
    if 'Redshift' in pars['model-primary']:
        if pars['plot-astrophysical']: plt.plot(m1_array, pars['wrappers']['m1w'].pdf(m1_array, np.median(samps_dict_astrophysical['z'])), c = '#153B60', label = 'Injected z-median')
    else:
        plt.plot(m1_array, pars['wrappers']['m1w'].pdf(m1_array), c = '#153B60', label = 'Injected')
    plt.title(title)
    plt.xlabel('m1')
    #if not pars['log10-PDF']:
        #plt.yscale('log')
        #plt.ylim(1e-5, 1)
    plt.legend()
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title = 'm1_detector_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    if not pars['log10-PDF']:
        if pars['plot-astrophysical']: plt.hist(samps_dict_astrophysical['m1d'], bins = nbins, color = colors[0], alpha = alpha, label = 'Astrophysical')
        plt.hist(                               samps_dict_observed[     'm1d'], bins = nbins, color = colors[1], alpha = alpha, label = 'Observed')
    else:
        if pars['plot-astrophysical']: plt.hist(np.log10(samps_dict_astrophysical['m1d']), bins = nbins, color = colors[0], alpha = alpha, label = 'Astrophysical')
        plt.hist(                               np.log10(samps_dict_observed[     'm1d']), bins = nbins, color = colors[1], alpha = alpha, label = 'Observed')
    plt.title(title)
    plt.xlabel('m1')
    plt.legend()
    plt.savefig('{}.pdf'.format(figname))
    plt.close()
    
    title = 'm1z_source_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    if not pars['log10-PDF']:
        if pars['plot-astrophysical']: plt.scatter(samps_dict_astrophysical['m1s'], samps_dict_astrophysical['z'], s = 0.1, c = colors[0], label = 'Astrophysical')
        plt.scatter(                               samps_dict_observed[     'm1s'], samps_dict_observed[     'z'], s = 4.0, c = colors[1], label = 'Observed', alpha = alpha)
    else:
        if pars['plot-astrophysical']: plt.scatter(np.log10(samps_dict_astrophysical['m1s']), samps_dict_astrophysical['z'], s = 0.1, c = colors[0], label = 'Astrophysical')
        plt.scatter(                               np.log10(samps_dict_observed[     'm1s']), samps_dict_observed[     'z'], s = 4.0, c = colors[1], label = 'Observed', alpha = alpha)
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('z')
    plt.legend()
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title = 'm1dL_detector_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    if not pars['log10-PDF']:
        if pars['plot-astrophysical']: plt.scatter(samps_dict_astrophysical['m1d'], samps_dict_astrophysical['dL'], s = 0.1, c = colors[0], label = 'Astrophysical')
        plt.scatter(                               samps_dict_observed[     'm1d'], samps_dict_observed[     'dL'], s = 4.0, c = colors[1], label = 'Observed', alpha = alpha)
    else:
        if pars['plot-astrophysical']: plt.scatter(np.log10(samps_dict_astrophysical['m1d']), samps_dict_astrophysical['dL'], s = 0.1, c = colors[0], label = 'Astrophysical')
        plt.scatter(                               np.log10(samps_dict_observed[     'm1d']), samps_dict_observed[     'dL'], s = 4.0, c = colors[1], label = 'Observed', alpha = alpha)
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('dL')
    plt.legend()
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title = 'z_distribution'
    figname = os.path.join(pars['output'], 'plots', title)
    if pars['plot-astrophysical']: plt.hist(samps_dict_astrophysical['z'], density = 1, bins = nbins, color = colors[0], alpha = alpha, label = 'Astrophysical')
    plt.hist(                               samps_dict_observed[     'z'], density = 1, bins = nbins, color = colors[1], alpha = alpha, label = 'Observed'     )
    plt.xlabel('$z$')
    plt.ylabel('$p(z)$')
    plt.legend()
    plt.tight_layout()
    plt.savefig('{}.pdf'.format(figname), transparent = True)
    plt.close()

    return 0

def plot_injected_distribution(pars, x_array, wrapper, title, redshift = False, rate_evolution = 0, q_samps = 0, z_samps = 0):

    if redshift:
        N_z = 10
        zx  = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
        zy        = np.linspace(zx[0], zx[-1], N_z)
        _, z_grid = np.meshgrid(zx, zy)

        _, ax  = plt.subplots(1, 2, figsize=(10, 5))
        colors = sns.color_palette('RdBu_r', N_z)

        for zi, z_array in enumerate(z_grid):
            pdf = wrapper.pdf(x_array, z_array)
            z = z_array[0]
            ax[0].plot(x_array, pdf+z, color = colors[zi])
            ax[1].plot(x_array, pdf,   color = colors[zi])

        ax[0].set_xlabel(r'$m_1\ [M_{\odot}]$')
        ax[1].set_xlabel(r'$m_1\ [M_{\odot}]$')
        ax[0].set_ylabel('$z$')
        ax[1].set_ylabel('$p(m_1)$')
        ax[1].set_xlim(0, 100)
        ax[1].set_ylim(1e-5, 1)
        ax[1].set_yscale('log')

        plt.subplots_adjust(wspace = 0.13)
        plt.tight_layout()
        figname = os.path.join(pars['output'], 'plots', title)
        plt.savefig('{}.pdf'.format(figname), transparent = True)
        plt.close()
    
    # FIXME: Should relax this condition.
    else:
        if not rate_evolution:
            figname = os.path.join(pars['output'], 'plots', title)
            pdf = wrapper.pdf(x_array)
            if not pars['log10-PDF']:
                plt.hist(q_samps, density = 1, bins = 40, color = '#0771AB', alpha = 0.5)
            else:
                plt.hist(np.log10(q_samps), density = 1, bins = 40, color = '#0771AB', alpha = 0.5)
            plt.plot(x_array, pdf, c = '#153B60', label = 'Injected')
            plt.xlabel('$q$')
            plt.ylabel('$p(q)$')
            plt.legend()
            plt.tight_layout()
            plt.savefig('{}.pdf'.format(figname), transparent = True)
            plt.close()

        else:
            if not pars['use-icarogw-sim-inj']:
                if not 'RedshiftProbability' in pars['model-rate']:
                    figname = os.path.join(pars['output'], 'plots', 'rate_evolution')
                    pdf = wrapper.rate.log_evaluate(x_array)
                    plt.plot(x_array, pdf, c = '#153B60', label = 'Injected')
                    plt.xlabel('$z$')
                    plt.ylabel('$R(z)/R0$')
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig('{}.pdf'.format(figname), transparent = True)
                    plt.close()
                else:
                    figname = os.path.join(pars['output'], 'plots', 'redshift_distribution')
                    pdf = wrapper.rate.evaluate(x_array)
                    plt.hist(z_samps, density = 1, bins = 40, color = '#0771AB', alpha = 0.5)
                    plt.plot(x_array, pdf, c = '#153B60', label = 'Injected')
                    plt.xlabel('$z$')
                    plt.ylabel('$p(z)$')
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig('{}.pdf'.format(figname), transparent = True)
                    plt.close()
    return 0



# Execute the main function.
if __name__=='__main__': main()
