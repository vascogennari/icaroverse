import os, sys, shutil, configparser
import numpy as np, pickle, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from optparse import OptionParser
from tqdm import tqdm

# Internal imports
sys.path.append('../')
import initialise, icarogw_runner as icarorun
from snr_computation import DetectorNetwork
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
    parser.add_option('--config-file', type='string', metavar = 'config_file', default = None)
    (opts, _) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    # Initialise input parameters dictionary.
    input_pars = initialise.InitialiseOptions(Config)

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

    # Initialise the model wrappers.
    tmp = icarorun.Wrappers(input_pars)
    m1w, m2w, rw, cw, ref_cosmo = tmp.return_Wrappers()
    input_pars['wrappers'] = {'m1w': m1w, 'm2w': m2w, 'rw': rw, 'cw': cw, 'ref-cosmo': ref_cosmo}

    # Save and print true parameters.
    population_parameters = input_pars['wrappers']['m1w'].population_parameters + input_pars['wrappers']['m2w'].population_parameters + input_pars['wrappers']['rw'].population_parameters + input_pars['wrappers']['cw'].population_parameters
    print('\n * Using the following population parameters.\n')
    print('\t{}'.format({key: input_pars['truths'][key] for key in population_parameters}))
    save_truths(input_pars['output'], {key: input_pars['truths'][key] for key in population_parameters})

    # Postprocess existing data.
    if input_pars['postprocessing']:
        try:
            dict_source_frame = pd.read_pickle(os.path.join(input_pars['output'], '{}_source_frame.pickle'  ).format(input_pars['run-type']))
            dict_detect_frame = pd.read_pickle(os.path.join(input_pars['output'], '{}_detector_frame.pickle').format(input_pars['run-type']))
            print(    '\n * Reading existing {}.\n'.format(             input_pars['run-type']))
        except: print('\n * Existing {} not found. Exiting...\n'.format(input_pars['run-type']))

    # Generate new data.
    else:
        print('\n * Generating new {}.'.format(input_pars['run-type']))
        save_settings(input_pars['output'], input_pars)
        
        # Generate either a synthetic population or a set of injections for selection effects.
        if   input_pars['run-type'] == 'population': dict_source_frame, dict_detect_frame = generate_population(input_pars)
        elif input_pars['run-type'] == 'injections': dict_source_frame, dict_detect_frame = generate_injections(input_pars)
        else: raise ValueError('Invalid type of run. Please either select "population" or "injections".')

        with open(os.path.join(input_pars['output'], '{}_source_frame.pickle'  ).format(input_pars['run-type']), 'wb') as handle:
            pickle.dump(dict_source_frame,   handle, protocol = pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(input_pars['output'], '{}_detector_frame.pickle').format(input_pars['run-type']), 'wb') as handle:
            pickle.dump(dict_detect_frame, handle, protocol = pickle.HIGHEST_PROTOCOL)

    # Plots section.
    plot_population(input_pars, dict_source_frame, dict_detect_frame)
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
    # Use the full waveform to compute the SNR.
    if   pars['SNR-method'] == 'full-waveform':
        SNR, idx_detected = compute_SNR_full_waveform( pars, m1d, m2d, dL)
    # Use the approximate waveform to compute the SNR.
    elif pars['SNR-method'] == 'proxy-waveform':
        SNR, idx_detected = compute_SNR_proxy_waveform(pars, m1s, m2s, zs)
    # Use the flat PSD to compute the SNR.
    elif pars['SNR-method'] == 'flat-PSD':
        SNR, idx_detected = compute_SNR_flat_PSD(      pars, zs)
    else:
        raise ValueError('Unknows method to compute the SNR. Exiting...')

    # Save the number of detected events.
    print('\n * Number of detections: {}\n'.format(len(idx_detected)), flush = True)
    with open(os.path.join(pars['output'], 'number_detected_events.txt'), 'w') as f: f.write('{}'.format(len(idx_detected)))
    
    # Select the detected events.
    m1d_det = m1d[idx_detected]
    m2d_det = m2d[idx_detected]
    dL_det  = dL[ idx_detected]
    snr     = SNR[idx_detected]
    
    # Save the population.
    samps_detector_dict = {'m1d': m1d_det, 'm2d': m2d_det, 'dL': dL_det, 'snr': snr}
    samps_source_dict   = {'m1s': m1s,     'm2s': m2s,     'z' : zs,     'snr': snr}
    
    return samps_source_dict, samps_detector_dict


def generate_injections(pars):
    '''
        Generate injections for selection effects.

        The injection values are extracted in the source frame, but the hierarchical likelihood is expressed in the detector frame.
        We convert the injections from source to detector frame, and include this transformation in the prior weight.
    '''

    # Initialize the dictionaries.
    samps_detector_dict = {'m1d':[], 'm2d':[], 'dL':[], 'prior':[], 'snr':[]}
    samps_source_dict   = {'m1s':[], 'm2s':[], 'z' :[]}

    c = 0
    inj_number_tmp = 0
    # Generate the injections.
    with tqdm(total = pars['injections-number']) as pbar:
        # Loop on banks of injections until the set number of detected injections is reached.
        while inj_number_tmp < pars['injections-number']:

            # Get distribution samples.
            if not pars['use-icarogw-sim-inj']:
                m1s, m2s, zs, m1d, m2d, dL, prior = get_distribution_samples(pars)

            else: # Use the icarogw.simulation implementation.
                # Draw the masses. Available options: PowerLaw, PowerLawPeak, MultiPeak.
                m1s, m2s, pdf_m = icarosim.generate_mass_inj(Nsamp = int(pars['injections-number-bank']), mass_model = pars['icarogw-sim-dict']['mass-model'], dic_param = pars['truths'])

                # Draw the luminosity distance.
                if   pars['icarogw-sim-dict']['draw-dL'] == 'uniform-dL':
                    dL, pdf_dL = icarosim.generate_dL_inj_uniform(  Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-dict']['z-max'])
                elif pars['icarogw-sim-dict']['draw-dL'] == 'uniform-z':
                    dL, pdf_dL = icarosim.generate_dL_inj_z_uniform(Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-dict']['z-max'])
                elif pars['icarogw-sim-dict']['draw-dL'] == 'uniform-volume':
                    dL, pdf_dL = icarosim.generate_dL_inj(          Nsamp = int(pars['injections-number-bank']), zmax = pars['icarogw-sim-dict']['z-max'])
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
            # Use the full waveform to compute the SNR.
            if   pars['SNR-method'] == 'full-waveform':
                SNR, idx_detected = compute_SNR_full_waveform( pars, m1d, m2d, dL)
            # Use the approximate waveform to compute the SNR.
            elif pars['SNR-method'] == 'proxy-waveform':
                SNR, idx_detected = compute_SNR_proxy_waveform(pars, m1s, m2s, zs)
            # Use the flat PSD to compute the SNR.
            elif pars['SNR-method'] == 'flat-PSD':
                SNR, idx_detected = compute_SNR_flat_PSD(      pars, zs)
            else:
                raise ValueError('Unknows method to compute the SNR. Exiting...')
            
            number_detected_chuck = len(idx_detected)
            inj_number_tmp = inj_number_tmp + number_detected_chuck
            c = c+1

            samps_source_dict['m1s'].append(m1s)
            samps_source_dict['m2s'].append(m2s)
            samps_source_dict['z'  ].append(zs )

            samps_detector_dict['m1d'  ].append(m1d[  idx_detected])
            samps_detector_dict['m2d'  ].append(m2d[  idx_detected])
            samps_detector_dict['dL'   ].append(dL[   idx_detected])
            samps_detector_dict['prior'].append(prior[idx_detected])
            samps_detector_dict['snr'  ].append(SNR[  idx_detected])

            if inj_number_tmp > pars['injections-number']: pbar.update(pars['injections-number'])
            else:               pbar.update(number_detected_chuck)

    # Clean the dictionary.
    for key in samps_source_dict.keys():   samps_source_dict[  key] = np.concatenate(samps_source_dict[  key])
    for key in samps_detector_dict.keys(): samps_detector_dict[key] = np.concatenate(samps_detector_dict[key])

    # Save the number of detected events.
    print('\n * Generated {} injections.'.format(pars['injections-number-bank'] * c))
    with open(os.path.join(pars['output'], 'injections_number.txt'), 'w') as f:
        f.write('Generated: {}\n'.format(int(pars['injections-number-bank'] * c)))
        f.write('Detected:  {}'.format(len(samps_detector_dict['m1d'])))

    return samps_source_dict, samps_detector_dict



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

def save_settings(path, dictionary):

    with open(os.path.join(path, 'population_settings.txt'), 'w') as f:
        for key in dictionary.keys():
            max_len = len(max(dictionary.keys(), key = len))
            f.write('{}  {}\n'.format(key.ljust(max_len), dictionary[key]))

def read_settings(path):

    with open(os.path.join(path, 'population_settings.txt'), 'r') as f:
        res = dict([line.strip().split('  ', 1) for line in f])
    for key in res.keys():
        if key == 'positive-peak' or key == 'low-smoothing' or key == 'single-mass' or key == 'flat-PSD':
            res[key] = int(  res[key])
        if key == 'SNR-cut' or key == 'frequency-cut':
            res[key] = float(res[key])
        if key == 'model-primary' or key == 'model-secondary' or key == 'model-rate' or key == 'redshift-transition' or key == 'output':
            res[key] = res[key].replace(" ", "")

    return res

def build_filter_subsample(N_evs, N_subset):

    filt = np.full(N_evs, False, dtype = bool)
    idxs = np.random.choice(filt.shape[0], N_subset, replace = False)
    for idx in idxs: filt[idx] = True
    return filt

def get_distribution_samples(pars):
    '''
        Draw samples from the selected sitribution.
        Returns the source and detector frame samples.
    '''
    # Initialise the arrays.
    m1_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])
    m2_array = np.linspace(pars['bounds-m2'][0], pars['bounds-m2'][1], pars['N-points'])
    q_array  = np.linspace(pars['bounds-q'][0],  pars['bounds-q'][1],  pars['N-points'])
    z_array  = np.linspace(pars['bounds-z'][0],  pars['bounds-z'][1],  pars['N-points'])

    # Rate evolution.
    update_weights(pars['wrappers']['rw'], pars['truths'])
    # Convert from rate to probability distribution.
    tmp = pars['wrappers']['rw'].rate.evaluate(z_array) * pars['wrappers']['ref-cosmo'].dVc_by_dzdOmega_at_z(z_array) / (1+z_array)
    zs, pdf_z = rejection_sampling_1D(z_array, tmp, pars['events-number'])

    # Primary mass.
    update_weights(pars['wrappers']['m1w'], pars['truths'])
    if 'Redshift' in pars['model-primary']:
        m1s    = np.zeros(pars['events-number'])
        pdf_m1 = np.zeros(pars['events-number'])
        for i,z in tqdm(enumerate(zs),  total = len(zs)):
            tmp = pars['wrappers']['m1w'].pdf(m1_array, z)
            # For redshift evolving distributions, we use the redshift samples to draw the masses.
            m1s[i], pdf_m1[i] = rejection_sampling_1D(m1_array, tmp, 1)
        plot_injected_distribution(pars, m1_array, pars['wrappers']['m1w'], 'm1z_redshift', redshift = True)
    else:
        tmp = pars['wrappers']['m1w'].pdf(m1_array)
        m1s, pdf_m1 = rejection_sampling_1D(m1_array, tmp, pars['events-number'])

    # Secondary mass.
    if 'MassRatio' in pars['model-secondary']:
        update_weights(pars['wrappers']['m2w'], pars['truths'])
        tmp = pars['wrappers']['m2w'].pdf(q_array)
        qs, pdf_q = rejection_sampling_1D(q_array, tmp, pars['events-number'])
        m2s = qs * m1s # Compute m2 from q.
        pdf_m2 = pdf_q / m1s # Compute the Jacobian (m1s,q)->(m1s,m2s).
        pdf_m1m2 = pdf_m1 * pdf_m2
        plot_injected_distribution(pars, q_array, pars['wrappers']['m2w'], 'q_source_frame', q_samps = qs)
    else:
        if 'Redshift' in pars['model-primary']:
            raise ValueError('The conditional secondary with redshift evolution in the primary is not implemented. Exiting...')
        pars['wrappers']['m2w'] = icarowrap.m1m2_conditioned_lowpass_m2(pars['wrappers']['m2w']) # Condition the secondary on the primary.
        update_weights(pars['wrappers']['m2w'], pars['truths'])
        m1s, m2s = pars['wrappers']['m2w'].prior.sample(pars['events-number'])
        pdf_m1m2 = pars['wrappers']['m2w'].prior.pdf(m1s, m2s)

    # Get detector frame quantities.
    m1d = m1s * (1 + zs)
    m2d = m2s * (1 + zs)
    dL  = pars['wrappers']['ref-cosmo'].z2dl(zs)

    # Transform the prior from source to detector frame.
    prior = (pdf_m1m2 * pdf_z) / ((1 + zs)**2 * pars['wrappers']['ref-cosmo'].ddl_by_dz_at_z(zs))

    return m1s, m2s, zs, m1d, m2d, dL, prior

def compute_SNR_full_waveform(pars, m1d, m2d, dL):
    '''
        Compute the SNR with the full waveform.
        Returns the SNR and the indices of the detected events.
    '''
    print('\n * Computing the SNR with the full waveform.')
    SNR = np.zeros(pars['events-number'])
    detector_network = DetectorNetwork(
        observing_run = pars['snr-fw-observing-run'], 
        flow          = pars['snr-fw-f-low'        ], 
        delta_f       = pars['snr-fw-delta-f'      ],
        sample_rate   = pars['snr-fw-sampling-rate'],
        network       = pars['snr-fw-detectors'    ],
        psd_directory = pars['snr-fw-PSD-path'     ],
    )
    detector_network.load_psds()
    for i, (_m1, _m2, _dL) in tqdm(enumerate(zip(m1d, m2d, dL)), total=len(m1d)):
        SNR[i] = detector_network.hit_network(
            m1=_m1, m2=_m2, dL=_dL,
            t_gps       = np.random.uniform(1240215503.0, 1240215503.0+3e7), # GW190425 trigtime. FIXME: Improve this.
            approximant = pars['snr-fw-waveform'  ],
            precessing  = pars['snr-fw-precession'],
            snr_method  = pars['snr-fw-method'    ],
        )
    idx_detected = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-cut'])
    return SNR, idx_detected

def compute_SNR_proxy_waveform(pars, m1s, m2s, zs):
    '''
        Compute the SNR with the approximated waveform.
        Returns the SNR and the indices of the detected events.
    '''
    print('\n * Computing the SNR with the approximate waveform.')
    theta        = icarosim.rvs_theta(pars['events-number'], 0., 1.4, pars['snr-ap-theta-path']) # Average on extrinsic parameters.
    SNR, _, _    = icarosim.snr_samples(
        m1s, m2s, zs, theta = theta,
        numdet = pars['snr-ap-N-detectors'  ],
        rho_s  = pars['snr-ap-SNR-reference'],
        dL_s   = pars['snr-ap-dL-reference' ],
        Md_s   = pars['snr-ap-Mc-reference' ],
    )
    idx_detected = icarosim.snr_and_freq_cut(m1s, m2s, zs, SNR, snrthr = pars['SNR-cut'], fgw_cut = pars['frequency-cut'])
    return SNR, idx_detected

def compute_SNR_flat_PSD(pars, zs):
    '''
        Compute the SNR with the flat PSD.
        Returns the SNR and the indices of the detected events.
    '''
    print('\n * Computing the SNR with the flat PSD.')
    SNR          = icarosim.snr_samples_flat(zs)
    idx_detected = icarosim.snr_cut_flat(SNR, snrthr = pars['SNR-cut'])
    return SNR, idx_detected

def draw_samples_CDF_1D(x, PDF, N):
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

def rejection_sampling_1D(x, PDF, N):
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



######################
# Plotting functions #
######################

def plot_population(pars, dict_source_frame, dict_detect_frame):

    title = 'm1_source_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    m1_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])
    plt.hist(dict_source_frame['m1s'], density = 1, bins = 40, color = 'k', alpha = 0.5)
    update_weights(pars['wrappers']['m1w'], pars['truths'])
    plt.plot(m1_array, pars['wrappers']['m1w'].pdf(m1_array))
    plt.title(title)
    plt.xlabel('m1')
    plt.yscale('log')
    plt.ylim(1e-5, 1)
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title = 'm1_detector_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    plt.hist(dict_detect_frame['m1d'], bins = 40, color = 'k', alpha = 0.5)
    plt.title(title)
    plt.xlabel('m1')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()
    
    title = 'm1z_source_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    plt.scatter(dict_source_frame['m1s'], dict_source_frame['z'], s = 0.1, c = 'k')
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('z')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title = 'm1dL_detector_frame'
    figname = os.path.join(pars['output'], 'plots', title)
    plt.scatter(dict_detect_frame['m1d'], dict_detect_frame['dL'], s = 0.1, c = 'k')
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('dL')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    return 0

def plot_injected_distribution(pars, x_array, wrapper, title, redshift = False, q_samps = 0):

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

        ax[0].set_xlabel('$m_1\ [M_{\odot}]$')
        ax[1].set_xlabel('$m_1\ [M_{\odot}]$')
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
        figname = os.path.join(pars['output'], 'plots', title)
        pdf = wrapper.pdf(x_array)
        plt.hist(q_samps, density = 1, bins = 40, color = 'k', alpha = 0.5)
        plt.plot(x_array, pdf)
        plt.xlabel('$q$')
        plt.ylabel('$p(q)$')
        plt.tight_layout()
        plt.savefig('{}.pdf'.format(figname), transparent = True)
        plt.close()

    return 0



# Execute the main function.
if __name__=='__main__': main()
