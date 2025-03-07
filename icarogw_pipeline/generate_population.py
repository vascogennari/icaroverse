import argparse, numpy as np, pickle, pandas as pd, os
import matplotlib.pyplot as plt, seaborn as sns
from tqdm import tqdm
from distutils.dir_util import copy_tree

from snr_computation import DetectorNetwork

import icarogw.simulation as icarosim
import icarogw_runner as icarorun


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
        if key == 'snr-cut' or key == 'fgw-cut':
            res[key] = float(res[key])
        if key == 'model-primary' or key == 'model-secondary' or key == 'model-rate' or key == 'redshift-transition' or key == 'output':
            res[key] = res[key].replace(" ", "")

    return res

def build_filter_subsample(N_evs, N_subset):

    filt = np.full(N_evs, False, dtype = bool)
    idxs = np.random.choice(filt.shape[0], N_subset, replace = False)
    for idx in idxs: filt[idx] = True
    return filt


def true_population_PDF_source(pars, truths, plot_dir, Ndetgen, return_wrappers = False):
    '''
        Extract a set of samples from the specified probability distribution.

        The event values are sampled from the source frame distribution and
        then converted in the detector frame, assuming a reference cosmology.

        The information on the prior is not present here, and Jacobians for the
        source-detector conversion need to be handled in the running script.
    '''

    N = 100000
    m_array = np.linspace(0.,   200., N)
    q_array = np.linspace(0.1,   1.,  N)
    z_array = np.linspace(1e-6, 1.5,  N)
    
    # Initialize ICAROGW wrappers
    # Primary, secondary, rate evolution, reference cosmology
    tmp = icarorun.Wrappers(pars)
    m1w, m2w, rw, _, cosmo_ref = tmp.return_Wrappers()
    inj_wrappers               = {'m1w': m1w, 'm2w': m2w, 'm_array': m_array, 'q_array': q_array, 'z_array': z_array}
    if return_wrappers: return inj_wrappers

    population_parameters = m1w.population_parameters + m2w.population_parameters + rw.population_parameters
    save_truths(pars['output'], {key: truths[key] for key in population_parameters})

    # Rate evolution
    update_weights(rw, truths)
    pdf = rw.rate.evaluate(z_array) * cosmo_ref.dVc_by_dzdOmega_at_z(z_array) / (1+z_array)    # Convert from rate to probability distribution.
    zs  = np.random.choice(z_array, size = Ndetgen, p = pdf/pdf.sum(), replace = True)

    # Primary mass
    update_weights(m1w, truths)
    m1s = np.zeros(Ndetgen)
    for i,z in tqdm(enumerate(zs),  total = len(zs)):
        if 'Redshift' in pars['model-primary']: pdf = m1w.pdf(m_array, z)
        else:                                   pdf = m1w.pdf(m_array)
        m1s[i] = np.random.choice(m_array, size = 1, p = pdf/pdf.sum(), replace = True)
    plot_injected_distribution(m_array, z_array, m1w, truths, plot_dir, pars['model-primary'], redshift = True)

    # Mass ratio
    update_weights(m2w, truths)
    pdf = m2w.pdf(q_array)
    qs  = np.random.choice(q_array, size = Ndetgen, p = pdf/pdf.sum(), replace = True)
    m2s = qs * m1s
    plot_injected_distribution(q_array, z_array, m2w, truths, plot_dir, pars['model-secondary'], q_samps = qs)

    # Detector frame
    m1d = m1s * (1 + zs)
    m2d = m2s * (1 + zs)
    dL  = cosmo_ref.z2dl(zs)

    if not pars['snr-approx']:
        rho_true_det = np.zeros_like(m1d)
        detector_network = DetectorNetwork(
            observing_run = 'O3', 
            flow          = 16., 
            delta_f       = 1./128.,
            sample_rate   = 1024.,
            network       = ['H1', 'L1', 'V1'],
            psd_directory = pars['psd_dir'],
        )
        detector_network.load_psds()
        print('\n * Computing the SNR with the full waveform.')
        for i, (_m1, _m2, _dL) in tqdm(enumerate(zip(m1d, m2d, dL)), total=len(m1d)):
            rho_true_det[i] = detector_network.hit_network(
                m1=_m1, m2=_m2, dL=_dL,
                t_gps       = np.random.uniform(1240215503.0, 1240215503.0+3e7), #GW190425
                approximant = 'IMRPhenomXHM',
                precessing  = False,
                snr_method  = 'mf_fast',
            )
        idx_cut_det = icarosim.snr_cut_flat(rho_true_det, snrthr = pars['snr-cut'])

    elif not pars['flat-PSD']:
        # Average on extrinsic parameters.
        theta              = icarosim.rvs_theta(Ndetgen, 0., 1.4, os.path.join('/Users/tbertheas/Documents/icarogw_pipeline/data/simulations', 'Pw_three.dat'))
        rho_true_det, _, _ = icarosim.snr_samples(     m1s, m2s, zs, numdet = 3, rho_s = 9, dL_s = 1.5, Md_s = 25, theta = theta)
        idx_cut_det        = icarosim.snr_and_freq_cut(m1s, m2s, zs, rho_true_det, snrthr = pars['snr-cut'], fgw_cut = pars['fgw-cut'])
    # Simulate a detection with a flat PSD.
    else:
        rho_true_det = icarosim.snr_samples_flat(zs)
        idx_cut_det  = icarosim.snr_cut_flat(rho_true_det, snrthr = pars['snr-cut'])

    print('\n * Number of detections: {}\n'.format(len(idx_cut_det)), flush = True)
    with open(os.path.join(results_dir, 'number_detected_events.txt'), 'w') as f: f.write('{}'.format(len(idx_cut_det)))
    
    m1d_det = m1d[idx_cut_det]
    m2d_det = m2d[idx_cut_det]
    dL_det  = dL[ idx_cut_det]
    snr = rho_true_det[idx_cut_det]
    
    samps_detector_dict = {'m1d': m1d_det, 'm2d': m2d_det, 'dL': dL_det, 'snr':snr}
    samps_source_dict   = {'m1s': m1s,     'm2s': m2s,     'z':  zs,     'snr':snr}
    
    return samps_source_dict, samps_detector_dict, inj_wrappers


def plot_population(source_dict, detector_dict, plot_dir, inj_wrappers, truths, model_primary):

    title   = 'm1_source_frame'
    figname = os.path.join(plot_dir, title)
    plt.hist(source_dict['m1s'], density = 1, bins = 40, color = 'k', alpha = 0.5)
    update_weights(inj_wrappers['m1w'], truths)
    if 'Redshift' in model_primary:
        plt.plot(inj_wrappers['m_array'], inj_wrappers['m1w'].pdf(inj_wrappers['m_array'], inj_wrappers['z_array'][0] ))
        plt.plot(inj_wrappers['m_array'], inj_wrappers['m1w'].pdf(inj_wrappers['m_array'], inj_wrappers['z_array'][-1]))
    else:
        plt.plot(inj_wrappers['m_array'], inj_wrappers['m1w'].pdf(inj_wrappers['m_array']))
        plt.plot(inj_wrappers['m_array'], inj_wrappers['m1w'].pdf(inj_wrappers['m_array']))
    plt.title(title)
    plt.xlabel('m1')
    plt.yscale('log')
    plt.ylim(1e-5, 1)
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title   = 'm1_detector_frame'
    figname = os.path.join(plot_dir, title)
    plt.hist(detector_dict['m1d'], bins = 40, color = 'k', alpha = 0.5)
    plt.title(title)
    plt.xlabel('m1')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()
    
    title   = 'm1z_source_frame'
    figname = os.path.join(plot_dir, title)
    plt.scatter(source_dict['m1s'], source_dict['z'], s = 0.1, c = 'k')
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('z')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title   = 'm1dL_detector_frame'
    figname = os.path.join(plot_dir, title)
    plt.scatter(detector_dict['m1d'], detector_dict['dL'], s = 0.1, c = 'k')
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('dL')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    return 0

def plot_injected_distribution(m_array, zx, mw, truths, plot_dir, wrap_name, redshift = False, q_samps = 0):

    if redshift:
        N_z = 10
        zy        = np.linspace(zx[0], zx[-1], N_z)
        _, z_grid = np.meshgrid(zx, zy)

        _, ax = plt.subplots(1, 2, figsize=(10, 5))
        colors = sns.color_palette('RdBu_r', N_z)

        for zi, z_array in enumerate(z_grid):
            if 'Redshift' in wrap_name: pdf = mw.pdf(m_array, z_array)
            else:                       pdf = mw.pdf(m_array)
            z = z_array[0]
            ax[0].plot(m_array, pdf + z, color = colors[zi])
            if not (zi == len(z_grid)-1):
                ax[1].plot(m_array, pdf, color = colors[zi])
            else:
                ax[1].plot(m_array, pdf, color = colors[zi], label='$\\alpha_0={}, \\alpha_1={}, mmin_0={}, mmin_1={}, \\mu_0={}, \\mu_1={}, mix z_0={}, mix z_1={}$'.format(
                                                                    truths['alpha_z0'], truths['alpha_z1'], truths['mmin_z0'], truths['mmin_z1'], truths['mu_z0'], truths['mu_z1'], truths['mix_z0'], truths['mix_z1']))
        ax[0].set_xlabel('$m_1\ [M_{\odot}]$')
        ax[1].set_xlabel('$m_1\ [M_{\odot}]$')
        ax[0].set_ylabel('$z$')
        ax[1].set_ylabel('$p(m_1)$')
        ax[1].set_xlim(0, 100)
        ax[1].set_ylim(1e-5, 1)
        ax[1].set_yscale('log')

        plt.subplots_adjust(wspace = 0.13)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'injected_m1.pdf'), transparent = True)
        plt.close()
        np.savetxt(os.path.join(plot_dir, 'injected_distribution.txt'), np.stack([m_array, pdf], axis=1), delimiter="\t", header="m_array\tpdf")
    
    else:
        # FIXME: Add option with m2 instead of q.
        title   = 'q_distribution'
        figname = os.path.join(plot_dir, title)
        pdf = mw.pdf(m_array)
        plt.hist(q_samps, density = 1, bins = 40, color = 'k', alpha = 0.5)
        plt.plot(m_array, pdf)
        plt.xlabel('$q$')
        plt.ylabel('$p(q)$')
        plt.tight_layout()
        plt.savefig('{}.pdf'.format(figname), transparent = True)
        plt.close()

    return 0


if __name__=='__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-t',   '--additional-text',     type = str,   metavar = 'additional_text',     default = ''                               )
    parser.add_argument('-g',   '--generated-events',    type = int,   metavar = 'generated_events',    default = 1000                             )
    parser.add_argument('-m1',  '--model-primary',       type = str,   metavar = 'model_primary',       default = 'PowerLawRedshiftLinear'         )
    parser.add_argument('-m2',  '--model-secondary',     type = str,   metavar = 'model_secondary',     default = 'MassRatio-Gaussian'             )
    parser.add_argument('-r',   '--model-rate',          type = str,   metavar = 'model_rate',          default = 'MadauDickinson'                 )
    parser.add_argument('-tr',  '--redshift-transition', type = str,   metavar = 'redshift_transition', default = 'linear'                         )
    parser.add_argument('-snr', '--snr-cut',             type = float, metavar = 'snr_cut',             default = 12                               )
    parser.add_argument('-fgw', '--fgw-cut',             type = float, metavar = 'fgw_cut',             default = 15                               )
    parser.add_argument('-apx', '--snr-approx',          type = int,   metavar = 'snr_approx',          default = 1                                )

    args                = parser.parse_args()
    additional_text     = args.additional_text
    N_events            = args.generated_events

    generate_population = 1
    subset_events       = 0

    print(args.snr_approx, type(args.snr_approx))

    input_pars  = {
        # Model parameters
        'model-primary'        : args.model_primary,
        'model-secondary'      : args.model_secondary,
        'model-rate'           : args.model_rate,

        'redshift-transition'  : args.redshift_transition,
        'positive-gaussian-z0' : 0,
        'positive-gaussian-z'  : 0,
        'separate-gaussians-z0': 0,
        'separate-gaussians-z' : 0,
        'redshift-mixture'     : 0,
        'low-smoothing'        : 1  ,
        'single-mass'          : 0,

        'snr-cut'              : args.snr_cut,
        'fgw-cut'              : args.fgw_cut,
        'snr-approx'           : args.snr_approx,
        'flat-PSD'             : 0,
    }

    true_values = {
        # Cosmology
        'H0'          : 67.7,
        'Om0'         : 0.308,

        # Primary mass distribution
        'delta_m'     : 3.0,

        'alpha'       : 50.,
        'mmin'        : 10.,
        'mmax'        : 100.,
        'mu'          : 0.,
        'sigma'       : 0.,

        'mu_g'        : 35.,
        'sigma_g'     : 2.,
        'lambda_peak' : 0.04,

        'alpha_z0'    : 50.,
        'alpha_z1'    : 20.,

        'alpha_a'     : 40.,
        'alpha_b'     : 18.,
        'alpha_c'     : 7.,
        'mmin_a'      : 10.,
        'mmin_b'      : 16.,
        'mmin_c'      : 28.,
        'mmax_a'      : 100.,
        'mmax_b'      : 100.,
        'mmax_c'      : 100.,

        'delta_m_a'   : 0.5,
        'delta_m_b'   : 10.,
        'delta_m_c'   : 8.,

        'mmin_z0'     : 10.,
        'mmin_z1'     : 0.,
        'mmax_z0'     : 100.,
        'mmax_z1'     : 0.,

        'mu_z0'       : 20.,
        'mu_z1'       : -80.,     # <----- Gaussian evolution
        'sigma_z0'    : 15.,
        'sigma_z1'    : 20.,

        'mix_z0'      : 0.8,
        'mix_z1'      : 0.6,
        
        # Secondary mass distribution
        'mu_q'        : 0.7,
        'sigma_q'     : 0.1,
        'alpha_q'     : 2.4,

        # Rate evolution
        'gamma'       : 3.,     # <----- Rate evolution
        'kappa'       : 0.,
        'zp'          : 3.,
        'R0'          : 0.,

        # Double peak
        'mu_z0_a'     : 25.,
        'mu_z1_a'     : 0.,
        'sigma_z0_a'  : 15.,
        'sigma_z1_a'  : 0.,
        'mu_z0_b'     : 35.,
        'mu_z1_b'     : 0.,
        'sigma_z0_b'  : 2.5,
        'sigma_z1_b'  : 0.,

        'mix_alpha_z0': 0.85,
        'mix_alpha_z1': 0.85,
        'mix_beta_z0' : 0.05,
        'mix_beta_z1' : 0.05,

        'mix_alpha'   : 0.8,
        'mix_beta'    : 0.1,
    }

    filename = 'pop-{}_{}_{}_{}{}'.format(int(N_events), input_pars['model-primary'], input_pars['model-secondary'], input_pars['model-rate'], additional_text)
    icarogw_pipeline_dir = '/sps/virgo/USERS/tbertheas/icarogw_pipeline'
    base_dir = os.path.join(icarogw_pipeline_dir, 'data/simulations')
    psd_dir  = os.path.join(icarogw_pipeline_dir, 'data/psd')
    results_dir = os.path.join(base_dir,    'simulated_population/H0_prospects', filename)
    plot_dir    = os.path.join(results_dir, 'population_plots')
    if not os.path.exists(results_dir): os.makedirs(results_dir)
    if not os.path.exists(plot_dir   ): os.makedirs(plot_dir)
    input_pars['output'] = results_dir
    input_pars['psd_dir'] = psd_dir

    if not generate_population:
        try:
            sample_source_dict_inj   = pd.read_pickle(os.path.join(results_dir, 'events_source_dict_{}.pickle'.format(  filename)))
            sample_detector_dict_inj = pd.read_pickle(os.path.join(results_dir, 'events_detector_dict_{}.pickle'.format(filename)))
            print('\n * Reading existing population.\n')
        except: print('\n * Existing population not found. Exiting...\n')

        if not subset_events == 0:
            N_events = len(sample_detector_dict_inj['m1d'])
            if subset_events > N_events: raise ValueError('The number of subset events need to be smaller than initial population.')
            print('\n * Sampling a subset of {} events out of the initial {}.\n'.format(subset_events, N_events))
            filter = build_filter_subsample(N_events, subset_events)
            sample_detector_dict_inj = {key: sample_detector_dict_inj[key][filter] for key in sample_detector_dict_inj.keys()}
            results_dir_new = results_dir + '_subset-{}'.format(subset_events)
            plot_dir        = os.path.join(results_dir_new, 'population_plots')
            input_pars['output'] = results_dir_new
            copy_tree(results_dir, results_dir_new)
            with open(os.path.join(input_pars['output'], 'events_detector_dict_{}_subset-{}.pickle'.format(filename, subset_events)), 'wb') as handle:
                pickle.dump(sample_detector_dict_inj, handle, protocol = pickle.HIGHEST_PROTOCOL)
            
            true_values  = read_truths(  results_dir_new)
            input_pars   = read_settings(results_dir_new)
            inj_wrappers = true_population_PDF_source(input_pars, true_values, plot_dir, Ndetgen = N_events, return_wrappers = True)

    else:
        print('\n * Generating new population.')
        save_settings(results_dir, input_pars)
        sample_source_dict_inj, sample_detector_dict_inj, inj_wrappers = true_population_PDF_source(input_pars, true_values, plot_dir, Ndetgen = N_events)
        with open(os.path.join(results_dir, 'events_source_dict_{}.pickle'.format(  filename)), 'wb') as handle:
            pickle.dump(sample_source_dict_inj,   handle, protocol = pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(results_dir, 'events_detector_dict_{}.pickle'.format(filename)), 'wb') as handle:
            pickle.dump(sample_detector_dict_inj, handle, protocol = pickle.HIGHEST_PROTOCOL)

    plot_population(sample_source_dict_inj, sample_detector_dict_inj, plot_dir, inj_wrappers, true_values, input_pars['model-primary'])
    print(' * Finished.\n')
