import numpy as np, pickle, pandas as pd, os
from astropy.cosmology import FlatLambdaCDM
import icarogw, icarogw.simulation as icarosim
import matplotlib.pyplot as plt, seaborn as sns
from tqdm import tqdm
import icarogw_runner as icarorun


def update_weights(w, val_dict):

    filt = {key: val_dict[key] for key in w.population_parameters}
    w.update(**filt)

def save_truths(path, dictionary):
    
    with open(os.path.join(path, 'true_parameters.txt'), 'w') as f:
        for key in dictionary.keys():
            max_len = len(max(dictionary.keys(), key = len))
            f.write('{}  {}\n'.format(key.ljust(max_len), dictionary[key]))


def true_population_PDF_source(pars, truths, plot_dir, Ndetgen):
    '''
        Extract a set of samples from the specified probability distribution.

        The event values are sampled from the source frame distribution and
        then converted in the detector frame, assuming a reference cosmology.

        The information on the prior is not present here, and Jacobians for the
        source-detector conversion need to be handled in the running script.
    '''

    N = 100000
    m_array = np.linspace(0.,   200., N)
    q_array = np.linspace(0.1,  1.,   N)
    z_array = np.linspace(1e-5, 2,    N)
    
    # Initialize ICAROGW wrappers
    # Cosmology
    cosmo_ref = icarogw.cosmology.astropycosmology(zmax = 20.)
    cosmo_ref.build_cosmology(FlatLambdaCDM(H0 = truths['H0'], Om0 = truths['Om0']))

    # Primary, secondary, rate evolution
    tmp = icarorun.Wrappers(pars)
    m1w, m2w, rw, _ = tmp.return_Wrappers()

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
        pdf = m1w.pdf(m_array, z)
        m1s[i] = np.random.choice(m_array, size = 1, p = pdf/pdf.sum(), replace = True)
    plot_injected_distribution(m_array, z_array, m1w, truths, plot_dir)

    # Mass ratio
    update_weights(m2w, truths)
    pdf = m2w.pdf(q_array)
    qs  = np.random.choice(q_array, size = Ndetgen, p = pdf/pdf.sum(), replace = True)
    m2s = qs * m1s

    theta      = icarosim.rvs_theta(Ndetgen, 0., 1.4, 'Pw_three.dat')
    rand_theta = np.random.choice(theta, Ndetgen)
    thetadet   = rand_theta
    
    rho_true_det, _, _ = icarosim.snr_samples(     m1s, m2s, zs, numdet = 3, rho_s = 9, dL_s = 1.5, Md_s = 25, theta = thetadet)
    idx_cut_det        = icarosim.snr_and_freq_cut(m1s, m2s, zs, rho_true_det, snrthr = pars['snr-cut'], fgw_cut = pars['fgw-cut'])
    print('\n * Number of detections: {}\n'.format(len(idx_cut_det)), flush = True)
    with open(os.path.join(results_dir, 'number_detected_events.txt'), 'w') as f: f.write('{}'.format(len(idx_cut_det)))
    
    # Detector frame
    m1d = m1s * (1 + zs)
    m2d = m2s * (1 + zs)
    dL  = cosmo_ref.z2dl(zs)
    
    m1d_det = m1d[idx_cut_det]
    m2d_det = m2d[idx_cut_det]
    dL_det  = dL[ idx_cut_det]
    
    sampe_detector_dict = { 'm1d': m1d_det, 'm2d': m2d_det, 'dL': dL_det}
    samps_source_dict   = { 'm1s': m1s,     'm2s': m2s,     'z':  zs}
    
    return samps_source_dict, sampe_detector_dict


def plot_population(source_dict, detector_dict, plot_dir):

    title   = 'm1_source_frame'
    figname = os.path.join(plot_dir, title)
    plt.hist(source_dict['m1s'], bins = 50, color = 'k', alpha = 0.5)
    plt.title(title)
    plt.xlabel('m1')
    plt.yscale('log')
    plt.savefig('{}.pdf'.format(figname))
    plt.close()

    title   = 'm1_detector_frame'
    figname = os.path.join(plot_dir, title)
    plt.hist(detector_dict['m1d'], bins = 50, color = 'k', alpha = 0.5)
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

def plot_injected_distribution(m_array, zx, mw, truths, plot_dir):

    N_z = 10
    zy        = np.linspace(zx[0], zx[-1], N_z)
    _, z_grid = np.meshgrid(zx, zy)

    _, ax = plt.subplots(1, 2, figsize=(10, 5))
    colors = sns.color_palette('RdBu_r', N_z)

    for zi, z_array in enumerate(z_grid):
        pdf = mw.pdf(m_array, z_array)
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

    return 0


# MAIN
generate_population = 1
additional_text     = '_evolving'
N_events            = 252000

input_pars  = {
    # Model parameters
    'model-primary'         : 'PowerLawRedshiftLinear-GaussianRedshiftLinear',                     
    'model-secondary'       : 'MassRatio',
    'model-rate'            : 'PowerLaw',

    'redshift-transition'   : 'linear',
    'positive-peak'         : 0,
    'low-smoothing'         : 1,
    'priors'                : {},

    'snr-cut'               : 12.,
    'fgw-cut'               : 15.,
}

true_values = {
    # Cosmology
    'H0'          : 67.7,
    'Om0'         : 0.308,

    # Primary mass distribution
    'delta_m'     : 5.,

    'alpha'       : 3.,
    'mmin'        : 8.,
    'mmax'        : 100.,
    'mu'          : 30.,
    'sigma'       : 5.,

    'alpha_z0'    : 3.8,
    'alpha_z1'    : 0.,

    'mmin_z0'     : 7.,
    'mmin_z1'     : 00.,
    'mmax_z0'     : 150.,
    'mmax_z1'     : 0.,

    'mu_z0'       : 30.,
    'mu_z1'       : 20.,
    'sigma_z0'    : 6.,
    'sigma_z1'    : 0.,

    'mix_z0'      : 0.9,
    'mix_z1'      : 0.9,
    
    # Secondary mass distribution
    'mu_q'        : 0.8,
    'sigma_q'     : 0.15,

    # Rate evolution
    'gamma'       : 0.,
    'kappa'       : 3.,
    'zp'          : 2.,
    'R0'          : 20.,
}


filename = 'pop-{}_{}_{}_{}{}'.format(int(N_events), input_pars['model-primary'], input_pars['model-secondary'], input_pars['model-rate'], additional_text)
base_dir = '/Users/vgennari/Documents/work/code/python/icarogw/scripts/injection_campaign'
results_dir = os.path.join(base_dir,    'simulated_population', filename)
plot_dir    = os.path.join(results_dir, 'population_plots')
if not os.path.exists(results_dir): os.makedirs(results_dir)
if not os.path.exists(plot_dir   ): os.makedirs(plot_dir)
input_pars['output'] = results_dir

if not generate_population:
    try:
        sample_source_dict_inj   = pd.read_pickle(os.path.join(results_dir, 'sample_source_dict_inj_up3_50m.pickle'  ))
        sample_detector_dict_inj = pd.read_pickle(os.path.join(results_dir, 'sample_detector_dict_inj_up3_50m.pickle'))
        print('\n * Reading existing population.\n')
    except: print('\nExisting population not found. Exiting...\n')

else:
    print('\n * Generating new population.\n')
    sample_source_dict_inj, sample_detector_dict_inj = true_population_PDF_source(input_pars, true_values, plot_dir, Ndetgen = N_events)     #50000000
    with open(os.path.join(results_dir, 'events_source_dict_{}.pickle'.format(filename)  ), 'wb') as handle:
        pickle.dump(sample_source_dict_inj,   handle, protocol = pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(results_dir, 'events_detector_dict_{}.pickle'.format(filename)), 'wb') as handle:
        pickle.dump(sample_detector_dict_inj, handle, protocol = pickle.HIGHEST_PROTOCOL)

plot_population(sample_source_dict_inj, sample_detector_dict_inj, plot_dir)
print(' * Finished.\n')
