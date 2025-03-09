import ast


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                      : 'icarogw_population',
        'run-type'                    : 'population',
        'screen-output'               : 0,
        'postprocessing'              : 0,

        'events-number'               : 1000,
        'injections-number'           : 1000,
        'injections-number-bank'      : 100,
        'selection-effects-cut'       : 'snr',
        'SNR-cut'                     : 12.,

        'SNR-method'                  : 'full-waveform',
        'frequency-cut'               : 15,
        'flat-PSD'                    : 0,

        'snr-fw-detectors'            : ['H1', 'L1', 'V1'],
        'snr-fw-observing-run'        : 'O3',
        'snr-fw-sampling-rate'        : 1024.,
        'snr-fw-delta-f'              : 1./128.,
        'snr-fw-f-low'                : 16.,
        'snr-fw-waveform'             : 'IMRPhenomXHM',
        'snr-fw-method'               : 'mf_fast',
        'snr-fw-precession'           : 0,
        'snr-fw-PSD-path'             : '',

        'snr-app-N-detectors'         : 2,
        'snr-app-SNR-reference'       : 9.,
        'snr-app-dL-reference'        : 1.5,
        'snr-app-Mc-reference'        : 25.,
        'snr-app-theta-path'          : 'Pw_three.dat',
        
        'use-icarogw-sim-inj'         : 0,
        'icarogw-sim-mass-model'      : 'PowerLawPeak',
        'icarogw-sim-draw-dL'         : 'uniform-dL',
        'icarogw-sim-z-max'           : 1.5,

        # Model
        'model-primary'               : 'PowerLaw-Gaussian',                     
        'model-secondary'             : 'MassRatio-Gaussian',
        'model-rate'                  : 'PowerLaw',

        'redshift-transition'         : 'linear',
        'redshift-mixture'            : 1,
        'low-smoothing'               : 0,
        'single-mass'                 : 0,
        'truths'                      : {},

        # Plots
        'N-points'                    : 100000,
        'bounds-m1'                   : [0, 200],
        'bounds-m2'                   : [0, 200],
        'bounds-q'                    : [0.1, 1],
        'bounds-dL'                   : [0, 10000],
        'bounds-z'                    : [1e-6, 3.0],
    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if ('output' in key) or ('run-type' in key) or ('selection-effects-cut' in key) or ('SNR-method' in key) or ('snr-fw-observing-run' in key) or ('snr-fw-waveform' in key) or ('snr-fw-method' in key) or ('snr-fw-PSD-path' in key) or ('snr-app-theta-path' in key) or ('icarogw-sim-mass-model' in key) or ('icarogw-sim-draw-dL' in key):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if ('SNR-cut' in key) or ('frequency-cut' in key) or ('snr-fw-sampling-rate' in key) or ('snr-fw-delta-f' in key) or ('snr-fw-f-low' in key) or ('snr-app-SNR-reference' in key) or ('snr-app-dL-reference' in key) or ('snr-app-Mc-reference' in key) or ('icarogw-sim-z-max' in key):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if ('events-number' in key) or ('injections-number' in key) or ('injections-number-bank' in key) or ('snr-app-N-detectors' in key):
            try: input_pars[key] = Config.getint('input', key)
            except: pass
        if ('screen-output' in key) or ('postprocessing' in key) or ('flat-PSD' in key) or ('use-icarogw-sim-inj' in key) or ('snr-fw-precession' in key) or ('use-icarogw-sim-inj' in key):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass
        if ('snr-fw-detectors' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass

        # Model
        if ('model-primary' in key) or ('model-secondary' in key) or ('model-rate' in key) or ('redshift-transition' in key):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if ('redshift-mixture' in key) or ('low-smoothing' in key) or ('single-mass' in key):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if ('truths' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass

        # Plots
        if ('N-points' in key):
            try: input_pars[key] = Config.getint('plots', key)
            except: pass
        if ('bounds-m1' in key) or ('bounds-m2' in key) or ('bounds-q' in key) or ('bounds-dL' in key) or ('bounds-z' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass
    
    # Initialise the population true values.
    input_pars['all-truths'] = default_population()
    if not input_pars['truths'] == {}:
        for key in input_pars['truths']: input_pars['all-truths'][key] = input_pars['truths'][key]
    else:
        input_pars['truths'] = input_pars['all-truths']

    return input_pars


def default_population():

    pop = {
        # Cosmology
        'H0'            : 67.7,
        'Om0'           : 0.308,

        # Primary mass distribution
        'delta_m'       : 5.,
        'delta_m_a'     : 5.,
        'delta_m_b'     : 5.,
        'delta_m_c'     : 5.,

        'alpha'         : 4.,
        'alpha_z0'      : 4.,
        'alpha_z1'      : 0.,
        'mu_alpha'      : 4.,

        'alpha_a'       : 50.,
        'alpha_b'       : 10.,
        'alpha_c'       : 10.,
        'break_p'       : 0.5,

        'alpha_a_z0'    : 50.,
        'alpha_b_z0'    : 10.,
        'alpha_c_z0'    : 10.,
        'alpha_a_z1'    : 0.,
        'alpha_b_z1'    : 0.,
        'alpha_c_z1'    : 0.,
        'mmin_a_z0'     : 10.,
        'mmin_b_z0'     : 15.,
        'mmin_c_z0'     : 35.,
        'mmin_a_z1'     : 0.,
        'mmin_b_z1'     : 0.,
        'mmin_c_z1'     : 0.,
        'mmax_a_z0'     : 100.,
        'mmax_b_z0'     : 100.,
        'mmax_c_z0'     : 100.,
        'mmax_a_z1'     : 0.,
        'mmax_b_z1'     : 0.,
        'mmax_c_z1'     : 0.,

        'mmin'          : 8.,
        'mmin_z0'       : 8.,
        'mmin_z1'       : 0.,
        'mmax'          : 100.,
        'mmax_z0'       : 100.,
        'mmax_z1'       : 0.,

        'mu_zt'         : 0.5,
        'mu_delta_zt'   : 50.,
        'sigma_zt'      : 0.5,
        'sigma_delta_zt': 50.,

        'mmin_a'        : 10.,
        'mmin_b'        : 15.,
        'mmin_c'        : 35.,
        'mmax_a'        : 100.,
        'mmax_b'        : 100.,
        'mmax_c'        : 100.,

        'mu_g'          : 35.,
        'mu_z0'         : 35.,
        'mu_z1'         : 0.,
        'mu_z2'         : 0.,
        'sigma_g'       : 5.,
        'sigma_z0'      : 5.,
        'sigma_z1'      : 0.,
        'sigma_z2'      : 0.,
        'mmin_g'        : 8.,

        'mu_z0_a'       : 10.,
        'mu_z0_b'       : 15.,
        'mu_z0_c'       : 35.,
        'mu_z1_a'       : 0.,
        'mu_z1_b'       : 0.,
        'mu_z1_c'       : 0.,
        'sigma_z0_a'    : 5.,
        'sigma_z0_b'    : 5.,
        'sigma_z0_c'    : 5.,
        'sigma_z1_a'    : 0.,
        'sigma_z1_b'    : 0.,
        'sigma_z1_c'    : 0.,

        'lambda_peak'   : 0.1,
        'mix_z0'        : 0.9,
        'mix_z1'        : 0.9,
        'mix_alpha_z0'  : 0.9,
        'mix_alpha_z1'  : 0.9,
        'mix_beta_z0'   : 0.05,
        'mix_beta_z1'   : 0.,
        'mix'           : 0.9,
        'mix_alpha'     : 0.9,
        'mix_beta'      : 0.05,

        'amp'           : 0.1,
        'freq'          : 0.,
        'zt'            : 0.5,
        'delta_zt'      : 50.,

        # Secondary mass distribution
        'beta'          : 4.,
        'mu_q'          : 0.8,
        'sigma_q'       : 0.1,

        # Rate evolution
        'gamma'         : 0.,
        'kappa'         : 0.,
        'zp'            : 3.,
        'R0'            : 20.,
    }

    return pop


# FIXME: Add the helper.
usage = """

\nWelcome to the ICAROGW runner helper!\n

    # ----- #
    # input #
    # ----- #

        output                      Directory from which input samples are read. The path is relative to the 'samples' namib directory. Default: ''
        
        injections-number           Flag to read the evidence from the samples. Default: 0
        snr-cut                     Flag to plot the prior. Option currently implemented only in corner plot and not fully implemented. Default: 0
        selection-effects-cut       Flag to read the evidence from the samples. Default: 0

        O3-cosmology                List of modes for which the QNMs are computed. This option is used only when QNMs are computed from {Mf, af}. Default: [(2,2)]
        simulation                  Flag to convert the damping time in [ms] and scale amplitudes as [1e-21]. The option is used to compare samples from Damped Sinusoids with other models. Default: 0

    # ----- #
    # model #
    # ----- #

        model-primary               Options: {'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear'}
        model-secondary             Options: {'Mass2-PowerLaw', 'MassRatio-Gaussian', 'MassRatio-PowerLaw'}
        model-rate                  Options: {'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution', 'PowerLaw'}
        redshift-transition         Options: {'linear', 'sigmoid', 'sinusoid'}

        positive-peak               List of parameters bounds that are used in the plots. Default: []. Syntax: [[a, b], [c, d], ...]
        low-smoothing               List of the selected parameters ordering. Default: []
        priors                      List of the compared parameters ordering. Default: []
        scale-free                  List of the compared parameters ordering. Default: []

"""