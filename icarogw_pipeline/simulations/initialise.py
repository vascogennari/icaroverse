import ast


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                        : 'icarogw_population',
        'run-type'                      : 'population',
        'screen-output'                 : 0,
        'postprocessing'                : 0,

        'events-number'               : 1000,
        'injections-number'           : 1000,
        'injections-number-bank'      : 100,
        'selection-effects-cut'       : 'snr',
        'SNR-cut'                     : 12.,
        'estimate-events-number'      : 0,
        'R0'                          : 25.,
        'observation-time'            : 1.,

        'SNR-method'                    : 'bilby',
        'PSD-path'                      : '',

        'snr-bilby-observing-run'       : 'O3',
        'snr-bilby-reference-frequency' : 20.,
        'snr-bilby-sampling-frequency'  : 2048.,
        'snr-bilby-waveform'            : 'IMRPhenomXHM',

        'snr-pycbc-detectors'           : ['H1', 'L1', 'V1'],
        'snr-pycbc-observing-run'       : 'O3',
        'snr-pycbc-sampling-rate'       : 1024.,
        'snr-pycbc-delta-f'             : 1./128.,
        'snr-pycbc-f-low'               : 16.,
        'snr-pycbc-waveform'            : 'IMRPhenomXHM',
        'snr-pycbc-method'              : 'mf_fast',
        'snr-pycbc-precession'          : 0,

        'snr-proxy-N-detectors'         : 2,
        'snr-proxy-fgw_cut'             : 15,
        'snr-proxy-SNR-reference'       : 9.,
        'snr-proxy-dL-reference'        : 1.5,
        'snr-proxy-Mc-reference'        : 25.,
        'snr-proxy-theta-path'          : 'Pw_three.dat',
        
        'use-icarogw-sim-inj'           : 0,
        'icarogw-sim-mass-model'        : 'PowerLawPeak',
        'icarogw-sim-draw-dL'           : 'uniform-dL',
        'icarogw-sim-z-max'             : 1.5,

        'log10-PDF'                     : 0,
        'inverse-mass-ratio'            : 1,

        # Model
        'model-primary'                 : 'PowerLaw-Gaussian',                     
        'model-secondary'               : 'MassRatio-Gaussian',
        'model-rate'                    : 'PowerLaw',

        'redshift-transition'           : 'linear',
        'redshift-mixture'              : 1,
        'low-smoothing'                 : 0,
        'single-mass'                   : 0,
        'truths'                        : {},

        # Plots
        'N-points'                      : 100000,
        'bounds-m1'                     : [0, 200],
        'bounds-m2'                     : [0, 200],
        'bounds-q'                      : [0.1, 1],
        'bounds-dL'                     : [0, 10000],
        'bounds-z'                      : [1e-6, 3.0],
    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if (key == 'output') or (key == 'run-type') or (key == 'selection-effects-cut') or (key == 'SNR-method') or (key == 'PSD-path') or (key == 'snr-bilby-observing-run') or (key == 'snr-bilby-waveform') or (key == 'snr-pycbc-observing-run') or (key == 'snr-pycbc-waveform') or (key == 'snr-pycbc-method') or (key == 'snr-proxy-theta-path') or (key == 'icarogw-sim-mass-model') or (key == 'icarogw-sim-draw-dL'):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if (key == 'SNR-cut') or (key == 'frequency-cut') or (key == 'snr-bilby-reference-frequency') or (key == 'snr-bilby-sampling-frequency') or (key == 'snr-pycbc-sampling-rate') or (key == 'snr-pycbc-delta-f') or (key == 'snr-pycbc-f-low') or (key == 'snr-proxy-SNR-reference') or (key == 'snr-proxy-dL-reference') or (key == 'snr-proxy-Mc-reference') or (key == 'icarogw-sim-z-max') or (key == 'R0') or (key == 'observation-time'):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if (key == 'events-number') or (key == 'injections-number') or (key == 'injections-number-bank') or (key == 'snr-proxy-N-detectors'):
            try: input_pars[key] = Config.getint('input', key)
            except: pass
        if (key == 'screen-output') or (key == 'postprocessing') or (key == 'flat-PSD') or (key == 'use-icarogw-sim-inj') or (key == 'snr-pycbc-precession') or (key == 'use-icarogw-sim-inj') or (key == 'log10-PDF') or (key == 'inverse-mass-ratio') or (key == 'estimate-events-number'):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass
        if (key == 'snr-pycbc-detectors'):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass

        # Model
        if (key == 'model-primary') or (key == 'model-secondary') or (key == 'model-rate') or (key == 'redshift-transition'):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if (key == 'redshift-mixture') or (key == 'low-smoothing') or (key == 'single-mass'):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if (key == 'truths'):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass

        # Plots
        if (key == 'N-points'):
            try: input_pars[key] = Config.getint('plots', key)
            except: pass
        if (key == 'bounds-m1') or (key == 'bounds-m2') or (key == 'bounds-q') or (key == 'bounds-dL') or (key == 'bounds-z'):
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

        # LISA only
        'm_b'           : 6.,
        'delta'         : 0.1,
        'a_gamma'       : 1.2,
        'theta'         : 0.3,
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