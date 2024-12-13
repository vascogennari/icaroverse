import ast


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                      : 'icarogw_run',
        'screen-output'               : 0,

        # Wrappers
        'model-primary'               : 'PowerLaw-Gaussian',                     
        'model-secondary'             : 'MassRatio-Gaussian',
        'model-rate'                  : 'PowerLaw',

        'redshift-transition'         : 'linear',
        'positive-gaussian-z0'        : 0,
        'positive-gaussian-z'         : 0,
        'separate-gaussians-z0'       : 0,
        'separate-gaussians-z'        : 0,
        'redshift-mixture'            : 1,
        'low-smoothing'               : 0,
        'priors'                      : {},
        'scale-free'                  : 0,
        'single-mass'                 : 0,

        # Selection effects
        'injections-path'             : '',
        'injections-number'           : 1,
        'snr-cut'                     : 12.,
        'ifar-cut'                    : 4.,
        'selection-effects-cut'       : 'snr',

        # Data
        'O3-cosmology'                : 0,
        'simulation'                  : 1,
        'data-path'                   : '',
        'distance-prior-PE'           : 1,
        'remove-events'               : ['GW190412_053044'],
    
        # Likelihood
        'nparallel'                   : 1,
        'neffPE'                      : 1,
        'neffINJ'                     : None,

        # Sampler
        'sampler'                     : 'nessai',
        'nlive'                       : 200,
        'naccept'                     : 60,
        'npool'                       : 1,
        'print_method'                : 'interval-60',
        'sample'                      : 'acceptance-walk',
        'nwalkers'                    : 64,
        'nsteps'                      : 1000,
        'ntemps'                      : 10,
        'threads'                     : 10,

        # Plots
        'N-points'                    : 500,
        'N-z-slices'                  : 10,
        'N-z-slices-log'              : 5,
        'bounds-m1'                   : [0, 100],
        'bounds-m2'                   : [0, 100],
        'bounds-q'                    : [0, 1],
        'bounds-dL'                   : [0, 10000],
        'bounds-z'                    : [1e-5, 0.8],
        'true-values'                 : {},
        'selection-effects'           : 0,
        'plot-prior'                  : 1,
        'N-points-KDE'                : 500,
        'N-samps-prior'               : 1000,
    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if ('output' in key) or ('injections-path' in key) or ('selection-effects-cut' in key) or ('data-path' in key):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if ('injections-number' in key) or ('snr-cut' in key) or ('ifar-cut' in key):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if ('O3-cosmology' in key) or ('simulation' in key) or ('distance-prior-PE' in key):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass
        if ('remove-events' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass

        # Model
        if ('model-primary' in key) or ('model-secondary' in key) or ('model-rate' in key) or ('redshift-transition' in key):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if ('positive-gaussian-z0' in key) or ('positive-gaussian-z' in key) or ('separate-gaussians-z0' in key) or ('separate-gaussians-z' in key) or ('redshift-mixture' in key) or ('low-smoothing' in key) or ('scale-free' in key) or ('single-mass' in key):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if ('priors' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass

        # Sampler
        if ('sampler' in key):
            try: input_pars[key] = Config.get('sampler', key)
            except: pass
        if ('nparallel' in key) or ('neffPE' in key) or ('neffINJ' in key) or ('nlive' in key) or ('npool' in key) or ('nwalkers' in key) or ('nsteps' in key) or ('ntemps' in key):
            try: input_pars[key] = Config.getint('sampler', key)
            except: pass

        # Plots
        if ('N-points' in key) or ('N-z-slices' in key) or ('N-z-slices-log' in key) or ('N-points-KDE' in key) or ('N-samps-prior' in key):
            try: input_pars[key] = Config.getint('plots', key)
            except: pass
        if ('true-values' in key) or ('bounds-m1' in key) or ('bounds-m2' in key) or ('bounds-q' in key) or ('bounds-dL' in key) or ('bounds-z' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass
        if ('selection-effects' in key) or ('plot-prior' in key):
            try: input_pars[key] = Config.getboolean('plots', key)
            except: pass
    
    # Initialise the prior bounds.
    input_pars['all-priors'] = default_priors()
    if not input_pars['priors'] == {}:
        for key in input_pars['priors']: input_pars['all-priors'][key] = input_pars['priors'][key]
    
    return input_pars


def default_priors():
      
    prior = {
        # Cosmology
        'H0'          : 67.7,
        'Om0'         : 0.308,

        # Primary mass distribution
        'delta_m'     : [0.,  10. ],
        'delta_m_a'   : [0.,  30. ],
        'delta_m_b'   : [0.,  30. ],
        'delta_m_c'   : [0.,  30. ],

        'alpha'       : [-4., 120.],
        'alpha_z0'    : [-4., 120.],
        'alpha_z1'    : [-100.,100.],
        'mu_alpha'    : [0.,  100.],

        'alpha_a'     : [-4., 120.],
        'alpha_b'     : [-4.,  20.],
        'alpha_c'     : [-4.,  20.],
        'break_p'     : [0.,   1. ],

        'alpha_a_z0'  : [-4., 120.],
        'alpha_b_z0'  : [-4., 120.],
        'alpha_c_z0'  : [-4., 120.],
        'alpha_a_z1'  : [-100.,100.],
        'alpha_b_z1'  : [-100.,100.],
        'alpha_c_z1'  : [-100.,100.],
        'mmin_a_z0'   : [1. , 100.],
        'mmin_b_z0'   : [1. , 100.],
        'mmin_c_z0'   : [1. , 100.],
        'mmin_a_z1'   : [-100.,100.],
        'mmin_b_z1'   : [-100.,100.],
        'mmin_c_z1'   : [-100.,100.],
        'mmax_a_z0'   : [30., 200.],
        'mmax_b_z0'   : [30., 200.],
        'mmax_c_z0'   : [30., 200.],
        'mmax_a_z1'   : 0.,
        'mmax_b_z1'   : 0.,
        'mmax_c_z1'   : 0.,

        'mmin'        : [1. , 100.],
        'mmin_z0'     : [1. , 100.],
        'mmin_z1'     : [-100.,100.],
        'mmax'        : [30., 200.],
        'mmax_z0'     : [30., 200.],
        'mmax_z1'     : 0.,

        'mu_zt'       : [0. , 1.  ],
        'mu_delta_zt' : [1. , 100.],
        'sigma_zt'    : [0. , 1.  ],
        'sigma_delta_zt': [1. , 100.],

        'mmin_a'      : [1. , 100.],
        'mmin_b'      : [1. , 100.],
        'mmin_c'      : [1. , 100.],
        'mmax_a'      : [30., 200.],
        'mmax_b'      : [30., 200.],
        'mmax_c'      : [30., 200.],

        'mu_g'        : [20., 60. ],
        'mu_z0'       : [20., 60. ],
        'mu_z1'       : [-80., 80.],
        'mu_z2'       : [-80., 80.],
        'sigma_g'     : [1. , 30. ],
        'sigma_z0'    : [1. , 30. ],
        'sigma_z1'    : [0.,  20. ],
        'sigma_z2'    : [0.,  20. ],
        'mmin_g'      : [2. , 50. ],

        'mu_z0_a'     : [1.,  100.],
        'mu_z0_b'     : [1.,  100.],
        'mu_z0_c'     : [1.,  100.],
        'mu_z1_a'     : [-100.,100.],
        'mu_z1_b'     : [-100.,100.],
        'mu_z1_c'     : [-100.,100.],
        'sigma_z0_a'  : [1. , 50. ],
        'sigma_z0_b'  : [1. , 50. ],
        'sigma_z0_c'  : [1. , 50. ],
        'sigma_z1_a'  : [0.,  100.],
        'sigma_z1_b'  : [0.,  100.],
        'sigma_z1_c'  : [0.,  100.],

        'lambda_peak' : [0. , 1.  ],
        'mix_z0'      : [0. , 1.  ],
        'mix_z1'      : [0. , 1.  ],
        'mix_alpha_z0': [0. , 1.  ],
        'mix_alpha_z1': [0. , 1.  ],
        'mix_beta_z0' : [0. , 1.  ],
        'mix_beta_z1' : [0. , 1.  ],
        'mix'         : [0. , 1.  ],
        'mix_alpha'   : [0. , 1.  ],
        'mix_beta'    : [0. , 1.  ],

        'amp'         : [0. , 0.2 ],
        'freq'        : [-20., 20.],
        'zt'          : [0. , 1.  ],
        'delta_zt'    : [1. , 100.],

        # Secondary mass distribution
        'beta'        : [-20., 20.],
        'mu_q'        : [0.1, 1.  ],
        'sigma_q'     : [0.01, 0.9],
        'alpha_q'     : [-20., 20.],

        # Rate evolution
        'gamma'       : [-50., 30.],
        'kappa'       : [-20., 10.],
        'zp'          : [0. , 4.  ],
        'R0'          : [0. , 100.],
    }

    return prior


usage = """

\nWelcome to the ICAROGW runner helper!\n

    # ----- #
    # input #
    # ----- #

        output                      Path where the run output is saved. Default: 'icarogw_run'
        screen-output               Flag to deviate the standard output to screen. Default: '0'

        injections-path             Path of the injections used to evaluate selection effects. Default: ''
        injections-number           Number of generated injections used for selection effects. Default: 1
        selection-effects-cut       Type of cut to select events and injections. Options: 'snr', 'ifar'. Default: 'snr'
        snr-cut                     Value of signal-to-noise ratio used to selected the events and injections. Default: 12
        ifar-cut                    Value of inverse false-alarm-rate used to selected the events and injections. Default: 4

        data-path                   Path of the single event data. Default: ''.
        O3-cosmology                Option to process PE samples using O3 data from the GWTC-3 cosmology paper. Default: 0
        simulation                  Option to process PE samples using simulated events. Default: 1
        distance-prior-PE           Flag to re-weight the PE samples on the luminosity distance uniform in volume, accoringly to the GWTC-3 cosmology paper. Default: 1
        remove-events               List of events to be removed from the analysis. Default: ['GW190412_053044']

    # ----- #
    # model #
    # ----- #

        model-primary               Model distribution for the primary object. Options: 'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-Gaussian', 'PowerLaw-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftQuadratic', 'PowerLaw-GaussianRedshiftPowerLaw', 'PowerLaw-GaussianRedshiftSigmoid', 'PowerLawBroken-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear', 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_GaussianRedshiftLinear'. Default: 'PowerLaw-Gaussian'
        model-secondary             Model distribution for the secondary object. Options: 'Mass2-PowerLaw', 'MassRatio-PowerLaw', 'MassRatio-Gaussian'. Default: 'MassRatio-Gaussian'
        model-rate                  Model distribution for the rate evolution. Options: 'PowerLaw', 'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution'. Default: 'PowerLaw'
        
        low-smoothing               Flag to apply a smoothing function to the Powerlaw. Available in all primary mass models with a Powerlaw. Default: 0
        priors                      Dictionary of the prior bounds for the population parameters. Default values can be found in 'icarogw_pipeline.options.default_priors'
        scale-free                  Flag to use the scale-free likelihood fromulation. Default: 0
        single-mass                 Flag to use only one mass for the single-event parameters. Default: 0

        redshift-transition         Model function for the mixture redshift evolution. Only available for primary mass redshift evolving models. Options: 'linear', 'sigmoid', 'sinusoid'. Default: 'linear'
        redshift-mixture            Flag to make the mixture transition function stationary in redshift. Only available for primary mass redshift evolving models. Default: 1
        positive-gaussian-z0        Flag to constrain the Gaussian peak to be positive within 3 sigmas at redshift zero. Only available for primary mass redshift evolving models with a Gaussian. Default: 0
        positive-gaussian-z         Flag to constrain the Gaussian peak to be positive within 3 sigmas at all redshifts. Only available for primary mass redshift evolving models with a Gaussian. Default: 0
        separate-gaussians-z0       Flag to impose an ordering on the two Gaussian peaks at redshift zero.  Only available for primary mass redshift evolving models with two Gaussians. Default: 0
        separate-gaussians-z        Flag to impose an ordering on the two Gaussian peaks at all redshifts.  Only available for primary mass redshift evolving models with two Gaussians. Default: 0

    # ------- #
    # sampler #
    # ------- #

        nparallel                   Default: 1
        neffPE                      Number of effective PE samples per event contributing to the numerical evaluation of the likelihood. Default: 1
        neffINJ                     Number of effective injections contributing to the numerical evaluation of the likelihood. Default: None

        sampler                     Type of sampler to be used to draw samples from the likelihood. Options: 'dynesty', 'nessai', 'ptemcee'. Default: 'dynesty'
        nlive                       Number of live points used by the nested sampler. Option not available for the MCMC samplers. Default: 200
        npool                       Default: 1
        naccept                     Default: 60
        print_method                Default: 'interval-60'
        sample                      Default: 'acceptance-walk'
        nwalkers                    Default: 64
        nsteps                      Default: 1000
        ntemps                      Default: 10
        threads                     Default: 10

    # ----- #
    # plots #
    # ----- #

        N-points                    Default: 500,
        N-z-slices                  Default: 10,
        N-z-slices-log              Default: 5,
        bounds-m1                   Default: [0, 100],
        bounds-m2                   Default: [0, 100],
        bounds-q                    Default: [0, 1],
        bounds-dL                   Default: [0, 10000],
        bounds-z                    Default: [1e-5, 0.8],
        true-values                 Default: {},
        selection-effects           Default: 0,
        plot-prior                  Default: 1,
        N-points-KDE                Default: 500,
        N-samps-prior               Default: 1000,
"""