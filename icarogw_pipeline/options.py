import ast


def default_priors():
      
    prior = {
        # Cosmology
        'H0'          : 67.7,
        'Om0'         : 0.308,

        # Primary mass distribution
        'delta_m'     : [0.1, 10. ],

        'alpha'       : [1.,  6.  ],
        'alpha_z0'    : [1.,  6.  ],
        'alpha_z1'    : [-5., 20. ],

        'mmin'        : [2. , 20. ],
        'mmin_z0'     : [2. , 20. ],
        'mmin_z1'     : [-5., 50. ],
        'mmax'        : [65., 200.],
        'mmax_z0'     : [65., 200.],
        'mmax_z1'     : 0.,

        'mu_g'        : [20., 60. ],
        'mu_z0'       : [20., 60. ],
        'mu_z1'       : [-80., 80.],
        'sigma_g'     : [5. , 30. ],
        'sigma_z0'    : [5. , 30. ],
        'sigma_z1'    : [0.,  20. ],
        'mix_z0'      : [0. , 1.  ],
        'mix_z1'      : [0. , 1.  ],

        'lambda_peak' : [0. , 1.  ],
        
        # Secondary mass distribution
        'beta'        : [1.,  6.  ],
        'mu_q'        : [0.4, 1.  ],
        'sigma_q'     : [0.01, 0.9],

        # Rate evolution
        'gamma'       : [-50., 30. ],
        'kappa'       : [-6. , 6.  ],
        'zp'          : [0. , 4.   ],
        'R0'          : [0. , 100. ],
    }

    return prior


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                : 'icarogw_run',
        'screen-output'         : 0,

        # Wrappers
        'model-primary'         : 'PowerLaw-Gaussian',                     
        'model-secondary'       : 'MassRatio',
        'model-rate'            : 'PowerLaw',

        'redshift-transition'   : '',
        'positive-peak'         : 0,
        'low-smoothing'         : 0,
        'priors'                : {},
        'conditional-prior'     : 0,
        'scale-free'            : 0,
        'single-mass'           : 0,

        # Selection effects
        'injections-path'       : '',
        'injections-number'     : 1,
        'snr-cut'               : 12.,
        'ifar-cut'              : 4.,
        'selection-effects-cut' : 'snr',

        # Data
        'O3-cosmology'          : 0,
        'simulation'            : 1,
        'data-path'             : '',
    
        # Likelihood
        'nparallel'             : 1,
        'neffPE'                : 1,
        'neffINJ'               : None,

        # Sampler
        'sampler'               : 'dynesty',
        'nlive'                 : 200,
        'naccept'               : 60,
        'npool'                 : 1,
        'print_method'          : 'interval-60',
        'sample'                : 'acceptance-walk',

        # Plots
        'N-points'              : 1000,
        'N-z-slices'            : 10,
        'bounds-m1'             : [1, 100],
        'bounds-m2'             : [1, 80],
        'bounds-q'              : [0, 1],
        'bounds-dL'             : [0, 10000],
        'bounds-z'              : [1e-5, 0.8],
        'true-values'           : {},
        'selection-effects'     : 0,
        'N-points-KDE'          : 1000,
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
        if ('O3-cosmology' in key) or ('simulation' in key):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass

        # Model
        if ('model-primary' in key) or ('model-secondary' in key) or ('model-rate' in key) or ('redshift-transition' in key):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if ('positive-peak' in key) or ('low-smoothing' in key) or ('scale-free' in key) or ('conditional-prior' in key) or ('single-mass' in key):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if ('priors' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass

        # Sampler
        if ('sampler' in key):
            try: input_pars[key] = Config.get('sampler', key)
            except: pass
        if ('nparallel' in key) or ('neffPE' in key) or ('nlive' in key) or ('npool' in key):
            try: input_pars[key] = Config.getint('sampler', key)
            except: pass

        # Plots
        if ('N-points' in key) or ('N-z-slices' in key) or ('N-points-KDE' in key):
            try: input_pars[key] = Config.getfloat('plots', key)
            except: pass
        if ('bounds-m1' in key) or ('bounds-m2' in key) or ('bounds-q' in key) or ('bounds-dL' in key) or ('bounds-z' in key):
            try: input_pars[key] = Config.get('plots', key)
            except: pass
        if ('true-values' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass
        if ('selection-effects' in key):
            try: input_pars[key] = Config.getboolean('plots', key)
            except: pass
    
    # Initialise the prior bounds.
    input_pars['all-priors'] = default_priors()
    if not input_pars['priors'] == {}:
        for key in input_pars['priors']: input_pars['all-priors'][key] = input_pars['priors'][key]
    
    return input_pars


usage = """

\nWelcome to the ICAROGW runner helper!\n

    # ----- #
    # input #
    # ----- #

        output                      Directory from which input samples are read. The path is relative to the 'samples' namib directory. Default: ''
        
        injections-path             Flag to activate a hard comparison in case only two options are compared. This parameter is used only in violin and ridgeline plots. Default: 0
        injections-number           Flag to read the evidence from the samples. Default: 0
        snr-cut                     Flag to plot the prior. Option currently implemented only in corner plot and not fully implemented. Default: 0
        ifar-cut                    List of true values of the selected parameters. Currently implemented only with corner-sns=0. Default: []
        selection-effects-cut       Flag to read the evidence from the samples. Default: 0

        O3-cosmology                List of modes for which the QNMs are computed. This option is used only when QNMs are computed from {Mf, af}. Default: [(2,2)]
        simulation                  Flag to convert the damping time in [ms] and scale amplitudes as [1e-21]. The option is used to compare samples from Damped Sinusoids with other models. Default: 0
        data-path                   Flag to use pyRing fits to compute the QNMs. Default: 1

    # ----- #
    # model #
    # ----- #

        model-primary               Options: {'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear'}
        model-secondary             Options: {'PowerLaw', 'PowerLaw-Gaussian', 'MassRatio'}
        model-rate                  Options: {'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution', 'PowerLaw'}
        redshift-transition         Options: {'sigmoid', 'double-sigmoid', 'linear', 'linear-sinusoid'}

        positive-peak               List of parameters bounds that are used in the plots. Default: []. Syntax: [[a, b], [c, d], ...]
        low-smoothing               List of the selected parameters ordering. Default: []
        priors                      List of the compared parameters ordering. Default: []
        scale-free                  List of the compared parameters ordering. Default: []

    # ------- #
    # sampler #
    # ------- #

        nparallel                   Flag to save the imput samples filtered on the selected parameters. They are saved in 'output/reduced_posteriors'. Default: 0
        neffPE                      Flag to save the medians of the selected parameters. They are saved in 'output/output_medians'. Default: 0
        neffINJ                     Option to downsample the input samples, taking the corresponding percent value (downsampling=1, 100% of initial samples). Default: 1

        sampler                     Flag to read the evidence from the samples. Default: 0
        nlive                       Flag to read the evidence from the samples. Default: 0
        npool                       Flag to read the evidence from the samples. Default: 0

"""