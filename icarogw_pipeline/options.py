import ast

def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                      : 'icarogw_run',
        'screen-output'               : False,

        'injections-path'             : '',
        'injections-number'           : 1,
        'selection-effects-cut'       : 'snr',
        'snr-cut'                     : 12.,
        'ifar-cut'                    : 4.,

        'data-path'                   : '',
        'O3-cosmology'                : False,
        'simulation'                  : True,
        'remove-events'               : [],
        'inverse-mass-ratio'          : False,
        'PE-prior-distance'           : 'dL3',
        'PE-prior-masses'             : 'm1-m2',
        'true-data'                   : False,

        # Model
        'model-primary'               : 'PowerLaw-Gaussian',                     
        'model-secondary'             : 'MassRatio-Gaussian',
        'model-rate'                  : 'PowerLaw',
        'model-cosmology'             : 'FlatLambdaCDM',
        'model-bkg-cosmo'             : 'FlatLambdaCDM',

        'redshift-transition'         : 'linear',
        'redshift-mixture'            : True,
        'low-smoothing'               : False,
        'priors'                      : {},
        'scale-free'                  : False,
        'single-mass'                 : False,

        # Sampler
        'sampler'                     : 'dynesty',
        'neffPE'                      : 10,
        'neffINJ'                     : None,
        'loglike-var'                 : 0,

        'nlive'                       : 500,
        'print-method'                : 'interval-60',
        'sample'                      : 'acceptance-walk',
        'naccept'                     : 60,
        'queue-size'                  : 1,

        'nwalkers'                    : 64,
        'nsteps'                      : 1000,
        'ntemps'                      : 10,
        'threads'                     : 1,
        'nparallel'                   : 1,
        'npool'                       : 10,

        # Plots
        'N-points'                    : 500,
        'N-z-slices'                  : 5,
        'bounds-m1'                   : [0, 100],
        'bounds-m2'                   : [0, 100],
        'bounds-q'                    : [0, 1],
        'bounds-dL'                   : [0, 10000],
        'bounds-z'                    : [1e-5, 0.8],
        'm1-logscale'                 : True,
        'log10-PDF'                   : False,  
        'true-values'                 : {},
        'selection-effects'           : False,
        'plot-prior'                  : True,
        'N-samps-prior'               : 500,

        'estimate-observed-method'    : 'KDE',
        'estimate-observed-method-m1' : 'GMM',
        'KDE-bandwidth-scale'         : 3,
        'KDE-bandwidth-scale-m1'      : 8,
        'GMM-components'              : 6,
        'N-points-KDE-GMM'            : 500,

        'percentiles'                 : {'ll': 5, 'l': 16, 'm': 50, 'h': 84, 'hh': 95},
        'downsample-postprocessing'   : 1,
    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if (key == 'output') or (key == 'injections-path') or (key == 'selection-effects-cut') or (key == 'data-path') or (key == 'PE-prior-distance') or (key == 'PE-prior-masses'):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if (key == 'injections-number') or (key == 'snr-cut') or (key == 'ifar-cut'):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if (key == 'O3-cosmology') or (key == 'simulation') or (key == 'distance-prior-PE') or (key == 'screen-output') or (key == 'true-data'):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass
        if (key == 'remove-events'):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass

        # Model
        if (key == 'model-primary') or (key == 'model-secondary') or (key == 'model-rate') or (key == 'model-cosmology') or (key == 'model-bkg-cosmo') or (key == 'redshift-transition'):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if (key == 'redshift-mixture') or (key == 'low-smoothing') or (key == 'scale-free') or (key == 'single-mass') or (key == 'inverse-mass-ratio'):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if (key == 'priors'):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass

        # Sampler
        if (key == 'sampler') or (key == 'print-method'):
            try: input_pars[key] = Config.get('sampler', key)
            except: pass
        if (key == 'nparallel') or (key == 'neffPE') or (key == 'neffINJ') or (key == 'nlive') or (key == 'queue-size') or (key == 'nwalkers') or (key == 'nsteps') or (key == 'ntemps') or (key == 'threads') or (key == 'npool'):
            try: input_pars[key] = Config.getint('sampler', key)
            except: pass
        if (key == 'loglike-var'):
            try: input_pars[key] = Config.getfloat('sampler', key)
            except: pass

        # Plots
        if (key == 'estimate-observed-method') or (key == 'estimate-observed-method-m1'):
            try: input_pars[key] = Config.get('plots', key)
            except: pass
        if (key == 'N-points') or (key == 'N-z-slices') or (key == 'N-points-KDE-GMM') or (key == 'N-samps-prior') or (key == 'GMM-components'):
            try: input_pars[key] = Config.getint('plots', key)
            except: pass
        if (key == 'true-values') or (key == 'bounds-m1') or (key == 'bounds-m2') or (key == 'bounds-q') or (key == 'bounds-dL') or (key == 'bounds-z')  or (key == 'percentiles'):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass
        if (key == 'selection-effects') or (key == 'plot-prior') or (key == 'm1-logscale') or (key == 'log10-PDF'):
            try: input_pars[key] = Config.getboolean('plots', key)
            except: pass
        if (key == 'downsample-postprocessing') or (key == 'KDE-bandwidth-scale') or (key == 'KDE-bandwidth-scale-m1'):
            try: input_pars[key] = Config.getfloat('plots', key)
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
        'w0'          : -1.,
        'wa'          : 0.,
        'xi'          : 0.,
        'eps0'        : 0.,
        'Xi0'         : 1.,
        'n'           : 1.,
        'D'           : 4.,
        'Rc'          : 1., # Mpc: same units as dL 
        'cM'          : 0.,
        'alphalog_1'  : 0.,
        'alphalog_2'  : 0.,
        'alphalog_3'  : 0.,

        # Primary mass distribution
        'delta_m'       : [   0.  ,  10.  ],
        'delta_m_a'     : [   0.  ,  30.  ],
        'delta_m_b'     : [   0.  ,  30.  ],
        'delta_m_c'     : [   0.  ,  30.  ],
        'delta'         : [   0.  ,   0.15],

        'alpha'         : [  -4.  , 120.  ],
        'alpha_z0'      : [  -4.  , 120.  ],
        'alpha_z1'      : [-100.  , 100.  ],
        'mu_alpha'      : [   0.  , 100.  ],

        'alpha_a'       : [  -4.  , 120.  ],
        'alpha_b'       : [  -4.  ,  20.  ],
        'alpha_c'       : [  -4.  ,  20.  ],
        'break_p'       : [   0.  ,   1.  ],
        'm_b'           : [   5.  ,   7.  ],

        'alpha_a_z0'    : [  -4.  , 120.  ],
        'alpha_b_z0'    : [  -4.  , 120.  ],
        'alpha_c_z0'    : [  -4.  , 120.  ],
        'alpha_a_z1'    : [-100.  , 100.  ],
        'alpha_b_z1'    : [-100.  , 100.  ],
        'alpha_c_z1'    : [-100.  , 100.  ],
        'mmin_a_z0'     : [   1.  , 100.  ],
        'mmin_b_z0'     : [   1.  , 100.  ],
        'mmin_c_z0'     : [   1.  , 100.  ],
        'mmin_a_z1'     : [-100.  , 100.  ],
        'mmin_b_z1'     : [-100.  , 100.  ],
        'mmin_c_z1'     : [-100.  , 100.  ],
        'mmax_a_z0'     : [  30.  , 200.  ],
        'mmax_b_z0'     : [  30.  , 200.  ],
        'mmax_c_z0'     : [  30.  , 200.  ],
        'mmax_a_z1'     : 0.,
        'mmax_b_z1'     : 0.,
        'mmax_c_z1'     : 0.,

        'mmin'          : [   1.  , 100.  ],
        'mmin_z0'       : [   1.  , 100.  ],
        'mmin_z1'       : [-100.  , 100.  ],
        'mmax'          : [  30.  , 200.  ],
        'mmax_z0'       : [  30.  , 200.  ],
        'mmax_z1'       : 0.,

        'mu_zt'         : [   0.  ,   1.  ],
        'mu_delta_zt'   : [   1.  , 100.  ],
        'sigma_zt'      : [   0.  ,   1.  ],
        'sigma_delta_zt': [   1.  , 100.  ],

        'mmin_a'        : [   1.  , 100.  ],
        'mmin_b'        : [   1.  , 100.  ],
        'mmin_c'        : [   1.  , 100.  ],
        'mmax_a'        : [  30.  , 200.  ],
        'mmax_b'        : [  30.  , 200.  ],
        'mmax_c'        : [  30.  , 200.  ],

        'mu_g'          : [  20.  ,  60.  ],
        'mu_z0'         : [  20.  ,  60.  ],
        'mu_z1'         : [ -80.  ,  80.  ],
        'mu_z2'         : [ -80.  ,  80.  ],
        'sigma_g'       : [   1.  ,  30.  ],
        'sigma_z0'      : [   1.  ,  30.  ],
        'sigma_z1'      : [   0.  ,  20.  ],
        'sigma_z2'      : [   0.  ,  20.  ],
        'mmin_g'        : [   2.  ,  50.  ],

        'mu_z0_a'       : [   1.  , 100.  ],
        'mu_z0_b'       : [   1.  , 100.  ],
        'mu_z0_c'       : [   1.  , 100.  ],
        'mu_z1_a'       : [-100.  , 100.  ],
        'mu_z1_b'       : [-100.  , 100.  ],
        'mu_z1_c'       : [-100.  , 100.  ],
        'sigma_z0_a'    : [   1.  ,  50.  ],
        'sigma_z0_b'    : [   1.  ,  50.  ],
        'sigma_z0_c'    : [   1.  ,  50.  ],
        'sigma_z1_a'    : [   0.  , 100.  ],
        'sigma_z1_b'    : [   0.  , 100.  ],
        'sigma_z1_c'    : [   0.  , 100.  ],

        'lambda_peak'   : [   0.  ,   1.  ],
        'mix_z0'        : [   0.  ,   1.  ],
        'mix_z1'        : [   0.  ,   1.  ],
        'mix_alpha_z0'  : [   0.  ,   1.  ],
        'mix_alpha_z1'  : [   0.  ,   1.  ],
        'mix_beta_z0'   : [   0.  ,   1.  ],
        'mix_beta_z1'   : [   0.  ,   1.  ],
        'mix_alpha'     : [   0.  ,   1.  ],
        'mix_beta'      : [   0.  ,   1.  ],
        'mix'           : [   0.  ,   1.  ],

        'amp'           : [   0.  ,   0.2 ],
        'freq'          : [ -20.  ,  20.  ],
        'zt'            : [   0.  ,   1.  ],
        'delta_zt'      : [   1.  , 100.  ],

        'a_johnson'     : [   0.  ,   3.  ],
        'b_johnson'     : [   0.1 ,   5.  ],
        'loc_johnson'   : [   3.  ,   8.  ],
        'scale_johnson' : [   0.01,  10.  ],

        # Secondary mass distribution
        'beta'          : [ -20.  ,  20.  ],
        'mu_q'          : [   0.1 ,   1.  ],
        'sigma_q'       : [   0.01,   0.9 ],
        'alpha_q'       : [ -20.  ,  20.  ],

        'a_gamma'       : [   1.  ,  10.  ],
        'theta'         : [   0.01,   1.  ],

        # Rate evolution
        'gamma'         : [ -50.  ,  30.  ],
        'kappa'         : [ -20.  ,  10.  ],
        'zp'            : [   0.  ,   4.  ],
        'R0'            : [   0.  , 100.  ],
        'mu_r'          : [-100.  ,   0.  ],
        'sigma_r'       : [   1.  , 100.  ],
        'amp_r'         : [   1.  , 200.  ],

        'a'             : [   1.  ,   5.  ],
        'b'             : [   1.  ,  20.  ],
        'loc'           : [   0.  ,  10.  ],
        'scale'         : [   1.  , 150.  ],
    }

    return prior


usage = """

\nWelcome to the ICAROGW runner helper!\n

    # ----- #
    # input #
    # ----- #

        output                      [str  ]  Path where the run output is saved. Default: 'icarogw_run'.
        screen-output               [bool ]  Flag to deviate the standard output to screen. Default: '0'.

        injections-path             [str  ]  Path of the injections used to evaluate selection effects. Default: ''.
        injections-number           [int  ]  Number of generated injections used for selection effects. Default: 1.
        selection-effects-cut       [str  ]  Type of cut to select events and injections. Options: 'snr', 'ifar'. Default: 'snr'.
        snr-cut                     [float]  Value of signal-to-noise ratio used as detection threshold for events and injections. Default: 12.
        ifar-cut                    [float]  Value of inverse false-alarm-rate used as detection threshold for events and injections. Default: 4

        data-path                   [str  ]  Path of the single event data. Default: ''.
        O3-cosmology                [bool ]  Option to process PE samples using O3 data from the LVK GWTC-3 cosmology paper (https://arxiv.org/abs/2111.03604). Default: 0.
        simulation                  [bool ]  Option to process PE samples using simulated events. Default: 1.
        remove-events               [list ]  List of events to be removed from the analysis. Example: ['GW190412_053044', 'GW190521_030229']. Default: [].
        PE-prior-distance           [str  ]  Option to re-weight the PE samples on the luminosity distance prior used in the single event parameter estimation. Options: 'dL' (uniform in luminosity distance), 'dL3' (uniform in comoving volume). Default: 'dL3'.
        PE-prior-masses             [str  ]  Option to re-weight the PE samples on the mass prior used in the single event parameter estimation. Options: 'm1-m2' (uniform in component masses), 'Mc-q' (uniform in chirp mass and mass ratio). Default: 'm1-m2'.
        true-data                   [bool ]  Flag to only use the true values for the events in the analysis instead of full posteriors. This is equivalent to use one PE sample for each event. Default: 0.
    
    # ----- #
    # model #
    # ----- #

        model-primary               [str  ]  Model distribution for the primary object. Options: 'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-Gaussian', 'DoublePowerlaw', 'PowerLaw-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftQuadratic', 'PowerLaw-GaussianRedshiftPowerLaw', 'PowerLaw-GaussianRedshiftSigmoid', 'PowerLawBroken-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear', 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_GaussianRedshiftLinear'. Default: 'PowerLaw-Gaussian'.
        model-secondary             [str  ]  Model distribution for the secondary object. Options: 'Mass2-PowerLaw', 'MassRatio-PowerLaw', 'MassRatio-Gaussian', 'MassRatio-Gamma'. Default: 'MassRatio-Gaussian'.
        model-rate                  [str  ]  Model distribution for the rate evolution. Options: 'PowerLaw', 'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution', 'Gaussian'. Default: 'PowerLaw'.
        model-cosmology             [str  ]  Model for cosmology. Options: 'FlatLambdaCDM', 'FlatwCDM', 'Flatw0waCDM', 'wIDS_linDE', 'Xi0', 'eps0', 'extraD', 'cM', 'alphalog'. Default: 'FlatLambdaCDM'
        model-bkg-cosmo             [str  ]  Model for background cosmology if model-cosmology is a modified gravity model. Options: 'FlatLambdaCDM', 'FlatwCDM', 'Flatw0waCDM', 'wIDS_linDE'. Default: 'FlatLambdaCDM'

        redshift-transition         [str  ]  Model function for the redshift evolution of the mixture functions. The option only applies to primary mass redshift evolving models. Options: 'linear', 'sigmoid'. Default: 'linear'.
        redshift-mixture            [bool ]  Flag to allow for the mixture functions to evolve in redshift. If zero, the mixture functions are stationary in redshift. The option only applies to primary mass redshift evolving models. Default: 1.
        low-smoothing               [bool ]  Flag to apply a smoothing function to the Powerlaws minimum mass. The option only applies to the mass models including Powerlaws. Default: 0.
        priors                      [dict ]  Dictionary of the prior bounds for the population parameters. Default values are set in 'icarogw_pipeline.options.default_priors'.
        scale-free                  [bool ]  Flag to use the scale-free likelihood fromulation. This is equivant to marginalizing over the expected number of events assuming a Jeffrey prior. Default: 0.
        single-mass                 [bool ]  Flag to use only one mass for the single-event parameters. Default: 0.
        inverse-mass-ratio          [bool ]  Flag to use the inverse mass ratio as the secondary mass parameter, defined as q=m1/m2 with m1>m2. Default: 0.

    # ------- #
    # sampler #
    # ------- #

        sampler                     [str  ]  Sampler to be used to draw samples from the likelihood. The samplers are called from the Bilby package (https://pypi.org/project/bilby/). Options: 'dynesty', 'nessai', 'ptemcee'. Default: 'dynesty'.
        neffPE                      [float]  Threshold on the number of effective PE samples per event contributing to the numerical evaluation of the likelihood. The likelihood is set to zero for the population samples bringing neff PE lower than this threshold, to ensure numerical stability. Default: 10.
        neffINJ                     [float]  Threshold on the number of effective injections contributing to the numerical evaluation of the likelihood. The likelihood is set to zero for the population samples bringing neff INJ lower than this threshold, to ensure numerical stability. If 'None', the threshold is set to four times the number of events. Default: None.
        loglike-var                 [float]  Threshold on the log-likelihood variance for the numerical evaluation of the likelihood. The likelihood is set to zero for the population samples bringing log-likehood variance lower than this threshold, to ensure numerical stability. If set to zero, the option is deactivated, otherwise it uses the float passed. If not zero, the option superseed the Neff PE and injections. Default: 0.

        nlive                       [int  ]  Number of live points used by the nested sampler. Option only available for Nested Samplers. Default: 500.
        print-method                [str  ]  Method for printing the sampler output. Dynesty uses a tqdm bar by default, otherwise passing 'interval-$TIME' it prints to sdtout every $TIME seconds. Default: 'interval-60'.
        sample                      [str  ]  Methods to perform the MCMC evolution to find a new point with a nested sampler. Option only available for Nested Samplers. More information on the different methods can be found in the related Bilby documentation (https://bilby-dev.github.io/bilby/dynesty-guide.html). Options: 'act-walk', 'acceptance-walk', 'rwalk'. Default: 'acceptance-walk'.
        npool                       [int  ]  Number of parallel process to be executed (see dynesty documentation: https://dynesty.readthedocs.io/en/stable/quickstart.html#parallelization). If running on a cluster, must match the number of . Default: 1.
        # queue-size                  [int  ]  Number of parallel process to be executed (see dynesty documentation: https://dynesty.readthedocs.io/en/stable/quickstart.html#parallelization). It corresponds to the number of threads used. Default: 1.

        naccept                     [int  ]  The length of the MCMC chains during the run follows a Poisson distribution with mean naccept. Option only available for Nested Samplers and only applies to the sample method 'acceptance-walk'. Default: 60.
        nwalkers                    [int  ]  Number of parallel chains (walkers) running in the MCMC ensemble. Option only available for MCMC samplers. Default: 64.
        nsteps                      [int  ]  Number of steps taken by each walker in the MCMC samplers. Option only available for MCMC samplers. Default: 1000.
        ntemps                      [int  ]  Number of parallel-tempered chains of the MCMC sampler. Option only available for MCMC samplers. Default: 10.
        threads                     [int  ]  Number of CPU threads used for parallel computation. Option only available for MCMC samplers. Default: 1.
        nparallel                   [int  ]  Number of likelihood evaluations performed simultaneously. While 'threads' distributes MCMC steps across CPU threads, 'nparallel' parallelizes across multiple processes. Option only available for MCMC samplers. Default: 1.

    # ----- #
    # plots #
    # ----- #

        N-points                    [int  ]  Number of points used to evaluate the reconstructed distributions. Default: 500.
        N-z-slices                  [int  ]  Number of points in redshift for which the conditional primary mass distribution is evaluated. Default: 5.
        bounds-m1                   [list ]  Bounds used to reconstruct the primary mass distribution. Default: [0, 100].
        bounds-m2                   [list ]  Bounds used to reconstruct the secondary mass distribution. Default: [0, 100].
        bounds-q                    [list ]  Bounds used to reconstruct the mass ratio distribution. Default: [0, 1].
        bounds-dL                   [list ]  Bounds used to reconstruct the distribution in luminosity distance. Default: [0, 10000].
        bounds-z                    [list ]  Bounds used to reconstruct the distribution in redshift and the rate evolution. Default: [1e-5, 0.8].
        m1-logscale                 [bool ]  Option to plot the primary mass distribution in log-scale. Default: 1.
        true-values                 [bool ]  Option to plot the true values of the events population. Default: {}.
        selection-effects           [bool ]  Option to show on the primary mass plot the 90% CI of the injections used to compute selection effects, qulitatively corresponding to the detector sensitivity. Default: 0.
        plot-prior                  [bool ]  Option to plot reconstructed distributions from samples drawn from the analysis priors. Default: 1.
        N-samps-prior               [int  ]  Number of samples drawn from the priors to estimate the prior reconstructed distributions. Default: 500.

        estimate-observed-method    [str  ]  Method used to estimate the observed distribution from the injection samples. Options: 'KDE', 'GMM'. Default: 'KDE'.
        estimate-observed-method-m1 [str  ]  Method used to estimate the observed distribution from the injection samples for the primary mass. Options: 'KDE', 'GMM'. Default: 'GMM'.
        KDE-bandwidth-scale         [float]  Float to scale the bandwidth of the KDE estimator compared with the base Silverman's rule. Default: 3.
        KDE-bandwidth-scale-m1      [float]  Float to scale the bandwidth of the KDE estimator compared with the base Silverman's rule for the primary mass. Default: 8.
        GMM-components              [int  ]  Number of components used in the Gaussian Mixture Model estimator. Default: 6.
        N-points-KDE-GMM            [int  ]  Number of points used to evaluate the KDE and GMM estimators. Default: 500.

        percentiles                 [dict ]  Dictionary with the percentiles used to estimate the credible intervals. Default: {'ll': 5, 'l': 16, 'm': 50, 'h': 84, 'hh': 95}.
        downsample-postprocessing   [float]  Option to downsample the posterior samples used in the postprocessing, taking the corresponding percent value (downsampling=1, 100% of initial samples). Default: 1.
"""