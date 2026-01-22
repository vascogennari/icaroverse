import ast

def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                        : 'icarogw_population',
        'run-type'                      : 'population',
        'screen-output'                 : False,
        'postprocessing'                : False,
        'drawing-method'                : 'rejection-sampling',

        'events-number'                 : 1000,
        'estimate-events-number'        : False,
        'R0'                            : 25.,
        'observation-time'              : 1.,
        'seed'                          : 42,
        'save-strain'                   : False,

        'injections-number'             : 1000,
        'injections-number-bank'        : 100,
        'inverse-checkpoint-rate'       : 1,
        'parallel'                      : False,
        'n-processes'                   : 10,

        'selection-effects-cut'         : 'snr',
        'SNR-method'                    : 'bilby',
        'observing-run'                 : 'O3',
        'SNR-cut'                       : 12.,
        'SNR-soft-cut'                  : -1.,

        # SNR-options
        'snr-bilby-waveform'            : 'IMRPhenomXHM',
        'snr-bilby-precessing-wf'       : False,
        'snr-bilby-reference-frequency' : 20.,
        'snr-bilby-sampling-frequency'  : 2048.,

        'snr-proxy-N-detectors'         : 2,
        'snr-proxy-fgw-cut'             : 15.,
        'snr-proxy-SNR-reference'       : 9.,
        'snr-proxy-dL-reference'        : 1.5,
        'snr-proxy-Mc-reference'        : 25.,
        'snr-proxy-theta-path'          : 'Pw_three.dat',
        
        'use-icarogw-sim-inj'           : False,
        'icarogw-sim-mass-model'        : 'PowerLawPeak',
        'icarogw-sim-draw-dL'           : 'uniform-dL',
        'icarogw-sim-z-max'             : 1.5,

        # Model
        'model-primary'                 : 'PowerLaw-Gaussian',                     
        'model-secondary'               : 'MassRatio-Gaussian',
        'model-rate'                    : 'PowerLaw',
        'model-cosmology'               : 'FlatLambdaCDM',
        'model-bkg-cosmo'               : 'FlatLambdaCDM',

        'zmax'                          : 20.,
        'redshift-transition'           : 'linear',
        'redshift-mixture'              : True,
        'low-smoothing'                 : False,
        'single-mass'                   : False,
        'truths'                        : {},
        'zmax'                          : 20.,

        'log10-PDF'                     : False,
        'inverse-mass-ratio'            : False,

        'splines-number'                : 10,
        'splines-order'                 : 3,
        'spacing'                       : 'log',

        # Plots
        'N-points'                      : 10000,
        'bounds-m1'                     : [0, 200],
        'bounds-m2'                     : [0, 200],
        'bounds-q'                      : [0.1, 1],
        'bounds-dL'                     : [0, 100000],
        'bounds-z'                      : [1e-6, 3.0],
        'plot-astrophysical'            : False,
    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if (key == 'output') or (key == 'run-type') or (key == 'selection-effects-cut') or (key == 'SNR-method') or (key == 'observing-run') or (key == 'drawing-method'):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if (key == 'SNR-cut') or (key == 'SNR-soft-cut') or (key == 'frequency-cut') or (key == 'R0') or (key == 'observation-time'):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if (key == 'events-number') or (key == 'injections-number') or (key == 'injections-number-bank') or (key == 'inverse-checkpoint-rate') or (key == 'n-processes') or (key == 'seed'):
            try: input_pars[key] = Config.getint('input', key)
            except: pass
        if (key == 'screen-output') or (key == 'postprocessing') or (key == 'flat-PSD') or (key == 'log10-PDF') or (key == 'estimate-events-number') or (key == 'parallel') or (key == 'save-strain'):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass

        # SNR-options
        if (key == 'snr-bilby-waveform') or (key == 'snr-proxy-theta-path') or (key == 'icarogw-sim-mass-model') or (key == 'icarogw-sim-draw-dL') or (key == 'spacing'): 
            try: input_pars[key] = Config.get('snr-options', key)
            except: pass
        if (key == 'snr-bilby-reference-frequency') or (key == 'snr-bilby-sampling-frequency') or (key == 'snr-proxy-SNR-reference') or (key == 'snr-proxy-dL-reference') or (key == 'snr-proxy-Mc-reference') or (key == 'snr-proxy-fgw-cut') or (key == 'icarogw-sim-z-max'):
            try: input_pars[key] = Config.getfloat('snr-options', key)
            except: pass
        if (key == 'snr-proxy-N-detectors'):
            try: input_pars[key] = Config.getint('snr-options', key)
            except: pass
        if (key == 'use-icarogw-sim-inj') or (key == 'snr-bilby-precessing-wf'):
            try: input_pars[key] = Config.getboolean('snr-options', key)
            except: pass

        # Model
        if (key == 'model-primary') or (key == 'model-secondary') or (key == 'model-rate') or (key == 'model-cosmology') or (key == 'model-bkg-cosmo') or (key == 'redshift-transition'):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if (key == 'redshift-mixture') or (key == 'low-smoothing') or (key == 'single-mass') or (key == 'log10-PDF') or (key == 'inverse-mass-ratio'):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if (key == 'zmax'):
            try: input_pars[key] = Config.getfloat('model', key)
            except: pass
        if (key == 'truths'):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass
        if (key == 'splines-number') or (key == 'splines-order'):
            try: input_pars[key] = Config.getint('snr-options', key)
            except: pass

        # Plots
        if (key == 'N-points'):
            try: input_pars[key] = Config.getint('plots', key)
            except: pass
        if (key == 'plot-astrophysical'):
            try: input_pars[key] = Config.getboolean('plots', key)
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
        'delta_m'       : 5.,
        'delta_m_a'     : 5.,
        'delta_m_b'     : 5.,
        'delta_m_c'     : 5.,
        'delta'         : 0.1,

        'alpha'         : 4.,
        'alpha_z0'      : 4.,
        'alpha_z1'      : 0.,
        'mu_alpha'      : 4.,

        'alpha_a'       : 50.,
        'alpha_b'       : 10.,
        'alpha_c'       : 10.,
        'break_p'       : 0.5,
        'm_b'           : 6.,
        'm_b_z0'        : 6.,
        'm_b_z10'       : 6.,
        'm_b_zt'        : 3.,
        'm_b_delta_zt'  : 3.,

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

        'mmin'          : 10.,
        'mmin_z0'       : 10.,
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

        'mu_a_z0'       : 10.,
        'mu_b_z0'       : 15.,
        'mu_c_z0'       : 35.,
        'mu_a_z1'       : 0.,
        'mu_b_z1'       : 0.,
        'mu_c_z1'       : 0.,
        'sigma_a_z0'    : 5.,
        'sigma_b_z0'    : 5.,
        'sigma_c_z0'    : 5.,
        'sigma_a_z1'    : 0.,
        'sigma_b_z1'    : 0.,
        'sigma_c_z1'    : 0.,

        'lambda_peak'   : 0.1,
        'mix_z0'        : 0.9,
        'mix_z1'        : 0.9,
        'mix_alpha_z0'  : 0.9,
        'mix_alpha_z1'  : 0.9,
        'mix_beta_z0'   : 0.05,
        'mix_beta_z1'   : 0.,
        'mix_alpha'     : 0.9,
        'mix_beta'      : 0.05,
        'mix'           : 0.9,

        'amp'           : 0.1,
        'freq'          : 0.,
        'zt'            : 0.5,
        'delta_zt'      : 50.,

        'skew_j'        : 0.8,
        'sharp_j'       : 1.,
        'peak_j'        : 5.5,
        'scale_j'       : 0.3,
        'mmin_j'        : 3.,
        'mmax_j'        : 9.,

        'c1'            : 10.,
        'c2'            : 10.,
        'c3'            : 10.,
        'c4'            : 10.,
        'c5'            : 10.,
        'c6'            : 10.,
        'c7'            : 10.,
        'c8'            : 10.,
        'c9'            : 10.,
        'c10'           : 10.,
        'c11'           : 10.,
        'c12'           : 10.,
        'c13'           : 10.,
        'c14'           : 10.,
        'c15'           : 10.,
        'c16'           : 10.,

        # Secondary mass distribution
        'beta'          : 4.,
        'mu_q'          : 0.8,
        'sigma_q'       : 0.1,
        'alpha_q'       : 1.2,

        'low_b'         : 1.,
        'high_b'        : 10.,
        'start_b'       : 0.,
        'scale_b'       : 5.,

        'a_gamma'       : 5.,
        'theta'         : 0.5,

        # Rate evolution
        'gamma'         : 0.,
        'kappa'         : 0.,
        'zp'            : 3.,
        'R0'            : 20.,
        'mu_r'          : -20.,
        'sigma_r'       : 20.,
        'amp_r'         : 20.,

        'low_b_r'       : 2.,
        'high_b_r'      : 6.,
        'start_b_r'     : 0.25,
        'scale_b_r'     : 15.,

        'z_min'         : 0.,
        'z_max'         : 1.,
    }

    return pop


usage = """

\nWelcome to the ICAROGW simulation pipeline helper!\n

    # ----- #
    # input #
    # ----- #

        output                        [str  ]  Path where the output is saved. Default: 'icarogw_population'.
        run-type                      [str  ]  Type of simulation to run. Options: 'population', 'injections', 'noise'. Default: 'population'.
        screen-output                 [bool ]  Flag to deviate the standard output to screen. Default: '0'.
        postprocessing                [bool ]  Flag to only postprocess an existing simulation. Default: '0'.
        drawing-method                [str  ]  Method used to draw samples from the target population distribution. Options: 'rejection-sampling', 'inverse-transform', 'deterministic-inverse-transform'. Default: 'rejection-sampling'.

        events-number                 [int  ]  Number of generated events to draw from the astrophysical popoulation. Option used if the run-type is 'population'. Default: 1000.
        estimate-events-number        [bool ]  Flag to set the number of generated events directly from the population rate evolution. If activated, it overwrites 'events-number'. Option used if the run-type is 'population'. Default: '0'.
        R0                            [float]  Astrophysical rate of events today [Gpc^(-3)yr^(-1)]. Used to compute the number of generated events if 'estimate-events-number' is active. Option used if the run-type is 'population'. Default: 25.
        observation-time              [float]  Time of observation [yr]. Used to compute the number of generated events if 'estimate-events-number' is active. Is overwritten by official IGWN observing-run duration if a negative value is given. Option used if the run-type is 'population'. Default: 1.
        drawing-method                [str  ]  Method used to draw samples from the target population distribution. Options: 'rejection-sampling', 'inverse-transform', 'deterministic-inverse-transform'. Default: 'rejection-sampling'.
        seed                          [int  ]  Seed for the random number generator. Default: 42.
        save-strain                   [bool ]  Flag to save strain data when generating a population. Default: False

        injections-number             [int  ]  Number of detected injections to draw from the popoulation. Option used if the run-type is 'injections'. Default: 1000.
        injections-number-bank        [int  ]  Number of injections to draw from the popoulation until 'injections-number' is obtained. Option used if the run-type is 'injections'. Default: 100.
        inverse-checkpoint-rate       [int  ]  Injections checkpoint save inverse rate. In units of number of injections batches (see 'injections-number-bank' for number of injections per batch). Option used if the run-type is 'injections'. Default: 1.
        parallel                      [bool ]  Flag to parallelize the generation of injections, launching multiple batches in parallel processes (see 'n-processes' to set the maximum number of parallel processes). Default: '0'.
        n-processes                   [int  ]  Maximum number of parallel processes to generate injections. Please make sure that this matches the number of available cores on your machine. Default: 10

        selection-effects-cut         [str  ]  Method to evaluate events detectability. Options: 'snr'. Default: 'snr'.
        SNR-cut                       [float]  SNR threshold to label an event as detectable. Default: 12.
        SNR-soft-cut                  [float]  SNR threshold to label an event as worth saving (should in principle be less than SNR-cut). Option used for injections generation only (NB: all the events are kept in case of a population generation). Default: -1. (i.e. keeps everything)
        SNR-method                    [str  ]  Method to compute the SNR. Options: 'bilby', 'proxy', 'flat-PSD', 'lisabeta'. Default: 'bilby'.
        observing-run                 [str  ]  IGWN observing run. Further defines the detectors sensitivity for SNR computation (with Bilby) and PE (with Bilby), as well as observation-time (if a negative value is given). Options: 'O3', 'O4', 'O5'. Default: 'O3'.

    # ----------- #
    # snr-options #
    # ----------- #

        snr-bilby-waveform            [str  ]  Waveform model used to compute the SNR with Bilby. Default: 'IMRPhenomXHM'.
        snr-bilby-reference-frequency [float]  Reference frequency [Hz] used by Bilby to compute the duration of the signal injected in detectors (Bilby's time to merger is defined as the duration from the time the signal is at reference-frequency, to the time of the merger. See https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.utils.calculate_time_to_merger.html#bilby.gw.utils.calculate_time_to_merger). Default: 20.
        snr-bilby-sampling-frequency  [float]  Sampling rate [Hz] used to generate the waveform with Bilby. Default: 2048.

        snr-proxy-N-detectors         [int  ]  Number of detectors used to compute the SNR with the proxy method. Default: 2.
        snr-proxy-fgw-cut             [float]  Condition on detection additional to the SNR cut. An event is detected if its estimated frequency at ISCO is higher than this value. Default: 15.
        snr-proxy-SNR-reference       [float]  Reference SNR used to compute the approximate analytical SNR. Default: 9.
        snr-proxy-dL-reference        [float]  Reference luminosity distance that should give an approximate SNR of 'snr-proxy-SNR-reference'. Default: 1.5.
        snr-proxy-Mc-reference        [float]  Reference chirp mass that should give an approximate SNR of 'snr-proxy-SNR-reference'. Default: 25.
        snr-proxy-theta-path          [str  ]  Path to the file with average on extrinsic parameters. Default: 'Pw_three.dat'.

        use-icarogw-sim-inj           [bool ]  Flag to use the ICAROGW simulation script to generate the injections. Default: 0.
        icarogw-sim-mass-model        [str  ]  Mass model used to generate the injections with the ICAROGW simulation script. The only mass models available with this method are: 'PowerLaw', 'PowerLawPeak', 'MultiPeak'. Default: 'PowerLawPeak'.
        icarogw-sim-draw-dL           [str  ]  Method used to draw the luminosity distance samples. Options: 'uniform-dL' (uniform in luminosity distance), 'uniform-z' (uniform in redshift), 'uniform-volume' (uniform in comoving volume). Default: 'uniform-dL'.
        icarogw-sim-z-max             [float]  Maximum redshift used to draw the samples. Default: 1.5.

    # ----- #
    # model #
    # ----- #

        model-primary                 [str  ]  Model distribution for the primary object. Options: 'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-Gaussian', 'DoublePowerlaw', 'PowerLaw-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftQuadratic', 'PowerLaw-GaussianRedshiftPowerLaw', 'PowerLaw-GaussianRedshiftSigmoid', 'PowerLawBroken-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear', 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_GaussianRedshiftLinear'. Default: 'PowerLaw-Gaussian'.
        model-secondary               [str  ]  Model distribution for the secondary object. Options: 'Mass2-PowerLaw', 'MassRatio-PowerLaw', 'MassRatio-Gaussian', 'MassRatio-Gamma'. Default: 'MassRatio-Gaussian'.
        model-rate                    [str  ]  Model distribution for the rate evolution. Options: 'PowerLaw', 'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution', 'Gaussian'. Default: 'PowerLaw'.
        model-cosmology               [str  ]  Model for cosmology. Options: 'FlatLambdaCDM', 'FlatwCDM', 'Flatw0waCDM', 'wIDS_linDE', 'Xi0', 'eps0', 'extraD', 'cM', 'alphalog'. Default: 'FlatLambdaCDM'
        model-bkg-cosmo               [str  ]  Model for background cosmology if model-cosmology is a modified gravity model. Options: 'FlatLambdaCDM', 'FlatwCDM', 'Flatw0waCDM', 'wIDS_linDE'. Default: 'FlatLambdaCDM'
        zmax                          [float]  Maximum redshift up to which the cosmology wrappers are initialized. Default: 20.
        redshift-transition           [str  ]  Model function for the redshift evolution of the mixture functions. The option only applies to primary mass redshift evolving models. Options: 'linear', 'sigmoid'. Default: 'linear'.
        redshift-mixture              [bool ]  Flag to allow for the mixture functions to evolve in redshift. If zero, the mixture functions are stationary in redshift. The option only applies to primary mass redshift evolving models. Default: 1.
        low-smoothing                 [bool ]  Flag to apply a smoothing function to the Powerlaws minimum mass. The option only applies to the mass models including Powerlaws. Default: 0.
        single-mass                   [bool ]  Flag to use only one mass for the single-event parameters. Default: 0.
        truths                        [dict ]  Dictionary with the true values of the population parameters. Default: {}.
        zmax                          [float]  Maximum redshift up to which the cosmology wrappers are initialized. Default: 20.
        log10-PDF                     [bool ]  Flag to use distributions defined in log10 scale. Default: 0.
        inverse-mass-ratio            [bool ]  Flag to use the inverse mass ratio as the secondary mass parameter, defined as q=m1/m2 with m1>m2. Default: 0.

        splines-number                [int  ]  Number of splines used for the spline models. The option only applies to models including splines. Default: 10.
        splines-order                 [int  ]  Order of the splines used for the spline models. The option only applies to models including splines. Default: 3.
        spacing                       [str  ]  Spacing of the spline knots. Options: 'linear', 'log'. The option only applies to models including splines. Default: 'log'.

    # ----- #
    # plots #
    # ----- #

        N-points                      [float]  Number of points used to generate the plots. Default: 10000.
        bounds-m1                     [list ]  Bounds for the primary mass plots. Default: [0, 200].
        bounds-m2                     [list ]  Bounds for the secondary mass plots. Default: [0, 200].
        bounds-q                      [list ]  Bounds for the mass ratio plots. Default: [0.1, 1].
        bounds-dL                     [list ]  Bounds for the luminosity distance plots. Default: [0, 10000].
        bounds-z                      [list ]  Bounds for the redshift plots. Default: [1e-6, 3.0].
        plot-astrophysical            [bool ]  Flag to plot the astrophysical generated population. If False, then only the detected events are shown in the plots (NB: this option is relevant for population generation, or injections when all the astrophysical events are saved). Default: False.
"""