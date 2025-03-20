import ast

def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                        : 'icarogw_population',
        'run-type'                      : 'population',
        'screen-output'                 : 0,
        'postprocessing'                : 0,

        'events-number'                 : 1000,
        'injections-number'             : 1000,
        'injections-number-bank'        : 100,
        'selection-effects-cut'         : 'snr',
        'SNR-cut'                       : 12.,
        'estimate-events-number'        : 0,
        'R0'                            : 25.,
        'observation-time'              : 1.,

        'SNR-method'                    : 'bilby',
        'PSD-path'                      : '',

        # SNR-options
        'snr-bilby-observing-run'       : 'O3',
        'snr-bilby-waveform'            : 'IMRPhenomXHM',
        'snr-bilby-reference-frequency' : 20.,
        'snr-bilby-sampling-frequency'  : 2048.,

        'snr-pycbc-detectors'           : ['H1', 'L1', 'V1'],
        'snr-pycbc-observing-run'       : 'O3',
        'snr-pycbc-sampling-rate'       : 2048.,
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

        # Model
        'model-primary'                 : 'PowerLaw-Gaussian',                     
        'model-secondary'               : 'MassRatio-Gaussian',
        'model-rate'                    : 'PowerLaw',

        'redshift-transition'           : 'linear',
        'redshift-mixture'              : 1,
        'low-smoothing'                 : 0,
        'single-mass'                   : 0,
        'truths'                        : {},

        'log10-PDF'                     : 0,
        'inverse-mass-ratio'            : 1,

        # Plots
        'N-points'                      : 10000,
        'bounds-m1'                     : [0, 200],
        'bounds-m2'                     : [0, 200],
        'bounds-q'                      : [0.1, 1],
        'bounds-dL'                     : [0, 10000],
        'bounds-z'                      : [1e-6, 3.0],
    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if (key == 'output') or (key == 'run-type') or (key == 'selection-effects-cut') or (key == 'SNR-method') or (key == 'PSD-path'):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if (key == 'SNR-cut') or (key == 'frequency-cut') or (key == 'R0') or (key == 'observation-time'):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if (key == 'events-number') or (key == 'injections-number') or (key == 'injections-number-bank'):
            try: input_pars[key] = Config.getint('input', key)
            except: pass
        if (key == 'screen-output') or (key == 'postprocessing') or (key == 'flat-PSD') or (key == 'log10-PDF') or (key == 'inverse-mass-ratio') or (key == 'estimate-events-number'):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass

        # SNR-options
        if (key == 'snr-bilby-observing-run') or (key == 'snr-bilby-waveform') or (key == 'snr-pycbc-observing-run') or (key == 'snr-pycbc-waveform') or (key == 'snr-pycbc-method') or (key == 'snr-proxy-theta-path') or (key == 'icarogw-sim-mass-model') or (key == 'icarogw-sim-draw-dL'): 
            try: input_pars[key] = Config.get('snr-options', key)
            except: pass
        if (key == 'snr-bilby-reference-frequency') or (key == 'snr-bilby-sampling-frequency') or (key == 'snr-pycbc-sampling-rate') or (key == 'snr-pycbc-delta-f') or (key == 'snr-pycbc-f-low') or (key == 'snr-proxy-SNR-reference') or (key == 'snr-proxy-dL-reference') or (key == 'snr-proxy-Mc-reference') or (key == 'icarogw-sim-z-max'):
            try: input_pars[key] = Config.getfloat('snr-options', key)
            except: pass
        if (key == 'snr-proxy-N-detectors'):
            try: input_pars[key] = Config.getint('snr-options', key)
            except: pass
        if (key == 'snr-pycbc-detectors'):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass
        if (key == 'use-icarogw-sim-inj') or (key == 'snr-pycbc-precession'):
            try: input_pars[key] = Config.getboolean('snr-options', key)
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


usage = """

\nWelcome to the ICAROGW runner helper!\n

    # ----- #
    # input #
    # ----- #

        output                        Path where the output is saved. Default: 'icarogw_population'.
        run-type                      Type of simulation to run. Options: 'population', 'injections'. Default: 'population'.
        screen-output                 Flag to deviate the standard output to screen. Default: '0'.
        postprocessing                Flag to only postprocess an existing simulation. Default: '0'.
        
        events-number                 Number of generated events to draw from the astrophysical popoulation. Option used if the run-type is 'population'. Default: 1000.
        injections-number             Number of detected injections to draw from the popoulation. Option used if the run-type is 'injections'. Default: 1000.
        injections-number-bank        Number of injections to draw from the popoulation until 'injections-number' is obtained. Option used if the run-type is 'injections'. Default: 100.
        selection-effects-cut         Method to evaluate events detectability. Options: 'snr'. Default: 'snr'.
        SNR-cut                       SNR threshold to consider an event detectable. Default: 12.
        estimate-events-number        Flag to set the number of generated events directly from the population rate evolution. If activated, it overwrites 'events-number'. Option used if the run-type is 'population'. Default: '0'.
        R0                            Astrophysical rate of events today [Gpc^(-3)yr^(-1)]. Used to compute the number of generated events if 'estimate-events-number' is active. Default: 25.
        observation-time              Time of observation [yr]. Used to compute the number of generated events if 'estimate-events-number' is active. Default: 1.

        SNR-method                    Method to compute the SNR. Options: 'bilby', 'pycbc', 'proxy', 'flat-PSD', 'lisabeta'. Default: 'bilby'.
        PSD-path                      Path to the PSD file used to compute the SNR. This is only used if SNR-method is 'pycbc'. Default: ''.

    # ----------- #
    # snr-options #
    # ----------- #

        snr-bilby-observing-run       Detector sensitivity used to compute the SNR with Bilby. Options: 'O3', 'O4', 'O5'. Default: 'O3'.
        snr-bilby-waveform            Waveform model used to compute the SNR with Bilby. Default: 'IMRPhenomXHM'.
        snr-bilby-reference-frequency Frequency at which the binary parameters are defined when Bilby generates the waveforms. Default: 20.
        snr-bilby-sampling-frequency  Sampling rate used to generate the waveform with Bilby. Default: 2048.
        
        snr-pycbc-detectors           List of detectors used to compute the SNR with PyCBC. Options: 'H1', 'L1', 'V1', 'K1'. Default: ['H1', 'L1', 'V1'].
        snr-pycbc-observing-run       Detector sensitivity used to compute the SNR with PyCBC. Options: 'O3', 'O4', 'O5'. Default: 'O3'.
        snr-pycbc-sampling-rate       Sampling rate used to generate the waveform with PyCBC. Default: 2048.
        snr-pycbc-delta-f             Frequency resolution used to generate the waveform with PyCBC. Default: 1/128.
        snr-pycbc-f-low               Lower frequency cut-off used to generate the waveform with PyCBC. Default: 16.
        snr-pycbc-waveform            Waveform model used to compute the SNR with PyCBC. Default: 'IMRPhenomXHM'.
        snr-pycbc-method              Method used to compute the SNR with PyCBC. Options: 'opt' (optimal SNR), 'mf_fast' (matched-filter SNR with noise contribution estimated from a unit Gaussian). Default: 'mf_fast'.
        snr-pycbc-precession          Flag to include precession in the waveform model. Default: 0.

        snr-proxy-N-detectors         Number of detectors used to compute the SNR with the proxy method. Default: 2.
        snr-proxy-fgw_cut             Condition on detection additional to the SNR cut. An event is detected if its estimated frequency at ISCO is higher than this value. Default: 15.
        snr-proxy-SNR-reference       Reference SNR used to compute the approximate analytical SNR. Default: 9.
        snr-proxy-dL-reference        Reference luminosity distance that should give an approximate SNR of 'snr-proxy-SNR-reference'. Default: 1.5.
        snr-proxy-Mc-reference        Reference chirp mass that should give an approximate SNR of 'snr-proxy-SNR-reference'. Default: 25.
        snr-proxy-theta-path          Path to the file with average on extrinsic parameters. Default: 'Pw_three.dat'.

        use-icarogw-sim-inj           Flag to use the ICAROGW simulation script to generate the injections. Default: 0.
        icarogw-sim-mass-model        Mass model used to generate the injections with the ICAROGW simulation script. The only mass models available with this method are: 'PowerLaw', 'PowerLawPeak', 'MultiPeak'. Default: 'PowerLawPeak'.
        icarogw-sim-draw-dL           Method used to draw the luminosity distance samples. Options: 'uniform-dL' (uniform in luminosity distance), 'uniform-z' (uniform in redshift), 'uniform-volume' (uniform in comoving volume). Default: 'uniform-dL'.
        icarogw-sim-z-max             Maximum redshift used to draw the samples. Default: 1.5.

    # ----- #
    # model #
    # ----- #

        model-primary                 Model distribution for the primary object. Options: 'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-PowerLaw', 'PowerLaw-PowerLaw-Gaussian', 'DoublePowerlaw', 'PowerLaw-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftQuadratic', 'PowerLaw-GaussianRedshiftPowerLaw', 'PowerLaw-GaussianRedshiftSigmoid', 'PowerLawBroken-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear', 'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear', 'GaussianRedshiftLinear-GaussianRedshiftLinear-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear', 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_GaussianRedshiftLinear'. Default: 'PowerLaw-Gaussian'.
        model-secondary               Model distribution for the secondary object. Options: 'Mass2-PowerLaw', 'MassRatio-PowerLaw', 'MassRatio-Gaussian', 'MassRatio-Gamma'. Default: 'MassRatio-Gaussian'.
        model-rate                    Model distribution for the rate evolution. Options: 'PowerLaw', 'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution', 'Gaussian'. Default: 'PowerLaw'.
        
        redshift-transition           Model function for the redshift evolution of the mixture functions. The option only applies to primary mass redshift evolving models. Options: 'linear', 'sigmoid'. Default: 'linear'.
        redshift-mixture              Flag to allow for the mixture functions to evolve in redshift. If zero, the mixture functions are stationary in redshift. The option only applies to primary mass redshift evolving models. Default: 1.
        low-smoothing                 Flag to apply a smoothing function to the Powerlaws minimum mass. The option only applies to the mass models including Powerlaws. Default: 0.
        single-mass                   Flag to use only one mass for the single-event parameters. Default: 0.
        truths                        Dictionary with the true values of the population parameters. Default: {}.

        log10-PDF                     Flag to use distributions defined in log10 scale. Default: 0.
        inverse-mass-ratio            Flag to use the inverse mass ratio as the secondary mass parameter, defined as q=m1/m2 with m1>m2. Default: 0.

    # ----- #
    # plots #
    # ----- #

        N-points                      Number of points used to generate the plots. Default: 10000.
        bounds-m1                     Bounds for the primary mass plots. Default: [0, 200].
        bounds-m2                     Bounds for the secondary mass plots. Default: [0, 200].
        bounds-q                      Bounds for the mass ratio plots. Default: [0.1, 1].
        bounds-dL                     Bounds for the luminosity distance plots. Default: [0, 10000].
        bounds-z                      Bounds for the redshift plots. Default: [1e-6, 3.0].
"""