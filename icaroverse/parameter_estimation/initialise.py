import ast, numpy as np


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                       : '',
        'screen-output'                : 0,
        'PSD-path'                     : '',
        'event-from-simulation'        : 0,
        'plot-waveform'                : False,
        'plot-skymap'                  : False,

        # model
        'event-parameters'             : {},
        'strain-file'                  : '',
        'priors'                       : {},
        'mass-parameters-uniform-prior': 'm1-m2',
        'mass-parameters-sampled'      : 'm1-m2_custom',
        'observing-run'                : 'O3',
        'waveform'                     : 'IMRPhenomXHM',
        'precession'                   : False,
        'reference-frequency'          : 20.,
        'sampling-frequency'           : 2048.,
        'phase-marginalization'        : False,
        'time-marginalization'         : False,
        'distance-marginalization'     : False,
        'calibration-marginalization'  : False,

        # Sampler
        'sampler'                      : 'dynesty',
        'print-method'                 : 'interval-60',

        'nlive'                        : 500,
        'npool'                        : 10,
        'sample'                       : 'acceptance-walk',
        'naccept'                      : 60,
        'queue-size'                   : 1,
        'nwalkers'                     : 64,
        'nsteps'                       : 1000,
        'ntemps'                       : 10,
        'threads'                      : 1,
        'nparallel'                    : 1,

    }

    # Read options from config file.
    for key in input_pars.keys():

        # Input
        if (key == 'output') or (key == 'PSD-path'):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if (key == 'screen-output') or (key == 'event-from-simulation') or (key == 'plot-waveform') or (key == 'plot-skymap'):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass

        # Model
        if (key == 'event-parameters') or (key == 'priors'):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass
        if (key == 'observing-run') or (key == 'waveform') or (key == 'mass-parameters-sampled') or (key == 'strain-file') or (key == 'mass-parameters-uniform-prior'): 
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if (key == 'precession') or (key == 'phase-marginalization') or (key == 'time-marginalization') or (key == 'distance-marginalization') or (key == 'calibration-marginalization'): 
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if (key == 'reference-frequency') or (key == 'sampling-frequency'):
            try: input_pars[key] = Config.getfloat('model', key)
            except: pass

        # Sampler
        if (key == 'sampler') or (key == 'print-method'):
            try: input_pars[key] = Config.get('sampler', key)
            except: pass
        if (key == 'nparallel') or (key == 'nlive') or (key == 'queue-size') or (key == 'nwalkers') or (key == 'nsteps') or (key == 'ntemps') or (key == 'npool'):
            try: input_pars[key] = Config.getint('sampler', key)
            except: pass


    # Set precession flag based on the waveform.
    input_pars = set_precession(input_pars)

    # Initialise the event parameters.
    input_pars['event-parameters-default'] = default_event_parameters()
    if not input_pars['event-parameters'] == {}:
        for key in input_pars['event-parameters']: input_pars['event-parameters-default'][key] = input_pars['event-parameters'][key]
    input_pars['event-parameters'] = input_pars['event-parameters-default']

    del input_pars['event-parameters-default']

    # Initialise the PE priors.
    input_pars['priors-default'] = default_PE_priors(precession = input_pars['precession'])
    if not input_pars['priors'] == {}:
        for key in input_pars['priors']: input_pars['priors-default'][key] = input_pars['priors'][key]
    input_pars['priors'] = input_pars['priors-default']

    del input_pars['priors-default']

    return input_pars


def default_event_parameters():

    dict = {
        'seed'               : 8848,
        'mass_1'             : 35.0,
        'mass_2'             : 30.0,
        'a_1'                : 0.0,
        'a_2'                : 0.0,
        'tilt_1'             : 0.0,
        'tilt_2'             : 0.0,
        'phi_jl'             : 0.0,
        'phi_12'             : 0.0,
        'luminosity_distance': 2000.0,
        'theta_jn'           : 0.0,
        'psi'                : 0.0,
        'phase'              : 0.0,
        'ra'                 : 0.0,
        'dec'                : 0.0,
        'geocent_time'       : 1126259462.4,
    }

    return dict

def default_PE_priors(precession = False):

    dict = {
        'mass_1'             : [1., 300.],
        'mass_2'             : [1., 300.],
        'chirp_mass'         : [1., 300.],
        'mass_ratio'         : [0.05, 1.0],
        'luminosity_distance': [0.,1.6e4],
        'theta_jn'           : [0., np.pi],
        'psi'                : [0., 2*np.pi],
        'phase'              : [0., 2*np.pi],
        'ra'                 : [0., 2*np.pi],
        'dec'                : [-np.pi/2., np.pi/2.],
    }

    if precession:
        dict['a_1']          = [0., 0.99],
        dict['a_2']          = [0., 0.99],
        dict['tilt_1']       = [0., np.pi]
        dict['tilt_2']       = [0., np.pi]
        dict['phi_jl']       = [0., 2*np.pi]
        dict['phi_12']       = [0., 2*np.pi]
    else:
        dict['chi_1']        = [-0.99, 0.99]
        dict['chi_2']        = [-0.99, 0.99]

    return dict

def set_precession(input_pars):

    # Full list of known approximants
    approx_nonprecessing = [
        "IMRPhenomD", "IMRPhenomXAS", "IMRPhenomXHM", "SEOBNRv4", "SEOBNRv4HM", "TaylorF2", "EOBNRv2", "SEOBNRv2"]
    approx_precessing = [
        "IMRPhenomPv2", "IMRPhenomPv3", "IMRPhenomPv3HM", "SEOBNRv4PHM", "IMRPhenomTPHM", "SEOBNRv3", "SEOBNRv4", "IMRPhenomPv2", "IMRPhenomXPHM"]

    if   input_pars['waveform'] in approx_precessing   : input_pars['precession'] = True
    elif input_pars['waveform'] in approx_nonprecessing: input_pars['precession'] = False
    else:
        raise ValueError(f"Waveform {input_pars['waveform']} not recognized. Please use one of the following: {approx_nonprecessing + approx_precessing}")

    return input_pars


usage = """

\nWelcome to the Bilby pipeline helper!\n

    # ----- #
    # input #
    # ----- #

        output                        [str  ] Path where the event PE output is saved. Default: 'parameter_estimation'.
        screen-output                 [bool ] Flag to deviate the standard output to screen. Default: 0.
        PSD-path                      [str  ] Path to the directory where PSD files are stored. PSD files should be named of the form '{ifo}_{observing_run}.txt'. Default: ''.
        event-from-simulation         [bool ] True if the event was created with icaroverse simulation pipeline. If True, accepts m1d, m2d, dL entries in event-parameters dict, instead of the mass_1, mass_2, luminosity_distance normally required for bilby's waveform generator. Default: False.
        plot-waveform                 [bool ] If True, plots the waveform in each detector at the end of the sampling procedure. Default: False
        plot-skymap                   [bool ] If True, plots the skymap at the end of the sampling procedure. Default: False

    # ----- #
    # model #
    # ----- #

        event-parameters              [dict ] Dictionary containing the event parameters. Default: {}.
        priors                        [dict ] Prior ranges to use for event parameters inference. Default: {}.
        mass_parameters_uniform_prior [str  ] Mass parameters for which the prior is Uniform. Options: 'm1-m2' (uniform in componenet masses), 'Mc-q' (uniform in chirp mass and mass ratio). Default: 'm1-m2'.
        mass_parameters_sampled       [str  ] Mass parameters that are sampled over, with which bilby prior. 'm1-m2_bilby' uses bilby.gw.prio.BBHPriorDict ; 'm1-m2_custom' uses bilby.core.prior.PriorDict with custom constraint on mass_ratio ; 'Mc-q' uses bilby.core.prior.PriorDict. Options: 'm1-m2_bilby', 'm1-m2_custom', 'Mc-q'. Default: 'm1-m2_custom'.
        observing-run                 [str  ] Detector sensitivity used to compute the SNR with Bilby. Options: 'O3', 'O4', 'O5'. Default: 'O3'.
        waveform                      [str  ] Waveform model used to compute the SNR with Bilby. Default: 'IMRPhenomXHM'.
        precession                    [bool ] Flag to turn on if the waveform model is a precessing one. Default: False.
        reference-frequency           [float] Frequency at which the binary parameters are defined when Bilby generates the waveforms. Default: 20.
        sampling-frequency            [float] Sampling rate used to generate the waveform with Bilby. Default: 2048.
        phase-marginalization         [bool ] Flag to marginalize over the waveform's phase. Default: False
        time-marginalization          [bool ] Flag to marginalize over the waveform's time. Default: False
        distance-marginalization      [bool ] Flag to marginalize over the waveform's distance. Default: False
        calibration-marginalization   [bool ] Flag to marginalize over the waveform's calibration parameters. Default: False

    # ------- #
    # sampler #
    # ------- #

        sampler                       [str  ] Sampler to be used to draw samples from the likelihood. The samplers are called from the Bilby package (https://pypi.org/project/bilby/). Options: 'dynesty', 'nessai', 'ptemcee'. Default: 'dynesty'.

        nlive                         [int  ] Number of live points used by the nested sampler. Option only available for Nested Samplers. Default: 500.
        print-method                  [str  ] Method for printing the sampler output. Dynesty uses a tqdm bar by default, otherwise passing 'interval-$TIME' it prints to sdtout every $TIME seconds. Default: 'interval-60'.
        sample                        [str  ] Methods to perform the MCMC evolution to find a new point with a nested sampler. Option only available for Nested Samplers. More information on the different methods can be found in the related Bilby documentation (https://bilby-dev.github.io/bilby/dynesty-guide.html). Options: 'act-walk', 'acceptance-walk', 'rwalk'. Default: 'acceptance-walk'.
        naccept                       [int  ] The length of the MCMC chains during the run follows a Poisson distribution with mean naccept. Option only available for Nested Samplers and only applies to the sample method 'acceptance-walk'. Default: 60.
        queue-size                    [int  ] Number of parallel process to be executed (see dynesty documentation: https://dynesty.readthedocs.io/en/stable/quickstart.html#parallelization). It corresponds to the number of threads used or the size of the multiprocessing pool. Default: 1.
        npool                         [int  ] Number of available CPUs for sampling parallelization (see bilby documentation: https://lscsoft.docs.ligo.org/bilby/api/bilby.core.sampler.run_sampler.html#bilby.core.sampler.run_sampler). Default: 1.

        nwalkers                      [int  ] Number of parallel chains (walkers) running in the MCMC ensemble. Option only available for MCMC samplers. Default: 64.
        nsteps                        [int  ] Number of steps taken by each walker in the MCMC samplers. Option only available for MCMC samplers. Default: 1000.
        ntemps                        [int  ] Number of parallel-tempered chains of the MCMC sampler. Option only available for MCMC samplers. Default: 10.
        nparallel                     [int  ] Number of likelihood evaluations performed simultaneously. While 'threads' distributes MCMC steps across CPU threads, 'nparallel' parallelizes across multiple processes. Option only available for MCMC samplers. Default: 1.

"""