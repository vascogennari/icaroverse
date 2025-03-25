import ast


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'                  : '',
        'screen-output'           : 0,
        'PSD-path'                : '',

        # model
        'event-parameters'        : {},
        'priors'                  : {},
        'observing-run'           : 'O3',
        'waveform'                : 'IMRPhenomXHM',
        'reference-frequency'     : 20.,
        'sampling-frequency'      : 2048.,

        # Sampler
        'sampler'                 : 'dynesty',
        'print-method'            : 'interval-60',

        'nlive'                   : 500,

        'sample'                  : 'acceptance-walk',
        'naccept'                 : 60,
        'queue-size'              : 1,
        'nwalkers'                : 64,
        'nsteps'                  : 1000,
        'ntemps'                  : 10,
        'threads'                 : 1,
        'nparallel'               : 1,

        # Plots
    }

    # Read options from config file.
    for key in input_pars.keys():

        # FIXME: add all the options here.
        pass

        # Input
        if (key == 'output') or (key == 'PSD-path'):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if (key == 'screen-output'):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass

        # Model
        if (key == 'event-parameters') or (key == 'priors'):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass
        if (key == 'observing-run') or (key == 'waveform'): 
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if (key == 'reference-frequency') or (key == 'sampling-frequency'):
            try: input_pars[key] = Config.getfloat('model', key)
            except: pass

        # Sampler
        if (key == 'sampler') or (key == 'print-method'):
            try: input_pars[key] = Config.get('sampler', key)
            except: pass
        if (key == 'nparallel') or (key == 'nlive') or (key == 'queue-size') or (key == 'nwalkers') or (key == 'nsteps') or (key == 'ntemps') or (key == 'threads'):
            try: input_pars[key] = Config.getint('sampler', key)
            except: pass


        # Plots

    
    # Initialise the event parameters.
    input_pars['event-parameters-default'] = default_event_parameters()
    if not input_pars['event-parameters'] == {}:
        for key in input_pars['event-parameters']: input_pars['event-parameters-default'][key] = input_pars['event-parameters'][key]
    else:
        input_pars['event-parameters'] = input_pars['event-parameters-default']

    # Initialise the PE priors.
    input_pars['priors-default'] = default_PE_priors()
    if not input_pars['priors'] == {}:
        for key in input_pars['priors']: input_pars['priors-default'][key] = input_pars['priors'][key]
    else:
        input_pars['priors'] = input_pars['priors-default']

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
        'phi_12'             : 0.0,
        'phi_jl'             : 0.0,
        'luminosity_distance': 2000.0,
        'theta_jn'           : 0.0,
        'psi'                : 0.0,
        'phase'              : 0.0,
        'geocent_time'       : 1126259642.0,
        'ra'                 : 0.0,
        'dec'                : 0.0,
    }

    return dict

# FIXME decide what default priors to choose
def default_PE_priors():

    dict = {
        'mass_1'             : [],
        'mass_2'             : [],
        'a_1'                : [],
        'a_2'                : [],
        'tilt_1'             : [],
        'tilt_2'             : [],
        'phi_12'             : [],
        'phi_jl'             : [],
        'luminosity_distance': [],
        'theta_jn'           : [],
        'psi'                : [],
        'phase'              : [],
        'geocent_time'       : [],
        'ra'                 : [],
        'dec'                : [],
    }

    return dict


usage = """

\nWelcome to the Bilby pipeline helper!\n

    # ----- #
    # input #
    # ----- #

        output                      Path where the event PE output is saved. Default: 'parameter_estimation'.
        screen-output               Flag to deviate the standard output to screen. Default: 0.

        PSD-path                    Path to the directory where PSD files are stored. PSD files should be named of the form '{ifo}_{observing_run}.txt'. Default: ''.

    # ----- #
    # model #
    # ----- #

        event-parameters            Dictionary containing the event parameters. Default: {}.
        priors                      Prior ranges to use for event parameters inference. Default: {}.
        observing-run               Detector sensitivity used to compute the SNR with Bilby. Options: 'O3', 'O4', 'O5'. Default: 'O3'.
        waveform                    Waveform model used to compute the SNR with Bilby. Default: 'IMRPhenomXHM'.
        reference-frequency         Frequency at which the binary parameters are defined when Bilby generates the waveforms. Default: 20.
        sampling-frequency          Sampling rate used to generate the waveform with Bilby. Default: 2048.

    # ------- #
    # sampler #
    # ------- #

        sampler                     Sampler to be used to draw samples from the likelihood. The samplers are called from the Bilby package (https://pypi.org/project/bilby/). Options: 'dynesty', 'nessai', 'ptemcee'. Default: 'dynesty'.

        nlive                       Number of live points used by the nested sampler. Option only available for Nested Samplers. Default: 500.
        print-method                Method for printing the sampler output. Dynesty uses a tqdm bar by default, otherwise passing 'interval-$TIME' it prints to sdtout every $TIME seconds. Default: 'interval-60'.
        sample                      Methods to perform the MCMC evolution to find a new point with a nested sampler. Option only available for Nested Samplers. More information on the different methods can be found in the related Bilby documentation (https://bilby-dev.github.io/bilby/dynesty-guide.html). Options: 'act-walk', 'acceptance-walk', 'rwalk'. Default: 'acceptance-walk'.
        naccept                     The length of the MCMC chains during the run follows a Poisson distribution with mean naccept. Option only available for Nested Samplers and only applies to the sample method 'acceptance-walk'. Default: 60.
        queue-size                  Number of parallel process to be executed (see dynesty documentation: https://dynesty.readthedocs.io/en/stable/quickstart.html#parallelization). It corresponds to the number of threads used. Default: 1.

        nwalkers                    Number of parallel chains (walkers) running in the MCMC ensemble. Option only available for MCMC samplers. Default: 64.
        nsteps                      Number of steps taken by each walker in the MCMC samplers. Option only available for MCMC samplers. Default: 1000.
        ntemps                      Number of parallel-tempered chains of the MCMC sampler. Option only available for MCMC samplers. Default: 10.
        threads                     Number of CPU threads used for parallel computation. Option only available for MCMC samplers. Default: 1.
        nparallel                   Number of likelihood evaluations performed simultaneously. While 'threads' distributes MCMC steps across CPU threads, 'nparallel' parallelizes across multiple processes. Option only available for MCMC samplers. Default: 1.

    # ----- #
    # plots #
    # ----- #

"""