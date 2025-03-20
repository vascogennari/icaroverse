import ast


def InitialiseOptions(Config):

    # Dictionary with the default options.
    input_pars = {

        # Input
        'output'             : 'output',
        'duration'           : 4.,
        'sampling-rate'      : 1024.,
        'approximant'        : 'IMRPhenomXPHM',
        'reference-frequency': 50.,
        'minimum-frequency'  : 20.,
        'detectors'          : ['H1', 'L1', 'V1'],

        # Sampler
        'sampler'            : 'nessai',
        'nlive'              : 500,
        'naccept'            : 60,
        'queue-size'         : 1,
        'print-method'       : 'interval-60',

        # Plots
    }

    # Read options from config file.
    for key in input_pars.keys():

        # FIXME: add all the options here.
        pass

        # Input


        # Sampler


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


    # ------- #
    # sampler #
    # ------- #


    # ----- #
    # plots #
    # ----- #

"""