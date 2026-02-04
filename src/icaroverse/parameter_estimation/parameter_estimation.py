import os, sys, shutil, configparser, pickle, bilby
from optparse import OptionParser

# Internal imports
from . import initialise
from .generate_pe_configs import chirp_mass
from ..generate_events import snr_computation as snr
from ..hierarchical_inference import hierarchical_inference

latex_labels = {
    'mass_1': r"$m_1$", 
    'mass_2': r"$m_2$", 
    'chirp_mass': r"$\mathcal{M}_c$", 
    'mass_ratio': r"$q$", 
    'luminosity_distance': r"$d_L$", 
    'ra': r"$\alpha$", 
    'psi': r"$\psi$", 
    'phase': r"$\phi$", 
    'dec': r"$\delta$", 
    'theta_jn': r"$\theta_{JN}$", 
    'a_1': r"$a_1$", 
    'a_2': r"$a_2$", 
    'tilt_1': r"$\theta_1$", 
    'tilt_2': r"$\theta_2$", 
    'phi_12': r"$\phi_{12}$", 
    'phi_jl': r"$\phi_{JL}$", 
    'chi_1': r"$\chi_1$", 
    'chi_2': r"$\chi_2$", 
    'geocent_time': r"$t_{geocent}$", 
}


def set_mass_constraint(params):
    converted_params = params.copy()
    converted_params['mass_ratio'] = params['mass_2'] / params['mass_1']
    return converted_params


def initialise_prior(dict_in, dict_out, input_pars, trigtime=None, precession=False):

    # Masses 
    valid_mass_prior_options = False
    if input_pars['mass-parameters-uniform-prior'] == 'm1-m2':

        if   type(dict_in['mass_1']) == list:  dict_out['mass_1'] = bilby.core.prior.Uniform(dict_in['mass_1'][0], dict_in['mass_1'][1], latex_label=latex_labels['mass_1'])
        elif type(dict_in['mass_1']) == float: dict_out['mass_1'] = dict_in['mass_1']
        else:  raise ValueError(f"Unknown type for prior on mass_1: {dict_in['mass_1']}. Please provide either a 2-list for prior bounds, or a float to fix mass_1.")

        if   input_pars['mass-parameters-sampled'] == 'm1-m2_bilby':
            valid_mass_prior_options = True

            if   type(dict_in['mass_2']) == list:  
                dict_out['mass_2'] = bilby.core.prior.ConditionalUniform(
                    condition_func = bilby.gw.prior.secondary_mass_condition_function,
                    minimum        = dict_in['mass_2'][0], 
                    maximum        = dict_in['mass_2'][1],
                    latex_label=latex_labels['mass_2'],
                )
            elif type(dict_in['mass_2']) == float: dict_out['mass_2'] = dict_in['mass_2']
            else: raise ValueError(f"Unknown type for prior on mass_2: {dict_in['mass_2']}. Please provide either a 2-list for prior bounds, or a float to fix mass_2.")
            # Clean priors on other mass parameters that could have been initialized if dict_out was a bilby.gw.priors.BBHPriorDict object.
            snr.clean_dict(dict_out, ['chirp_mass', 'mass_ratio'])
            
        elif input_pars['mass-parameters-sampled'] == 'm1-m2_custom':
            valid_mass_prior_options = True

            if   type(dict_in['mass_2']) == list:  dict_out['mass_2'] = bilby.core.prior.Uniform(dict_in['mass_2'][0], dict_in['mass_2'][1], latex_label=latex_labels['mass_2'])
            elif type(dict_in['mass_2']) == float: dict_out['mass_2'] = dict_in['mass_2']
            else: raise ValueError(f"Unknown type for prior on mass_2: {dict_in['mass_2']}. Please provide either a 2-list for prior bounds, or a float to fix mass_2.")

            if   type(dict_in['mass_ratio']) == list and len(dict_in['mass_ratio']) == 2:  dict_out['mass_ratio'] = bilby.core.prior.Constraint(minimum = dict_in['mass_ratio'][0], maximum = dict_in['mass_ratio'][1])
            else: raise ValueError(f"Unknown type for constraints bounds on mass_ratio: {dict_in['mass_ratio']}. Please provide a 2-list for prior bounds.")
            # Clean priors on other mass parameters that could have been initialized if dict_out was a bilby.gw.priors.BBHPriorDict object.
            snr.clean_dict(dict_out, ['chirp_mass'])
        
        elif input_pars['mass-parameters-sampled'] == 'Mc-q':
            valid_mass_prior_options = True

            if   type(dict_in['chirp_mass']) == list:  dict_out['chirp_mass'] = bilby.gw.prior.UniformInComponentsChirpMass(minimum=dict_in['chirp_mass'][0], maximum=dict_in['chirp_mass'][1], latex_label=latex_labels['chirp_mass'])
            elif type(dict_in['chirp_mass']) == float: dict_out['chirp_mass'] = dict_in['chirp_mass']
            else:  raise ValueError(f"Unknown type for prior on chirp_mass: {dict_in['chirp_mass']}. Please provide either a 2-list for prior bounds, or a float to fix chirp_mass.")

            if   type(dict_in['mass_ratio']) == list:  dict_out['mass_ratio'] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=dict_in['mass_ratio'][0], maximum=dict_in['mass_ratio'][1], latex_label=latex_labels['mass_ratio'])
            elif type(dict_in['mass_ratio']) == float: dict_out['mass_ratio'] = dict_in['mass_ratio']
            else:  raise ValueError(f"Unknown type for prior on mass_ratio: {dict_in['mass_ratio']}. Please provide either a 2-list for prior bounds, or a float to fix mass_ratio.")
            # Clean priors on other mass parameters that could have been initialized if dict_out was a bilby.gw.priors.BBHPriorDict object.
            snr.clean_dict(dict_out, ['mass_1', 'mass_2'])

        else:
            pass

    elif input_pars['mass-parameters-uniform-prior'] == 'Mc-q' and input_pars['mass-parameters-sampled'] == 'Mc-q':
        valid_mass_prior_options = True

        if   type(dict_in['chirp_mass']) == list:  dict_out['chirp_mass'] = bilby.core.prior.Uniform(dict_in['chirp_mass'][0], dict_in['chirp_mass'][1], latex_label=latex_labels['chirp_mass'])
        elif type(dict_in['chirp_mass']) == float: dict_out['chirp_mass'] = dict_in['chirp_mass']
        else:  raise ValueError(f"Unknown type for prior on chirp_mass: {dict_in['chirp_mass']}. Please provide either a 2-list for prior bounds, or a float to fix chirp_mass.")

        if   type(dict_in['mass_ratio']) == list:  dict_out['mass_ratio'] = bilby.core.prior.Uniform(dict_in['mass_ratio'][0], dict_in['mass_ratio'][1], latex_label=latex_labels['mass_ratio'])
        elif type(dict_in['mass_ratio']) == float: dict_out['mass_ratio'] = dict_in['mass_ratio']
        else:  raise ValueError(f"Unknown type for prior on mass_ratio: {dict_in['mass_ratio']}. Please provide either a 2-list for prior bounds, or a float to fix mass_ratio.")

    if not valid_mass_prior_options:
        raise ValueError("Unknown combination of mass_parameters_uniform_prior and mass_parameters_sampled. Possible combinations:\n\tmass_parameters_uniform_prior = 'm1-m2', mass_parameters_sampled = 'm1-m2_bilby' \n\tmass_parameters_uniform_prior = 'm1-m2', mass_parameters_sampled = 'm1-m2_custom' \n\tmass_parameters_uniform_prior = 'm1-m2', mass_parameters_sampled = 'Mc-q' \n\tmass_parameters_uniform_prior = 'Mc-q',  mass_parameters_sampled = 'Mc-q' ")
    else:
        pass

    # Luminosity distance
    dict_in_dL = dict_in['luminosity_distance']
    if   type(dict_in_dL) == list and len(dict_in_dL) > 2:
        if dict_in_dL[-1] == "PowerLaw":
            dict_out['luminosity_distance'] = bilby.core.prior.PowerLaw(2., dict_in_dL[0], dict_in_dL[1], latex_label=latex_labels['luminosity_distance'])
        elif dict_in_dL[-1] == "UniformSourceFrame":
            dict_out['luminosity_distance'] = bilby.gw.prior.UniformSourceFrame(dict_in_dL[0], dict_in_dL[1], name="luminosity_distance", latex_label=latex_labels['luminosity_distance'])
        elif dict_in_dL[-1] == "UniformComovingVolume":
            dict_out['luminosity_distance'] = bilby.gw.prior.UniformComovingVolume(dict_in_dL[0], dict_in_dL[1], name="luminosity_distance", latex_label=latex_labels['luminosity_distance'])
        elif dict_in_dL[-1] == "Uniform":
            dict_out['luminosity_distance'] = bilby.gw.prior.Uniform(dict_in_dL[0], dict_in_dL[1], latex_label=latex_labels['luminosity_distance'])
        else:
            raise ValueError(f"Unknown prior on luminosity_distance: {dict_in_dL[-1]}. Please choose between [Uniform, PowerLaw, UniformSourceFrame, UniformComovingVolume].")
    elif   type(dict_in_dL) == list and len(dict_in_dL) == 2:  
        dict_out['luminosity_distance'] = bilby.core.prior.Uniform(dict_in_dL[0], dict_in_dL[1], latex_label=latex_labels['luminosity_distance'])
    elif type(dict_in_dL) == float: 
        dict_out['luminosity_distance'] = dict_in_dL
    else:  
        raise ValueError(f"Unknown type for prior on luminosity_distance: {dict_in_dL}. Please provide either a 2-list for prior bounds, or a float to fix luminosity_distance.")

    # Inclination angle
    if   type(dict_in['theta_jn']) == list:  dict_out['theta_jn'] = bilby.core.prior.analytical.Sine(minimum = dict_in['theta_jn'][0], maximum = dict_in['theta_jn'][1], latex_label=latex_labels['theta_jn'])
    elif type(dict_in['theta_jn']) == float: dict_out['theta_jn'] = dict_in['theta_jn']
    else: raise ValueError(f"Unknown type for prior on theta_jn: {dict_in['theta_jn']}. Please provide either a 2-list for prior bounds, or a float to fix theta_jn.")

    # Sky localization & phase & polarization
    for par in ['psi', 'phase', 'ra']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1], boundary='periodic', latex_label=latex_labels[par])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")

    if   type(dict_in['dec']) == list:  dict_out['dec'] = bilby.core.prior.analytical.Cosine(minimum = dict_in['dec'][0], maximum = dict_in['dec'][1], latex_label=latex_labels['dec'])
    elif type(dict_in['dec']) == float: dict_out['dec'] = dict_in['dec']
    else: raise ValueError(f"Unknown type for prior on dec: {dict_in['dec']}. Please provide either a 2-list for prior bounds, or a float to fix dec.")

    # Spins
    if precession:
        for par in ['a_1', 'a_2']:
            if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1], latex_label=latex_labels[par])
            elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
            else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")
        
        for par in ['tilt_1', 'tilt_2']:
            if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.analytical.Sine(minimum = dict_in[par][0], maximum = dict_in[par][1], latex_label=latex_labels[par])
            elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
            else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")
        
        for par in ['phi_12', 'phi_jl']:
            if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1], boundary='periodic', latex_label=latex_labels[par])
            elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
            else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")
        
        snr.clean_dict(dict_out, ['chi_1', 'chi_2'])
    
    else:
        for par in ['chi_1', 'chi_2']:
            if   type(dict_in[par]) == list:  
                dict_out[par] = bilby.gw.prior.AlignedSpin(
                    a_prior = bilby.core.prior.Uniform(minimum = 0.,  maximum = max(map(abs, dict_in[par]))), 
                    z_prior = bilby.core.prior.Uniform(minimum = -1., maximum = 1.),
                    minimum = dict_in[par][0], maximum = dict_in[par][1],
                    latex_label=latex_labels[par],
                )
            elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
            else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")
        
        snr.clean_dict(dict_out, ['a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])

    # Geocentric time
    if 'geocent_time' in dict_in and type(dict_in['geocent_time']) == float: 
        dict_out['geocent_time'] = dict_in['geocent_time']
    elif 'geocent_time' in dict_in: 
        raise ValueError(f"Unknown type for prior on dec: {dict_in['geocent_time']}. Please provide either a 2-list for prior bounds, or a float to fix dec.")
    else: 
        dict_out['geocent_time'] = bilby.core.prior.Uniform(trigtime - 0.1, trigtime + 0.1, latex_label=latex_labels['geocent_time'])

    print('\n * Using the following priors.\n')
    for key in dict_out.keys():
        max_len = len(max(dict_out.keys(), key = len))
        print('\t{}  {}'.format(key.ljust(max_len), dict_out[key]))
    
    if input_pars['mass-parameters-sampled'] == 'm1-m2_custom':
        dict_out = bilby.core.prior.PriorDict(dict_out, conversion_function=set_mass_constraint)
    elif input_pars['mass-parameters-sampled'] == 'Mc-q_custom' or input_pars['mass-parameters-uniform-prior'] == 'Mc-q':
        dict_out = bilby.core.prior.PriorDict(dict_out)

    return dict_out


######################
# Main function call #
######################

def main():

    # ------------------------------------------------ #
    # Read config file and initialise input parameters #
    # ------------------------------------------------ #

    parser = OptionParser(initialise.usage)
    parser.add_option(      '--config-file', type='string', metavar = 'config_file', default = None)
    parser.add_option('-n', '--n-processes', type='int',    metavar = 'n_processes', default = -1, help="Set the number of processes for parallelized injections generation from command line, if this should match some external structure (e.g. number of CPUs allocated to the simulation on a computing cluster job.)")
    (opts, _) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    # Initialise input parameters dictionary.
    input_pars = initialise.InitialiseOptions(Config)

    # Set output directory.
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])

    # Copy config file to output.
    try:    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))
    except: pass # Config file already copied.

    # Set the number of parallel processes according to command line if provided to match hardware structure
    if opts.n_processes > 0:
        print(f"\n * Number of processes set via command-line option: n_processes = {opts.n_processes} \n")
        input_pars['npool'] = opts.n_processes

    # Deviate stdout and stderr to file.
    if not input_pars['screen-output']:
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_bilby.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_bilby.txt'), 'w')
    else: pass
    print('\n\n Starting  b i l b y  pipeline\n\n')

    # Print run parameters.
    print(' * I will be running with the following parameters.\n')
    hierarchical_inference.print_dictionary(input_pars)

    # Read the event parameters. (We just need to rename masses and distance to match the waveform generator's syntax)
    injection_parameters = input_pars['event-parameters']

    if input_pars['event-from-simulation']:
        injection_parameters['mass_1'             ] = injection_parameters['m1d']
        injection_parameters['mass_2'             ] = injection_parameters['m2d']
        injection_parameters['luminosity_distance'] = injection_parameters['dL' ]
    
    strain_record = {}
    if input_pars['strain-file'] != "":
        with open(input_pars['strain-file'], 'rb') as f:
            strain_record = pickle.load(f)
        strain_data = strain_record.pop('strain')

    # Initialise the Bilby class.
    BilbyClass = snr.BilbyDetectionPipeline(
        observing_run       = input_pars['observing-run'], 
        reference_frequency = input_pars['reference-frequency'], 
        sampling_frequency  = input_pars['sampling-frequency'], 
        approximant         = input_pars['waveform'],
        precessing_apx      = input_pars['precession'],
        **strain_record # Contains technical data such as duration and start_time
    )
    
    # Set the events and interferometers within Bilby.
    if input_pars['strain-file'] != "":
        BilbyClass.set_event_dict(injection_parameters, flag_set_duration_and_start_time=False)
        BilbyClass.set_ifos_list()
        BilbyClass.set_waveform_generator()
        BilbyClass.set_strain_data_from_arrays(strain_data)
    else:
        BilbyClass.set_event_dict(injection_parameters, flag_set_duration_and_start_time=True)
        BilbyClass.set_ifos_list()
        BilbyClass.set_waveform_generator()
        BilbyClass.set_strain_data_from_psd()
    
    BilbyClass.inject_signal()

    # Initialise priors.
    if   input_pars['mass-parameters-sampled'] == 'm1-m2_bilby':
        priors = bilby.gw.prior.BBHPriorDict()
        priors = initialise_prior(
            input_pars['priors'], 
            priors, 
            input_pars,
            trigtime   = injection_parameters['geocent_time'], 
            precession = input_pars['precession']
        )
    elif input_pars['mass-parameters-sampled'] == 'm1-m2_custom' or input_pars['mass-parameters-sampled'] == 'Mc-q':
        priors = {}
        priors = initialise_prior(
            input_pars['priors'], 
            priors, 
            input_pars,
            trigtime   = injection_parameters['geocent_time'], 
            precession = input_pars['precession']
        )

    # Initialise the likelihood.
    optional_likelihood_arguments = {}
    if input_pars["sampler"] == "nessai" and (input_pars['time-marginalization'] or input_pars['distance-marginalization']):
        optional_likelihood_arguments["distance_marginalization_lookup_table"] = True
        optional_likelihood_arguments["jitter_time"] = False
    else:
        pass

    likelihood = bilby.gw.GravitationalWaveTransient(
        interferometers             = BilbyClass.ifos_list,
        waveform_generator          = BilbyClass.waveform_generator,
        priors                      = priors,
        phase_marginalization       = input_pars['phase-marginalization'],
        time_marginalization        = input_pars['time-marginalization'],
        distance_marginalization    = input_pars['distance-marginalization'],
        calibration_marginalization = input_pars['calibration-marginalization'],
        **optional_likelihood_arguments
    )

    # Run the sampler.
    print("\n * Starting the sampler...")
    result = bilby.run_sampler(
        likelihood           = likelihood,
        priors               = priors,
        outdir               = os.path.join(input_pars['output'], "sampler"),
        sampler              = input_pars['sampler'],
        nlive                = input_pars['nlive'],
        npool                = input_pars['npool'],
        naccept              = input_pars['naccept'],
        sample               = input_pars['sample'],
        injection_parameters = BilbyClass.projected_event_dict,
        conversion_function  = bilby.gw.conversion.generate_all_bbh_parameters,
        result_class         = bilby.gw.result.CBCResult,
    )
    print("\t...finished sampling !\n")

    # Plot the inferred waveform superposed on the actual data.
    print("\n * Plotting the waveforms in each detector...")
    result.plot_waveform_posterior(n_samples = 1000)
    print("\t...done !\n")

    # Make a corner plot.
    print("\n * Producing the corner plot...")
    if (input_pars['mass-parameters-uniform-prior'] == 'Mc-q') or (input_pars['mass-parameters-uniform-prior'] == 'm1-m2' and input_pars['mass-parameters-sampled'] == 'Mc-q'):
        pars_to_plot = ['mass_1', 'mass_2'] + [key for key in result.search_parameter_keys if key not in {'chirp_mass', 'mass_ratio'}] + ['chirp_mass', 'mass_ratio']
        BilbyClass.projected_event_dict['chirp_mass'] = chirp_mass(m1=BilbyClass.projected_event_dict['mass_1'], m2=BilbyClass.projected_event_dict['mass_2'])
        BilbyClass.projected_event_dict['mass_ratio'] = BilbyClass.projected_event_dict['mass_2'] / BilbyClass.projected_event_dict['mass_1']
    elif input_pars['mass-parameters-uniform-prior'] == 'm1-m2':
        pars_to_plot = result.search_parameter_keys
    result.plot_corner(
        # priors = True,
        truths = {key: BilbyClass.projected_event_dict[key] for key in pars_to_plot},
    )

    print("\t...done !\n")
