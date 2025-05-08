import os, sys, shutil, configparser
import bilby, numpy as np
from optparse import OptionParser

# Internal imports
import initialise
sys.path.append('../simulations')
import snr_computation as snr
sys.path.append('../')
import icarogw_runner as icarorun

latex_labels = {
    'mass_1': r"$m_1$", 
    'mass_2': r"$m_2$", 
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


def initialise_prior(dict_in, dict_out, trigtime=None, precession=False):

    for par in ['mass_1', 'luminosity_distance']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1], latex_label=latex_labels[par])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else:  raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")

    for par in ['ra', 'psi', 'phase']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1], boundary='periodic', latex_label=latex_labels[par])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")

    for par in ['dec']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.analytical.Cosine(minimum = dict_in[par][0], maximum = dict_in[par][1], latex_label=latex_labels[par])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")
        
    for par in ['theta_jn']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.analytical.Sine(minimum = dict_in[par][0], maximum = dict_in[par][1], latex_label=latex_labels[par])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else: raise ValueError(f"Unknown type for prior on {par}: {dict_in[par]}. Please provide either a 2-list for prior bounds, or a float to fix {par}.")

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

    dict_out['mass_2'] = bilby.core.prior.ConditionalUniform(
        condition_func = bilby.gw.prior.secondary_mass_condition_function,
        minimum        = dict_in['mass_2'][0], 
        maximum        = dict_in['mass_2'][1],
        latex_label=latex_labels['mass_2'],
    )
    
    dict_out['geocent_time'] = bilby.core.prior.Uniform(trigtime - 0.1, trigtime + 0.1, latex_label=latex_labels['geocent_time'])

    snr.clean_dict(dict_out, ['chirp_mass', 'mass_ratio'])

    print('\n * Using the following priors.\n')
    for key in dict_out.keys():
        max_len = len(max(dict_out.keys(), key = len))
        print('\t{}  {}'.format(key.ljust(max_len), dict_out[key]))

    return dict_out


######################
# Main function call #
######################

def main():

    # ------------------------------------------------ #
    # Read config file and initialise input parameters #
    # ------------------------------------------------ #

    parser = OptionParser(initialise.usage)
    parser.add_option('--config-file', type='string', metavar = 'config_file', default = None)
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
    if not os.path.exists(os.path.join(input_pars['output'], 'plots')): os.makedirs(os.path.join(input_pars['output'], 'plots'))

    # Copy config file to output.
    try:    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))
    except: pass # Config file already copied.

    # Deviate stdout and stderr to file.
    if not input_pars['screen-output']:
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_bilby.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_bilby.txt'), 'w')
    else: pass
    print('\n\n Starting  b i l b y  pipeline\n\n')

    # Print run parameters.
    print(' * I will be running with the following parameters.\n')
    icarorun.print_dictionary(input_pars)

    # Read the event parameters. (We just need to rename masses and distance to match the waveform generator's syntax)
    injection_parameters = input_pars['event-parameters']

    if input_pars['event-from-simulation']:
        injection_parameters['mass_1'             ] = injection_parameters['m1d']
        injection_parameters['mass_2'             ] = injection_parameters['m2d']
        injection_parameters['luminosity_distance'] = injection_parameters['dL' ]

    # Initialise the Bilby class.
    BilbyClass = snr.BilbyDetectionPipeline(
        psd_dir             = input_pars['PSD-path'], 
        observing_run       = input_pars['observing-run'], 
        reference_frequency = input_pars['reference-frequency'], 
        sampling_frequency  = input_pars['sampling-frequency'], 
        approximant         = input_pars['waveform'],
        precessing_apx      = input_pars['precession']
    )
    
    # Set the events and interferometers within Bilby.
    BilbyClass.set_event_dict(injection_parameters)
    BilbyClass.set_ifos_and_inject_signal()

    # Initialise priors.
    priors = bilby.core.prior.BBHPriorDict()
    priors = initialise_prior(
        input_pars['priors'], 
        priors, 
        trigtime = injection_parameters['geocent_time'], 
        precession=input_pars['precession']
    )

    # Initialise the likelihood.
    likelihood = bilby.gw.GravitationalWaveTransient(
        interferometers       = BilbyClass.ifos_list,
        waveform_generator    = BilbyClass.waveform_generator,
        priors                = priors,
        phase_marginalization = True,
    )

    # Run the sampler.
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

    # Plot the inferred waveform superposed on the actual data.
    result.plot_waveform_posterior(n_samples = 1000)

    # Make a corner plot.
    result.plot_corner(
        priors = True,
        truths = {key: BilbyClass.projected_event_dict[key] for key in result.search_parameter_keys}
    )



# Execute the main function.
if __name__=='__main__': main()