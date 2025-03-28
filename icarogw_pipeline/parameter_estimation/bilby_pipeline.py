import os, sys, shutil, configparser
import bilby, numpy as np
from optparse import OptionParser

# Internal imports
import initialise
sys.path.append('../simulations')
import snr_computation as snr
sys.path.append('../')
import icarogw_runner as icarorun


def initialise_prior(dict_in, dict_out, trigtime = None):

    for par in ['mass_1', 'luminosity_distance', 'a_1', 'a_2']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else:
            raise ValueError('Unknown type for prior on {}'.format(dict_in[par]))
    
    for par in 'phi_12', 'phi_jl', 'ra', 'psi', 'phase':
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1], boundary='periodic')
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else:
            raise ValueError('Unknown type for prior on {}'.format(dict_in[par]))

    for par in 'tilt_1', 'tilt_2':
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.analytical.Sine(minimum = dict_in[par][0], maximum = dict_in[par][1])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else:
            raise ValueError('Unknown type for prior on {}'.format(dict_in[par]))
        
    for par in ['dec', 'theta_jn']:
        if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.analytical.Cosine(minimum = dict_in[par][0], maximum = dict_in[par][1])
        elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
        else:
            raise ValueError('Unknown type for prior on {}'.format(dict_in[par]))
        
    dict_out['mass_2'] = bilby.core.prior.ConditionalUniform(
        condition_func = bilby.gw.prior.secondary_mass_condition_function,
        minimum = dict_in['mass_2'][0], maximum = dict_in['mass_2'][1])
    
    dict_out['geocent_time'] = bilby.core.prior.Uniform(trigtime - 0.1, trigtime + 0.1)

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
    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))

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
        approximant         = input_pars['waveform']
    )
    
    # Set the events and interferometers within Bilby.
    BilbyClass.set_event_dict(injection_parameters)
    BilbyClass.set_ifos_and_inject_signal()

    # Initialise priors.
    #priors = bilby.gw.prior.BBHPriorDict()
    priors = bilby.core.prior.PriorDict()
    priors = initialise_prior(input_pars['priors'], priors, trigtime = injection_parameters['geocent_time'])

    # Initialise the likelihood.
    likelihood = bilby.gw.GravitationalWaveTransient(
        interferometers    = BilbyClass.ifos_list,
        waveform_generator = BilbyClass.waveform_generator,
        priors             = priors,
    )

    # Run the sampler.
    result = bilby.run_sampler(
        likelihood           = likelihood,
        priors               = priors,
        sampler              = input_pars['sampler'],
        nlive                = input_pars['nlive'],
        outdir               = input_pars['output'],
        n_pool               = 1,
        injection_parameters = injection_parameters,
        conversion_function  = bilby.gw.conversion.generate_all_bbh_parameters,
        result_class         = bilby.gw.result.CBCResult,
    )

    # Plot the inferred waveform superposed on the actual data.
    result.plot_waveform_posterior(n_samples = 1000)

    # Make a corner plot.
    result.plot_corner()



# Execute the main function.
if __name__=='__main__': main()