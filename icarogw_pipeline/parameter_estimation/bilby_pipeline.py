import os, sys, shutil, configparser
import bilby, numpy as np
from optparse import OptionParser

# Internal imports
import initialise
sys.path.append('../simulations')
import snr_computation as snr
sys.path.append('../')
import icarogw_runner as icarorun



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

    # --------------------------------- #
    # Generate population or injections #
    # --------------------------------- #

    # Read the event parameters.
    injection_parameters = {
        'mass_1'             : input_pars['event-parameters']['mass_1'],
        'mass_2'             : input_pars['event-parameters']['mass_2'],
        'a_1'                : input_pars['event-parameters']['a_1'],
        'a_2'                : input_pars['event-parameters']['a_2'],
        'tilt_1'             : input_pars['event-parameters']['tilt_1'],
        'tilt_2'             : input_pars['event-parameters']['tilt_2'],
        'phi_12'             : input_pars['event-parameters']['phi_12'],
        'phi_jl'             : input_pars['event-parameters']['phi_jl'],
        'luminosity_distance': input_pars['event-parameters']['luminosity_distance'],
        'theta_jn'           : input_pars['event-parameters']['theta_jn'],
        'psi'                : input_pars['event-parameters']['psi'],
        'phase'              : input_pars['event-parameters']['phase'],
        'geocent_time'       : input_pars['event-parameters']['geocent_time'],
        'ra'                 : input_pars['event-parameters']['ra'],
        'dec'                : input_pars['event-parameters']['dec'],
        'ifos_on'            : input_pars['event-parameters']['detectors'],
        'seed'               : input_pars['event-parameters']['seed'],
    }

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
    # FIXME: Need to be implemented. We can use the implemenetation in the runner: https://github.com/vascogennari/icarogw_pipeline/blob/main/icarogw_pipeline/icarogw_runner.py#L427
    # Default priors are here: input_pars['priors']
    priors = bilby.gw.prior.BBHPriorDict()

    # FIXME: Do we need this? Not necessarily, especially since we are not sure H1 is the firrst
    time_delay = BilbyClass.ifos_list[0].time_delay_from_geocenter(
        injection_parameters['ra'],
        injection_parameters['dec'],
        injection_parameters['geocent_time'],
    )
    priors['H1_time'] = bilby.core.prior.Uniform(
        minimum     = injection_parameters['geocent_time'] + time_delay - 0.1,
        maximum     = injection_parameters['geocent_time'] + time_delay + 0.1,
        name        = 'H1_time',
        latex_label = '$t_H$',
        unit        = '$s$',
    )

    # Initialise the likelihood.
    likelihood = bilby.gw.GravitationalWaveTransient(
        interferometers          = BilbyClass.ifos_list,
        waveform_generator       = BilbyClass.waveform_generator,
        priors                   = priors,
        distance_marginalization = True,
        phase_marginalization    = False,
        time_marginalization     = False,
        reference_frame          = 'H1L1',
        time_reference           = 'H1',
    )

    # Run the sampler.
    result = bilby.run_sampler(
        likelihood           = likelihood,
        priors               = priors,
        sampler              = input_pars['sampler'],
        nlive                = input_pars['nlive'],
        npool                = input_pars['npool'],
        outdir               = input_pars['output'],
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