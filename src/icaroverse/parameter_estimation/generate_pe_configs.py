import os
from argparse import ArgumentParser
from tqdm import tqdm

import pickle, pandas as pd, json
import numpy as np

# Internal imports
from . import initialise


# ----------------------------------------------------------------------------------------- #
# Read a population folder, and generate the corresponding subfolders & config files for PE #
# ----------------------------------------------------------------------------------------- #

def chirp_mass(m1, m2):
    """Computes the binary chirp mass from individual components masses."""
    return np.power(m1*m2, 3/5) * np.power(m1+m2, -1/5)

def primary_mass_expected_std(mc, q, snr):
    """
    Compute an estimate of the uncertainty on m1 based on SNR, 
    estimates of uncertainties on Mc and q from [Farr et al. (2016), arXiv:1508.05336] (fig.2)
    and the error propagation formula (https://en.wikipedia.org/wiki/Propagation_of_uncertainty)
    
    Parameters
    ----------
    mc: float or array-like 
        chirp mass
    q: float or array-like
        mass ratio
    snr: float or array-like
        signal to noise ratio
    
    Return
    ------
    m1_std: float or array-like
        error estimate
    """
    return (mc / snr) * np.power(q**3 * (1+q)**4, -1/5) * np.sqrt( 1e-4 * np.power(1+q, 2) + (1.21/25) * np.power(3 + 2*q, 2) * snr )

def linear_tailored_prior_component_masses_upper_bound(m1d):
    """
    Linear adaptive prior (upper bound)
    Interpolates between:
    prior = [1,  60] Msun for m1d = 10 Msun event
        and
    prior = [1, 300] Msun for m1d = 90 Msun event

    Parameters
    ----------
    mc: float or array-like 
        primary mass (detector frame)
    
    Return
    ------
    m1d_prior_upper_bound: float or array-like
        prior upper bound
    """
    return 60. + 3. * (m1d - 10.)

def linear_tailored_prior_luminosity_distance_upper_bound(dL):
    """
    Linear adaptive prior (upper bound)
    Interpolates between:
    prior = [0, 1e4] Mpc for dL = 0 Mpc event
        and
    prior = [0, 6e4] Mpc for dL = 2e4 Mpc event

    Parameters
    ----------
    dL: float or array-like 
        luminosity distance
    
    Return
    ------
    dL_prior_upper_bound: float or array-like
        prior upper bound
    """
    return 2.5 * dL + 1e4

def chirp_mass_prior_width(mc, snr):
    """
    Parameters
    ----------
    mc: float or array-like 
        chirp mass (detector frame)
    
    Return
    ------
    width_Mc_prior: float or array-like
        prior width
    """
    return mc**(5/3)/15 * (12/snr)

def true_value_projected_spins(par, event_parameters, precession):
    """
    Return the true value of the parameter from the true values dict event_parameters which contains full polar coordinate 3D spins. If precession is False, then par can be chi_1 or chi_2 and the function returns the projected component of the spin along the z-axis: chi_i = a_i * cos(tilt_i)
    """
    if not precession and par == 'chi_1':
        res = event_parameters['a_1'] * np.cos(event_parameters['tilt_1'])
    elif not precession and par == 'chi_2':
        res = event_parameters['a_2'] * np.cos(event_parameters['tilt_2'])
    else:
        res = event_parameters[par]
    return float(res)




config_template = """[input]
output                        = {output}
screen-output                 = {screen_output}
plot-waveform                 = {plot_waveform}
plot-skymap                   = {plot_skymap}


[model]
event-parameters              = {event_parameters}
strain-file                   = {strain_file}
priors                        = {priors}
mass_parameters_uniform_prior = {mass_parameters_uniform_prior}
mass_parameters_sampled       = {mass_parameters_sampled}
observing-run                 = {observing_run}
waveform                      = {waveform}
precession                    = {precession}
reference-frequency           = {reference_frequency}
sampling-frequency            = {sampling_frequency}
{marginalization_settings}

[sampler]
sampler                       = {sampler}
print-method                  = {print_method}
{sampler_dependent_config}
"""


def main():

    # parser = OptionParser(options.usage)
    parser = ArgumentParser(usage=initialise.usage)
    parser.add_argument('-d', '--pop_dir',                       type=str, default = None             )
    # model options
    parser.add_argument('-m', '--mass_parameters_uniform_prior', type=str, default = 'm1-m2',           help="Which mass parameters have flat priors. Options: m1-m2, Mc-q")
    parser.add_argument(      '--mass_parameters_sampled',       type=str, default = 'm1-m2_custom',    help="Which mass parameters are sampled over. Options: m1-m2_bilby, m1-m2_custom, Mc-q")
    parser.add_argument(      '--distance_prior',                type=str, default = 'Uniform',         help="Distance prior. Options: Uniform, PowerLaw, UniformSourceFrame, UniformComovingVolume")
    parser.add_argument('-t', '--tailor_priors',                 type=str, default = 'linear',          help="Options: linear, farr16, chirp")
    parser.add_argument('-k', '--nsigma',                        type=int, default = 7,                 help="Width of the tailored prior ranges, in units of estimate of m1_std (prior: m1 +/- nsigma * m1_std).")
    parser.add_argument(      '--marginalization',               type=str, default = None,              nargs="+", help="Parameters over which to marginalize the likelihood. Available: phase, time, distance, calibration (can stack multiple options). Default: None")
    parser.add_argument('-r', '--use_recorded_strain',                     default = False,             action="store_true", help="Flag to use saved strain data from detection pipeline, if there is any.")
    parser.add_argument('-f', '--fixed_parameters',              type=str, default = None,              nargs="+", help="Parameters to fix when running PE. Meta-options : all-but-mass-dist, all-but-mass-dist-incl. Default: None")
    parser.add_argument(      '--plot_waveform',                           default = False,             action="store_true", help="Flag to plot waveform when sampling ends.")
    parser.add_argument(      '--plot_skymap',                             default = False,             action="store_true", help="Flag to plot skymap when sampling ends.")
    # sampler
    parser.add_argument('-s', '--sampler',                       type=str, default = 'dynesty'        )
    parser.add_argument(      '--print_method',                  type=str, default = 'interval-60'    )
    # nested samplers options
    parser.add_argument('-l', '--nlive',                         type=int, default = 1000             )
    parser.add_argument('-n', '--npool',                         type=int, default = 10               )
    parser.add_argument('-a', '--naccept',                       type=int, default = 60               )
    parser.add_argument(      '--sample',                        type=str, default = 'acceptance-walk')
    # MCMC samplers options
    parser.add_argument('-q', '--queue_size',                    type=int, default = 1                )
    parser.add_argument('-w', '--nwalkers',                      type=int, default = 64               )
    parser.add_argument(      '--nsteps',                        type=int, default = 1000             )
    parser.add_argument(      '--ntemps',                        type=int, default = 10               )
    parser.add_argument(      '--threads',                       type=int, default = 1                )
    parser.add_argument(      '--nparallel',                     type=int, default = 1                )
    args = parser.parse_args()

    samplers_types = {
        'dynesty': 'nested', 
        'nessai':  'nested', 
        'ptemcee': 'MCMC',
    }

    # Retrieve population and transform it into a list of single event dict
    data_filepath = os.path.join(args.pop_dir, "population_observed.pickle")
    if os.path.exists(data_filepath):
        with open(data_filepath, 'rb') as f: data = pickle.load(f)
        data['ifos_on'] = data['ifos_on'].tolist()
        events_list = pd.DataFrame(data).to_dict('records')
    else:
        raise FileNotFoundError(f"Unknown file {data_filepath}. Please make sure the population directory you passed does contain a 'population_observed.pickle' file.")
    
    strain_records_filepath = os.path.join(args.pop_dir, "strain_records.pickle")
    if os.path.exists(strain_records_filepath):
        with open(strain_records_filepath, 'rb') as f: 
            strain_records = pickle.load(f)
            there_is_strain = True
    else:
        strain_records = [None for _ in events_list]
        there_is_strain = False
        print(f"No strain data found. I will generate config files for parameter estimation runs where strain data will be generated from scratch.")

    # Create 'parameter_estimation' subdirectory in args.pop_dir
    parameter_estimation_subdir = os.path.join(args.pop_dir, "parameter_estimation")
    if not os.path.exists(parameter_estimation_subdir): os.makedirs(parameter_estimation_subdir)

    # Fetch a 'bilby_settings' file where are stored observing-run, waveform, reference-frequency, sampling-frequency from the population generation
    try: 
        print('\n * Reading the population generation settings.\n')
        with open(os.path.join(args.pop_dir, "analysis_settings.json")) as f: population_generation_pars = json.load(f)
    except FileNotFoundError: 
        raise FileNotFoundError("Please make sure the population directory does contain an 'analysis_settings.json' file")
    
    if population_generation_pars['SNR-method'] == 'bilby':
        observing_run       = population_generation_pars['observing-run']
        waveform            = population_generation_pars['snr-bilby-waveform']
        precession          = population_generation_pars['snr-bilby-precessing-wf']
        reference_frequency = population_generation_pars['snr-bilby-reference-frequency']
        sampling_frequency  = population_generation_pars['snr-bilby-sampling-frequency']
    else:
        raise ValueError("Please make sure the simulated population you are trying to work with has been processed with bilby within the present pipeline (for SNR computation)")

    # Fixing some of the parameters for PE.
    if args.fixed_parameters is not None:
        parameters = ['luminosity_distance', 'ra', 'psi', 'phase', 'dec', 'theta_jn', 'geocent_time']
        parameters += ['mass_1', 'mass_2']*(args.mass_parameters_uniform_prior == "m1-m2") + ['chirp_mass', 'mass_ratio']*(args.mass_parameters_uniform_prior == "Mc-q")
        parameters += ['a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']*(precession) + ['chi_1', 'chi_2']*(not precession)
        
        if "all-but-mass-dist" in args.fixed_parameters:
            parameters_to_fix = [par for par in parameters if par not in {'luminosity_distance', 'mass_1', 'mass_2', 'chirp_mass', 'mass_ratio'}]
        elif "all-but-mass-dist-incl" in args.fixed_parameters:
            parameters_to_fix = [par for par in parameters if par not in {'luminosity_distance', 'mass_1', 'mass_2', 'chirp_mass', 'mass_ratio', 'theta_jn'}]
        else:
            parameters_to_fix = [par for par in args.fixed_parameters if par in parameters]
    else:
        parameters_to_fix = []
    print(f"\n * Fixing following parameters: {parameters_to_fix}\n")

    # Marginalization settings
    if args.marginalization is None:
        marginalization_settings = ""
        print(f"\n * Marginalizing the likelihood over: Nothing.\n")
    else:
        available_settings = ["phase", "time", "distance", "calibration"]
        marginalization_settings = "\n".join([f"{ms+'-marginalization'     :<29} = True" for ms in args.marginalization if ms in available_settings])
        print(f"\n * Marginalizing the likelihood over: {', '.join([ms for ms in args.marginalization if ms in available_settings])}.\n")

    print('\n * Generating config files.\n')
    # Checks for validity of parameters
    mass_priors_possible_combinations = [
        ('m1-m2', 'm1-m2_bilby'), 
        ('m1-m2', 'm1-m2_custom'), 
        ('m1-m2', 'Mc-q'), 
        ('Mc-q', 'Mc-q'), 
    ]
    if (args.mass_parameters_uniform_prior, args.mass_parameters_sampled) not in mass_priors_possible_combinations:
        raise ValueError("Unknown combination of mass_parameters_uniform_prior and mass_parameters_sampled. Possible combinations:\n\tmass_parameters_uniform_prior = 'm1-m2', mass_parameters_sampled = 'm1-m2_bilby' \n\tmass_parameters_uniform_prior = 'm1-m2', mass_parameters_sampled = 'm1-m2_custom' \n\tmass_parameters_uniform_prior = 'm1-m2', mass_parameters_sampled = 'Mc-q' \n\tmass_parameters_uniform_prior = 'Mc-q',  mass_parameters_sampled = 'Mc-q' ")
    else:
        pass

    for i, (event_parameters, strain_record) in tqdm(enumerate(zip(events_list, strain_records)), total=len(events_list)):

        # Remove unused IFOs from list, which appear as str with only white spaces.
        event_parameters['ifos_on'] = [ifo for ifo in event_parameters['ifos_on'] if ifo.strip()]

        # Use bilby conventions for the event parameters.
        event_parameters['mass_1']              = event_parameters.pop('m1d')
        event_parameters['mass_2']              = event_parameters.pop('m2d')
        event_parameters['luminosity_distance'] = event_parameters.pop('dL' )

        priors = {}
        # Masses
        if args.mass_parameters_sampled == 'Mc-q' and args.tailor_priors == 'chirp':
            mc  = chirp_mass(event_parameters['mass_1'], event_parameters['mass_2'])
            snr = event_parameters['snr']
            mc_prior_width = chirp_mass_prior_width(mc, snr)
            mc_prior_min, mc_prior_max = max(1., mc - mc_prior_width), mc + mc_prior_width
            priors['chirp_mass'] = [mc_prior_min, mc_prior_max]

        elif args.mass_parameters_sampled in ['m1-m2_custom', 'm1-m2_bilby'] and args.tailor_priors == 'farr16':
            mc  = chirp_mass(event_parameters['mass_1'], event_parameters['mass_2'])
            q   = event_parameters['mass_2']/event_parameters['mass_1']
            snr = event_parameters['snr']
            m1_std = primary_mass_expected_std(mc, q, snr)
            m1_prior_min, m1_prior_max = max(1., event_parameters['mass_1'] - args.nsigma * m1_std), event_parameters['mass_1'] + args.nsigma * m1_std
            priors['mass_1'] = [m1_prior_min, m1_prior_max]
            priors['mass_2'] = [1.0         , m1_prior_max]

        elif args.mass_parameters_sampled in ['m1-m2_custom', 'm1-m2_bilby'] and args.tailor_priors == 'linear':
            m1d_prior_upper_bound = linear_tailored_prior_component_masses_upper_bound(event_parameters['mass_1'])
            priors['mass_1'] = [1.0 , m1d_prior_upper_bound]
            priors['mass_2'] = [1.0 , m1d_prior_upper_bound] # This is actually redundant with q < 1 constraint, we keep it for consistency
        
        else:
            raise ValueError("Unknown combination ({}, {}) of mass_parameters_sampled and tailor_priors. Possible combinations:\n\tmass_parameters_sampled = 'Mc-q' and tailor_priors = 'chirp' \n\tmass_parameters_sampled = 'm1-m2' and tailor_priors = 'farr16' \n\tmass_parameters_sampled = 'm1-m2' and tailor_priors = 'linear' ".format(args.mass_parameters_sampled, args.tailor_priors))

        if i+1 == len(events_list): print(f"\n * Using following settings for mass parameters:\n\tSampling over {args.mass_parameters_sampled.split("_")[0]}, \n\tPriors are uniform over {args.mass_parameters_uniform_prior}, \n\tPrior bounds are tailored according to the {args.tailor_priors} option.\n")

        # Distance
        dL_prior_upper_bound = linear_tailored_prior_luminosity_distance_upper_bound(event_parameters['luminosity_distance'])
        if args.distance_prior in ["Uniform", "PowerLaw", "UniformSourceFrame", "UniformComovingVolume"]:
            if i+1 == len(events_list): print(f"\n * Using a {args.distance_prior} prior for the luminosity_distance.\n")
            priors['luminosity_distance'] = [0., dL_prior_upper_bound, args.distance_prior] # Mpc
        else:
            raise ValueError("Unknown distance prior. Please choose between:\n\tUniform\n\tPowerLaw\n\tUniformSourceFrame\n\tUniformComovingVolume.")

        # Fixing the parameters the user asked to fix for PE
        for par in parameters_to_fix:
            priors[par] = true_value_projected_spins(par, event_parameters, precession)

        # Create a sub directory for the results of each event
        event_subdir = os.path.join(parameter_estimation_subdir, f"event_{i:04d}")
        if not os.path.exists(event_subdir): os.makedirs(event_subdir)
        
        if args.use_recorded_strain and there_is_strain:
            strain_record_filepath = os.path.join(event_subdir, f"strain_record_{os.path.basename(args.pop_dir)}_event_{i:04d}.pickle")
            with open(os.path.join(event_subdir, f"strain_record_{os.path.basename(args.pop_dir)}_event_{i:04d}.pickle"), 'wb') as f:
                pickle.dump(strain_record, f, protocol = pickle.HIGHEST_PROTOCOL)
        elif args.use_recorded_strain:
            print('\n * No strain data found. Carrying on without.\n')
            strain_record_filepath = ""
        else:
            strain_record_filepath = ""

        if   (args.sampler in samplers_types) and (samplers_types[args.sampler] == 'nested'):
            sampler_dependent_config = "\n".join([
                f"{'nlive'     :<29} = {args.nlive     }",
                f"{'npool'     :<29} = {args.npool     }",
                f"{'naccept'   :<29} = {args.naccept   }",
                f"{'sample'    :<29} = {args.sample    }",
            ])
        elif (args.sampler in samplers_types) and (samplers_types[args.sampler] == 'MCMC'):
            sampler_dependent_config = "\n".join([
                f"{'sample'    :<29} = {args.sample    }",
                f"{'naccept'   :<29} = {args.naccept   }",
                f"{'queue-size':<29} = {args.queue_size}",
                f"{'nwalkers'  :<29} = {args.nwalkers  }",
                f"{'nsteps'    :<29} = {args.nsteps    }",
                f"{'ntemps'    :<29} = {args.ntemps    }",
                f"{'threads'   :<29} = {args.threads   }",
                f"{'nparallel' :<29} = {args.nparallel }",
            ])
        else:
            raise ValueError("Unknown sampler. Please choose from the available samplers\n\t{}".format("\n\t".join(samplers_types.keys())))
        
        event_config_content = config_template.format(
            output                        = event_subdir,
            screen_output                 = False,
            event_parameters              = event_parameters,
            plot_waveform                 = args.plot_waveform,
            plot_skymap                   = args.plot_skymap,
            strain_file                   = strain_record_filepath,
            priors                        = priors,
            mass_parameters_sampled       = args.mass_parameters_sampled,
            mass_parameters_uniform_prior = args.mass_parameters_uniform_prior,
            observing_run                 = observing_run,
            waveform                      = waveform,
            precession                    = precession,
            reference_frequency           = reference_frequency,
            sampling_frequency            = sampling_frequency,
            marginalization_settings      = marginalization_settings,
            sampler                       = args.sampler,
            print_method                  = args.print_method,
            sampler_dependent_config      = sampler_dependent_config,
        )

        with open(os.path.join(event_subdir, f"config_PE_{os.path.basename(args.pop_dir)}_event_{i:04d}.ini"), 'w') as f:
            f.write(event_config_content)

    print('\n * Finished.\n')
