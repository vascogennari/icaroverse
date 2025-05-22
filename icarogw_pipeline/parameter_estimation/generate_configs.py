import os, sys, configparser, shutil
from optparse import OptionParser
from tqdm import tqdm

import pickle, pandas as pd, json
import numpy as np

# Internal imports
import initialise


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


config_template = """[input]
output                = {output}
PSD-path              = {psd_dir}
screen-output         = {screen_output}

[model]
event-parameters      = {event_parameters}
priors                = {priors}
priors-dict           = {priors_dict}
observing-run         = {observing_run}
waveform              = {waveform}
precession            = {precession}
reference-frequency   = {reference_frequency}
sampling-frequency    = {sampling_frequency}
phase-marginalization = {phase_marginalization}

[sampler]
sampler               = {sampler}
print-method          = {print_method}
{sampler_dependent_config}
"""


def main():

    # parser = OptionParser(options.usage)
    parser = OptionParser(usage=initialise.usage)
    parser.add_option('-d', '--pop-dir',           type='string', metavar = 'pop_dir',      default = None             )
    # model options
    parser.add_option('-p', '--priors-dict',       type='string', metavar = 'priors_dict',  default = 'bilby',           help="Options: bilby, custom")
    parser.add_option('-t', '--tailor_priors', action="store_true", dest='tailor_priors',   default = False,             help="Flag to choose the mass priors adapted to true values.")
    parser.add_option('-k', '--nsigma',            type='int',    metavar = 'nsigma',       default = 7,                 help="Width of the tailored prior ranges, in units of estimate of m1_std (prior: m1 +/- nsigma * m1_std).")
    parser.add_option(      '--phase', action="store_true", dest='phase_marginalization',   default = False,             help="Flag to turn on phase marginalization in the PE corresponding to the generated config files.")
    # sampler
    parser.add_option('-s', '--sampler',           type='string', metavar = 'sampler',      default = 'dynesty'        )
    parser.add_option('-m', '--print-method',      type='string', metavar = 'print_method', default = 'interval-60'    )
    # nested samplers options
    parser.add_option('-l', '--nlive',             type='int',    metavar = 'nlive',        default = 1000             )
    parser.add_option('-n', '--npool',             type='int',    metavar = 'npool',        default = 10               )
    parser.add_option('-a', '--naccept',           type='int',    metavar = 'naccept',      default = 60               )
    parser.add_option(      '--sample',            type='string', metavar = 'sample',       default = 'acceptance-walk')
    # MCMC samplers options
    parser.add_option('-q', '--queue-size',        type='int',    metavar = 'queue_size',   default = 1                )
    parser.add_option('-w', '--nwalkers',          type='int',    metavar = 'nwalkers',     default = 64               )
    parser.add_option(      '--nsteps',            type='int',    metavar = 'nsteps',       default = 1000             )
    parser.add_option(      '--ntemps',            type='int',    metavar = 'ntemps',       default = 10               )
    parser.add_option('-r', '--threads',           type='int',    metavar = 'threads',      default = 1                )
    parser.add_option(      '--nparallel',         type='int',    metavar = 'nparallel',    default = 1                )
    (opts, _) = parser.parse_args()

    samplers_types = {
        'dynesty': 'nested', 
        'nessai':  'nested', 
        'ptemcee': 'MCMC',
    }

    # Retrieve population and transform it into a list of single event dict
    data_filepath = os.path.join(opts.pop_dir, "population_observed.pickle")
    if os.path.exists(data_filepath):
        with open(data_filepath, 'rb') as f: data = pickle.load(f)
        data['ifos_on'] = data['ifos_on'].tolist()
        events_list = pd.DataFrame(data).to_dict('records')
    else:
        raise FileNotFoundError(f"Unknown file {data_filepath}. Please make sure the population directory you passed does contain a 'population_observed.pickle' file.")

    # Create 'parameter_estimation' subdirectory in opts.pop_dir
    parameter_estimation_subdir = os.path.join(opts.pop_dir, "parameter_estimation")
    if not os.path.exists(parameter_estimation_subdir): os.makedirs(parameter_estimation_subdir)

    # Fetch a 'bilby_settings' file where are stored observing-run, waveform, reference-frequency, sampling-frequency from the population generation
    try: 
        print('\n * Reading the population generation settings.\n')
        with open(os.path.join(opts.pop_dir, "analysis_settings.json")) as f: population_generation_pars = json.load(f)
    except FileNotFoundError: 
        raise FileNotFoundError("Please make sure the population directory does contain an 'analysis_settings.json' file")
    
    if population_generation_pars['SNR-method'] == 'bilby':
        psd_dir             = population_generation_pars['PSD-path']
        observing_run       = population_generation_pars['observing-run']
        waveform            = population_generation_pars['snr-bilby-waveform']
        precession          = population_generation_pars['snr-bilby-precessing-wf']
        reference_frequency = population_generation_pars['snr-bilby-reference-frequency']
        sampling_frequency  = population_generation_pars['snr-bilby-sampling-frequency']
    else:
        raise ValueError("Please make sure the simulated population you are trying to work with has been processed with bilby within the present pipeline (for SNR computation)")

    print('\n * Generating config files.\n')
    for i, event_parameters in tqdm(enumerate(events_list)):

        # Use bilby conventions for the event parameters.
        event_parameters['mass_1']              = event_parameters.pop('m1d')
        event_parameters['mass_2']              = event_parameters.pop('m2d')
        event_parameters['luminosity_distance'] = event_parameters.pop('dL' )

        # Remove unused IFOs from list, which appear as str with only white spaces.
        event_parameters['ifos_on'] = [ifo for ifo in event_parameters['ifos_on'] if ifo.strip()]

        priors = {}
        if opts.tailor_priors:
            mc  = chirp_mass(event_parameters['mass_1'], event_parameters['mass_2'])
            q   = event_parameters['mass_2']/event_parameters['mass_1']
            snr = event_parameters['snr']
            m1_std = primary_mass_expected_std(mc, q, snr)
            if   observing_run == 'O3': m1_prior_abs_max = 200.
            elif observing_run == 'O4': m1_prior_abs_max = 300.
            elif observing_run == 'O5': m1_prior_abs_max = 400.
            m1_prior_min, m1_prior_max = max(1., event_parameters['mass_1'] - opts.nsigma * m1_std), min(m1_prior_abs_max, event_parameters['mass_1'] + opts.nsigma * m1_std)
            priors['mass_1'] = [m1_prior_min, m1_prior_max]
            priors['mass_2'] = [m1_prior_min, m1_prior_max]
        if observing_run == 'O4':
            priors['luminosity_distance'] = [0., 1.6e4] # up to redshift ~2 with Planck15 cosmo
        if observing_run == 'O5':
            priors['luminosity_distance'] = [0., 3e4] # up to redshift ~3.37 with Planck15 cosmo

        # Create a sub directory for the results of each event
        event_subdir = os.path.join(parameter_estimation_subdir, f"event_{i:04d}")
        if not os.path.exists(event_subdir): os.makedirs(event_subdir)
        
        if   (opts.sampler in samplers_types) and (samplers_types[opts.sampler] == 'nested'):
            sampler_dependent_config = "\n".join([
                f"{'nlive'     :<21} = {opts.nlive     }",
                f"{'npool'     :<21} = {opts.npool     }",
                f"{'naccept'   :<21} = {opts.naccept   }",
                f"{'sample'    :<21} = {opts.sample    }",
            ])
        elif (opts.sampler in samplers_types) and (samplers_types[opts.sampler] == 'MCMC'):
            sampler_dependent_config = "\n".join([
                f"{'sample'    :<21} = {opts.sample    }",
                f"{'naccept'   :<21} = {opts.naccept   }",
                f"{'queue-size':<21} = {opts.queue_size}",
                f"{'nwalkers'  :<21} = {opts.nwalkers  }",
                f"{'nsteps'    :<21} = {opts.nsteps    }",
                f"{'ntemps'    :<21} = {opts.ntemps    }",
                f"{'threads'   :<21} = {opts.threads   }",
                f"{'nparallel' :<21} = {opts.nparallel }",
            ])
        else:
            raise ValueError("Unknown sampler. Please choose from the available samplers\n\t{}".format("\n\t".join(samplers_types.keys())))
        
        event_config_content = config_template.format(
            output                   = event_subdir,
            psd_dir                  = psd_dir,
            screen_output            = False,
            event_parameters         = event_parameters,
            priors                   = priors,
            priors_dict              = opts.priors_dict,
            observing_run            = observing_run,
            waveform                 = waveform,
            precession               = precession,
            reference_frequency      = reference_frequency,
            sampling_frequency       = sampling_frequency,
            phase_marginalization    = opts.phase_marginalization,
            sampler                  = opts.sampler,
            print_method             = opts.print_method,
            sampler_dependent_config = sampler_dependent_config,
        )

        with open(os.path.join(event_subdir, f"config_PE_{os.path.basename(opts.pop_dir)}_event_{i:04d}.ini"), 'w') as f:
            f.write(event_config_content)

    print('\n * Finished.\n')



# Execute the main function.
if __name__=='__main__': main()
