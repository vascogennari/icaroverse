import os, sys, configparser, shutil, time
from optparse import OptionParser

import pickle, h5py, pandas as pd, json
import numpy as np
import icarogw, bilby, astropy

# Internal imports
import options, icarogw_postprocessing



def get_wrapper(wrap_name, input_wrapper = None, order = None):

    print('\t{}'.format(wrap_name))
    wrap = getattr(icarogw.wrappers, wrap_name)
    if order == None:
        if not input_wrapper == None:
            return wrap(input_wrapper)
        else:
            return wrap()
    else:
        return wrap(order = order)

def print_dictionary(dictionary):
      
    for key in dictionary.keys():
        max_len = len(max(dictionary.keys(), key = len))
        if not key == 'all-priors': print('\t{}  {}'.format(key.ljust(max_len), dictionary[key]))

def check_effective_number_injections(pars, likelihood, n_events):

    def single_likelihood_eval():
        count = time.time()
        _     = likelihood.log_likelihood()
        print('\n\tA single likelihood evaluation takes {0:.5f} [s].'.format(time.time() - count))

    if pars['simulation'] and (not pars['true-values'] == {}):
        likelihood.parameters = {key: pars['true-values'][key] for key in likelihood.rate_model.population_parameters}
        single_likelihood_eval()
        N_eff_inj = likelihood.injections.effective_injections_number()
        stability = N_eff_inj / (4 * n_events)
        print('\n\tThe effective number of injections for the injected model is {0:.1f}. N_eff_inj/4*N_events is {1:.1f}.'.format(N_eff_inj, stability))
        if stability < 1: print('\n\tWARNING: The number of injections is not enough to ensure numerical stability in the computation of selection effects in the likelihood. Please consider using a larger set of injections.')


class Wrappers:

    def __init__(self, pars):

        self.pars = pars
        print('\n * Initialising the wrappers.\n')

    def PrimaryMass(self, pars):

        # Stationary PowerLaw.
        if not ('PowerLawRedshift' in pars['model-primary'] or 'GaussianRedshift-order-' in pars['model-primary']):

            if   pars['model-primary'] == 'PowerLaw' or pars['model-primary'] == 'PowerLaw-GaussianRedshiftLinear': 
                                                                       w = get_wrapper('massprior_PowerLaw')
            elif pars['model-primary'] == 'PowerLaw-Gaussian':         w = get_wrapper('massprior_PowerLawPeak')

            elif not pars['model-primary'] == 'PowerLawRedshiftLinear-GaussianRedshiftLinear':
                raise ValueError('Unknown model for the primary mass {}. See the help for the available models.'.format(pars['model-primary']))

            if pars['low-smoothing']:                                  w = get_wrapper('lowSmoothedwrapper',                                  input_wrapper = w)

            # Evolving Gaussian.
            if 'GaussianRedshift' in pars['model-primary']:
                if pars['redshift-transition'] == '':                  w = get_wrapper('mixed_mass_redshift_evolving',                        input_wrapper = w)
                if pars['redshift-transition'] == 'sigmoid':           w = get_wrapper('mixed_mass_redshift_evolving_sigmoid',                input_wrapper = w)
                if pars['redshift-transition'] == 'double-sigmoid':    w = get_wrapper('double_mixed_mass_redshift_evolving_sigmoid',         input_wrapper = w)
                if pars['redshift-transition'] == 'linear':            w = get_wrapper('double_mixed_mass_redshift_evolving_linear_prior',    input_wrapper = w)
                if pars['redshift-transition'] == 'linear-sinusoid':   w = get_wrapper('double_mixed_mass_redshift_evolving_linear_sinusoid', input_wrapper = w)

        # Evolving PowerLaw and Gaussian.
        else:
            if   pars['model-primary'] == 'PowerLawRedshiftLinear-GaussianRedshiftLinear':
                                                                       w = get_wrapper('PowerLawLinear_GaussianLinear_TransitionLinear')
            elif 'GaussianRedshift-order-' in pars['model-primary']:
                                                                       order = int(pars['model-primary'].split('GaussianRedshift-order-')[-1])
                                                                       w = get_wrapper('GaussianEvolving', order = order)
            else:
                raise ValueError('Unknown model for the primary mass {}. Please consult the available models.'.format(pars['model-primary']))          

        # Only for stationary PowerLaw and Gaussian.
        # FIXME: This model should be substituted with a conditional Bilby prior.
        if pars['positive-peak']: w = icarogw.wrappers.massprior_PowerLawPeakPositive(w)
        return w

    def SecondaryMass(self, pars, m1w = None):

        if not pars['single-mass']:
            if   pars['model-secondary'] == 'PowerLaw-Gaussian':           w = get_wrapper('m1m2_conditioned',                                    input_wrapper = m1w)
            elif pars['model-secondary'] == 'PowerLaw':                    w = get_wrapper('m1m2_conditioned',                                    input_wrapper = m1w)
            elif pars['model-secondary'] == 'MassRatio':                   w = get_wrapper('mass_ratio_prior_Gaussian')
            else:
                raise ValueError('Unknown model for the secondary mass {}. Please consult the available models.'.format(pars['model-secondary']))
        else:
            print('\t * Skipping secondary mass wrapper')
            w = None
        return w

    def RateEvolution(self, pars):

        if   pars['model-rate'] == 'MadauDickinson':                   w = get_wrapper('rateevolution_Madau')
        elif pars['model-rate'] == 'BetaDistribution':                 w = get_wrapper('rateevolution_beta')
        elif pars['model-rate'] == 'BetaDistribution-Line':            w = get_wrapper('rateevolution_beta_line')
        elif pars['model-rate'] == 'MadauDickinson-GammaDistribution': w = get_wrapper('rateevolution_Madau_gamma')
        elif pars['model-rate'] == 'PowerLaw':                         w = get_wrapper('rateevolution_PowerLaw')
        else:
            raise ValueError('Unknown model for the rate evolution {}. Please consult the available models.'.format(pars['model-rate']))   
        return w
    
    def Cosmology(self):

        w = icarogw.wrappers.FlatLambdaCDM_wrap(zmax = 20.)
        return w

    def ReferenceCosmology(self):
          
        w = icarogw.cosmology.astropycosmology(10)
        w.build_cosmology(astropy.cosmology.FlatLambdaCDM(H0 = 67.7, Om0 = 0.308))
        return w

    def return_Wrappers(self):

        self.Wrapper_PrimaryMass   = self.PrimaryMass(self.pars)
        self.Wrapper_SecondaryMass = self.SecondaryMass(self.pars, self.Wrapper_PrimaryMass)
        self.Wrapper_RateEvolution = self.RateEvolution(self.pars)
        self.Wrapper_Cosmology     = self.Cosmology()
        self.Wrapper_RefCosmology  = self.ReferenceCosmology()

        return self.Wrapper_PrimaryMass, self.Wrapper_SecondaryMass, self.Wrapper_RateEvolution, self.Wrapper_Cosmology, self.Wrapper_RefCosmology



class Rate():
      
    def __init__(self, pars, m1w, m2w, rw, cw):

        if not pars['single-mass']:
            if  (not 'Redshift' in pars['model-primary']) and (not 'MassRatio' in pars['model-secondary']):
                self.w = icarogw.rates.CBC_vanilla_rate(             cw,      m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_vanilla_rate'))
            elif not  'Redshift' in pars['model-primary'] and      'MassRatio' in pars['model-secondary']:
                self.w = icarogw.rates.CBC_rate_m1_q(                cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m1_q'))
            elif      'Redshift' in pars['model-primary'] and  not 'MassRatio' in pars['model-secondary']:
                self.w = icarogw.rates.CBC_rate_m1_given_redshift_m2(cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m1_given_redshift_m2'))
            elif      'Redshift' in pars['model-primary'] and     'MassRatio' in pars['model-secondary']:
                self.w = icarogw.rates.CBC_rate_m1_given_redshift_q( cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m1_given_redshift_q'))
        else:
            self.w = icarogw.rates.CBC_rate_m_given_redshift(        cw, m1w,      rw, scale_free = pars['scale-free'])

        print('\n * Population parameters.\n')
        print('\t{}'.format('[%s]' % ', '.join(map(str, self.w.population_parameters))))
        pars['population-parameters'] = self.w.population_parameters
    
    def return_Rate(self):
        return self.w



class SelectionEffects:
        
    def __init__(self, pars, ref_cosmo):

        # Load file with the injections used to compute selection effects.
        print('\n * Loading injections for selection effects.\n\n\t{}'.format(pars['injections-path']))
        try:
            if   '.hdf5'   in pars['injections-path']:
                data_inj = h5py.File(pars['injections-path'])
            elif '.pickle' in pars['injections-path']:
                with open(pars['injections-path'], 'rb') as f: data_inj = pickle.load(f)
            else:
                raise ValueError('Unknown format for the file containing the injections for selection effects. Please make sure the file is correct:\n{}'.format(pars['injections-path']))
        except:
            raise ValueError('Could not open the file containing the injections for selection effects. Please verify that the path is correct:\n{}'.format(pars['injections-path']))

        # O3 Cosmology paper injections
        if   pars['O3-cosmology']:

            obs_time = (28519200 / 86400) / 365
            pars['injections-number'] = data_inj.attrs['total_generated']
            prior  = icarogw.cupy_pal.np2cp(data_inj['injections/mass1_source_mass2_source_sampling_pdf'][()] * data_inj['injections/redshift_sampling_pdf'][()])
            # Converting the injections from source to detector frame, we need to correct the injections prior by the Jacobian of the transformation (m_1s, m_2s, z) --> (m_1d, m_2d, d_L).
            prior *= icarogw.conversions.source2detector_jacobian(icarogw.cupy_pal.np2cp(data_inj['injections/redshift'][()]), ref_cosmo)

            tmp = np.vstack([data_inj['injections'][key] for key in ['ifar_cwb', 'ifar_gstlal', 'ifar_mbta', 'ifar_pycbc_bbh', 'ifar_pycbc_hyperbank']])
            ifarmax = np.max(tmp, axis = 0)
            mass_1 = 'injections/mass1'
            mass_2 = 'injections/mass2'
            dist   = 'injections/distance'

            inj_dict = {
                'mass_1':              data_inj[mass_1][()],
                'mass_2':              data_inj[mass_2][()],
                'luminosity_distance': data_inj[dist][()]}
            
            # If using the mass ratio, correct the prior with the Jacobian m2->q.
            if pars['model-secondary'] == 'MassRatio':
                prior *= data_inj[mass_1][()]   # (m_1, m_2) --> (m_1, q)
                inj_dict['mass_ratio'] = inj_dict.pop('mass_2') / data_inj[mass_1][()]

        # Internal simulations
        elif pars['simulation']:

            prior = data_inj['prior']
            obs_time = 1
            mass_1 = 'm1d'
            mass_2 = 'm2d'
            dist   = 'dL'

            inj_dict = {
                'mass_1':              data_inj[mass_1],
                'mass_2':              data_inj[mass_2],
                'luminosity_distance': data_inj[dist]}
            
            if not pars['single-mass']:
                # If using the mass ratio, correct the prior with the Jacobian m2->q.
                if pars['model-secondary'] == 'MassRatio':
                    prior *= data_inj[mass_1]
                    inj_dict['mass_ratio'] = inj_dict.pop('mass_2') / data_inj[mass_1]
            else:
                # If only using one mass, remove the Jacobian contribution from the secondary.
                # This operation depends on the injection prior used to generate the injections.
                prior *= (1 + ref_cosmo.dl2z(data_inj[dist]))
                inj_dict.pop('mass_2')
        else:
            raise ValueError('Unknown option to compute selection effects.')

        self.injections = icarogw.injections.injections(inj_dict, prior = prior, ntotal = pars['injections-number'], Tobs = obs_time)
        if   pars['selection-effects-cut'] == 'snr' : self.injections.update_cut(data_inj['snr'] >= pars['snr-cut' ])
        elif pars['selection-effects-cut'] == 'ifar': self.injections.update_cut(ifarmax         >= pars['ifar-cut'])
        else:
            raise ValueError('Unknown option to compute the selection effects cut.')

        print('\n\tUsing {} injections out of {} to compute selection effects.\n'.format(len(self.injections.injections_data_original['mass_1']), len(inj_dict['mass_1'])))

    def return_SelectionEffects(self):
        return self.injections



class Data:
       
    def __init__(self, pars, ref_cosmo):
        
        print('\n * Loading data.\n\n\t{}'.format(pars['data-path']))
        if not pars['distance-prior-PE']:              print('\n\tUsing a flat prior for PE samples on the luminosity distance.')
        if not pars['single-mass']:
            if pars['model-secondary'] == 'MassRatio': print('\n\tCorrecting the PE samples prior for mass ratio.')
        else:                                          print('\n\tCorrecting the PE samples prior for mass ratio.')
              
        # O3 Cosmology paper injections
        if   pars['O3-cosmology']:

            sys.path.append(pars['data-path'])
            from analyses_dictionaries import BBHs_O3_IFAR_4

            samps_dict = {}
            for ev in list(BBHs_O3_IFAR_4.keys()):
                
                # Skip the two low mass ratio events in Rinaldi+. This is to avoid problems in assuming a Gaussian distirbution in the mass ratio.
                if ('190412' in ev) or ('190917' in ev): continue
                else:                                    print('    {}'.format(ev))

                tmp = h5py.File(BBHs_O3_IFAR_4[ev]['PE'])
                data_evs = tmp[BBHs_O3_IFAR_4[ev]['PE_waveform']]['posterior_samples']

                pos_dict  = {
                    'mass_1'             : data_evs['mass_1'][()],
                    'mass_2'             : data_evs['mass_2'][()],
                    'luminosity_distance': data_evs['luminosity_distance'][()]}

                # This assumes that the luminosity distance prior for the single events PE is uniform in volume, thus p(d_L)=d_L^2.
                # If not, set the prior to one.
                if pars['distance-prior-PE']: prior = np.power(   data_evs['luminosity_distance'][()], 2.)
                else:                         prior = np.ones(len(data_evs['luminosity_distance'][()]))

                # If using the mass ratio, correct the prior with the Jacobian m2->q.
                if pars['model-secondary'] == 'MassRatio':
                    pos_dict['mass_ratio'] = pos_dict.pop('mass_2') / data_evs['mass_1'][()]
                    prior *= data_evs['mass_1'][()]
                
                samps_dict[ev] = icarogw.posterior_samples.posterior_samples(pos_dict, prior = prior)

        # Internal simulations
        elif pars['simulation']:

            if '.pickle' in pars['data-path']:
                with open(pars['data-path'], 'rb') as f: data_evs = pickle.load(f)
            else:
                raise ValueError('Unknown format for the file containing the single events samples. Please make sure the file is correct:\n{}'.format(pars['data-path']))

            samps_dict = {}
            for i in range(len(data_evs['m1d'])):
                pos_dict = {
                    'mass_1':              np.array([data_evs['m1d'][i]]),
                    'mass_2':              np.array([data_evs['m2d'][i]]),
                    'luminosity_distance': np.array([data_evs['dL'][i]])}

                # This assumes that the luminosity distance prior for the single events PE is uniform in volume, thus p(d_L)=d_L^2.
                # If not, set the prior to one.
                if pars['distance-prior-PE']: prior = np.array([data_evs['dL'][i]**2])
                else:                         prior = np.array([1.])

                if not pars['single-mass']:
                    # If using the mass ratio, correct the prior with the Jacobian m2->q.
                    if pars['model-secondary'] == 'MassRatio':
                        pos_dict['mass_ratio'] = pos_dict.pop('mass_2') / np.array([data_evs['m1d'][i]])
                        #prior *= np.array([data_evs['m1d'][i]])    # FIXME: Include this option for future simulations with PE samples.
                else:
                    pos_dict.pop('mass_2')

                samps_dict['{}'.format(i)] = icarogw.posterior_samples.posterior_samples(pos_dict, prior = prior)
        else:
            raise ValueError('Unknown option to process single events data.')
        
        self.data = icarogw.posterior_samples.posterior_samples_catalog(samps_dict)
        print('\n\tUsing a population of {} events.'.format(self.data.n_ev))

    def return_Data(self):
        return self.data



class LikelihoodPrior:

    def __init__(self, pars, data, injections, wrapper):
          
        self.pars       = pars
        self.data       = data
        self.injections = injections
        self.wrapper    = wrapper

    def Likelihood(self, data, injections, wrapper):
        
        res = icarogw.likelihood.hierarchical_likelihood(
                        data, injections, wrapper,
                        nparallel = self.pars['nparallel'],
                        neffPE    = self.pars['neffPE'],
                        neffINJ   = self.pars['neffINJ'])
        return res

    def Prior(self, pars, w):
    
        def conditional_prior(prior):
            prior['mean_three_sigmas'] = prior['mu_z0'] - 3 * prior['sigma_z0']
            return prior

        def initialise_prior(dict_in, dict_out, w):
              
            for par in w.population_parameters:
                if   type(dict_in[par]) == list:  dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1])
                elif type(dict_in[par]) == float: dict_out[par] = dict_in[par]
                else:
                    raise ValueError('Unknown type for prior on {}'.format(dict_in[par]))
            
            print('\n * Using the following priors.\n')
            print_dictionary({key: dict_in[key] for key in dict_out.keys()})

            return dict_out

        # FIXME: Add option to deal with additional conditions on the prior.
        if pars['conditional-prior']: prior = bilby.core.prior.PriorDict(conversion_function = conditional_prior)
        else:                         prior = bilby.core.prior.PriorDict()
        prior = initialise_prior(pars['all-priors'], prior, w)

        if pars['conditional-prior']:
              prior['mean_three_sigmas'] = bilby.prior.Constraint(0., 1000.)
              print('\n * Adding conditional priors.\n')
              print_dictionary({'mean_three_sigmas': [0., 1000.]})

        return prior
        
    def return_LikelihoodPrior(self):

        self.likelihood = self.Likelihood(self.data, self.injections, self.wrapper)
        self.prior      = self.Prior(self.pars, self.wrapper)

        return self.likelihood, self.prior



def main():

    # ------------------------------------------------ #
    # Read config file and initialise input parameters #
    # ------------------------------------------------ #

    parser = OptionParser(options.usage)
    parser.add_option('--config-file', type='string', metavar = 'config_file', default = None)
    (opts, _) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    # Initialise input parameters dictionary.
    input_pars = options.InitialiseOptions(Config)

    # Set output directory.
    # FIXME: Add option to control that only one of the two between 'O3-cosmology' and 'simulation' is active.
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])

    # Copy config file to output.
    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))

    # Deviate stdout and stderr to file.
    if not input_pars['screen-output']:
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_icarogw.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_icarogw.txt'), 'w')
    else: pass
    print('\n\n Starting  i c a r o g w  runner\n\n')

    # Print run parameters.
    print(' * I will be running with the following parameters.\n')
    print_dictionary(input_pars)

    # ------------------------------------------------------------------------------ #
    # Initialise the ICAROGW run: model, data, selection effects, likelihood, priors #
    # ------------------------------------------------------------------------------ #

    # Initialise the model wrappers.
    tmp = Wrappers(input_pars)
    m1w, m2w, rw, cw, ref_cosmo = tmp.return_Wrappers()

    tmp = Rate(input_pars, m1w, m2w, rw, cw)
    wrapper = tmp.return_Rate()

    # Read injections for selection effects.
    tmp = SelectionEffects(input_pars, ref_cosmo)
    injections = tmp.return_SelectionEffects()

    # Read events data.
    tmp = Data(input_pars, ref_cosmo)
    data = tmp.return_Data()

    # Initialise hierarchical likelihood and set the priors.
    tmp = LikelihoodPrior(input_pars, data, injections, wrapper)
    likelihood, prior = tmp.return_LikelihoodPrior()

    # Control the effective number of injections on the injected model.
    check_effective_number_injections(input_pars, likelihood, data.n_ev)

    # ----------------------------------------------- #
    # Start the sampler and run hierarchical analysis #
    # ----------------------------------------------- #

    print('\n * Running hierarchical analysis with this settings.\n')
    print_dictionary({key: input_pars[key] for key in ['sampler', 'nlive', 'naccept', 'npool', 'print_method', 'sample']})

    # Start Bilby sampler.
    print('\n * Starting the sampler.\n')
    hierarchical = bilby.run_sampler(
            likelihood, prior,
            sampler      = input_pars['sampler'],
            nlive        = input_pars['nlive'],
            naccept      = input_pars['naccept'],
            npool        = input_pars['npool'],
            print_method = input_pars['print_method'],
            sample       = input_pars['sample'],
            outdir       = os.path.join(input_pars['output'], 'sampler'),
    )
    if input_pars['true-values'] == {}: hierarchical.plot_corner()
    else:                               hierarchical.plot_corner(truth = {key: input_pars['true-values'][key] for key in hierarchical.search_parameter_keys})
    
    # Get the samples.
    samp_path = os.path.join(input_pars['output'], 'sampler', 'label_result.json')
    with open(samp_path) as f:
        tmp = json.load(f)
        df  = pd.DataFrame(tmp['posterior']['content'])
        priors_dict = tmp['priors']

    # Save the evidence.
    with open('{}/log_evidence.txt'.format(input_pars['output']), 'w') as f:
        f.write('{}\n'.format('# log_Z_base_e\tlog_Z_err\tmax_log_L'))
        f.write('{}\t{}\t\t{}'.format(round(tmp['log_evidence'], 2), round(tmp['log_evidence_err'], 2), round(max(df['log_likelihood']), 2)))

    # ----------------------------------- #
    # Plots production and postprocessing #
    # ----------------------------------- #

    print('\n * Producing plots.')
    # FIXME: Should define a new variable instead of overwriting the output path.
    input_pars['output'] = os.path.join(input_pars['output'], 'plots')
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])

    tmp = icarogw_postprocessing.Plots(input_pars, df, m1w, m2w, rw, cw, ref_cosmo, wrapper, priors_dict, injections)
    tmp.ProducePlots()

    print('\n * Finished.\n')

if __name__=='__main__':
    main()
