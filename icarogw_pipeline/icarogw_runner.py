import os, sys, configparser, ast, shutil
from optparse import OptionParser
from options import usage

import pickle, h5py, pandas as pd
import numpy as np
import icarogw, bilby
import icarogw_postprocessing as icaroproc



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


def default_priors():
      
    prior = {
        # Cosmology
        'H0'          : 67.7,
        'Om0'         : 0.308,

        # Primary mass distribution
        'delta_m'     : [0.1, 10. ],

        'alpha'       : [1.,  6.  ],
        'alpha_z0'    : [1.,  6.  ],
        'alpha_z1'    : [-5., 20. ],

        'mmin'        : [2. , 20. ],
        'mmin_z0'     : [2. , 20. ],
        'mmin_z1'     : [-5., 50. ],
        'mmax'        : [65., 200.],
        'mmax_z0'     : [65., 200.],
        'mmax_z1'     : 0.,

        'mu_g'        : [20., 60. ],
        'mu_z0'       : [20., 60. ],
        'mu_z1'       : [-80., 80.],
        'sigma_g'     : [5. , 30. ],
        'sigma_z0'    : [5. , 30. ],
        'sigma_z1'    : [0.,  20. ],
        'mix_z0'      : [0. , 1.  ],
        'mix_z1'      : [0. , 1.  ],

        'lambda_peak' : [0. , 1.  ],
        
        # Secondary mass distribution
        'beta'        : [1.,  6.  ],
        'mu_q'        : [0.4, 1.  ],
        'sigma_q'     : [0.01, 0.9],

        # Rate evolution
        'gamma'       : [-50., 30. ],
        'kappa'       : [-6. , 6.  ],
        'zp'          : [0. , 4.   ],
        'R0'          : [0. , 100. ],
    }

    return prior
    



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

        if   pars['model-secondary'] == 'PowerLaw-Gaussian':           w = get_wrapper('m1m2_conditioned',                                    input_wrapper = m1w)
        elif pars['model-secondary'] == 'PowerLaw':                    w = get_wrapper('m1m2_conditioned',                                    input_wrapper = m1w)
        elif pars['model-secondary'] == 'MassRatio':                   w = get_wrapper('mass_ratio_prior_Gaussian')
        else:
            raise ValueError('Unknown model for the secondary mass {}. Please consult the available models.'.format(pars['model-secondary']))
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
    
    def return_Wrappers(self):

        self.Wrapper_PrimaryMass   = self.PrimaryMass(self.pars)
        self.Wrapper_SecondaryMass = self.SecondaryMass(self.pars, self.Wrapper_PrimaryMass)
        self.Wrapper_RateEvolution = self.RateEvolution(self.pars)
        self.Wrapper_Cosmology     = self.Cosmology()

        return self.Wrapper_PrimaryMass, self.Wrapper_SecondaryMass, self.Wrapper_RateEvolution, self.Wrapper_Cosmology



class Rate():
      
    def __init__(self, pars, m1w, m2w, rw, cw):

        if  (not 'Redshift' in pars['model-primary']) and (not    'MassRatio' in pars['model-secondary']):
            self.w = icarogw.rates.CBC_vanilla_rate(             cw,      m2w, rw, scale_free = pars['scale-free'])
            print('\t{}'.format('CBC_vanilla_rate'))
        elif not  'Redshift' in pars['model-primary'] and     'MassRatio' in pars['model-secondary']:
            self.w = icarogw.rates.CBC_rate_m1_q(                cw, m1w, m2w, rw, scale_free = pars['scale-free'])
            print('\t{}'.format('CBC_rate_m1_q'))
        elif      'Redshift' in pars['model-primary'] and  not 'MassRatio' in pars['model-secondary']:
            self.w = icarogw.rates.CBC_rate_m1_given_redshift_m2(cw, m1w, m2w, rw, scale_free = pars['scale-free'])
            print('\t{}'.format('CBC_rate_m1_given_redshift_m2'))
        elif      'Redshift' in pars['model-primary'] and     'MassRatio' in pars['model-secondary']:
            self.w = icarogw.rates.CBC_rate_m1_given_redshift_q( cw, m1w, m2w, rw, scale_free = pars['scale-free'])
            print('\t{}'.format('CBC_rate_m1_given_redshift_q'))

        print('\n * Population parameters.\n')
        print('\t{}'.format('[%s]' % ', '.join(map(str, self.w.population_parameters))))
    
    def return_Rate(self):
        return self.w



class SelectionEffects:
        
    def __init__(self, pars):

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

            from astropy.cosmology import FlatLambdaCDM
            obs_time = (28519200 / 86400) / 365
            cosmo_ref = icarogw.cosmology.astropycosmology(10)
            cosmo_ref.build_cosmology(FlatLambdaCDM(H0 = 67.7, Om0 = 0.308))
            pars['injections-number'] = data_inj.attrs['total_generated']
            prior  = icarogw.cupy_pal.np2cp(data_inj['injections/mass1_source_mass2_source_sampling_pdf'][()] * data_inj['injections/redshift_sampling_pdf'][()])
            prior *= icarogw.conversions.source2detector_jacobian(icarogw.cupy_pal.np2cp(data_inj['injections/redshift'][()]), cosmo_ref)

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
                prior *= data_inj[mass_1][()]
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
            
            # If using the mass ratio, correct the prior with the Jacobian m2->q.
            if pars['model-secondary'] == 'MassRatio':
                prior *= data_inj[mass_1]
                inj_dict['mass_ratio'] = inj_dict.pop('mass_2') / data_inj[mass_1]
        else:
            raise ValueError('Unknown option to compute selection effects.')

        self.injections = icarogw.injections.injections(inj_dict, prior = prior, ntotal = pars['injections-number'], Tobs = obs_time)
        if   pars['selection-effects-cut'] == 'snr' : self.injections.update_cut(data_inj['snr'] >= pars['snr-cut' ])
        elif pars['selection-effects-cut'] == 'ifar': self.injections.update_cut(ifarmax         >= pars['ifar-cut'])
        else:
            raise ValueError('Unknown option to compute the selection effects cut.')

    def return_SelectionEffects(self):
        return self.injections



class Data:
       
    def __init__(self, pars):
        
        print('\n * Loading data.\n\n\t{}\n'.format(pars['data-path']))
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
                    'luminosity_distance': data_evs['luminosity_distance'][()]
                }

                # If using the mass ratio, correct the prior with the Jacobian m2->q.
                if pars['model-secondary'] == 'MassRatio':
                    pos_dict['mass_ratio'] = pos_dict.pop('mass_2') / data_evs['mass_1'][()]

                # FIXME: What is that due to?
                pp_remove      = np.power(data_evs['luminosity_distance'][()], 2.) * data_evs['mass_1'][()]
                samps_dict[ev] = icarogw.posterior_samples.posterior_samples(pos_dict, prior = pp_remove)

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
                    'mass_ratio':          np.array([data_evs['m2d'][i]]) / np.array([data_evs['m1d'][i]]),
                    'luminosity_distance': np.array([data_evs['dL'][i]])}
                
                pp_remove = np.array([data_evs['dL'][i]**2 * data_evs['m1d'][i]])
                samps_dict['{}'.format(i)] = icarogw.posterior_samples.posterior_samples(pos_dict, prior = pp_remove)
        else:
            raise ValueError('Unknown option to process single events data.')
        
        self.data = icarogw.posterior_samples.posterior_samples_catalog(samps_dict)

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

        # FIXME: Add option to del with additional conditions on the prior
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

    parser = OptionParser(usage)
    parser.add_option('--config-file', type='string', metavar = 'config_file', default = None)
    (opts, args) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    # INPUT PARAMETERS
    input_pars = {

        # Input
        'output'                : 'icarogw_run',
        'screen-output'         : 0,

        # Wrappers
        'model-primary'         : 'PowerLaw-Gaussian',                     
        'model-secondary'       : 'MassRatio',
        'model-rate'            : 'PowerLaw',

        'redshift-transition'   : 'linear',
        'positive-peak'         : 0,
        'low-smoothing'         : 0,
        'priors'                : {},
        'conditional-prior'     : 0,
        'scale-free'            : 0,

        # Selection effects
        'injections-path'       : '',
        'injections-number'     : 1,
        'snr-cut'               : 12.,
        'ifar-cut'              : 4.,
        'selection-effects-cut' : 'snr',

        # Data
        'O3-cosmology'          : 0,
        'simulation'            : 1,
        'data-path'             : '',
    
        # Likelihood
        'nparallel'             : 1,
        'neffPE'                : 1,
        'neffINJ'               : None,

        # Sampler
        'sampler'               : 'dynesty',
        'nlive'                 : 200,
        'naccept'               : 60,
        'npool'                 : 1,
        'print_method'          : 'interval-60',
        'sample'                : 'acceptance-walk',

        # Plots
        'N-points'              : 1000,
        'N-z-slices'            : 10,
        'bounds-m1'             : [1, 100],
        'bounds-z'              : [1e-5, 0.8],
        'true-values'           : {},
    }

    # Read options from config file
    for key in input_pars.keys():

        # Input
        if ('output' in key) or ('injections-path' in key) or ('selection-effects-cut' in key) or ('data-path' in key):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if ('injections-number' in key) or ('snr-cut' in key) or ('ifar-cut' in key):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if ('O3-cosmology' in key) or ('simulation' in key):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass

        # Model
        if ('model-primary' in key) or ('model-secondary' in key) or ('model-rate' in key) or ('redshift-transition' in key):
            try: input_pars[key] = Config.get('model', key)
            except: pass
        if ('positive-peak' in key) or ('low-smoothing' in key) or ('scale-free' in key) or ('conditional-prior' in key):
            try: input_pars[key] = Config.getboolean('model', key)
            except: pass
        if ('priors' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('model', key))
            except: pass

        # Sampler
        if ('sampler' in key):
            try: input_pars[key] = Config.get('sampler', key)
            except: pass
        if ('nparallel' in key) or ('neffPE' in key) or ('nlive' in key) or ('npool' in key):
            try: input_pars[key] = Config.getint('sampler', key)
            except: pass

        # Plots
        if ('N-points' in key) or ('N-z-slices' in key):
            try: input_pars[key] = Config.getfloat('plots', key)
            except: pass
        if ('bounds-m1' in key) or ('bounds-z' in key):
            try: input_pars[key] = Config.get('plots', key)
            except: pass
        if ('true-values' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass

    # Set output directory
    # FIXME: Add option to control that only one of the two between 'O3-cosmology' and 'simulation' is active.
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])

    # Save the config file to output
    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))

    # Deviate stdout and stderr to file
    if not input_pars['screen-output']:
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_icarogw.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_icarogw.txt'), 'w')
    else: pass
    print('\n\n Starting  i c a r o g w  runner\n\n')

    # Initialise the prior
    input_pars['all-priors'] = default_priors()
    if not input_pars['priors'] == {}:
        for key in input_pars['priors']: input_pars['all-priors'][key] = input_pars['priors'][key]

    # Print run parameters
    print(' * I will be running with the following parameters.\n')
    print_dictionary(input_pars)

    # Initialise the wrappers
    tmp = Wrappers(input_pars)
    m1w, m2w, rw, cw = tmp.return_Wrappers()

    tmp = Rate(input_pars, m1w, m2w, rw, cw)
    wrapper = tmp.return_Rate()

    tmp = SelectionEffects(input_pars)
    injections = tmp.return_SelectionEffects()

    tmp = Data(input_pars)
    data = tmp.return_Data()

    tmp = LikelihoodPrior(input_pars, data, injections, wrapper)
    likelihood, prior = tmp.return_LikelihoodPrior()

    # Run the hierarchical analysis
    print('\n * Running hierarchical analysis with this settings.\n')
    print_dictionary({key: input_pars[key] for key in ['sampler', 'nlive', 'naccept', 'npool', 'print_method', 'sample']})

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
    hierarchical.plot_corner()

    # Get the samples
    df = pd.DataFrame(hierarchical.samples, columns = hierarchical.search_parameter_keys)
    priors_dict = hierarchical.priors

    # Make plots
    print(' * Producing plots.\n')
    input_pars['output'] = os.path.join(input_pars['output'], 'plots')
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])

    # Primary mass
    if not 'Redshift' in input_pars['model-primary']:
        curves, plot_dict = icaroproc.PrimaryMassFunction(df, m1w, priors_dict, input_pars)
        if input_pars['true-values'] == {}:
            icaroproc.plot_curves(curves, plot_dict, logscale = True)
        else:
            curve_true, _ = icaroproc.PrimaryMassFunction(pd.DataFrame(input_pars['true-values'], index = [0]), m1w, priors_dict, input_pars)
            icaroproc.plot_curves(curves, plot_dict, truth = curve_true, logscale = True)
    else:
        curves, plot_dict = icaroproc.PrimaryMassFunction(df, m1w, priors_dict, input_pars)
        if input_pars['true-values'] == {}:
            icaroproc.plot_curves_evolving(curves, plot_dict)
        else:
            curve_true, _ = icaroproc.PrimaryMassFunction(pd.DataFrame(input_pars['true-values'], index = [0]), m1w, priors_dict, input_pars)
            icaroproc.plot_curves_evolving_long(curves, plot_dict, truth = curve_true)

    # Secondary mass
    curves, plot_dict = icaroproc.SecondaryMassFunction(df, m2w, priors_dict, input_pars)
    if input_pars['true-values'] == {}:
        icaroproc.plot_curves(curves, plot_dict, figsize = (8,8))
    else:
        curve_true, _ = icaroproc.SecondaryMassFunction(pd.DataFrame(input_pars['true-values'], index = [0]), m2w, priors_dict, input_pars)
        icaroproc.plot_curves(curves, plot_dict, figsize = (8,8), truth = curve_true[0])

    # Rate evolution
    curves, plot_dict = icaroproc.RateEvolutionFunction(df, rw, priors_dict, input_pars)
    if input_pars['true-values'] == {}:
        icaroproc.plot_curves(curves, plot_dict)
    else:
        curve_true, _ = icaroproc.RateEvolutionFunction(pd.DataFrame(input_pars['true-values'], index = [0]), rw, priors_dict, input_pars)
        icaroproc.plot_curves(curves, plot_dict, truth = curve_true[0])

    print('\n * Finished.\n')

if __name__=='__main__':
    main()
