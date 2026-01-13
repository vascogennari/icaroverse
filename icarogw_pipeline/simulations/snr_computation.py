import numpy as np
import os, logging
import time

import bilby
from astropy.time import Time


# -------------------------------------------------------------------------- #
#                            bilby implementation                            #
# -------------------------------------------------------------------------- #

# Utils
def clean_dict(d, keys):
    """Remove dictionary entries"""
    for key in keys:
        if key in d: d.pop(key)

def bring_to_next_power_of_2(x):
    if x <= 0: 
        raise ValueError("x must be positive")
    else:
        return np.power(2, np.ceil(np.log2(x)))

def chieff_from_wf_params(event_pars_dict):
    """
    Compute the effective spin parameter from the masses and spins
    
    Parameters
    ----------
    event_pars_dict: dict of float or array-like
        Dictionary containing masses and spins with following keys:
        (mass_1, mass_2, a_1, a_2, tilt_1, tilt_2)
    """
    mass_1, mass_2 = event_pars_dict['mass_1'], event_pars_dict['mass_2']
    a_1, a_2       = event_pars_dict['a_1']   , event_pars_dict['a_2']
    tilt_1, tilt_2 = event_pars_dict['tilt_1'], event_pars_dict['tilt_2']
    chi_1, chi_2 = a_1 * np.cos(tilt_1), a_2 * np.cos(tilt_2)
    chi_eff = (mass_1 * chi_1 + mass_2 * chi_2) / (mass_1 + mass_2)
    return chi_eff

def tell_me_ifos_on(run='opt'):
    """
    Function to return the IFOs on. Based on IGWN reported and expected duty cycles.

    * For O3, values are based on S.Mastrogiovanni's scripts for MICEcat MDC.

    * For O4, total runtime spans May 24th 2023 to October 7th 2025.
      Two commisionning breaks are taken into account (that contribute to the [] output),
      V1 starts O4 after 1st break, K1 starts after second break.
      The present function only implements duty cycles of the periods outside the commisionning breaks.
      (ie. to be consistent in a simulation the total observation time should be set to the total O4 duration, 
      minus the duration of the two commisionning breaks. NB: this is implemented in run_to_time_window())
      (see https://observing.docs.ligo.org/plan/ and https://wiki.ligo.org/LSC/JRPComm/ObsRun4)

    * For O5, 70% duty cycle is assumed for each detector independently.

    Parameters
    ----------
    run: str
        observing run, available: O3, O4, O5, opt (NB: opt yields 100% full network)

    Returns
    -------
    ifos: list
        list of IFOs on, random subset of [H1, L1, V1, K1] based on run's duty cycle
    """
    if run == 'opt':
        ifos = ['H1', 'L1', 'V1', 'K1']
    elif run == 'O3':
        lucky = np.random.rand()
        if   lucky < 0.50: ifos = ['H1', 'L1', 'V1']
        elif lucky < 0.64: ifos = ['H1', 'L1', '  ']
        elif lucky < 0.78: ifos = ['H1', '  ', 'V1']
        elif lucky < 0.92: ifos = ['  ', 'L1', 'V1']
        elif lucky < 0.94: ifos = ['H1', '  ', '  ']
        elif lucky < 0.96: ifos = ['  ', 'L1', '  ']
        elif lucky < 0.98: ifos = ['  ', '  ', 'V1']
        else:              ifos = ['  ', '  ', '  ']
    elif run == 'O4':
        lucky = np.random.rand()
        if   lucky < 0.0371: ifos = ['L1', 'H1', 'V1', 'K1']
        elif lucky < 0.2438: ifos = ['L1', 'H1', 'V1', '  ']
        elif lucky < 0.2858: ifos = ['L1', 'H1', '  ', 'K1']
        elif lucky < 0.5200: ifos = ['L1', 'H1', '  ', '  ']
        elif lucky < 0.5400: ifos = ['L1', '  ', 'V1', 'K1']
        elif lucky < 0.6513: ifos = ['L1', '  ', 'V1', '  ']
        elif lucky < 0.6739: ifos = ['L1', '  ', '  ', 'K1']
        elif lucky < 0.8000: ifos = ['L1', '  ', '  ', '  ']
        elif lucky < 0.8093: ifos = ['  ', 'H1', 'V1', 'K1']
        elif lucky < 0.8610: ifos = ['  ', 'H1', 'V1', '  ']
        elif lucky < 0.8715: ifos = ['  ', 'H1', '  ', 'K1']
        elif lucky < 0.9300: ifos = ['  ', 'H1', '  ', '  ']
        elif lucky < 0.9350: ifos = ['  ', '  ', 'V1', 'K1']
        elif lucky < 0.9628: ifos = ['  ', '  ', 'V1', '  ']
        elif lucky < 0.9685: ifos = ['  ', '  ', '  ', 'K1']
        else:                ifos = ['  ', '  ', '  ', '  ']
    elif run == 'O5':
        lucky = np.random.rand()
        if   lucky < 0.2744: ifos = ['L1', 'H1', 'V1', 'K1']
        elif lucky < 0.3920: ifos = ['L1', 'H1', 'V1', '  ']
        elif lucky < 0.5096: ifos = ['L1', 'H1', '  ', 'K1']
        elif lucky < 0.5600: ifos = ['L1', 'H1', '  ', '  ']
        elif lucky < 0.6776: ifos = ['L1', '  ', 'V1', 'K1']
        elif lucky < 0.7280: ifos = ['L1', '  ', 'V1', '  ']
        elif lucky < 0.7784: ifos = ['L1', '  ', '  ', 'K1']
        elif lucky < 0.8000: ifos = ['L1', '  ', '  ', '  ']
        elif lucky < 0.8686: ifos = ['  ', 'H1', 'V1', 'K1']
        elif lucky < 0.8980: ifos = ['  ', 'H1', 'V1', '  ']
        elif lucky < 0.9274: ifos = ['  ', 'H1', '  ', 'K1']
        elif lucky < 0.9400: ifos = ['  ', 'H1', '  ', '  ']
        elif lucky < 0.9694: ifos = ['  ', '  ', 'V1', 'K1']
        elif lucky < 0.9820: ifos = ['  ', '  ', 'V1', '  ']
        elif lucky < 0.9946: ifos = ['  ', '  ', '  ', 'K1']
        else:                ifos = ['  ', '  ', '  ', '  ']
    else:
        raise AssertionError(f"{run} is not an available option for tell_me_ifos_on(). Please choose from 'O3', 'O4', 'O5', 'opt'.")
    return ifos

def run_to_time_window(observing_run='O3'):
    """
    Computes the gps time for the start and end of the observing run.

    For O4, the starting date is artificially forwarded in time 
    by the duration fo the two commmisionning breaks
    (see https://observing.docs.ligo.org/plan/ and https://wiki.ligo.org/LSC/JRPComm/ObsRun4)
    
    Parameters
    ----------
    observing_run: str
        Choose from O3, O4, O5
    
    Returns
    -------
    (start, end): 2-tuple of floats
        starting and ending gps time of the observing_run 
    """
    if observing_run == 'O3':
        return (Time('2019-04-01').to_value('gps'), Time('2020-03-27').to_value('gps'))
    elif observing_run == 'O4':
        return (Time('2023-10-20').to_value('gps'), Time('2025-10-07').to_value('gps'))
    elif observing_run == 'O5':
        return (Time('2028-01-01').to_value('gps'), Time('2030-12-31').to_value('gps'))
    else:
        raise AssertionError(f"{observing_run} is not an available option for run_to_time_window(). Please choose from 'O3', 'O4', 'O5'.")

def there_is_fully_parametrised_spins(event_dict):
    """
    Check if event_dict contains all six parameters for the binary spins
    in (a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl) parametrization.

    Parameters
    ----------
    event_dict: dict
        dictionary containing CBC parameters

    Returns
    -------
    res: bool
        True if all six spin parameters are given in the event_dict
    """
    return (
        ('a_1'    in event_dict) and 
        ('a_2'    in event_dict) and 
        ('tilt_1' in event_dict) and 
        ('tilt_2' in event_dict) and 
        ('phi_12' in event_dict) and 
        ('phi_jl' in event_dict) 
    )



# Class to compute SNR & draw missing single event parameters
class BilbyDetectionPipeline():

    def __init__(self, psd_dir, observing_run, reference_frequency=20., sampling_frequency=2048., approximant='IMRPhenomXPHM', precessing_apx=True, duration=None, start_time=None):
        """
        A class that initialises a simulated event with some single event parameters, 
        drawing the missing ones, and injecting the signal 
        generated from the desired approximant in some detector network.

        Parameters
        ----------
        psd_dir: str
            absolute path of the folder where the PSD files are stored.
        observing_run: str
            LVK observing run. Used to load the correct PSD for bilby detector objects. Choose from O3, O4, O5.
        reference_frequency: float
            lower bound of the frequency range used to compute signals inner products in frequency domain
        sampling_frequency: float
            twice the higher bound of the frequency range used to compute signals inner products in frequency domain 
            (the latter often being the Nyquist frequency)
        approximant: str
            waveform approximant (default: IMRPhenomXPHM)
        """
        self.psd_dir = psd_dir
        self.observing_run = observing_run
        self.load_psd_from_file()
        self.reference_frequency = reference_frequency
        self.sampling_frequency = sampling_frequency
        self.approximant = approximant
        self.precessing_apx = precessing_apx
        self.duration = duration
        self.start_time = start_time

    def load_psd_from_file(self):
        self.all_ifos_available_psd_dict = {
            ifo: bilby.gw.detector.PowerSpectralDensity(
                asd_file = os.path.join(self.psd_dir, f"{ifo}_{self.observing_run}.txt")
            )
            for ifo in ['H1', 'L1', 'V1'] + ['K1']*(self.observing_run in {'O4', 'O5'})
        }

    def set_event_dict(self, init_dict, flag_set_duration_and_start_time=False):
        """
        Stocks all single event parameters in a dictionary eventually containing 
        - geocent_time
        - mass_1
        - mass_2
        - luminosity_distance
        - dec
        - ra
        - theta_jn
        - psi
        - phase
        - a_1
        - a_2
        - tilt_1
        - tilt_2
        - phi_12
        - phi_jl
        - ifos_on
        
        Parameter
        ---------
        init_dict: dict
            dictionary with single event parameters, 
            at least (mass_1, mass_2, luminosity_distance)
        """
        self.event_dict = init_dict
        if 'geocent_time' not in self.event_dict: 
            self.draw_geocent_time()
        if not there_is_fully_parametrised_spins(self.event_dict): 
            self.draw_spins()
        if ('ra' not in self.event_dict) or ('dec' not in self.event_dict):
            self.draw_skyloc()
        if 'theta_jn' not in self.event_dict: 
            self.draw_inclination()
        if 'psi' not in self.event_dict: 
            self.draw_polarization()
        if 'phase' not in self.event_dict: 
            self.draw_phase()
        if 'ifos_on' not in self.event_dict: 
            self.draw_ifos_on()

        if flag_set_duration_and_start_time:
            self.set_duration_and_start_time()

    def draw_geocent_time(self):
        """
        Draw GPS time uniformly in the time window associated to self.observing_run
        """
        t_gps_start, t_gps_end = run_to_time_window(self.observing_run)
        self.event_dict['geocent_time'] = np.random.uniform(t_gps_start, t_gps_end)

    def draw_spins(self):
        """
        Draw spins, with the (a_1, tilt_1, a_2, tilt_2, phi_12, phi_jl) parametrisation.
        Spin vectors distribution is uniform in magnitude and istropic in orientation (see https://zenodo.org/records/7890437)
        (see eg. http://corysimon.github.io/articles/uniformdistn-on-sphere/ for unit 2-sphere sampling)

        * Kerr parameters a are drawn such that a is uniform in [0, 1]
        * cos(tilt) is drawn uniformly from [-1, 1] 
          (So that spin orientation is sampled uniformly on the unit sphere), 
        * phi_12 and phi_jl are drawn uniformly from [0, 2*pi]
        """
        self.event_dict['a_1']    = np.random.uniform(0., 1.)
        self.event_dict['a_2']    = np.random.uniform(0., 1.)
        self.event_dict['tilt_1'] = np.arccos(np.random.uniform(-1., 1.))
        self.event_dict['tilt_2'] = np.arccos(np.random.uniform(-1., 1.))
        self.event_dict['phi_12'] = np.random.uniform(0., 2*np.pi)
        self.event_dict['phi_jl'] = np.random.uniform(0., 2*np.pi)

    def draw_skyloc(self):
        """
        Draw the sky localization uniformly on the unit sphere
        (see http://corysimon.github.io/articles/uniformdistn-on-sphere/)
        
        * rignt ascention is drawn uniformly from [0, 2*pi]
        * cos(pi/2 - declination) is drawn randomly in [-1, 1]
        """
        self.event_dict['ra']           = np.random.uniform(0., 2*np.pi)
        self.event_dict['dec']          = np.pi/2 - np.arccos(np.random.uniform(-1., 1.))

    def draw_inclination(self):
        """
        Draw the inclination so that binary orientation is uniform on unit sphere
        
        * cos(theta_jn) is drawn randomly in [-1, 1]
        """
        self.event_dict['theta_jn']     = np.arccos(np.random.uniform(-1., 1.))

    def draw_polarization(self):
        """
        Draw polarization angle uniformly in [0, 2*pi]
        """
        self.event_dict['psi']          = np.random.uniform(0., 2*np.pi)

    def draw_phase(self):
        """
        Draw reference phare uniformly in [0, 2*pi]
        """
        self.event_dict['phase']        = np.random.uniform(0., 2*np.pi)

    def draw_ifos_on(self):
        """
        Draw ifos_on based on detectors duty cycles.
        """
        self.event_dict['ifos_on']      = tell_me_ifos_on(run = self.observing_run)
    
    def set_duration_and_start_time(self):
        """Title says all"""
        time_to_merger = bilby.gw.utils.calculate_time_to_merger(
            self.reference_frequency,
            self.event_dict['mass_1'],
            self.event_dict['mass_2'],
            chi = chieff_from_wf_params(self.event_dict),
            safety = 1.1,
        )
        # The duration of signal considered is time_to_merger + 4s
        self.duration = bring_to_next_power_of_2(int(np.ceil(time_to_merger)) + 4)
        # The start time is set to be at self.duration +1s before the geocent_time of the event.
        self.start_time = self.event_dict['geocent_time'] - (self.duration - 1.)


    def check_all_parameters_present(self):
        """
        Checks that all the expected CBC parameters are present in self.event_dict
        """
        if not (
            ('geocent_time' in self.event_dict) and
            ('mass_1' in self.event_dict) and
            ('mass_2' in self.event_dict) and
            ('luminosity_distance' in self.event_dict) and
            ('dec' in self.event_dict) and
            ('ra' in self.event_dict) and
            ('theta_jn' in self.event_dict) and
            ('psi' in self.event_dict) and
            ('phase' in self.event_dict) and
            ('a_1' in self.event_dict) and
            ('a_2' in self.event_dict) and
            ('tilt_1' in self.event_dict) and
            ('tilt_2' in self.event_dict) and
            ('phi_12' in self.event_dict) and
            ('phi_jl' in self.event_dict) and
            ('ifos_on' in self.event_dict)
        ):
            raise KeyError("Some single event parameters are missing for SNR computation. Please consider reloading the event dict with the set_event_dict() method")

    def set_random_seed(self):
        """
        Sets numpy's RNG seed to the one tied to the event  in self.event_dict.
        If no seed is present in self.event_dict, draws one.
        """
        if 'seed' not in self.event_dict:
            self.event_dict['seed'] = int(time.time()) + np.random.randint(1000000)
        np.random.seed(self.event_dict['seed'])

    def set_ifos_list(self):
        """
        Initialize bilby's interferometers list 
        from the (already loaded) available PSDs
        stored in self.all_ifos_available_psd_dict.
        """
        ifos_on = [ifo for ifo in self.event_dict['ifos_on'] if ifo != '  ']
        self.ifos_list = bilby.gw.detector.InterferometerList(ifos_on)
        for i, ifo in enumerate(ifos_on):
            self.ifos_list[i].power_spectral_density = self.all_ifos_available_psd_dict[ifo]

    def set_frequency_mask(self):
        """
        Manually sets the frequency_mask attribute of bilby's Interferometers objects
        based on the frequency_array, the minimum and maximum frequencies of each Interferometer.
        """
        for ifo in self.ifos_list:
            ifo.frequency_mask = ((ifo.strain_data.minimum_frequency < self.ifos_list.frequency_array) & 
                                  (self.ifos_list.frequency_array < ifo.strain_data.maximum_frequency))

    def set_waveform_generator(self):
        # Setup the Bilby waveform
        waveform_arguments = dict(
            waveform_approximant = self.approximant,
            reference_frequency  = self.reference_frequency
        )
        self.waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
            sampling_frequency            = self.sampling_frequency, 
            duration                      = self.duration,
            frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole,
            parameter_conversion          = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments            = waveform_arguments
        )
    
    def set_strain_data_from_psd(self):
        """
        Set the strain data in all the detectors of the network from the loaded PSDs.
        """
        # Setup ifos strain data from psd
        self.set_random_seed()
        self.ifos_list.set_strain_data_from_power_spectral_densities(
            sampling_frequency = self.sampling_frequency, 
            duration           = self.duration,
            start_time         = self.start_time,
        )

    def set_strain_data_from_arrays(self, strain_data_list):
        """
        Set strain data in all of the detectors of the network from given arrays.

        Note that the strain_data_list should contain the same number of strain arrays as detectors in the network (i.e len(strain_data_list) == len(self.ifos_list))
        and to ensure consistent simulation:
        > strain_data_list[i] should be the strain to put in the interferometer self.ifos_list[i].
        > duration and start time should match the ones that would be computed from the event contained in self.event_dict().

        """
        # Setup ifos strain data from psd
        for ifo, strain_data in zip(self.ifos_list, strain_data_list):
            ifo.set_strain_data_from_frequency_domain_strain(
                strain_data, 
                sampling_frequency = self.sampling_frequency, 
                duration           = self.duration, 
                start_time         = self.start_time,
            )

    def inject_signal(self):
        """
        Inject GW signal into each detector of the network.
        """

        self.check_all_parameters_present()
        # update the frequency mask (it doesn't update by itself : brute force way of solving the issue)
        self.set_frequency_mask()

        # Inject signals
        if self.precessing_apx:
            # keep 3D spins fro precessing approximant.
            self.projected_event_dict = self.event_dict
            self.ifos_list.inject_signal(
                waveform_generator = self.waveform_generator,
                parameters = self.event_dict
            )
        else:
            # project spins along the z-axis in case of a non precessing waveform approximant
            self.projected_event_dict = self.event_dict.copy()
            self.projected_event_dict['chi_1'] = self.event_dict['a_1'] * np.cos(self.event_dict['tilt_1'])
            self.projected_event_dict['chi_2'] = self.event_dict['a_2'] * np.cos(self.event_dict['tilt_2'])
            clean_dict(self.projected_event_dict, ['a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])
            self.ifos_list.inject_signal(
                waveform_generator = self.waveform_generator,
                parameters = self.projected_event_dict
            )

    def compute_matched_filter_SNR(self):
        """
        Compute the matched filter SNR for the event with parameters stored in self.event_dict.
        """
        # Compute network SNR (with interferometers on)
        matched_filter_SNR = np.sqrt(np.sum([
            np.real(self.ifos_list.meta_data[ifo]['matched_filter_SNR'])**2. 
            for ifo in self.event_dict['ifos_on']
            if ifo != '  '
        ]))
        self.event_dict['matched_filter_SNR'] = matched_filter_SNR

        return self.event_dict


# -------------------------------------------------------------------------- #
#                          lisabeta  implementation                          #
# -------------------------------------------------------------------------- #

try:
    import lisabeta.lisa.lisa as lisa

    def SNR_lisabeta(m1d, q, dL, Tobs = 4.):

        SNR = []
        N = len(m1d)

        # Randomize the remaining parameters
        t0     = np.random.uniform(0., Tobs, N)
        inc    = np.random.uniform(0., np.pi, N)
        phi    = np.random.uniform(0., 2 * np.pi, N)
        lambd  = np.random.uniform(0., 2 * np.pi, N)
        beta   = np.random.uniform(np.pi/2, np.pi/2, N)
        psi    = np.random.uniform(0., 2 * np.pi, N)

        # Draw spin orientations unifrormly on a unit sphere
        # Keep only the aligned components for the SNR computation
        a_1    = np.random.uniform(0., 1., N)
        a_2    = np.random.uniform(0., 1., N)
        tilt_1 = np.arccos(np.random.uniform(-1., 1., N))
        tilt_2 = np.arccos(np.random.uniform(-1., 1., N))
        chi1 = a_1 * np.cos(tilt_1)
        chi2 = a_2 * np.cos(tilt_2)

        i = 0
        for mi, qi, di in zip(m1d, q, dL):
            
            M = mi + mi / qi

            params = {
                'M':      M,         # Total *redshifted* mass M=m1+m2 [M_odot]
                'q':      qi,        # Mass ratio q=m1/m2
                'dist':   di,        # Luminosity distance [Mpc]
                'chi1':   chi1[i],   # Dimensionless spin 1 component along orbital momentum
                'chi2':   chi2[i],   # Dimensionless spin 2 component along orbital momentum
                'inc':    inc[i],    # Inclination, observer's colatitude in source-frame
                'phi':    phi[i],    # Phase, observer's longitude in source-frame
                'lambda': lambd[i],  # Longitude in the sky
                'beta':   beta[i],   # Latitude in the sky
                'psi':    psi[i],    # Polarization angle

                'Deltat': 0.,        # Time shift of coalescence, s -- coalescence is at t0*yr + Deltat*s, t0 in waveform_params
                'Lframe': True       # Flag indicating whether angles and Deltat pertain to the L-frame or SSB-frame
            }
            
            waveform_params = {
                
                # Frequency range
                'minf':                  1e-10,
                'maxf':                  0.5,
                
                't0':                    t0[i],  # Reference epoch of coalescence, yr -- coalescence is at t0*yr + Deltat*s, Deltat in params
                'timetomerger_max':      Tobs,   # Always cut signals timetomerger_max*yr before merger -- to avoid needlessly long signals using minf
                'DeltatL_cut':           None,   # Option to cut the signal pre-merger -- must be in L-frame
                
                # Further options to cut signals
                'fend':                  None,
                'tmin':                  None,
                'tmax':                  None,

                # Options for the time and phase alignment -- development/testing
                'phiref':                0.0,
                'fref_for_phiref':       0.0,
                'tref':                  0.0,
                'fref_for_tref':         0.0,
                'toffset':               0.0,
                'force_phiref_fref':     True,
                
                # TDI channels to generate
                'TDI':                   'TDIAET',

                # Internal accuracy params
                'acc':                   1e-4,
                'order_fresnel_stencil': 0,

                # Waveform approximant and set of harmonics to use
                'approximant':           'IMRPhenomXHM',
                'modes':                 None,

                # LISA response options
                'LISAconst':             'Proposal',
                'responseapprox':        'full',
                'frozenLISA':            False,
                'TDIrescaled':           True,
                
                # Noise options -- can also be given as a numpy array for interpolation
                'LISAnoise': {
                    'InstrumentalNoise':       'SciRDv1',
                    'WDbackground':            True,
                    'WDduration' :             Tobs,
                    'lowf_add_pm_noise_f0':    0.0,
                    'lowf_add_pm_noise_alpha': 2.0
                }
            }
            
            # Generate the waveform
            tdi_signal = lisa.GenerateLISATDISignal_SMBH(params, **waveform_params)
            SNR.append(tdi_signal['SNR'])
            i += 1
        
        return np.array(SNR)

except (ModuleNotFoundError, ImportError):
    logging.warning("Failed to import lisabeta. (See https://pypi.org/project/lisabeta/)")

def cut_SNR(snr, snr_thr = 10):
    '''
    Apply a cut in snr : snr > snrthr
    Returns the indices of each event that fits the criterion.
    '''
    indices = np.where((snr >= snr_thr))[0]

    return indices


# -------------------------------------------------------------------------- #
#                            pycbc implementation                            #
# -------------------------------------------------------------------------- #

try:
    import pycbc.psd # PYCBC MUST BE IMPORTED BEFORE ICAROGW
    from pycbc.waveform import get_fd_waveform
    from pycbc.filter import get_cutoff_indices
    from pycbc.detector import Detector

    class Event():

        def __init__(self, flow, fhigh, delta_f, m1, m2, dL):
            """
            A class to store sampling frequencies and source parameters,
            and to generate waveform.
            Events are labeled by their (m1, m2, dL)
            
            Parameters
            ----------
            flow: float
                low frequency cut for waveform generation [Hz]
            fhigh: float
                high frequency cut for waveform generation [Hz]
            delta_f: float
                frequency step for waveform generation [Hz]
            m1: float
                primary mass (detector frame) [solar mass]
            m2: float
                secondary mass (detector frame) [solar mass]
            dL: float
                luminosity distance [Mpc]
            """
            self.wf_params = {}
            self.set_event_params(m1, m2, dL)
            self.set_sampling_frequencies_params(flow, fhigh, delta_f)

        def set_sampling_frequencies_params(self, flow, fhigh, delta_f):
            """
            Fill a dict that will be passed to pycbc to generate the waveform,
            with the sampling frequencies parameters

            Parameters
            ----------
            flow: float
                low frequency cut for waveform generation [Hz]
            fhigh: float
                high frequency cut for waveform generation [Hz]
            delta_f: float
                frequency step for waveform generation [Hz]
            """
            self.wf_params['f_lower'] = flow
            self.wf_params['delta_f'] = delta_f
            self.wf_params['f_final'] = fhigh

        def set_event_params(self, m1, m2, dL):
            """
            Fill a dict that will be passed to pycbc to generate the waveform,
            with the event masses and distance

            Parameters
            ----------
            m1: float
                primary mass (detector frame) [solar mass]
            m2: float
                secondary mass (detector frame) [solar mass]
            dL: float
                luminosity distance [Mpc]
            """
            self.wf_params['mass1'] = m1
            self.wf_params['mass2'] = m2
            self.wf_params['distance'] = dL
        
        def set_projection_params(self, t_gps=0):
            self.projection_params = {
                'polarization_type'   : 'tensor', 
                't_gps'               : t_gps,
            }
            self.draw_skyloc_polarization()

        def draw_skyloc_polarization(self):
            """
            Draw sky localization uniformly on the sphere (see http://corysimon.github.io/articles/uniformdistn-on-sphere/)
            - right ascension is uniform in [0, 2*pi], 
            - cos(pi/2 - declination) is uniform in [-1, 1]

            Draw polarization angle uniformly in [0, 2*pi]
            """
            self.projection_params['right_ascension'] = 2 *np.pi/2 * np.random.rand(), 
            self.projection_params['declination']     = np.pi/2 - np.arccos(2*np.random.rand() - 1), 
            self.projection_params['polarization']    = 2 *np.pi/2 * np.random.rand(), 


        def draw_spins(self):
            """
            Draw spins of each binary component uniformly on the sphere (see http://corysimon.github.io/articles/uniformdistn-on-sphere/). 
            - azimuthal angle is uniform in [0, 2*pi], 
            - cos(colatitude) is uniform in [-1, 1]

            If the chosen appoximant is a non-precessing one, transverse components are set to 0
            """
            theta_s1, phi_s1 = np.arccos(2*np.random.rand() - 1), 2 * np.pi * np.random.rand()
            theta_s2, phi_s2 = np.arccos(2*np.random.rand() - 1), 2 * np.pi * np.random.rand()
            self.wf_params['spin1x'] = np.sin(theta_s1) * np.cos(phi_s1) * self.precessing
            self.wf_params['spin1y'] = np.sin(theta_s1) * np.cos(phi_s1) * self.precessing
            self.wf_params['spin1z'] = np.cos(theta_s1)
            self.wf_params['spin2x'] = np.sin(theta_s2) * np.cos(phi_s2) * self.precessing
            self.wf_params['spin2y'] = np.sin(theta_s2) * np.cos(phi_s2) * self.precessing
            self.wf_params['spin2z'] = np.cos(theta_s2)

        def draw_inclination(self, precessing=True):
            """
            Inclination is drawn such that the binary orientation is uniform 
            on the unit sphere (so cos(inclination) is uniform in [-1, 1])
            """
            incl = np.arccos(2*np.random.rand() - 1)
            self.wf_params['inclination'] = incl

        def generate_waveform(self, approximant='IMRPhenomXHM', precessing=False):
            """
            Computes the h_+ and h_x polarizations in FD using the specified approximant
            and the stored source parameters (see pycbc.waveform.get_fd_waveform)

            Parameters
            ----------
            approximant: str
                waveform model (see pycbc for available models. default='IMRPhenomXHM')
            precessing: bool
                True if approximant labels a precessing waveform model, False otherwise
            """
            self.precessing = precessing
            self.draw_inclination()
            self.draw_spins()
            self.wf_params['approximant'] = approximant

            self.hptilde, self.hctilde = get_fd_waveform(**self.wf_params)


    class Detector_custom(Detector):

        def __init__(self, name, flow, delta_f, sample_rate, psd_directory):
            super().__init__(name)
            self.delta_f = delta_f
            self.flow = flow
            self.sample_rate = sample_rate

            self.fhigh = self.sample_rate
            self.flen = int(self.sample_rate / self.delta_f)

            self.psd_directory = psd_directory

        def load_psd(self, observing_run):
            """
            Load the PSD of detector with id name, for observing period observing_run

            Parameters
            ----------
            observing_run: str
                Observing period ('O4', 'O5')
            """
            filename = f"{self.name}_{observing_run}.txt"
            filepath = os.path.join(self.psd_directory, filename)
            self.psd = pycbc.psd.from_txt(
                filepath,
                self.flen,
                self.delta_f,
                self.flow,
                is_asd_file=True,
            )
        
        def compute_optimal_snr(self, event):
            """
            The given event's waveform is put against the detector's PSD
            and the (optimal) SNR is computed. The event object must have generated waveforms.
            Geometric parameters (sky localization and polarization angle) are drawn at random:
            
            Parameters
            ----------
            event: Event() object
                The event of which to compute the SNR, with generated waveforms
            projection_params: dict
                Dictionnary containing right_ascension, declination, polarization, polarization_type, t_gps
            
            Return
            ------
            snr: float
                SNR of the event on the detector
            """
            fp, fc = self.antenna_pattern(**event.projection_params)
            htilde = fp * event.hptilde + fc * event.hctilde

            # Retrieve the index range in which to compute the SNR
            N = (len(htilde) - 1) * 2
            kmin, kmax = get_cutoff_indices(self.flow, self.fhigh, self.delta_f, N)

            # Computing the suaqre of the SNR. Note the factor 4 because 
            # of the definition of the psd which is the single-sided one.
            snr_sqr = 4 * ((htilde[kmin:kmax] / self.psd[kmin:kmax]).inner( htilde[kmin:kmax] ) * htilde.delta_f).real
            
            # Store and return the SNR (sqrt of above result)
            snr = np.sqrt(snr_sqr)
            return snr


    class DetectorNetwork():

        def __init__(self, observing_run, flow, delta_f, sample_rate, psd_directory, network=['H1', 'L1', 'V1', 'K1']):
            """
            A class to store information about a network of detector.

            Init parameters
            ---------------
            observing_run: str
                The observing period of which to use the PSD ()
            flow: float
                low frequency cut for waveform generation [Hz]
            delta_f: float
                frequency step for waveform generation [Hz]
            sample_rate: float
                Rate at which the signal would be sampled in the detectors [Hz]
            network: list of str
                List of 2-characters detectors' names 
                (default full 4 detectors network: ['H1', 'L1', 'V1', 'K1'])
            """

            self.flow = flow
            self.delta_f = delta_f
            self.sample_rate = sample_rate

            self.fhigh = self.sample_rate
            self.flen = int(self.sample_rate / self.delta_f)

            self.psd_directory = psd_directory

            self.network = {name: Detector_custom(name, self.flow, self.delta_f, self.sample_rate, self.psd_directory) for name in network}
            self.snrs = {name: 0. for name in network}
            self.run = observing_run

        def load_psds(self):
            """
            Load the PSD of detectors, for observing period self.run
            """
            for det in self.network:
                self.network[det].load_psd(self.run)

        def hit_network(self, m1, m2, dL, t_gps, approximant='IMRPhenomXHM', precessing=False, snr_method='opt'):
            """
            Generate Event() object with source parameter (m1, m2, dL). 
            Give the event as an imput to each detector of the network.
            - the sky localisation is drawn uniformly on the unit sphere, 
            - the polarization angle is drawn uniformly in [0, 2*pi]

            Parameters
            ----------
            m1: float
                primary mass (detector frame) [solar mass]
            m2: float
                secondary mass (detector frame) [solar mass]
            dL: float
                luminosity distance [Mpc]
            t_gps: float
                The time at which the event hits the network, LVK time convention [s]
            approximant: str
                waveform model (see pycbc for available models. default='IMRPhenomXHM')
            precessing: bool
                True if approximant labels a precessing waveform model, False otherwise (default: False)

            Returns
            -------
            snr_network: float
                SNR of the event in the network of detector 
                (SNR_network^2 = SNR_det1^2 + SNR_det2^2 + ...) 
            """

            if snr_method == 'opt' or snr_method == 'mf_fast':

                self.event = Event(
                    flow=self.flow,
                    fhigh=self.fhigh,
                    delta_f=self.delta_f,
                    m1=m1,
                    m2=m2,
                    dL=dL,
                )
                # the other single event parameters are initialized in the following line:
                self.event.generate_waveform(approximant=approximant, precessing=precessing)

                # the skylocaliszation and polarization are initilized in the following line
                self.event.set_projection_params(t_gps=t_gps)

                for det in self.network:
                    self.snrs[det] = self.network[det].compute_optimal_snr(self.event)
                    if snr_method == 'mf_fast': self.snrs[det] += np.random.normal()
                
                snr_network = np.sqrt( sum( self.snrs[det]**2 for det in self.network ) )
                return snr_network
            
            elif snr_method == 'mf_full':
                ################################## YET TO IMPLEMENT
                return 0.
            
            else:
                raise AssertionError(f"Unknown snr computation method: {snr_method}")

except (ModuleNotFoundError, ImportError, NameError):
    logging.warning(" Failed to import PyCBC. (See https://pycbc.org/)")

