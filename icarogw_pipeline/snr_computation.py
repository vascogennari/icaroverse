import numpy as np
import os
import time

import pycbc.psd # PYCBC MUST BE IMPORTED BEFORE ICAROGW
from pycbc.waveform import get_fd_waveform
from pycbc.filter import get_cutoff_indices
from pycbc.detector import Detector


# -------------------------------------------------------------------------- #
#                            pycbc implementation                            #
# -------------------------------------------------------------------------- #


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


# -------------------------------------------------------------------------- #
#                            bilby implementation                            #
# -------------------------------------------------------------------------- #


import bilby
from gwpy.time import to_gps

# Utils
def chieff_from_wf_params(event_pars_dict):
    
    mass_1, mass_2 = event_pars_dict['mass_1'], event_pars_dict['mass_2']
    a_1, a_2       = event_pars_dict['a_1']   , event_pars_dict['a_2']
    tilt_1, tilt_2 = event_pars_dict['tilt_1'], event_pars_dict['tilt_2']
    q = mass_2 / mass_1
    chi_eff = (mass_1 * a_1 * np.cos(tilt_1) + mass_2 * a_2 * np.cos(tilt_2)) / (mass_1 + mass_2)
    return chi_eff

def tell_me_ifos_on(run='opt'):
    """
    Function to return the IFOs on. Based on IGWN reported and expected duty cycles.

    > For O3, values are based on S. Mastrogiovanni's scripts for MICEcat MDC.

    > For O4, total runtime spans May 24th 2023 to October 7th 2025.
      Two commisionning breaks are taken into account (that contribute to the [] output),
      V1 starts O4 after 1st break, K1 starts after second break.
      (see https://observing.docs.ligo.org/plan/ and https://wiki.ligo.org/LSC/JRPComm/ObsRun4)

    > For O5, 70% duty cycle is assumed for each detector independently.

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
        elif lucky < 0.64: ifos = ['H1', 'L1'      ]
        elif lucky < 0.78: ifos = ['H1',       'V1']
        elif lucky < 0.92: ifos = [      'L1', 'V1']
        elif lucky < 0.94: ifos = ['H1'            ]
        elif lucky < 0.96: ifos = [      'L1'      ]
        elif lucky < 0.98: ifos = [            'V1']
        else:              ifos = [                ]
    elif run == 'O4':
        lucky = np.random.rand()
        if   lucky < 0.0174: ifos = ['L1', 'H1', 'V1', 'K1']
        elif lucky < 0.1385: ifos = ['L1', 'H1', 'V1'      ]
        elif lucky < 0.1660: ifos = ['L1', 'H1',       'K1']
        elif lucky < 0.3566: ifos = ['L1', 'H1'            ]
        elif lucky < 0.3716: ifos = ['L1',       'V1', 'K1']
        elif lucky < 0.4754: ifos = ['L1',       'V1'      ]
        elif lucky < 0.4990: ifos = ['L1',             'K1']
        elif lucky < 0.6625: ifos = ['L1'                  ]
        elif lucky < 0.6714: ifos = [      'H1', 'V1', 'K1']
        elif lucky < 0.7331: ifos = [      'H1', 'V1'      ]
        elif lucky < 0.7471: ifos = [      'H1',       'K1']
        elif lucky < 0.8442: ifos = [      'H1'            ]
        elif lucky < 0.8518: ifos = [            'V1', 'K1']
        elif lucky < 0.9047: ifos = [            'V1'      ]
        elif lucky < 0.9167: ifos = [                  'K1']
        else:                ifos = [                      ]
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

def run_to_time_window(run='O3'):
    if run == 'O3':
        return (to_gps('April   01 2019'), to_gps('March    27 2020'))
    elif run == 'O4':
        return (to_gps('May     24 2023'), to_gps('October  07 2025'))
    elif run == 'O5':
        return (to_gps('January 01 2028'), to_gps('December 01 2030'))
    else:
        raise AssertionError(f"{run} is not an available option for run_to_time_window(). Please choose from 'O3', 'O4', 'O5'.")

def there_is_fully_parametrised_spins(event_dict):
    """Check if event_dict contains all six parameters for the binary spins"""
    return (
        ('a_1'    in event_dict) and 
        ('a_2'    in event_dict) and 
        ('tilt_1' in event_dict) and 
        ('tilt_2' in event_dict) and 
        ('phi_12' in event_dict) and 
        ('phi_jl' in event_dict) 
    )


# Class to compute SNR & draw missing single event parameters
class BilbySNR():

    def __init__(self, psd_dir, observing_run):
        self.psd_dir = psd_dir
        self.observing_run = observing_run
        self.load_psd_from_file()

    def load_psd_from_file(self):
        self.all_ifos_available_psd_dict = {
            ifo: bilby.gw.detector.PowerSpectralDensity(
                asd_file = os.path.join(self.psd_directory, f"{ifo}_{self.observing_run}.txt")
            )
            for ifo in ['H1', 'L1', 'V1', 'K1']
        }

    def set_event_dict(self, init_dict):
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
            dictionary with single event parameters, at least (mass_1, mass_2, luminosity_distance)"""
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

    def draw_geocent_time(self):
        """Draw GPS time uniformly in the time window associated to self.observing_run"""
        t_gps_start, t_gps_end = run_to_time_window(self.observing_run)
        return np.random.uniform(t_gps_start, t_gps_end)

    def draw_spins(self):
        """
        Draw spins, with the a_1, tilt_1, a_2, tilt_2, phi_12, phi_jl parametrisation.

        > dimensionless spin parameters a are drawn uniformly from [0, 1]
        > cos(tilt) is drawn uniformly from [-1, 1] 
          (So that spin orientation is sampled uniformly on the unit sphere, 
          see http://corysimon.github.io/articles/uniformdistn-on-sphere/)
        > phi_12 and phi_jl are drawn uniformly from [0, 2*pi]
        """
        self.event_dict['a_1']      = np.random.uniform(0., 1.)
        self.event_dict['tilt_1']   = np.arccos(np.random.uniform(-1., 1.))
        self.event_dict['a_2']      = np.random.uniform(0., 1.)
        self.event_dict['tilt_2']   = np.arccos(np.random.uniform(-1., 1.))
        self.event_dict['phi_12']   = np.random.uniform(0., 2*np.pi)
        self.event_dict['phi_jl']   = np.random.uniform(0., 2*np.pi)

    def draw_skyloc(self):
        """
        Draw the sky localization uniformly on the unit sphere
        (see http://corysimon.github.io/articles/uniformdistn-on-sphere/)
        
        > rignt ascention is drawn uniformly from [0, 2*pi]
        > cos(pi/2 - declination) is drawn randomly in [-1, 1]
        """
        self.event_dict['ra']       = np.random.uniform(0., 2*np.pi)
        self.event_dict['dec']      = np.pi/2 - np.arccos(np.random.uniform(-1., 1.))

    def draw_inclination(self):
        """
        Draw the inclination so that binary orientation is uniform on unit sphere
        
        > cos(theta_jn) is drawn randomly in [-1, 1]
        """
        self.event_dict['theta_jn'] = np.arccos(np.random.uniform(-1., 1.))

    def draw_polarization(self):
        """Draw polarization angle uniformly in [0, 2*pi]"""
        self.event_dict['psi']      = np.random.uniform(0., 2*np.pi)

    def draw_phase(self):
        """Draw reference phare uniformly in [0, 2*pi]"""
        self.event_dict['phase']    = np.random.uniform(0., 2*np.pi)

    def draw_ifos_on(self):
        """
        Draw ifos_on based on detectors duty cycles.
        """
        self.event_dict['ifos_on']  = tell_me_ifos_on('opt')#run = self.observing_run)
    
    def set_ifos_list(self):
        """
        Initialize bilby's interferometers list.
        """
        self.ifos_list = bilby.gw.detector.InterferometerList(self.event_dict['ifos_on'])
        for i, ifo in enumerate(self.event_dict['ifos_on']):
            self.ifos_list[i].power_spectral_density = self.all_ifos_available_psd_dict[ifo]


    def compute_snr(self, reference_frequency=20., sampling_frequency=2048., approx='IMRPhenomXPHM'):

        # Calculate the time between the moment the binary's frequency enters the detector's band (default: 20 Hz) and the merger
        time_to_merger = bilby.gw.utils.calculate_time_to_merger(
            reference_frequency,
            self.event_dict['mass_1'],
            self.event_dict['mass_2'],
            chi = chieff_from_wf_params(self.event_dict),
            safety = 1.1,
        )
        # The duration of signal considered is time_to_merger + 4s
        duration = int(np.ceil(time_to_merger)) + 4

        # Setup the Bilby waveform
        waveform_arguments = dict(
            waveform_approximant = approx,
            reference_frequency  = reference_frequency
        )
        self.waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
            sampling_frequency            = sampling_frequency, 
            duration                      = duration,
            frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole,
            parameter_conversion          = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments            = waveform_arguments
        )

        # Setup ifos list & respective psd
        self.set_ifos_list()

        # Setup ifos strain data from psd
        self.set_random_seed()
        self.ifos.set_strain_data_from_power_spectral_densities(
            sampling_frequency = sampling_frequency, 
            duration           = duration,
            start_time         = self.event_dict['geocent_time'] - (duration - 1.)
        )

        # update the frequency mask (it doesn't update by itself : brute force way of solving the issue)
        self.set_frequency_mask()

        # Inject signals
        self.ifos.inject_signal(
            waveform_generator = self.waveform_generator,
            parameters = self.event_dict
        )

        # Compute network SNR (with interferometers on)
        matched_filter_SNR = np.sqrt(np.sum([
            np.real(self.ifos.meta_data[ifo]['matched_filter_SNR'])**2. 
            for ifo in self.event_dict['ifos_on']
        ]))
        self.event_dict['matched_filter_SNR'] = matched_filter_SNR


    def set_frequency_mask(self):
        for ifo in self.ifos:
            ifo.frequency_mask = ((ifo.strain_data.minimum_frequency < self.ifos.frequency_array) & 
                                  (self.ifos.frequency_array < ifo.strain_data.maximum_frequency))
    
    def set_random_seed(self):
        if 'seed' not in self.event_dict:
            self.event_dict['seed'] = int(time.time()) + np.random.randint(1000000)
        np.random.seed(self.event_dict['seed'])

