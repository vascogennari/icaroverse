import numpy as np

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

    def draw_spins(self):
        """
        Draw spins of each binary component uniformly on the sphere. 
        - azimuthal angle is uniform in [0, 2*pi], 
        - cos(colatitude) is uniform in [-1, 1])

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

    def __init__(self, name, flow, delta_f, sample_rate):
        super().__init__(name)
        self.delta_f = delta_f
        self.flow = flow
        self.sample_rate = sample_rate

        self.fhigh = self.sample_rate
        self.flen = int(self.sample_rate / self.delta_f)

    def load_psd(self, observing_run):
        """
        Load the PSD of detector with id name, for observing period observing_run

        Parameters
        ----------
        observing_run: str
            Observing period ('O4', 'O5')
        """
        filename = f"/Users/tbertheas/Documents/icarogw_pipeline/data/psd/{self.name}_{observing_run}.txt"
        self.psd = pycbc.psd.from_txt(
            filename,
            self.flen,
            self.delta_f,
            self.flow,
            is_asd_file=True,
        )
    
    def hit_detector(self, event, projection_params):
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
        fp, fc = self.antenna_pattern(**projection_params)
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

    def __init__(self, observing_run, flow, delta_f, sample_rate, network=['H1', 'L1', 'V1', 'K1']):
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
        # self.network = {
        #     'H1': Detector_custom('H1'), # LIGO-Handford
        #     'L1': Detector_custom('L1'), # LIGO-Livingstone
        #     'V1': Detector_custom('V1'), # Virgo
        #     'K1': Detector_custom('K1'), # KAGRA
        # }
        # self.snrs = {'H1': 0., 'L1': 0., 'V1': 0., 'K1': 0.}
        self.flow = flow
        self.delta_f = delta_f
        self.sample_rate = sample_rate

        self.fhigh = self.sample_rate
        self.flen = int(self.sample_rate / self.delta_f)

        self.network = {name: Detector_custom(name, self.flow, self.delta_f, self.sample_rate) for name in network}
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
            SNR for of the event in the network of detector 
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
            self.event.generate_waveform(approximant=approximant, precessing=precessing)

            projection_params = {
                'right_ascension'     : 2 *np.pi/2 * np.random.rand(), 
                'declination'         : np.pi/2 - np.arccos(2*np.random.rand() - 1), 
                'polarization'        : 2 *np.pi/2 * np.random.rand(), 
                'polarization_type'   : 'tensor', 
                't_gps'               : t_gps,
            }

            for det in self.network:
                self.snrs[det] = self.network[det].hit_detector(self.event, projection_params)
                if snr_method == 'mf_fast': self.snrs[det] += np.random.normal()
            
            snr_network = np.sqrt( sum( self.snrs[det]**2 for det in self.network ) )
            return snr_network
        
        elif snr_method == 'mf_full':
            # YET TO IMPLEMENT
            return 0.
        
        else:
            raise AssertionError(f"Unknown snr computation method: {snr_method}")
