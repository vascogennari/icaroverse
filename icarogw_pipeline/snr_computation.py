import numpy as np

import pycbc.psd # PYCBC MUST BE IMPORTED BEFORE ICAROGW
from pycbc.waveform import get_fd_waveform
from pycbc.filter import get_cutoff_indices
from pycbc.detector import Detector


# /!\ S O O N   D E P R E C A T E D /!\
class SignaltoNoiseRatio():
    """
    A class to compute SNR from CBC waveform and detector's PSD in FD
    """

    def __init__(self, flow, fhigh, delta_f, psd_str='aLIGOZeroDetHighPower'):
        """
        Initialize the class: sampling frequencies and PSD
        """
        self.flow = flow
        self.fhigh = fhigh
        self.delta_f = delta_f
        self.generate_psd_from_pycbc(psd_str=psd_str)
        self.dets = {ifo: Detector(ifo) for ifo in ['H1', 'L1', 'V1', 'K1']}
        self.opt_snr_per_det = {ifo: 0. for ifo in ['H1', 'L1', 'V1', 'K1']}
        pass

    def generate_psd_from_pycbc(self, psd_str):
        sample_rate = 1024
        flen = int(sample_rate / self.delta_f) + 1
        self.psd = pycbc.psd.from_string(psd_str, flen, self.delta_f, self.flow)

    def set_wf_params(self, m1, m2, dL):
        """
        Events are labeled by their (m1, m2, dL)
        Other event parameters (intrinsic and extrinsinc) are drawn randomly:
        - spins are drawn randomly on the unit sphere 
          (azimuthal angle is uniform in [0, 2*pi], 
          cos(colatitude) is uniform in [-1, 1])
        - inclination is drawn such that the binary orientation is
          uniform on the unit sphere (so cos(inclination) is uniform in [-1, 1])
        
        Parameters
        ----------
        m1: float
            primary mass (detector frame) [solar mass]
        m2: float
            secondary mass (detector frame) [solar mass]
        dL: float
            luminosity distance [Mpc]
        """
        theta_s1, phi_s1 = np.arccos(2*np.random.rand() - 1), 2 * np.pi * np.random.rand()
        theta_s2, phi_s2 = np.arccos(2*np.random.rand() - 1), 2 * np.pi * np.random.rand()
        incl = np.arccos(2*np.random.rand() - 1)
        self.wf_params = {
            'mass1'         : m1,
            'mass2'         : m2,
            # 'spin1x'        : np.sin(theta_s1) * np.cos(phi_s1),
            # 'spin1y'        : np.sin(theta_s1) * np.cos(phi_s1),
            'spin1z'        : np.cos(theta_s1),
            # 'spin2x'        : np.sin(theta_s2) * np.cos(phi_s2),
            # 'spin2y'        : np.sin(theta_s2) * np.cos(phi_s2),
            'spin2z'        : np.cos(theta_s2),
            'distance'      : dL,
            'inclination'   : incl,
            'f_lower'       : self.flow,
            'delta_f'       : self.delta_f
        }

    def generate_waveform(self, approximant='IMRPhenomXHM'):
        """
        Computes the h_+ and h_x polarizations in FD using the specified approximant
        and the stored source parameters (see pycbc.waveform.get_fd_waveform)

        Parameters
        ----------
        approximant: str
            waveform model (see pycbc for available models.
            default='IMRPhenomXHM')
        """
        try:
            self.wf_params['approximant'] = approximant
        except AttributeError:
            self.set_wf_params(2**(1./5)*30, 2**(1./5)*30, 1500)
            self.wf_params['approximant'] = approximant
        self.hptilde, self.hctilde = get_fd_waveform(**self.wf_params)

    def compute_optimal_snr_single_det(self, det_ifo, projection_params):
        """
        Computes the optimal SNR (see e.g. Maggiore's GW Vol. 1, Eq 7.51)
        with stored waveform and psd, performing the computation in the range 
        [self.flow, self.fhigh].

        Parameters
        ----------
        det_ifo: str
            detector 2 characters id ['H1', 'L1', 'V1', 'K1']
        projection_params: dict
            keys : 
        """
        # Project the GW prolarisations onto the detector response
        fp, fc = self.dets[det_ifo].antenna_pattern(**projection_params)
        htilde = fp * self.hptilde + fc * self.hctilde

        # Retrieve the index range in which to compute the SNR
        N = (len(htilde) - 1) * 2
        kmin, kmax = get_cutoff_indices(self.flow, self.fhigh, self.delta_f, N)

        # Computing the suaqre of the SNR. Note the factor 4 because 
        # of the definition of the psd which is the single-sided one.
        snr_sqr = 4 * ((htilde[kmin:kmax] / self.psd[kmin:kmax]).inner( htilde[kmin:kmax] ) * htilde.delta_f).real
        
        # Store and return the SNR (sqrt of above result)
        snr = np.sqrt(snr_sqr)
        self.opt_snr_per_det[det_ifo] = snr

    def compute_optimal_snr_multiple_dets(self, dets_ifos, t_gps):
        """
        Computes separately the SNR for each given detector. 
        The projection parameters are computed once : 
        - the sky localisation is drawn uniformly on the unit sphere, 
        - the polarization angle is drawn uniformly in [0, 2*pi]
        """
        projection_params = {
            'right_ascension'     : 2 *np.pi/2 * np.random.rand(), 
            'declination'         : np.pi/2 - np.arccos(2*np.random.rand() - 1), 
            'polarization'        : 2 *np.pi/2 * np.random.rand(), 
            'polarization_type'   : 'tensor', 
            't_gps'               : t_gps,
        }
        for ifo in dets_ifos:
            self.compute_optimal_snr_single_det(ifo, projection_params)

    def compute_optimal_snr(self, dets_ifos, t_gps):
        """
        Computes the SNR for any given network of LVK detector
        """
        self.compute_optimal_snr_multiple_dets(dets_ifos, t_gps)
        self.opt_snr = np.sqrt( sum(self.opt_snr_per_det[ifo]**2 for ifo in dets_ifos) )
        return self.opt_snr


class Event():

    def __init__(self, flow, fhigh, delta_f, m1, m2, dL):
        """
        Initialize the class: sampling frequencies and source parameters

        Events are labeled by their (m1, m2, dL)
        Other event parameters (intrinsic and extrinsinc) are drawn randomly:
        - spins are drawn randomly on the unit sphere 
          (azimuthal angle is uniform in [0, 2*pi], 
          cos(colatitude) is uniform in [-1, 1])
        - inclination is drawn such that the binary orientation is
          uniform on the unit sphere (so cos(inclination) is uniform in [-1, 1])
        
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
        self.flow = flow
        self.fhigh = fhigh
        self.delta_f = delta_f
        self.m1 = m1
        self.m2 = m2
        self.dL = dL
        self.set_wf_params

    def set_wf_params(self):
        theta_s1, phi_s1 = np.arccos(2*np.random.rand() - 1), 2 * np.pi * np.random.rand()
        theta_s2, phi_s2 = np.arccos(2*np.random.rand() - 1), 2 * np.pi * np.random.rand()
        incl = np.arccos(2*np.random.rand() - 1)
        self.wf_params = {
            'mass1'         : self.m1,
            'mass2'         : self.m2,
            # 'spin1x'        : np.sin(theta_s1) * np.cos(phi_s1),
            # 'spin1y'        : np.sin(theta_s1) * np.cos(phi_s1),
            'spin1z'        : np.cos(theta_s1),
            # 'spin2x'        : np.sin(theta_s2) * np.cos(phi_s2),
            # 'spin2y'        : np.sin(theta_s2) * np.cos(phi_s2),
            'spin2z'        : np.cos(theta_s2),
            'distance'      : self.dL,
            'inclination'   : incl,
            'f_lower'       : self.flow,
            'delta_f'       : self.delta_f
        }

    def generate_waveform(self, approximant='IMRPhenomXHM', precessing=False):
        """
        Computes the h_+ and h_x polarizations in FD using the specified approximant
        and the stored source parameters (see pycbc.waveform.get_fd_waveform)

        Parameters
        ----------
        approximant: str
            waveform model (see pycbc for available models.
            default='IMRPhenomXHM')
        """
        try:
            self.wf_params['approximant'] = approximant
        except AttributeError:
            self.set_wf_params(2**(1./5)*30, 2**(1./5)*30, 1500)
            self.wf_params['approximant'] = approximant
        self.hptilde, self.hctilde = get_fd_waveform(**self.wf_params)


class DetectorNetwork():
    """
    A class to store information about a network of detector.
    For now hard-coded as ['H1', 'L1', 'V1', 'K1']

    Init parameters
    ---------------
    observing_run: str
        The observing period of hwich to use the PSD ()
    """

    def __init__(self, observing_run):
        """"""

        pass

