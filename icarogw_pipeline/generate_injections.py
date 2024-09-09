import numpy as np, pickle, os
import icarogw.simulation as icarosim
import matplotlib.pyplot as plt
from tqdm import tqdm


def compute_theta(m1, path):
     
    theta = icarosim.rvs_theta(Nsamp = len(m1), a = 0, b = 1.4, path = path)
    return theta

def save_truths(path, dictionary):
    
    with open(os.path.join(path, 'true_parameters.txt'), 'w') as f:
        for key in dictionary.keys():
            max_len = len(max(dictionary.keys(), key = len))
            f.write('{}  {}\n'.format(key.ljust(max_len), dictionary[key]))

def plot_injections(true_param, plot_dir):

    title   = 'm1_detector_frame'
    figname = os.path.join(plot_dir, title)
    plt.hist(true_param['m1d'], bins = 50, color = 'k', alpha = 0.5)
    plt.title(title)
    plt.xlabel('m1')
    plt.savefig('{}.pdf'.format(figname))
    plt.savefig('{}.png'.format(figname))
    plt.close()

    title   = 'm1dL_detector_frame'
    figname = os.path.join(plot_dir, title)
    plt.scatter(true_param['m1d'], true_param['dL'], marker = '.', s = 0.1, c = 'k')
    plt.title(title)
    plt.xlabel('m1')
    plt.ylabel('dL')
    plt.savefig('{}.png'.format(figname))
    plt.close()

    return 0


# Injections set 
dic_param ={
    'alpha':       1.,
    'beta':        1.,
    'mmin':        1.,
    'mmax':        300.,
    'delta_m':     1.,
    'mu_g':        30.,
    'sigma_g':     10.,
    'lambda_peak': 0.8,
}

m_model   = 'PowerLawPeak'
Ninj      = 100000  #200000
Ndet_inj  = 0
N_ext     = 10e4
snr_thr   = 12
fgw_cut   = 15
unif_dist = 0
unif_z    = 0
zmax      = 1.5
flat_PSD  = 0
additional_text = ''

print('\n * Generating injections for selection effects.\n')

inj_name = 'inj_{}_N{}_SNR{}_fGW{}{}'.format(m_model, int(Ninj), snr_thr, fgw_cut, additional_text)

base_dir = '/Users/vgennari/Documents/work/code/python/icarogw/data/simulations'
results_dir = os.path.join(base_dir,    'injections_selection_effects', inj_name)
plot_dir    = os.path.join(results_dir, 'injections_plots')
if not os.path.exists(results_dir): os.makedirs(results_dir)
if not os.path.exists(plot_dir   ): os.makedirs(plot_dir)

save_truths(results_dir, dic_param)

true_param = {'m1d': [],'m2d': [],'dL': [],'prior': [],'snr': []}

c = 0
with tqdm(total = Ninj) as pbar:
    while Ndet_inj < Ninj:

        m1s_inj, m2s_inj, pdf_m          = icarosim.generate_mass_inj(        Nsamp = int(N_ext), mass_model = m_model, dic_param = dic_param)
        if not unif_z:
            if unif_dist: dL_inj, pdf_dL = icarosim.generate_dL_inj_uniform(  Nsamp = int(N_ext), zmax = zmax)
            else:         dL_inj, pdf_dL = icarosim.generate_dL_inj(          Nsamp = int(N_ext), zmax = zmax)
        else:
                          dL_inj, pdf_dL = icarosim.generate_dL_inj_z_uniform(Nsamp = int(N_ext), zmax = zmax)
        z_inj                            = icarosim.dl_to_z(dL_inj)

        # The injection values are extracted in the source frame, but the output values are transformed in the detector frame,
        # because the hierarchical likelihood is expressed in the detector frame, and so the input injections for selection effects.
        # This transformation needs to be tracked in the injections prior, by including the Jacobian (m1s, m2s) --> (m1d, m2d).
        # Note that we directly extract from the luminosity distance, thus there is no ddL/dz contribution.
        # FIXME: Implement a more general framework to account for different parameters: mass ratio or only one mass.
        jacobian                     = (1 + z_inj)**(-2)
        prior_inj                    = pdf_m * pdf_dL * jacobian
        if not flat_PSD:
            theta                    = compute_theta(m1s_inj, os.path.join(base_dir, 'Pw_three.dat'))
            snr_inj, _, _            = icarosim.snr_samples(m1s_inj, m2s_inj, z_inj, numdet = 3, rho_s = 9, dL_s = 1.5, Md_s = 25, theta = theta)
            idx_detected_inj         = icarosim.snr_and_freq_cut(m1s_inj, m2s_inj, z_inj, snr_inj, snrthr = snr_thr, fgw_cut = fgw_cut)
        else:
            snr_inj                  = icarosim.snr_samples_flat(z_inj)
            idx_detected_inj         = icarosim.snr_cut_flat(snr_inj, snrthr = snr_thr)
        
        m1d_inj = m1s_inj[idx_detected_inj] * (1 + z_inj[idx_detected_inj])
        m2d_inj = m2s_inj[idx_detected_inj] * (1 + z_inj[idx_detected_inj])
        
        Ndet     = len(idx_detected_inj)
        Ndet_inj = Ndet_inj + Ndet
        c = c+1
        
        true_param['m1d'].append(  m1d_inj)
        true_param['m2d'].append(  m2d_inj)
        true_param['dL'].append(   dL_inj[idx_detected_inj])
        true_param['prior'].append(prior_inj[idx_detected_inj])
        true_param['snr'].append(  snr_inj[idx_detected_inj])

        if Ndet_inj > Ninj: pbar.update(Ninj)
        else:               pbar.update(Ndet)

for keys in true_param.keys():
    true_param[keys] = np.concatenate(true_param[keys], axis = 0)

with open(os.path.join(results_dir, '{}.pickle'.format(inj_name)), 'wb') as handle:
    pickle.dump(true_param, handle, protocol = pickle.HIGHEST_PROTOCOL)

plot_injections(true_param, plot_dir)
print('\n * Generated {} injections.'.format(N_ext * c))
with open(os.path.join(results_dir, 'number_injections.txt'), 'w') as f:
    f.write('Generated: {}\n'.format(int(N_ext * c)))
    f.write('Detected:  {}'.format(  len(true_param['m1d'])))

print('\n * Finished.\n')
