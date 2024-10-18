import icarogw, bilby
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import tqdm, seaborn as sns
from scipy.stats import kde


# from matplotlib import rcParams
# from distutils.spawn import find_executable

#if find_executable('latex'): rcParams["text.usetex"] = True
# rcParams["xtick.labelsize"] = 18
# rcParams["ytick.labelsize"] = 13
# rcParams["xtick.direction"] = "in"
# rcParams["ytick.direction"] = "in"
# rcParams["legend.fontsize"] = 15
# rcParams["legend.frameon"]  = False
# rcParams["legend.loc"]      = "best"
# rcParams["axes.labelsize"]  = 24
# rcParams["axes.grid"]       = True
# rcParams["grid.alpha"]      = 0.6
# rcParams["grid.linestyle"]  = "dotted"
# rcParams["lines.linewidth"] = 0.7


def get_plot_parameters(pars, x_array, x_min, x_max, figname, color, x_label, y_label, label, colors = None, z_grid = None, y_label_R = None, y_label_L = None):

    plot_dict = {
        'x'         : x_array,
        'output'    : pars['output-plots'],
        'figname'   : figname,
        'color'     : color,

        'label'     : label,
        'x_min'     : x_min,
        'x_max'     : x_max,

        'z_grid'    : z_grid,
        'x_label'   : x_label,
        'y_label'   : y_label,

        # Only for redshift slices plots
        'colors'    : colors,
        'y_label_R' : y_label_R,
        'y_label_L' : y_label_L,
    }

    return plot_dict


def selection_effects_countour_level(x, y, ref_cosmo):

    y  = ref_cosmo.dl2z(y)
    x /= (1+y)

    k = kde.gaussian_kde([x, y])

    # Evaluate the KDE on a grid
    N = 100
    x_grid, y_grid = np.meshgrid(np.linspace(x.min(), x.max(), N), np.linspace(y.min(), y.max(), N))
    z = k(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)

    # Compute the credible interval
    sorted_z = np.sort(z.ravel())
    cumulative_sum = np.cumsum(sorted_z)
    cumulative_sum /= cumulative_sum[-1]  # Normalize to make it a CDF

    # Find the contour level for the credible interval
    credible_level = 0.1
    contour_level = sorted_z[np.searchsorted(cumulative_sum, credible_level)]

    return x_grid, y_grid, z, contour_level


def add_curves_to_dict(dictionary, x, y, label):

    dictionary[label] = {}
    dictionary[label]['x'] = x
    dictionary[label]['y'] = y



# -------------------------------- #
# Class with all the type of plots #
# -------------------------------- #

class PlotDistributions:

    def __init__(self):
        return 0

    def plot_curves(curves, pl_dct, logscale = False, figsize = (10,5), truth = np.array([0]), curves_prior = 0):

        _, ax = plt.subplots(figsize = figsize)

        percentiles = [50, 5, 16, 84, 95]
        if hasattr(curves_prior, "__len__"):
            MF = {}
            for perc in percentiles: MF[perc] = np.percentile(curves_prior, perc, axis = 0)
            ax.fill_between(pl_dct['x'], MF[5] , MF[95], color = '#AB7C41', alpha = 0.15)
            ax.fill_between(pl_dct['x'], MF[16], MF[84], color = '#AB7C41', alpha = 0.25, label = '$\mathrm{Prior}$')

        MF = {}
        for perc in percentiles: MF[perc] = np.percentile(curves, perc, axis = 0)
        ax.fill_between(pl_dct['x'], MF[5] , MF[95],   color = pl_dct['color'], alpha = 0.25)
        ax.fill_between(pl_dct['x'], MF[16], MF[84],   color = pl_dct['color'], alpha = 0.5, label = pl_dct['label'])
        ax.plot(        pl_dct['x'], MF[50], lw = 0.7, color = pl_dct['color'])

        if not (truth.all() == 0 and len(truth) == 1):
            ax.plot(pl_dct['x'], truth, lw = 0.3, color = '#494949')

        if logscale:
            plt.yscale('log')
            plt.ylim(1e-5, 1)

        plt.xlim( pl_dct['x_min'], pl_dct['x_max'])
        plt.xlabel(pl_dct['x_label'])
        plt.ylabel(pl_dct['y_label'])
        if pl_dct['figname'] == 'RateEvolutionDistribution_Probability': plt.ylim(-5, 7)

        plt.grid(linestyle='dotted')
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
        plt.close()


    def plot_curves_evolving(curves, pl_dct, ref_cosmo, truth = {}, curves_prior = {}, selection_effects = {}):

        _, ax = plt.subplots(figsize = (5, 9))

        if not selection_effects == {}:
            m1_grid, z_grid, height, contour_level = selection_effects_countour_level(selection_effects['mass_1'], selection_effects['luminosity_distance'], ref_cosmo)
            ax.contourf(m1_grid, z_grid, height, levels = [contour_level, height.max()], colors = '#C4B692', alpha = 0.2)

        if bool(curves_prior):
            for zi, z_array in enumerate(pl_dct['z_grid']):
                z = z_array[0]
                ax.fill_between(pl_dct['x'], curves_prior[zi][5] +z, curves_prior[zi][95]+z, color = '#AB7C41', alpha = 0.05)
                ax.fill_between(pl_dct['x'], curves_prior[zi][16]+z, curves_prior[zi][84]+z, color = '#AB7C41', alpha = 0.15)
    
        for zi, z_array in enumerate(pl_dct['z_grid']):
            z = z_array[0]
            ax.fill_between(pl_dct['x'], curves[zi][5] +z, curves[zi][95]+z, color = pl_dct['colors'][zi], alpha = 0.25)
            ax.fill_between(pl_dct['x'], curves[zi][16]+z, curves[zi][84]+z, color = pl_dct['colors'][zi], alpha = 0.5)
            #ax.plot(        pl_dct['x'], curves[zi][50]+z, lw = 0.7,         color = pl_dct['colors'][zi])
        
            if not truth == {}:
                ax.plot(    pl_dct['x'], truth[zi][50]+z,  lw = 0.5,         color = '#494949')

        ax.set_xlim( pl_dct['x_min'], pl_dct['x_max'])
        ax.set_xlabel(pl_dct['x_label'])
        ax.set_ylabel(pl_dct['y_label_L'])

        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
        plt.close()


    def plot_curves_evolving_log(curves, pl_dct, truth = {}, curves_prior = 0):

        nz = np.shape(pl_dct['z_grid'])[0]
        fig, ax = plt.subplots(nz, 1, figsize=(5, 2 * nz), sharex = True, constrained_layout = True)
        fig.set_constrained_layout_pads(hspace = 0.0, h_pad = 0.0) 

        if bool(curves_prior):
            for zi, z_array in enumerate(pl_dct['z_grid']):
                z = z_array[0]
                zi_inv = nz-1 - zi
                ax[zi_inv].fill_between(pl_dct['x'], curves_prior[zi][5] , curves_prior[zi][95], color = '#AB7C41', alpha = 0.05)
                ax[zi_inv].fill_between(pl_dct['x'], curves_prior[zi][16], curves_prior[zi][84], color = '#AB7C41', alpha = 0.15, label = '$\mathrm{Prior}$')       
    
        for zi, z_array in enumerate(pl_dct['z_grid']):
            z = z_array[0]
            zi_inv = nz-1 - zi
            ax[zi_inv].fill_between(pl_dct['x'], curves[zi][5] , curves[zi][95], color = pl_dct['colors'][zi], alpha = 0.25)
            ax[zi_inv].fill_between(pl_dct['x'], curves[zi][16], curves[zi][84], color = pl_dct['colors'][zi], alpha = 0.5)
            ax[zi_inv].plot(        pl_dct['x'], curves[zi][50], lw = 0.7,       color = pl_dct['colors'][zi], label = '$z={}$'.format(round(z, 2)))
        
            if not truth == {}:
                ax[zi_inv].plot(    pl_dct['x'], truth[zi][50],  lw = 0.3,       color = '#494949')
        
            ax[zi_inv].set_xlim(-3, 90)
            ax[zi_inv].set_yscale('log')
            ax[zi_inv].set_ylim(1e-5, 0.5)
            ax[zi_inv].set_ylabel(pl_dct['y_label_R'])
            ax[zi_inv].grid(linestyle = 'dotted', linewidth = 0.3)
            ax[zi_inv].legend(loc = 'best')
        ax[-1].set_xlabel(pl_dct['x_label'])

        fig.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname'] + '_logscale'), transparent = True)
        plt.close()


    def plot_curves_evolving_MassRedshift_Joint(curves, pl_dct):

        fig, ax = plt.subplots(1, 2, figsize=(10, 5))       
    
        for zi, z_array in enumerate(pl_dct['z_grid'][::-1]):
            z = z_array[0]
            ax[1].plot(        pl_dct['x'], curves[len(pl_dct['z_grid'])-1-zi][50],      lw = 1,   color = pl_dct['colors'][len(pl_dct['z_grid'])-1-zi])
            ax[0].fill_between(pl_dct['x'], curves[len(pl_dct['z_grid'])-1-zi][50]+z, z, lw = 0.7, color = pl_dct['colors'][len(pl_dct['z_grid'])-1-zi], alpha = 0.7)

        ax[0].set_xlim(0, 70)
        ax[0].set_xlabel(pl_dct['x_label'])
        ax[0].set_ylabel(pl_dct['y_label_L'])

        ax[1].set_xlim(0, 100)
        ax[1].set_xlabel(pl_dct['x_label'])
        ax[1].set_ylabel(pl_dct['y_label_R'])
        ax[1].set_yscale('log')
        ax[1].set_ylim(1e-5, 1)

        plt.legend()
        plt.tight_layout()
        fig.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
        plt.close()


# ----------------------------------------------- #
# Class to compute the reconstructed ditributions #
# ----------------------------------------------- #

class ReconstructDistributions:

    def __init__(self):
        return 0
    
    def PrimaryMassFunction(df, w, p_dct, pars, prior = False):

        mass_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])

        if 'Redshift' in pars['model-primary']:
            if prior:
                tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N-samps-prior']) for key in w.population_parameters}
                df  = pd.DataFrame(tmp)

            pdf = np.empty(shape = (pars['N-points'], pars['N-points']))
            curves = np.empty(shape = (len(df), pars['N-points']))

            percentiles = [50, 5, 16, 84, 95]
            curves_z = {zi: {pi: np.empty(pars['N-points']) for pi in percentiles} for zi in range(pars['N-z-slices'])}

            zx         = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
            zy         = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-z-slices'])
            _, z_grid  = np.meshgrid(zx, zy)

            colors = sns.color_palette('blend:#0A4F8A,#9F0C0C', pars['N-z-slices'])   # 'RdBu_r'
            zi = 0

            for z_array in tqdm.tqdm(z_grid, desc = 'Reconstructing primary distribution'):
                for idx, samp in df.iterrows():

                    samp_filt = {key: samp[key] for key in w.population_parameters}
                    w.update(**samp_filt)
                    pdf = w.pdf(mass_array, z_array)
                    if not (np.isnan(pdf).any()): curves[idx] = pdf
                    else: pass

                for perc in percentiles: curves_z[zi][perc] = np.percentile(curves, perc, axis = 0)
                zi += 1

            colors = sns.color_palette('blend:#0A4F8A,#9F0C0C', pars['N-z-slices'])
            plot_dict = get_plot_parameters(pars, mass_array, pars['bounds-m1'][0], pars['bounds-m1'][1], 'PrimaryMassDistribution', '#000000', '$m_1\ [M_{\odot}]$', '$z$', pars['model-primary'], colors = colors, z_grid = z_grid, y_label_L = '$z$', y_label_R = '$p(m_1)$')

            return curves_z, plot_dict

        else:
            if prior:
                tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N-samps-prior']) for key in w.population_parameters}
                df  = pd.DataFrame(tmp)

            pdf = np.empty(shape = (pars['N-points'], pars['N-points']))
            curves = np.empty(shape = (len(df), pars['N-points']))

            for idx, samp in df.iterrows():

                samp_filt = {key: samp[key] for key in w.population_parameters}
                w.update(**samp_filt)
                pdf = w.pdf(mass_array)
                curves[idx] = pdf

                plot_dict = get_plot_parameters(pars, mass_array, pars['bounds-m1'][0], pars['bounds-m1'][1], 'PrimaryMassDistribution', '#890C0A', '$m_1\ [M_{\odot}]$', '$p(m_1)$', pars['model-primary'])

            return curves, plot_dict


    def RemoveSelectionEffects(df, pars, rate_w, ref_cosmo, injections):

        N_samps     = len(df.index)
        # Number of samples to be extracted from the reconstructed
        # distribution of each PE sample, to compute the KDE.
        N_samps_KDE = pars['N-points-KDE']
        
        # Initialise arrays.
        mass_array  = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1] * (1+pars['bounds-z'][1]), pars['N-points'])
        mass2_array = np.linspace(pars['bounds-m2'][0], pars['bounds-m2'][1] * (1+pars['bounds-z'][1]), pars['N-points'])
        q_array     = np.linspace(pars['bounds-q'][0],  pars['bounds-q'][1],                            pars['N-points'])
        dL_array    = np.linspace(pars['bounds-dL'][0], pars['bounds-dL'][1],                           pars['N-points'])
        z_array     = np.linspace(pars['bounds-z'][0],  pars['bounds-z'][1],                            pars['N-points'])
        zy          = np.linspace(pars['bounds-z'][0],  pars['bounds-z'][1],                            pars['N-z-slices'])
        _, z_grid   = np.meshgrid(z_array, zy)

        if 'MassRatio' in pars['model-secondary']:
            pars['bounds-m2'] = pars['bounds-q']
            m2_array = q_array
        else:
            m2_array = mass2_array

        m1d         = np.zeros([N_samps, N_samps_KDE])
        m2d         = np.zeros([N_samps, N_samps_KDE])
        dL          = np.zeros([N_samps, N_samps_KDE])
        curves_m1d  = np.zeros([N_samps, pars['N-points']])
        curves_m2d  = np.zeros([N_samps, pars['N-points']])
        curves_dL   = np.zeros([N_samps, pars['N-points']])

        m1s_PDF     = {zi: np.zeros([N_samps, pars['N-points']]) for zi in range(pars['N-z-slices'])}
        curves_m2s  = np.zeros(     [N_samps, pars['N-points']])
        curves_z    = np.zeros(     [N_samps, pars['N-points']])
        z_array_kde = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])

        # Remove selection effects from the PE samples.
        for idx, samp in tqdm.tqdm(df.iterrows(), total = len(df.index), desc = 'Removing selection effects'):
            samp_filt = {key: samp[key] for key in rate_w.population_parameters}
            rate_w.update(**samp_filt)

            injections.update_weights(rate_w)
            tmp = injections.return_reweighted_injections(Nsamp = N_samps_KDE, replace = True)
            m1d[idx,:] = tmp['mass_1']
            if not pars['single-mass']:
                if   'MassRatio' in pars['model-secondary']: m2d[idx,:] = tmp['mass_ratio']
                elif 'PowerLaw'  in pars['model-secondary']: m2d[idx,:] = tmp['mass_2']
            dL[idx,:]  = tmp['luminosity_distance']

        # Compute KDE of the distribution for each PE sample.
        for i in tqdm.tqdm(range(N_samps), desc = 'Computing KDE detector frame distributions'):
            tmp          = kde.gaussian_kde(m1d[i,:])
            curves_m1d[i,:] = tmp.evaluate(mass_array)
            if not pars['single-mass']:
                tmp          = kde.gaussian_kde(m2d[i,:])
                curves_m2d[i,:] = tmp.evaluate(m2_array)
            tmp          = kde.gaussian_kde(dL[i,:])
            curves_dL[i,:]  = tmp.evaluate(dL_array)

        # Get the source frame distribution.
        m1s, m2s, zs = icarogw.conversions.detector2source(m1d, m2d, dL, ref_cosmo)

        # Redshift binning to plot the distribution p(m1|z) on redshift slices.
        m1s_z_binned = [{zi: np.empty([]) for zi in range(pars['N-z-slices'])} for i in range(N_samps)]
        for i in range(N_samps):
            z_binned = np.digitize(zs[i], zy.tolist())
            indices_dict = {value: np.where(z_binned == value)[0].tolist() for value in np.unique(z_binned)}
            for zi, bin in enumerate(indices_dict.keys()):
                m1s_z_binned[i][zi] = m1s[i][indices_dict[bin]]

        # Compute KDE of the distribution for each PE sample and redshift bin.
        for zi in range(pars['N-z-slices']):
            for i in range(N_samps):
                try:
                    tmp = kde.gaussian_kde(m1s_z_binned[i][zi])
                    m1s_PDF[zi][i,:] = tmp.evaluate(mass_array)
                except:
                    m1s_PDF[zi][i,:] = np.zeros(len(mass_array))

        for i in tqdm.tqdm(range(N_samps), desc = 'Computing KDE source frame distributions'):
            if not pars['single-mass']:
                tmp             = kde.gaussian_kde(m2s[i,:])
                curves_m2s[i,:] = tmp.evaluate(m2_array)
            tmp                 = kde.gaussian_kde(zs[i,:])
            curves_z[i,:]       = tmp.evaluate(z_array_kde)

        # Get confidence bundles.
        percentiles = [50, 5, 16, 84, 95]
        curves_z_m1s = {zi: {pi: np.empty([]) for pi in percentiles} for zi in range(pars['N-z-slices'])}
        for zi in m1s_PDF.keys():
            for perc in percentiles: curves_z_m1s[zi][perc] = np.percentile(m1s_PDF[zi], perc, axis = 0)

        # Get plot paramters
        plots_inputs = {
            'curves-m1d':   curves_m1d,   'curves-m2d': curves_m2d, 'curves-dL': curves_dL,
            'curves-z-m1s': curves_z_m1s, 'curves-m2s': curves_m2s, 'curves-z' : curves_z,
        }
        colors = sns.color_palette('blend:#0A4F8A,#9F0C0C', pars['N-z-slices'])
        plots_inputs['plot-dict-m1d'] = get_plot_parameters(pars, mass_array, pars['bounds-m1'][0], pars['bounds-m1'][1] * (1+pars['bounds-z'][1]), 'PrimaryMassDistribution_DetectorFrame',         '#0A4F8A', '$m_1\ [M_{\odot}]$', '$p(m_1)$', pars['model-primary'])
        plots_inputs['plot-dict-m2d'] = get_plot_parameters(pars, m2_array,   pars['bounds-m2'][0], pars['bounds-m2'][1],                           'SecondaryMassDistribution_DetectorFrame',       '#6A8820', '$m_2\ [M_{\odot}]$', '$p(m_2)$', pars['model-secondary'])
        plots_inputs['plot-dict-dL']  = get_plot_parameters(pars, dL_array,   pars['bounds-dL'][0], pars['bounds-dL'][1],                           'LuminosityDistranceDistribution_DetectorFrame', '#7E375B', '$d_L\ [Mpc]$',       '$p(d_L)$', pars['model-rate'])
        plots_inputs['plot-dict-m1s'] = get_plot_parameters(pars, mass_array, pars['bounds-m1'][0], pars['bounds-m1'][1],                           'PrimaryMassDistribution_NoSelectionEffects',    '#000000', '$m_1\ [M_{\odot}]$', '$z$'     , pars['model-primary'], colors = colors, z_grid = z_grid)
        plots_inputs['plot-dict-m2s'] = get_plot_parameters(pars, m2_array,   pars['bounds-m2'][0], pars['bounds-m2'][1],                           'SecondaryMassDistribution_NoSelectionEffects',  '#1F5623', '$m_2\ [M_{\odot}]$', '$p(m_2)$', pars['model-secondary'])
        plots_inputs['plot-dict-z']   = get_plot_parameters(pars, z_array,    pars['bounds-z'][0],  pars['bounds-z'][1],                            'RedshiftDistribution_NoSelectionEffects',       '#86042A', '$z$',                '$p(z)$'  , pars['model-rate'])

        return plots_inputs


    def SecondaryMassFunction(df, w, p_dct, pars, prior = False):

        if   'MassRatio' in pars['model-secondary']: bound = 'bounds-q'
        elif 'PowerLaw'  in pars['model-secondary']: bound = 'bounds-m2'
        else:
            raise ValueError('Unknown option for the secondary mass plot. Current implementation accounts for MassRatio and PowerLaw.')
        
        if prior:
            tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N-samps-prior']) for key in w.population_parameters}
            df  = pd.DataFrame(tmp)

        m_array = np.linspace(pars[bound][0], pars[bound][1], pars['N-points'])
        curves  = np.empty(shape = (len(df), pars['N-points']))
        pdf     = np.empty(shape = (pars['N-points']))

        for idx, samp in df.iterrows():
            samp_filt = {key: samp[key] for key in w.population_parameters}
            w.update(**samp_filt)
            if   'MassRatio' in pars['model-secondary']: pdf = w.pdf(m_array)
            elif 'PowerLaw'  in pars['model-secondary']: pdf = np.exp(w.prior.pdf2._log_pdf(m_array))
            curves[idx] = pdf

        plot_dict = get_plot_parameters(pars, m_array, pars[bound][0], pars[bound][1], 'SecondaryMassDistribution', '#2A4D00', '$m_2\ [M_{\odot}]$', '$p(m_2)$', pars['model-secondary'])

        return curves, plot_dict


    def RateEvolutionFunction(df, w, p_dct, pars, prior = False):

        if prior:
            tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N-samps-prior']) for key in w.population_parameters}
            df  = pd.DataFrame(tmp)

        z_array = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
        curves  = np.empty(shape = (len(df), pars['N-points']))

        for idx, samp in df.iterrows():

            samp_filt = {key: samp[key] for key in w.population_parameters}
            w.update(**samp_filt)        
            func = w.rate.log_evaluate(z_array)
            curves[idx] = func

        plot_dict = get_plot_parameters(pars, z_array, pars['bounds-z'][0], pars['bounds-z'][1], 'RateEvolutionFunction', '#825310', '$z$', '$ln[\Psi(z)/R_0]$', pars['model-rate'])

        return curves, plot_dict


    def RateEvolutionFunctionProbability(df, rw, cw, pars):

        z_array = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
        curves  = np.empty(shape = (len(df), pars['N-points']))

        for idx, samp in df.iterrows():

            samp_filt = {key: samp[key] for key in rw.population_parameters}
            rw.update(**samp_filt)
            func = np.exp(rw.rate.log_evaluate(z_array))
            if not pars['scale-free']: curves[idx] = func * samp.R0
            else:                      curves[idx] = func

            # Comoving volume and redshift (1/(1+z)*dV/dz)
            samp_filt = {key: samp[key] for key in cw.population_parameters}
            cw.update(**samp_filt)
            func = cw.cosmology.dVc_by_dzdOmega_at_z(z_array) * 4*np.pi / (1+z_array)
            curves[idx] *= func
            curves[idx] = np.log(curves[idx])
            curves[idx] -= 7

        plot_dict = get_plot_parameters(pars, z_array, pars['bounds-z'][0], pars['bounds-z'][1], 'RateEvolutionDistribution_Probability', '#AC9512', '$z$', '$\propto ln[p(z)]$', pars['model-rate'])

        return curves, plot_dict


    def RedshiftTransitionFunction(df, p_dct, pars, prior = False):

        z_array = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])

        if prior:
            if   pars['redshift-transition'] == 'sigmoid': tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N-samps-prior']) for key in ['zt', 'delta_zt', 'mix_z0']}
            elif pars['redshift-transition'] == 'linear':  tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N-samps-prior']) for key in ['mix_z0', 'mix_z1']}
            df  = pd.DataFrame(tmp)
            
        curves  = np.empty(shape = (len(df), pars['N-points']))

        for idx, samp in df.iterrows():
            if   pars['redshift-transition'] == 'sigmoid':         curves[idx] = icarogw.priors._mixed_sigmoid_function(        z_array, samp['zt'],     samp['delta_zt'], samp['mix_z0'])
            elif pars['redshift-transition'] == 'double-sigmoid':  curves[idx] = icarogw.priors._mixed_double_sigmoid_function( z_array, samp['zt'],     samp['delta_zt'], samp['mix_z0'], samp['mix_z1'])
            elif pars['redshift-transition'] == 'linear':          curves[idx] = icarogw.priors._mixed_linear_function(         z_array, samp['mix_z0'], samp['mix_z1'])
            elif pars['redshift-transition'] == 'linear-sinusoid': curves[idx] = icarogw.priors._mixed_linear_sinusoid_function(z_array, samp['mix_z0'], samp['mix_z1'],   samp['amp'],    samp['freq'])
        
        plot_dict = get_plot_parameters(pars, z_array, pars['bounds-z'][0], pars['bounds-z'][1], 'RedshiftTransitionFunction', '#212121', '$z$', '$\\sigma(z)$', pars['model-primary'])
        
        return curves, plot_dict


    # FIXME: Fix the 2D joint implementation.
    def MassRedshift_Joint(df, pars):

        raise ValueError('This plot is currently unavailable. Please fix it.')

        mass_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])

        m = initialize_PrimaryMass(df, pars['model'], pars)
        r = initialize_Rate(df, pars['rate'])
        c = icarogw.wrappers.FlatLambdaCDM_wrap(zmax = 20.)

        pdf = np.empty(shape = (pars['N-points'], pars['N-points']))
        percentiles = [50, 5, 16, 84, 95]
        curves = {zi: {pi: np.empty(shape = (pars['N-points'], pars['N-points'])) for pi in percentiles} for zi in range(pars['N-z-slices'])}
        p_grid = np.empty(shape = (len(df), pars['N-z-slices'], pars['N-points']))

        zx = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
        zy = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-z-slices'])
        _, zY = np.meshgrid(zx, zy)
        X,  Y  = np.meshgrid(mass_array, zy)

        for idx, samp in tqdm.tqdm(df.iterrows(), total = len(df), desc = 'Computing 2D mass-redshift joint distribution'):

            # Compute conditional 2D distribution p(m1|z)
            for zyi in range(pars['N-z-slices']):
                samp_filt = {key: samp[key] for key in m.population_parameters}
                m.update(**samp_filt)
                pdf = m.pdf(mass_array, zY[zyi])
                if not (np.isnan(pdf).any()): p_grid[idx][zyi] = pdf
                else: pass

            # Evolution rate (\psi)
            samp_filt = {key: samp[key] for key in r.population_parameters}
            r.update(**samp_filt)        
            func = np.exp(r.rate.log_evaluate(zy))
            tmp = func * samp.R0

            # Comoving volume and redshift (1/(1+z)*dV/dz)
            samp_filt = {key: samp[key] for key in c.population_parameters}
            c.update(**samp_filt)    
            func = c.cosmology.dVc_by_dzdOmega_at_z(zy) * 4*np.pi / (1+zy)
            tmp *= func

            # Multiply the conditional p(m1|z) by p(z) to get the joint 2D distribution
            for mi in range(pars['N-points']): p_grid[idx].T[mi] *= tmp

            plt.contourf(X, Y, p_grid[idx], alpha = 0.1, cmap = 'Blues')
        plt.xlim(0, 100)
        plt.ylim(0, pars['bounds-z'][1])
        plt.xlabel('$m_1\ [M_{\odot}]$')
        plt.ylabel('$z$')
        plt.savefig('{}/{}.pdf'.format(pars['output'], 'MassRedshift_Joint', transparent = True))

        for zi,zyi in enumerate(zy):
            for perc in percentiles:
                curves[zi][perc] = np.percentile(np.transpose(p_grid, (1,0,2))[zi], perc, axis = 0)
                if perc == 50: max = np.max(curves[zi][perc])
                curves[zi][perc] /= max * 10

        colors = sns.color_palette('blend:#0A4F8A,#9F0C0C', pars['N-z-slices'])
        plot_dict = get_plot_parameters(pars, mass_array, pars['bounds-m1'][0], pars['bounds-m1'][1], 'MassRedshift_JointDistribution_rescaled', '#000000', '$m_1\ [M_{\odot}]$', '$z$', pars['model-primary'], colors = colors, z_grid = zY, y_label_L = '$z$', y_label_R = '$p(m_1)$')

        return curves, plot_dict


# ------------------------------- #
# Class to generate all the plots #
# ------------------------------- #

class Plots:

    def __init__(self, pars, df, m1w, m2w, rw, cw, ref_cosmo, rate_w, priors, injections):

        self.pars      = pars
        self.df        = df
        self.m1w       = m1w
        self.m2w       = m2w
        self.rw        = rw
        self.cw        = cw
        self.rate_w    = rate_w
        self.ref_cosmo = ref_cosmo
        self.priors    = priors
        self.inj       = injections

        self.distributions = ReconstructDistributions
        self.plots         = PlotDistributions

        self.curves_dict = {}

    def PrimaryMass(self):

        if not 'Redshift' in self.pars['model-primary']:
            curves_prior = 0
            if self.pars['plot-prior']:
                curves_prior, _ = self.distributions.PrimaryMassFunction(self.df, self.m1w, self.priors, self.pars, prior = True)
            curves, plot_dict   = self.distributions.PrimaryMassFunction(self.df, self.m1w, self.priors, self.pars)
            add_curves_to_dict(self.curves_dict, plot_dict['x'], curves, plot_dict['figname'])
            if self.pars['true-values'] == {}:
                self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior, logscale = True)
            else:
                curve_true, _ = self.distributions.PrimaryMassFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.m1w, self.priors, self.pars)
                self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior, truth = curve_true, logscale = True)
        else:
            curves_prior = 0
            if self.pars['plot-prior']:
                curves_prior, _ = self.distributions.PrimaryMassFunction(self.df, self.m1w, self.priors, self.pars, prior = True)
            curves, plot_dict   = self.distributions.PrimaryMassFunction(self.df, self.m1w, self.priors, self.pars)
            # FIXME: Save evolving curves.
            add_curves_to_dict(self.curves_dict, plot_dict['x'], curves, plot_dict['figname'])
            if self.pars['true-values'] == {}:
                self.plots.plot_curves_evolving_log(curves, plot_dict, curves_prior = curves_prior)
                if not self.pars['selection-effects']: self.plots.plot_curves_evolving(curves, plot_dict, self.ref_cosmo)
                else:                                  self.plots.plot_curves_evolving(curves, plot_dict, self.ref_cosmo, selection_effects = self.inj.injections_data)                    
            else:
                curve_true, _ = self.distributions.PrimaryMassFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.m1w, self.priors, self.pars)
                self.plots.plot_curves_evolving_log( curves, plot_dict, curves_prior = curves_prior, truth = curve_true)
                if not self.pars['selection-effects']: self.plots.plot_curves_evolving(curves, plot_dict, self.ref_cosmo, truth = curve_true)
                else:                                  self.plots.plot_curves_evolving(curves, plot_dict, self.ref_cosmo, truth = curve_true, selection_effects = self.inj.injections_data)

    def SecondaryMass(self):

        curves_prior = 0
        if self.pars['plot-prior']:
            curves_prior, _ = self.distributions.SecondaryMassFunction(self.df, self.m2w, self.priors, self.pars, prior = True)
        curves, plot_dict   = self.distributions.SecondaryMassFunction(self.df, self.m2w, self.priors, self.pars)
        add_curves_to_dict(self.curves_dict, plot_dict['x'], curves, plot_dict['figname'])
        if self.pars['true-values'] == {}:
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior,)
        else:
            curve_true, _ = self.distributions.SecondaryMassFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.m2w, self.priors, self.pars)
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior, truth = curve_true[0])

    def RateEvolution(self):

        curves_prior = 0
        if self.pars['plot-prior']:
            curves_prior, _ = self.distributions.RateEvolutionFunction(self.df, self.rw, self.priors, self.pars, prior = True)
        curves, plot_dict   = self.distributions.RateEvolutionFunction(self.df, self.rw, self.priors, self.pars)
        add_curves_to_dict(self.curves_dict, plot_dict['x'], curves, plot_dict['figname'])
        if self.pars['true-values'] == {}:
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior)
        else:
            curve_true, _ = self.distributions.RateEvolutionFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.rw, self.priors, self.pars)
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior, truth = curve_true[0])

    def RateEvolutionProbability(self):

        curves_prior = 0
        if self.pars['plot-prior']:
            curves_prior, _ = self.distributions.RateEvolutionFunctionProbability(self.df, self.rw, self.cw, self.pars, prior = True)
        curves, plot_dict = self.distributions.RateEvolutionFunctionProbability(self.df, self.rw, self.cw, self.pars)
        add_curves_to_dict(self.curves_dict, plot_dict['x'], curves, plot_dict['figname'])
        if self.pars['true-values'] == {}:
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior)
        else:
            curve_true, _ = self.distributions.RateEvolutionFunctionProbability(pd.DataFrame(self.pars['true-values'], index = [0]), self.rw, self.cw, self.pars)
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior, truth = curve_true[0])
    
    def RedshiftTransition(self):

        curves_prior = 0
        if self.pars['plot-prior']:
            curves_prior, _ = self.distributions.RedshiftTransitionFunction(self.df, self.priors, self.pars, prior = True)
        curves, plot_dict = self.distributions.RedshiftTransitionFunction(self.df, self.priors, self.pars)
        add_curves_to_dict(self.curves_dict, plot_dict['x'], curves, plot_dict['figname'])
        if self.pars['true-values'] == {}:
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior)
        else:
            curve_true, _ = self.distributions.RedshiftTransitionFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.priors, self.pars)
            self.plots.plot_curves(curves, plot_dict, curves_prior = curves_prior, truth = curve_true[0])

    def NoSelectionEffects(self):

        plots_inputs   = self.distributions.RemoveSelectionEffects(self.df, self.pars, self.rate_w, self.ref_cosmo, self.inj)
        
        add_curves_to_dict(    self.curves_dict, plots_inputs['plot-dict-m1d']['x'], plots_inputs['curves-m1d'],   plots_inputs['plot-dict-m1d']['figname'])
        if not self.pars['single-mass']:
            add_curves_to_dict(self.curves_dict, plots_inputs['plot-dict-m2d']['x'], plots_inputs['curves-m2d'],   plots_inputs['plot-dict-m2d']['figname'])
        add_curves_to_dict(    self.curves_dict, plots_inputs['plot-dict-dL'][ 'x'], plots_inputs['curves-dL'],    plots_inputs['plot-dict-dL'][ 'figname'])
        add_curves_to_dict(    self.curves_dict, plots_inputs['plot-dict-m1s']['x'], plots_inputs['curves-z-m1s'], plots_inputs['plot-dict-m1s']['figname'])
        if not self.pars['single-mass']:
            add_curves_to_dict(self.curves_dict, plots_inputs['plot-dict-m2s']['x'], plots_inputs['curves-m2s'],   plots_inputs['plot-dict-m2s']['figname'])
        add_curves_to_dict(    self.curves_dict, plots_inputs['plot-dict-z'][  'x'], plots_inputs['curves-z'],     plots_inputs['plot-dict-z'][  'figname'])

        if self.pars['true-values'] == {}:
            self.plots.plot_curves(              plots_inputs['curves-m1d'],   plots_inputs['plot-dict-m1d'])
            if not self.pars['single-mass']:
                self.plots.plot_curves(          plots_inputs['curves-m2d'],   plots_inputs['plot-dict-m2d'])
            self.plots.plot_curves(              plots_inputs['curves-dL'],    plots_inputs['plot-dict-dL'])
            self.plots.plot_curves_evolving(     plots_inputs['curves-z-m1s'], plots_inputs['plot-dict-m1s'], self.ref_cosmo)
            if not self.pars['single-mass']:
                self.plots.plot_curves(          plots_inputs['curves-m2s'],   plots_inputs['plot-dict-m2s'])
            self.plots.plot_curves(              plots_inputs['curves-z'],     plots_inputs['plot-dict-z'])

        else:
            inputs_true = self.distributions.RemoveSelectionEffects(pd.DataFrame(self.pars['true-values'], index = [0]), self.pars, self.rate_w, self.ref_cosmo, self.inj)
            self.plots.plot_curves(              plots_inputs['curves-m1d'],   plots_inputs['plot-dict-m1d'],                 truth = inputs_true['curves-m1d'][0])
            if not self.pars['single-mass']:
                self.plots.plot_curves(          plots_inputs['curves-m2d'],   plots_inputs['plot-dict-m2d'],                 truth = inputs_true['curves-m2d'][0])
            self.plots.plot_curves(              plots_inputs['curves-dL'],    plots_inputs['plot-dict-dL'],                  truth = inputs_true['curves-dL'][0])
            self.plots.plot_curves_evolving(     plots_inputs['curves-z-m1s'], plots_inputs['plot-dict-m1s'], self.ref_cosmo, truth = inputs_true['curves-z-m1s'])
            if not self.pars['single-mass']:
                self.plots.plot_curves(          plots_inputs['curves-m2s'],   plots_inputs['plot-dict-m2s'],                 truth = inputs_true['curves-m2s'][0])
            self.plots.plot_curves(              plots_inputs['curves-z'],     plots_inputs['plot-dict-z'],                   truth = inputs_true['curves-z'][0])

    # Call the class functions to generate the plots.
    def ProducePlots(self):

        try:    self.PrimaryMass()
        except: pass
        try:    self.SecondaryMass()
        except: pass
        try:    self.RateEvolution()
        except: pass
        try:    self.RateEvolutionProbability()
        except: pass
        try:    self.NoSelectionEffects()
        except: pass
        try:    self.RedshiftTransition()
        except: pass

    def return_curves(self):
        return self.curves_dict
