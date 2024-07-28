import os, icarogw, bilby
import json, pandas as pd, numpy as np
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



def selection_effects_countour_level(x, y):

    from astropy.cosmology import FlatLambdaCDM

    cosmo_ref = icarogw.cosmology.astropycosmology(zmax = 20.)
    cosmo_ref.build_cosmology(FlatLambdaCDM(H0 = 67.7, Om0 = 0.308))

    y  = cosmo_ref.dl2z(y)
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
# ------------------------------------------------------------- #

# ------------------------------------------------------------- #
def plot_curves(curves, pl_dct, logscale = False, figsize = (10,5), truth = np.array([0]), curves_prior = 0):

    _, ax = plt.subplots(figsize = figsize)

    percentiles = [50, 5, 16, 84, 95]

    if not curves_prior == 0:
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

    plt.grid(linestyle='dotted')
    plt.legend()
    plt.tight_layout()
    plt.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
    plt.close()


def plot_curves_evolving_long(curves, pl_dct, truth = {}, curves_prior = {}, selection_effects = {}):

    _, ax = plt.subplots(figsize = (5, 9))

    if not selection_effects == {}:
        m1_grid, z_grid, height, contour_level = selection_effects_countour_level(selection_effects['mass_1'], selection_effects['luminosity_distance'])
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

    ax.set_xlim(0, 70)
    ax.set_xlabel(pl_dct['x_label'])
    ax.set_ylabel(pl_dct['y_label_a'])

    plt.legend()
    plt.tight_layout()
    plt.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
    plt.close()


def plot_curves_evolving(curves, pl_dct, truth = {}, curves_prior = {}):

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))

    if bool(curves_prior):
        for zi, z_array in enumerate(pl_dct['z_grid']):
            z = z_array[0]
            ax[0].fill_between(pl_dct['x'], curves_prior[zi][5] +z, curves_prior[zi][95]+z, color = '#AB7C41', alpha = 0.05)
            ax[0].fill_between(pl_dct['x'], curves_prior[zi][16]+z, curves_prior[zi][84]+z, color = '#AB7C41', alpha = 0.15)
            if zi == 0:
                ax[1].fill_between(pl_dct['x'], curves_prior[zi][5],  curves_prior[zi][95], color = '#AB7C41', alpha = 0.05)
                ax[1].fill_between(pl_dct['x'], curves_prior[zi][16], curves_prior[zi][84], color = '#AB7C41', alpha = 0.15, label = '$\mathrm{Prior}$')
            if zi == len(pl_dct['z_grid'])-1:
                ax[1].fill_between(pl_dct['x'], curves_prior[zi][5],  curves_prior[zi][95], color = '#AB7C41', alpha = 0.05)
                ax[1].fill_between(pl_dct['x'], curves_prior[zi][16], curves_prior[zi][84], color = '#AB7C41', alpha = 0.15)           
 
    for zi, z_array in enumerate(pl_dct['z_grid']):
        z = z_array[0]
        ax[0].fill_between(pl_dct['x'], curves[zi][5] +z, curves[zi][95]+z, color = pl_dct['colors'][zi], alpha = 0.25)
        ax[0].fill_between(pl_dct['x'], curves[zi][16]+z, curves[zi][84]+z, color = pl_dct['colors'][zi], alpha = 0.5)
        ax[0].plot(        pl_dct['x'], curves[zi][50]+z, lw = 0.7,         color = pl_dct['colors'][zi])
        ax[1].plot(        pl_dct['x'], curves[zi][50],   lw = 1,           color = pl_dct['colors'][zi])
    
        if not truth == {}:
            ax[0].plot(    pl_dct['x'], truth[zi][50]+z,  lw = 0.3,         color = '#494949')
            ax[1].plot(    pl_dct['x'], truth[zi][50],    lw = 0.3,         color = '#494949')
    
    # ------------------------------------------------------------------ #
    # White plot slides/poster
    # for zi, z_array in enumerate(pl_dct['z_grid']):
    #     z = z_array[0]
    #     #ax[0].fill_between(pl_dct['x'], curves[zi][16]+z, curves[zi][84]+z, color = '#E5E3C8', alpha = 0.2)
    #     ax[0].plot(        pl_dct['x'], curves[zi][50]+z, lw = 1.5,         color = pl_dct['colors'][zi])
    #     ax[1].plot(        pl_dct['x'], curves[zi][50],   lw = 1,           color = pl_dct['colors'][zi])

    # ax[0].grid(False)
    # #ax[1].grid(False)
    # import matplotlib.colors as clr
    # cmap = clr.LinearSegmentedColormap.from_list('custom blue', ['#E5E3C8','#914E63'], N=256)
    # sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = 0, vmax = 1))
    # cb = plt.colorbar(sm, aspect = 40)
    # cb.outline.set_edgecolor('#E5E3C8')
    # cb.ax.yaxis.set_tick_params(color = '#E5E3C8')
    # cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    # plt.setp(cbar_yticks, color = '#E5E3C8')

    # for a in ax:
    #     a.tick_params(color = '#E5E3C8', labelcolor = '#E5E3C8')
    #     for spine in a.spines.values():
    #         spine.set_edgecolor('#E5E3C8')
    # ax[1].tick_params(color = '#E5E3C8', labelcolor = '#E5E3C8')
    # ax[1].patch.set_facecolor('#E5E3C8')
    # ------------------------------------------------------------------ #

    ax[0].set_xlim(0, 70)
    ax[0].set_xlabel(pl_dct['x_label'])
    ax[0].set_ylabel(pl_dct['y_label_a'])

    ax[1].set_xlim(-3, 90)
    ax[1].set_xlabel(pl_dct['x_label'])
    ax[1].set_ylabel(pl_dct['y_label_b'])
    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-5, 0.5)

    plt.legend()
    plt.tight_layout()
    fig.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
    plt.close()


def plot_curves_evolving_MassRedshift_Joint(curves, pl_dct):

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))       
 
    for zi, z_array in enumerate(pl_dct['z_grid'][::-1]):
        z = z_array[0]
        ax[1].plot(        pl_dct['x'], curves[len(pl_dct['z_grid'])-1-zi][50],      lw = 1,   color = pl_dct['colors'][len(pl_dct['z_grid'])-1-zi])
        ax[0].fill_between(pl_dct['x'], curves[len(pl_dct['z_grid'])-1-zi][50]+z, z, lw = 0.7, color = pl_dct['colors'][len(pl_dct['z_grid'])-1-zi], alpha = 0.7)

    ax[0].set_xlim(0, 70)
    ax[0].set_xlabel(pl_dct['x_label'])
    ax[0].set_ylabel(pl_dct['y_label_a'])

    ax[1].set_xlim(0, 100)
    ax[1].set_xlabel(pl_dct['x_label'])
    ax[1].set_ylabel(pl_dct['y_label_b'])
    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-5, 1)

    plt.legend()
    plt.tight_layout()
    fig.savefig('{}/{}.pdf'.format(pl_dct['output'], pl_dct['figname']), transparent = True)
    plt.close()
# ------------------------------------------------------------- #

# ------------------------------------------------------------- #
def PrimaryMassFunction(df, w, p_dct, pars, prior = False):

    mass_array = np.linspace(pars['bounds-m1'][0], pars['bounds-m1'][1], pars['N-points'])

    if 'Redshift' in pars['model-primary']:
        if prior:
            tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N_samp_prior']) for key in w.population_parameters}
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

        for z_array in tqdm.tqdm(z_grid):
            for idx, samp in df.iterrows():

                samp_filt = {key: samp[key] for key in w.population_parameters}
                w.update(**samp_filt)
                pdf = w.pdf(mass_array, z_array)
                if not (np.isnan(pdf).any()): curves[idx] = pdf
                else: pass

            for perc in percentiles: curves_z[zi][perc] = np.percentile(curves, perc, axis = 0)
            zi += 1

        plot_dict = {
            'x'         : mass_array,
            'output'    : pars['output'],
            'figname'   : 'PrimaryMassFunction',
            'colors'    : colors,
            'z_grid'    : z_grid,
            'x_label'   : '$m_1\ [M_{\odot}]$',
            'y_label_a' : '$z$',
            'y_label_b' : '$p(m_1)$',
            'label'     : pars['model-rate'],
            'x_min'     : pars['bounds-m1'][0],
            'x_max'     : pars['bounds-m1'][1],
        }

        return curves_z, plot_dict

    else:
        if prior:
            tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N_samp_prior']) for key in w.population_parameters}
            df  = pd.DataFrame(tmp)

        pdf = np.empty(shape = (pars['N-points'], pars['N-points']))
        curves = np.empty(shape = (len(df), pars['N-points']))

        for idx, samp in df.iterrows():

            samp_filt = {key: samp[key] for key in w.population_parameters}
            w.update(**samp_filt)
            pdf = w.pdf(mass_array)
            curves[idx] = pdf

            plot_dict = {
                'x'       : mass_array,
                'output'  : pars['output'],
                'figname' : 'PrimaryMassFunction',
                'color'   : '#890C0A',
                'x_label' : '$m_1\ [M_{\odot}]$',
                'y_label' : '$p(m_1)$',
                'label'   : pars['model-primary'],
                'x_min'   : pars['bounds-m1'][0],
                'x_max'   : pars['bounds-m1'][1],
            }

        return curves, plot_dict


def SecondaryMassFunction(df, w, p_dct, pars, prior = False):

    if prior:
        tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N_samp_prior']) for key in w.population_parameters}
        df  = pd.DataFrame(tmp)

    q_min = pars['all-priors']['mu_q'][0]
    q_max = pars['all-priors']['mu_q'][1]
    q_array = np.linspace(q_min, q_max, pars['N-points'])
    curves  = np.empty(shape = (len(df), pars['N-points']))
    pdf     = np.empty(shape = (pars['N-points']))

    for idx, samp in df.iterrows():

        samp_filt = {key: samp[key] for key in w.population_parameters}
        w.update(**samp_filt)
        pdf = w.pdf(q_array)
        curves[idx] = pdf

    plot_dict = {
        'x'       : q_array,
        'output'  : pars['output'],
        'figname' : 'SecondaryMassFunction',
        'color'   : '#890C0A',
        'x_label' : '$q$',
        'y_label' : '$p(q)$',
        'label'   : 'gaussian',
        'x_min'   : q_min,
        'x_max'   : q_max,
    }

    return curves, plot_dict


def RateEvolutionFunction(df, w, p_dct, pars, prior = False):

    if prior:
        tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N_samp_prior']) for key in w.population_parameters}
        df  = pd.DataFrame(tmp)

    z_array = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
    curves  = np.empty(shape = (len(df), pars['N-points']))

    for idx, samp in df.iterrows():

        samp_filt = {key: samp[key] for key in w.population_parameters}
        w.update(**samp_filt)        
        func = w.rate.log_evaluate(z_array)
        curves[idx] = func

    plot_dict = {
        'x'       : z_array,
        'output'  : pars['output'],
        'figname' : 'RateEvolutionFunction',
        'color'   : '#0A3689',
        'x_label' : '$z$',
        'y_label' : '$ln[\Psi(z)/R_0]$',
        'label'   : pars['model-rate'],
        'x_min'   : pars['bounds-z'][0],
        'x_max'   : pars['bounds-z'][1],
    }

    return curves, plot_dict


def RateEvolutionFunctionProb(df, w, pars):

    c = icarogw.wrappers.FlatLambdaCDM_wrap(zmax = 20.)

    z_array = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])
    curves  = np.empty(shape = (len(df), pars['N-points']))

    for idx, samp in df.iterrows():

        samp_filt = {key: samp[key] for key in w.population_parameters}
        w.update(**samp_filt)        
        func = np.exp(w.rate.log_evaluate(z_array))
        curves[idx] = func * samp.R0

        # Comoving volume and redshift (1/(1+z)*dV/dz)
        samp_filt = {key: samp[key] for key in c.population_parameters}
        c.update(**samp_filt)    
        func = c.cosmology.dVc_by_dzdOmega_at_z(z_array) * 4*np.pi / (1+z_array)
        curves[idx] *= func
        curves[idx] = np.log(curves[idx])
        curves[idx] -= 7

    plot_dict = {
        'x'       : z_array,
        'output'  : pars['output'],
        'figname' : 'RateEvolutionFunctionProb',
        'color'   : '#164B0C',
        'x_label' : '$z$',
        'y_label' : '$\propto ln[p(z)]$',
        'label'   : pars['model-rate'],
        'x_min'   : pars['bounds-z'][0],
        'x_max'   : pars['bounds-z'][1],
    }

    return curves, plot_dict


def TransitionFunction(df, p_dct, pars, prior = False):

    z_array = np.linspace(pars['bounds-z'][0], pars['bounds-z'][1], pars['N-points'])

    if prior:
        if   pars['transition'] == 'sigmoid': tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N_samp_prior']) for key in ['zt', 'delta_zt', 'mix_z0']}
        elif pars['transition'] == 'linear':  tmp = {key: bilby.prior.Uniform(p_dct[key]['kwargs']['minimum'], p_dct[key]['kwargs']['maximum']).sample(pars['N_samp_prior']) for key in ['mix_z0', 'mix_z1']}
        df  = pd.DataFrame(tmp)
        
    curves  = np.empty(shape = (len(df), pars['N-points']))

    for idx, samp in df.iterrows():
        if   pars['transition'] == 'sigmoid':         curves[idx] = icarogw.priors._mixed_sigmoid_function(z_array, samp['zt'], samp['delta_zt'], samp['mix_z0'])
        elif pars['transition'] == 'double-sigmoid':  curves[idx] = icarogw.priors._mixed_double_sigmoid_function(z_array, samp['zt'], samp['delta_zt'], samp['mix_z0'], samp['mix_z1'])
        elif pars['transition'] == 'linear':          curves[idx] = icarogw.priors._mixed_linear_function(z_array, samp['mix_z0'], samp['mix_z1'])
        elif pars['transition'] == 'linear-sinusoid': curves[idx] = icarogw.priors._mixed_linear_sinusoid_function(z_array, samp['mix_z0'], samp['mix_z1'], samp['amp'], samp['freq'])
    
    plot_dict = {
        'x'       : z_array,
        'output'  : pars['output'],
        'figname' : 'TransitionFunction',
        'color'   : 'k',
        'x_label' : '$z$',
        'y_label' : '$\\sigma(z)$',
        'label'   : pars['transition'],
        'x_min'   : pars['bounds-z'][0],
        'x_max'   : pars['bounds-z'][1],
    }
    
    return curves, plot_dict


def MassRedshift_Joint(df, pars):

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

    colors = sns.color_palette('RdBu_r', pars['N-z-slices'])

    for idx, samp in tqdm.tqdm(df.iterrows(), total = len(df)):

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

    plot_dict = {
        'x'         : mass_array,
        'output'    : pars['output'],
        'figname'   : 'MassRedshift_Joint_rescaled',
        'colors'    : colors,
        'z_grid'    : zY,
        'x_label'   : '$m_1\ [M_{\odot}]$',
        'y_label_a' : '$z$',
        'y_label_b' : '$p(m_1)$',
        'label'     : 'PL+G(z)',
        'x_min'     : pars['bounds-m1'][0],
        'x_max'     : pars['bounds-m1'][1],
    }

    return curves, plot_dict
# ------------------------------------------------------------- #

# ------------------------------------------------------------- #
class Plots:

    def __init__(self, pars, df, m1w, m2w, rw, priors, injections):

        self.pars   = pars
        self.df     = df
        self.m1w    = m1w
        self.m2w    = m2w
        self.rw     = rw
        self.priors = priors
        self.inj    = injections

    def PrimaryMass(self):

        if not 'Redshift' in self.pars['model-primary']:
            curves, plot_dict = PrimaryMassFunction(self.df, self.m1w, self.priors, self.pars)
            if self.pars['true-values'] == {}:
                plot_curves(curves, plot_dict, logscale = True)
            else:
                curve_true, _ = PrimaryMassFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.m1w, self.priors, self.pars)
                plot_curves(curves, plot_dict, truth = curve_true, logscale = True)
        else:
            curves, plot_dict = PrimaryMassFunction(self.df, self.m1w, self.priors, self.pars)
            if self.pars['true-values'] == {}:
                plot_curves_evolving(curves, plot_dict)
            else:
                curve_true, _ = PrimaryMassFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.m1w, self.priors, self.pars)
                if not self.pars['selection-effects']: plot_curves_evolving_long(curves, plot_dict, truth = curve_true)
                else:                                  plot_curves_evolving_long(curves, plot_dict, truth = curve_true, selection_effects = self.inj.injections_data)

    def SecondaryMass(self):

        curves, plot_dict = SecondaryMassFunction(self.df, self.m2w, self.priors, self.pars)
        if self.pars['true-values'] == {}:
            plot_curves(curves, plot_dict, figsize = (8,8))
        else:
            curve_true, _ = SecondaryMassFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.m2w, self.priors, self.pars)
            plot_curves(curves, plot_dict, figsize = (8,8), truth = curve_true[0])

    def RateEvolution(self):

        curves, plot_dict = RateEvolutionFunction(self.df, self.rw, self.priors, self.pars)
        if self.pars['true-values'] == {}:
            plot_curves(curves, plot_dict)
        else:
            curve_true, _ = RateEvolutionFunction(pd.DataFrame(self.pars['true-values'], index = [0]), self.rw, self.priors, self.pars)
            plot_curves(curves, plot_dict, truth = curve_true[0])

    # Produce the plots
    def ProducePlots(self):

        self.PrimaryMass()
        self.SecondaryMass()
        self.RateEvolution()
