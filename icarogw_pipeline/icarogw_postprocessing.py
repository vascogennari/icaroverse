import os, icarogw, bilby
import json, pandas as pd, numpy as np
import matplotlib.pyplot as plt
import tqdm, seaborn as sns

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

# ---------------------------------------------------------------- #
# -------------------- Wrappers initialization ------------------- #
# ---------------------------------------------------------------- #
def initialize_PrimaryMass(df, func, pars):

    if   func == 'PowerLaw':     w = icarogw.wrappers.massprior_PowerLaw()
    elif func == 'PowerLawPeak': w = icarogw.wrappers.massprior_PowerLawPeak()
    else:
        raise ValueError('Unknown function for the primary mass wrapper.')

    if pars['smoothing']:
        w = icarogw.wrappers.lowSmoothedwrapper(w)
        
    if pars['evolving'] and pars['transition'] == '':                w = icarogw.wrappers.mixed_mass_redshift_evolving(w)
    if pars['evolving'] and pars['transition'] == 'sigmoid':         w = icarogw.wrappers.mixed_mass_redshift_evolving_sigmoid(w)
    if pars['evolving'] and pars['transition'] == 'double-sigmoid':  w = icarogw.wrappers.double_mixed_mass_redshift_evolving_sigmoid(w)
    if pars['evolving'] and pars['transition'] == 'linear':          w = icarogw.wrappers.double_mixed_mass_redshift_evolving_linear(w)
    if pars['evolving'] and pars['transition'] == 'linear-sinusoid': w = icarogw.wrappers.double_mixed_mass_redshift_evolving_linear_sinusoid(w)

    # FIXME: Change the structure
    if pars['evolving'] and func == 'PowerLawZPeakZ':                w = icarogw.wrappers.PowerLawLinear_GaussianLinear_TransitionLinear()

    if pars['positive']:
        w = icarogw.wrappers.massprior_PowerLawPeakPositive(w)
    if not (set(w.population_parameters) <= set(df.keys())):
        raise ValueError('Cannot find the parameters for the selected primary mass function. Please make sure that you are using the correct one.')

    return w

def initialize_SecondaryMass(df, func):

    if func == 'MassRatio': w = icarogw.wrappers.mass_ratio_prior_Gaussian()
    else:
        raise ValueError('Unknown function for the secondary mass wrapper.')
    if not (set(w.population_parameters) <= set(df.keys())):
        raise ValueError('Cannot find the parameters for the selected secondary mass function. Please make sure that you are using the correct one.')

    return w

def initialize_Rate(df, func):

    if   func == 'MadauDickinson':      w = icarogw.wrappers.rateevolution_Madau()
    elif func == 'beta':                w = icarogw.wrappers.rateevolution_beta()
    elif func == 'betaline':            w = icarogw.wrappers.rateevolution_beta_line()
    elif func == 'MadauDickinsonGamma': w = icarogw.wrappers.rateevolution_Madau_gamma()
    elif func == 'PowerLaw':            w = icarogw.wrappers.rateevolution_PowerLaw()
    else:
        raise ValueError('Unknown function for the rate wrapper.')
    if not (set(w.population_parameters) <= set(df.keys())):
        raise ValueError('Cannot find the parameters for the selected rate. Please make sure that you are using the correct one.')

    return w
# ---------------------------------------------------------------- #
# ---------------------------------------------------------------- #

# ------------------------------------------------------------- #
# ----------------------- Generate plots ---------------------- #
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


def plot_curves_evolving_long(curves, pl_dct, truth = {}, curves_prior = {}):

    _, ax = plt.subplots(figsize = (5, 9))

    if bool(curves_prior):
        for zi, z_array in enumerate(pl_dct['z_grid']):
            z = z_array[0]
            ax.fill_between(pl_dct['x'], curves_prior[zi][5] +z, curves_prior[zi][95]+z, color = '#AB7C41', alpha = 0.05)
            ax.fill_between(pl_dct['x'], curves_prior[zi][16]+z, curves_prior[zi][84]+z, color = '#AB7C41', alpha = 0.15)
 
    for zi, z_array in enumerate(pl_dct['z_grid']):
        z = z_array[0]
        ax.fill_between(pl_dct['x'], curves[zi][5] +z, curves[zi][95]+z, color = pl_dct['colors'][zi], alpha = 0.25)
        ax.fill_between(pl_dct['x'], curves[zi][16]+z, curves[zi][84]+z, color = pl_dct['colors'][zi], alpha = 0.5)
        ax.plot(        pl_dct['x'], curves[zi][50]+z, lw = 0.7,         color = pl_dct['colors'][zi])
    
        if not truth == {}:
            ax.plot(    pl_dct['x'], truth[zi][50]+z,  lw = 0.3,         color = '#494949')

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

# ------------------------------------------------------------------ #
# ----------------------- Return level curves ---------------------- #
# ------------------------------------------------------------------ #
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

        colors = sns.color_palette('RdBu_r', pars['N-z-slices'])   # 'blend:#E5E3C8,#914E63'
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
# ---------------------------------------------------------------- #
# ---------------------------------------------------------------- #


# ---------------------------------------------------------------- #
# ------------------------- Main function ------------------------ #
# ---------------------------------------------------------------- #
def main():

    # INPUT PARAMETERS
    input_pars = {
        'model'       : 'PowerLaw',          # Options: [PowerLaw, PowerLawPeak, PowerLawZPeakZ]
        'rate'        : 'MadauDickinson',    # Options: [MadauDickinson, beta, betaline, MadauDickinsonGamma, PowerLaw]
        'transition'  : 'linear',            # Options: [sigmoid, double-sigmoid, linear, linear-sinusoid]

        'Primary'     : 1,
        'Secondary'   : 1,
        'Rate'        : 1,
        'RateProb'    : 0,
        'Transition'  : 1,
        '2D_joint'    : 0,

        'evolving'    : 1,              # Evolution already implies an evolving gaussian peak: you don't need to use PowerLawPeak as mass model
        'smoothing'   : 1,
        'positive'    : 0,
        'logscale'    : 1,
        'priors'      : 0,
        'cosmology'   : 0,

        'N-points'    : 1000,
        'N-z-slices'  : 10,
        'N_samp_prior': 1000,
        'downsample'  : -1,

        'bounds-m1'   : [1, 100],
        'bounds-z'    : [1e-5, 0.8],

        'filename'    : 'PLGz-m1_G-q_EXP16_md_double-mixed-linear_COND-PRIOR',
        'root-path'   : '/Users/vgennari/Documents/work/code/python/icarogw/results/evolution_study/sigmoid',
    }

    # Set output directory
    if input_pars['evolving']:
        if not input_pars['cosmology']: input_pars['outdir-path'] = os.path.join(input_pars['root-path'], 'evolving', 'population', input_pars['filename'])
        else:                           input_pars['outdir-path'] = os.path.join(input_pars['root-path'], 'evolving', 'cosmology',  input_pars['filename'])
        input_pars['output'] = os.path.join(input_pars['outdir-path'], 'postprocessing')
        input_pars['outdir-path'] = os.path.join(input_pars['outdir-path'], 'GWTC-3_FAR-025_PLGz-m1_G-q')
    else:
        if not input_pars['cosmology']: input_pars['outdir-path'] = os.path.join(input_pars['root-path'], 'non-evolving', 'population', input_pars['filename'])
        else:                           input_pars['outdir-path'] = os.path.join(input_pars['root-path'], 'non-evolving', 'cosmology',  input_pars['filename'])
        input_pars['output'] = os.path.join(input_pars['outdir-path'], 'postprocessing')
        input_pars['outdir-path'] = os.path.join(input_pars['outdir-path'], 'GWTC-3_FAR-025_PLG-m1_G-q')
    if not os.path.exists(input_pars['output']):
        os.makedirs(input_pars['output'])

    # Load samples as Data Frame
    samp_path = os.path.join(input_pars['outdir-path'], 'label_result.json')
    with open(samp_path) as f:
        tmp = json.load(f)
        df  = pd.DataFrame(tmp['posterior']['content'])
        priors_dict = tmp['priors']

    # Read and save the evidence
    with open('{}/log_evidence.txt'.format(input_pars['output']), 'w') as f:
        f.write('{}\n'.format('# log_Z_base_e\tlog_Z_err\tmax_log_L'))
        f.write('{}\t{}\t\t{}'.format(round(tmp['log_evidence'], 2), round(tmp['log_evidence_err'], 2), round(max(df['log_likelihood']), 2)))

    # Downsample the posteriors
    if not input_pars['downsample'] == -1:
        print('Total number of samples: {}'.format(len(df)))
        df = df.iloc[::input_pars['downsample']]
        df = df.reset_index()
        print('Reduced number of samples: {}'.format(len(df)))

    if input_pars['Primary']:
        print('\nPlotting primary mass function.')
        if not input_pars['evolving']:
            curves_prior = np.zeros(5)
            if input_pars['priors']:
                curves_prior, _ = PrimaryMassFunction(  df, priors_dict, input_pars, prior = True)
            curves, plot_dict   = PrimaryMassFunction(  df, priors_dict, input_pars)
            plot_curves(curves, plot_dict, curves_prior = curves_prior, logscale = True)
        else:
            curves_prior = {}
            if input_pars['priors']:
                curves_prior, _ = PrimaryMassFunction(  df, priors_dict, input_pars, prior = True)
            curves, plot_dict   = PrimaryMassFunction(  df, priors_dict, input_pars)
            plot_curves_evolving(curves, plot_dict, curves_prior = curves_prior)

    if input_pars['Secondary']:
        print('\nPlotting secondary mass function.')
        curves_prior = np.zeros(5)
        if input_pars['priors']:
            curves_prior, _ = SecondaryMassFunction(df, priors_dict, input_pars, prior = True)
        curves, plot_dict   = SecondaryMassFunction(df, priors_dict, input_pars)
        plot_curves(curves, plot_dict, figsize = (8,8), curves_prior = curves_prior)

    if input_pars['Transition']:
        print('\nPlotting transition function.')
        curves_prior = np.zeros(5)
        if input_pars['priors']:
            curves_prior, _ = TransitionFunction(df, priors_dict, input_pars, prior = True)
        curves, plot_dict   = TransitionFunction(df, priors_dict, input_pars)
        plot_curves(curves, plot_dict, curves_prior = curves_prior)

    if input_pars['Rate']:
        print('\nPlotting rate evolution function.')
        curves_prior = np.zeros(5)
        if input_pars['priors']:
            curves_prior, _ = RateEvolutionFunction(df, priors_dict, input_pars, prior = True)
        curves, plot_dict   = RateEvolutionFunction(df, priors_dict, input_pars)
        plot_curves(curves, plot_dict, curves_prior = curves_prior)

    if input_pars['RateProb']:
        curves_prior = np.zeros(5)
        print('\nPlotting rate evolution probability distribution.')
        curves, plot_dict = RateEvolutionFunctionProb(df, input_pars)
        plot_curves(curves, plot_dict, curves_prior = curves_prior)

    if input_pars['2D_joint']:
        print('\nPlotting 2D joint Mass-Redshift distribution.')
        curves, plot_dict = MassRedshift_Joint(df, input_pars)
        plot_curves_evolving_MassRedshift_Joint(curves, plot_dict)
        
    print('\nFinished.\n')

if __name__=='__main__':
    main()
# ---------------------------------------------------------------- #
# ---------------------------------------------------------------- #