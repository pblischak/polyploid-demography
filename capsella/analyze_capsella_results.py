#!/usr/bin/env python3

import dadi
import dadi.Godambe
import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy as np
from dadi import Numerics, Inference
from sys import exit

CB_color_cycle = [
    '#377eb8', '#ff7f00', '#4daf4a',
    '#f781bf', '#a65628', '#984ea3',
    '#999999', '#e41a1c', '#dede00'
]


def plot_1d_comp_multinom(
    model1, model2, data, fig_num=None,
    residual='Anscombe', plot_masked=False
):
    """
    Mulitnomial comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.

    This comparison is multinomial in that it rescales the model to optimally
    fit the data.
    """
    model1 = Inference.optimally_scaled_sfs(model1, data)
    model2 = Inference.optimally_scaled_sfs(model2, data)

    plot_1d_comp_Poisson(model1, model2, data, fig_num, residual,
                         plot_masked)

def plot_1d_comp_Poisson(
    model1, model2, data, fig_num=None,
    residual='Anscombe', plot_masked=False, show=True
):
    """
    Poisson comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.
    show: If True, execute pylab.show command to make sure plot displays.
    """
    if fig_num is None:
        f = pylab.gcf()
    else:
        f = pylab.figure(fig_num, figsize=(5,4))
    pylab.clf()

    if data.folded and not model1.folded:
        model1 = model1.fold()

    if data.folded and not model2.folded:
        model2 = model2.fold()

    masked_model1, masked_data = Numerics.intersect_masks(model1, data)
    masked_model2, masked_data = Numerics.intersect_masks(model2, data)

    ax = pylab.subplot(2,1,1)
    pylab.semilogy(masked_data, '-o', color=CB_color_cycle[0])
    pylab.semilogy(masked_model1, '-o', color=CB_color_cycle[1])
    pylab.semilogy(masked_model2, '-o', color=CB_color_cycle[2])

    if plot_masked:
        pylab.semilogy(masked_data.data, '--ob', mfc='w', zorder=-100)
        pylab.semilogy(masked_model1.data, '--or', mfc='w', zorder=-100)
        pylab.semilogy(masked_model2.data, '--og', mfc='w', zorder=-100)

    pylab.subplot(2,1,2, sharex = ax)
    if residual == 'Anscombe':
        resid1 = Inference.Anscombe_Poisson_residual(masked_model1, masked_data)
        resid2 = Inference.Anscombe_Poisson_residual(masked_model2, masked_data)
    elif residual == 'linear':
        resid1 = Inference.linear_Poisson_residual(masked_model1, masked_data)
        resid2 = Inference.linear_Poisson_residual(masked_model2, masked_data)
    else:
        raise ValueError("Unknown class of residual '%s'." % residual)
    pylab.plot(resid1, '-oc')
    pylab.plot(resid2, '-om')
    #pylab.ylim(-160,120)
    if plot_masked:
        pylab.plot(resid1.data, '--oc', mfc='w', zorder=-100)
        pylab.plot(resid2.data, '--om', mfc='w', zorder=-100)

    ax.set_xlim(0, data.shape[0]-1)
    ax.set_xticks([0, 5, 10, 15, 20])
    if show:
        pylab.show()

# Defining demographic model
def allotetraploid_bottleneck(params, ns, pts):
    """
    params = (nu0,nuBot,T1,T2,misid)
    ns = (n1,n2)

    Split into two populations of specifed size.

    nu0: Size of populations after split.
    nuBot: Size of populations after bottleneck.
    T1: Length of time between split and bottleneck.
    T2: Length of bottleneck.
    n1,n2: Sample sizes of resulting Spectrum
    misid: probability of ancestral allele misidentification.
    pts: Number of grid points to use in integration.
    """
    nu0,nuBot,T1,T2,misid = params
    new_ns = [int(ns[0]/2),int(ns[0]/2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx,xx), pop_ids=['CbpA','CbpB']
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D,misid)
    fs_1D = fs_2D.combine_pops([1,2])
    return fs_1D

def segtetraploid_bottleneck(params, ns, pts):
    """
    params = (nu0,nuBot,T1,T2,dij,misid)
    ns = (n1,n2)

    Split into two populations of specifed size.

    nu0: Size of populations after split.
    nuBot: Size of populations after bottleneck
    T1: Time in the past between split and polyploid formation
        (in units of 2*Na generations)
    T2: Time since polyploid formation and bottleneck
    n1,n2: Sample sizes of resulting Spectrum
    dij: Effective homoeologous exchange rate (2Nm)
    misid: probability of ancestral allele misidentification
    pts: Number of grid points to use in integration
    """
    nu0,nuBot,T1,T2,dij,misid = params
    new_ns = [int(ns[0]/2),int(ns[0]/2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(
        phi, xx, T2, nuBot, nuBot, m12=dij, m21=dij
    )

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx,xx), pop_ids=['CbpA','CbpB']
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D,misid)
    fs_1D = fs_2D.combine_pops([1,2])
    return fs_1D

if __name__ == "__main__":
    data = dadi.Spectrum.from_file(
        "Capsella_intergene_4fold_corr_4pop_4_DSFS.fs"
    )
    L = np.array(data).sum()
    data.pop_ids = ["CbpA","CbpB","grand","ori"]
    data = data.marginalize([2,3]).combine_pops([1,2])
    all_boot = [data.sample() for i in range(100)]
    pts_l = [50,60,70]
    func1 = allotetraploid_bottleneck
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = segtetraploid_bottleneck
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt1 = [
        100.000000, 0.2023164, 2.064304,
        0.5392525, 0.01660832
    ]
    popt2 = [
        100.00000000, 7.166250, 100.00000,
        18.3723521, 0.0009271943, 0.031397990
    ]

    model1 = func1_ex(popt1, data.sample_sizes, pts_l)
    model2 = func2_ex(popt2, data.sample_sizes, pts_l)

    plot_1d_comp_multinom(model1, model2, data, fig_num=1)

    llik1  = -651.0186
    theta1 = 5393.922
    llik2  = -512.4066
    theta2 = 167.3206

    mu = 7.0e-9 # mutation rate
    print(L)
    # L  = 773748 # sequence length
    g = 1       # generation time
    scaler = 4.0 * L * mu
    print(scaler)
    Nref1 = theta1 / scaler
    Nref2 = theta2 / scaler
    print(Nref1, Nref2, sep='\t')
    X_sq = 2*(llik2 - llik1)
    print(f"Likelihood ratio statistic = {X_sq}")
    p_val = dadi.Godambe.sum_chi2_ppf(X_sq, weights=(0.5,0.5))
    print(f"p_val = {p_val}")

    print("Converted params for the models:")
    print(f"$N_A$     & {Nref1}            & {Nref2}")
    print(f"$N_0$     & {Nref1 * popt1[0]} & {Nref2 * popt2[0]}")
    print(f"$N_bot$   & {Nref1 * popt1[1]} & {Nref2 * popt2[1]}")
    print(f"$T_1$     & {2 * Nref1 * popt1[2]} & {2 * Nref2 * popt2[2]}")
    print(f"$T_2$     & {2 * Nref1 * popt1[3]} & {2 * Nref2 * popt2[3]}")
    print(f"$e_ij$    & NA                     & {popt2[4] / (2 * Nref2)}")
    print(f"\\% MisID & {popt1[4]}             & {popt2[5]}")

    # Look at loglikelihood surfaces
    # nu0_grid = [1,10,25,50,100,150,200,250,500]
    # allotet_llik_vals = [0 for i in nu0_grid]
    # segtet_llik_vals = [0 for i in nu0_grid]
    # popt1_tmp = popt1
    # popt2_tmp = popt2
    # for i,nu in enumerate(nu0_grid):
    #     print(f"Testing nu0 value {nu}...")
    #     popt1_tmp[0] = nu
    #     popt2_tmp[0] = nu
    #     allotet_llik_vals[i] = (
    #         dadi.Inference.ll_multinom(
    #             func1_ex(popt1_tmp, data.sample_sizes, pts_l),
    #             data
    #         )
    #     )
    #     segtet_llik_vals[i] = (
    #         dadi.Inference.ll_multinom(
    #             func2_ex(popt2_tmp, data.sample_sizes, pts_l),
    #             data
    #         )
    #     )
    # plt.plot(nu0_grid, allotet_llik_vals, '-o')
    # plt.show()
    # plt.plot(nu0_grid, segtet_llik_vals, '-o')
    # plt.show()

