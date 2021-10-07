#!/usr/bin/env python3

import dadi
import dadi.Godambe
import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy as np
from dadi import Numerics, Inference


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
        f = pylab.figure(fig_num, figsize=(10,8))
    pylab.clf()

    if data.folded and not model1.folded:
        model1 = model1.fold()

    if data.folded and not model2.folded:
        model2 = model2.fold()

    masked_model1, masked_data = Numerics.intersect_masks(model1, data)
    masked_model2, masked_data = Numerics.intersect_masks(model2, data)

    ax = pylab.subplot(2,1,1)
    pylab.semilogy(masked_data, '-ob')
    pylab.semilogy(masked_model1, '-or')
    pylab.semilogy(masked_model2, '-og')

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
    data.pop_ids = ["CbpA","CbpB","grand","ori"]
    data = data.marginalize([2,3]).combine_pops([1,2])
    all_boot = [data.sample() for i in range(100)]
    pts_l = [50,60,70]
    func1 = allotetraploid_bottleneck
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = segtetraploid_bottleneck
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt1 = [
        16.13550897,   0.1953863,  1.9726986,
        0.507677479, 0.016506372
    ]
    popt2 = [
        8.0429719, 1.22627327, 16.8074473,
        2.5559702, 5.576213e-03, 0.030813141
    ]

    model1 = func1_ex(popt1, data.sample_sizes, pts_l)
    model2 = func2_ex(popt2, data.sample_sizes, pts_l)

    #plot_1d_comp_multinom(model1, model2, data, fig_num=1)

    llik1  = -651.3589
    theta1 = 5584.862
    llik2  = -515.0822
    theta2 = 978.8552

    mu = 7.0e-9 # mutation rate
    L  = 773748 # sequence length
    g = 1       # generation time
    scalar = 4.0 * L * mu
    Nref1 = theta1 / scalar
    Nref2 = theta2 / scalar
    X_sq = 2*(llik2 - llik1)
    print(f"Likelihood ratio statistic = {X_sq}")
    p_val = dadi.Godambe.sum_chi2_ppf(X_sq, weights=(0.5,0.5))
    print(f"p_val = {p_val}")

    """
    Uncertainties for allotetraploid bottleneck (model 1)
    """
    check_grid_sizes1 = True
    if check_grid_sizes1:
        print("\nChecking uncertainties with different grid sizes:")
        for e in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    	    u = dadi.Godambe.GIM_uncert(
                func1_ex, pts_l, all_boot, popt1,
                data, log=True, multinom=True, eps=e
            )
    	    print(f"{e} = {u}")
    
    eps=1e-3 # the default is 1e-2
    uncerts1, GIM1 = dadi.Godambe.GIM_uncert(
        func1_ex, pts_l, all_boot, popt1, data, log=True,
        multinom=True, return_GIM=True, eps=eps
    )
    vcov1 = np.linalg.inv(GIM1)

    # Do conversions here for propogation of uncertainty
    # (ie, we're multiplying everything by theta, which is also estimated).
    uncerts1_2 = [
        np.sqrt(vcov1[-1,-1] + vcov1[0,0] + 2*vcov1[0,-1]),
        np.sqrt(vcov1[-1,-1] + vcov1[1,1] + 2*vcov1[1,-1]),
        np.sqrt(vcov1[-1,-1] + vcov1[2,2] + 2*vcov1[2,-1]),
        np.sqrt(vcov1[-1,-1] + vcov1[3,3] + 2*vcov1[3,-1]),
        uncerts1[4],
        uncerts1[-1]
    ]
    
    log_params1 = [
        np.log(theta1) + np.log(popt1[0]) + np.log(1/scalar),   # N0
        np.log(theta1) + np.log(popt1[1]) + np.log(1/scalar),   # Nbot
        np.log(theta1) + np.log(popt1[2]) + np.log(2*g/scalar), # T1
        np.log(theta1) + np.log(popt1[3]) + np.log(2*g/scalar), # T2
    ]
    
    print(
        f"Parameter standard deviations from GIM:\n{uncerts1}\n"
    )
    print(
        f"Parameter standard deviations from error propagation:\n{uncerts1_2}\n"
    )
    print("Variance-Covariance Matrix:\n{}\n".format(vcov1))

    """
    With propogation of uncertainty
    """
    print("\nParameter estimates and 95% confidence intervals:")
    print(
        "Nref  = {} ({}--{})".format(
            Nref1, np.exp(np.log(theta1)-1.96*uncerts1_2[-1])/scalar,
            np.exp(np.log(theta1)+1.96*uncerts1_2[-1])/scalar
        )
    )
    print(
        "N0    = {} ({}--{})".format(
            popt1[0]*Nref1, np.exp(log_params1[0]-1.96*uncerts1_2[0]),
            np.exp(log_params1[0]+1.96*uncerts1_2[0])
        )
    )
    print(
        "Nbot  = {} ({}--{})".format(
            popt1[1]*Nref1, np.exp(log_params1[1]-1.96*uncerts1_2[1]),
            np.exp(log_params1[1]+1.96*uncerts1_2[1])
        )
    )
    print(
        "T1    = {} ({}--{})".format(
            popt1[2]*2*g*Nref1, np.exp(log_params1[2]-1.96*uncerts1_2[2]),
            np.exp(log_params1[2]+1.96*uncerts1_2[2])
        )
    )
    print(
        "T2    = {} ({}--{})".format(
            popt1[3]*2*g*Nref1, np.exp(log_params1[3]-1.96*uncerts1_2[3]),
            np.exp(log_params1[3]+1.96*uncerts1_2[3])
        )
    )
    print(
        "misid = {} ({}--{})".format(
            popt1[4], np.exp(np.log(popt1[4])-1.96*uncerts1_2[4]),
            np.exp(np.log(popt1[4])+1.96*uncerts1_2[4])
        )
    )

    """
    Uncertainties for segtetraploid bottleneck (model 2)
    """

    check_grid_sizes2 = True
    if check_grid_sizes2:
        print("\nChecking uncertainties with different grid sizes:")
        for e in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    	    u = dadi.Godambe.GIM_uncert(
                func2_ex, pts_l, all_boot, popt2,
                data, log=True, multinom=True, eps=e
            )
    	    print(f"{e} = {u}")
    
    eps=1e-3 # the default is 1e-2
    uncerts2, GIM2 = dadi.Godambe.GIM_uncert(
        func2_ex, pts_l, all_boot, popt2, data, log=True,
        multinom=True, return_GIM=True, eps=eps
    )
    vcov2 = np.linalg.inv(GIM2)

    # Do conversions here for propogation of uncertainty
    # (ie, we're multiplying everything by theta, which is also estimated).
    uncerts2_2 = [
        np.sqrt(vcov2[-1,-1] + vcov1[0,0] + 2*vcov1[0,-1]),
        np.sqrt(vcov2[-1,-1] + vcov1[1,1] + 2*vcov1[1,-1]),
        np.sqrt(vcov2[-1,-1] + vcov1[2,2] + 2*vcov1[2,-1]),
        np.sqrt(vcov2[-1,-1] + vcov1[3,3] + 2*vcov1[3,-1]),
        np.sqrt(vcov2[-1,-1] + vcov1[4,4] + 2*vcov1[4,-1]),
        uncerts2[5],
        uncerts2[-1]
    ]
    
    log_params2 = [
        np.log(theta2) + np.log(popt2[0]) + np.log(1/scalar),   # N0
        np.log(theta2) + np.log(popt2[1]) + np.log(1/scalar),   # Nbot
        np.log(theta2) + np.log(popt2[2]) + np.log(2*g/scalar), # T1
        np.log(theta2) + np.log(popt2[3]) + np.log(2*g/scalar), # T2
        np.log(theta2) + np.log(popt2[4]) + np.log(scalar/2),   # M
    ]
    
    print(
        f"Parameter standard deviations from GIM:\n{uncerts2}\n"
    )
    print(
        f"Parameter standard deviations from error propagation:\n{uncerts2_2}\n"
    )
    print("Variance-Covariance Matrix:\n{}\n".format(vcov2))

    """
    With propogation of uncertainty
    """
    print("\nParameter estimates and 95% confidence intervals:")
    print(
        "Nref  = {} ({}--{})".format(
            Nref2, np.exp(np.log(theta2)-1.96*uncerts2_2[-1])/scalar,
            np.exp(np.log(theta2)+1.96*uncerts2_2[-1])/scalar
        )
    )
    print(
        "N0    = {} ({}--{})".format(
            popt2[0]*Nref2, np.exp(log_params2[0]-1.96*uncerts2_2[0]),
            np.exp(log_params2[0]+1.96*uncerts2_2[0])
        )
    )
    print(
        "Nbot  = {} ({}--{})".format(
            popt2[1]*Nref2, np.exp(log_params2[1]-1.96*uncerts2_2[1]),
            np.exp(log_params2[1]+1.96*uncerts2_2[1])
        )
    )
    print(
        "T1    = {} ({}--{})".format(
            popt2[2]*2*g*Nref2, np.exp(log_params2[2]-1.96*uncerts2_2[2]),
            np.exp(log_params2[2]+1.96*uncerts2_2[2])
        )
    )
    print(
        "T2    = {} ({}--{})".format(
            popt2[3]*2*g*Nref2, np.exp(log_params2[3]-1.96*uncerts2_2[3]),
            np.exp(log_params2[3]+1.96*uncerts2_2[3])
        )
    )
    print(
        "dij   = {} ({}--{})".format(
            popt2[4]/(2*Nref2), np.exp(log_params2[4]-1.96*uncerts2_2[4]),
            np.exp(log_params2[4]+1.96*uncerts2_2[4])
        )
    )
    print(
        "misid = {} ({}--{})".format(
            popt2[5], np.exp(np.log(popt2[5])-1.96*uncerts2_2[5]),
            np.exp(np.log(popt2[5])+1.96*uncerts2_2[5])
        )
    )