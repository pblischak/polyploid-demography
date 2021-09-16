#!/usr/bin/env python3

import dadi
import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy
from dadi import Numerics, Inference

def plot_1d_comp_multinom(model1, model2, data, fig_num=None, residual='Anscombe',
                          plot_masked=False):
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

def plot_1d_comp_Poisson(model1, model2, data, fig_num=None, residual='Anscombe',
                         plot_masked=False, show=True):
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
    pts_l = [50,60,70]
    func1 = allotetraploid_bottleneck
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = segtetraploid_bottleneck
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt1 = [16.13550897, 0.1953863, 1.9726986, 0.507677479, 0.016506372]
    popt2 = [8.0429719, 1.22627327, 16.8074473, 2.5559702, 5.576213e-03, 0.030813141]

    model1 = func1_ex(popt1, data.sample_sizes, pts_l)
    model2 = func2_ex(popt2, data.sample_sizes, pts_l)

    plot_1d_comp_multinom(model1, model2, data, fig_num=1)
