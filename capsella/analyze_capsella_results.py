#!/usr/bin/env python3

import sys

import numpy as np
import pylab

import dadi
import dadi.Godambe
from dadi import Numerics, Inference

CB_color_cycle = [
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#f781bf",
    "#a65628",
    "#984ea3",
    "#999999",
    "#e41a1c",
    "#dede00",
]


def plot_1d_comp_Poisson(
    model1,
    model2,
    model3,
    data,
    residual="Anscombe",
    plot_masked=False,
    show=True,
):
    """
    Poisson comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.
    show: If True, execute pylab.show command to make sure plot displays.
    """
    pylab.clf()

    if data.folded and not model1.folded:
        model1 = model1.fold()

    if data.folded and not model2.folded:
        model2 = model2.fold()

    if data.folded and not model3.folded:
        model3 = model3.fold()

    masked_model1, masked_data = Numerics.intersect_masks(model1, data)
    masked_model2, masked_data = Numerics.intersect_masks(model2, data)
    masked_model3, masked_data = Numerics.intersect_masks(model3, data)

    ax = pylab.subplot(2, 1, 1)
    pylab.semilogy(masked_data, "-o", color=CB_color_cycle[0], label="data")
    pylab.semilogy(
        masked_model1, "-o", color=CB_color_cycle[1], label="allotet"
    )
    pylab.semilogy(masked_model2, "-o", color=CB_color_cycle[2], label="segtet")
    pylab.semilogy(
        masked_model3, "-o", color=CB_color_cycle[3], label="douglas"
    )

    if plot_masked:
        pylab.semilogy(masked_data.data, "--ob", mfc="w", zorder=-100)
        pylab.semilogy(masked_model1.data, "--or", mfc="w", zorder=-100)
        pylab.semilogy(masked_model2.data, "--og", mfc="w", zorder=-100)
        pylab.semilogy(masked_model3.data, "--og", mfc="w", zorder=-100)

    ax2 = pylab.subplot(2, 1, 2, sharex=ax)
    ax2.hlines([0], 0, 20, linestyles="dashed")
    if residual == "Anscombe":
        resid1 = Inference.Anscombe_Poisson_residual(masked_model1, masked_data)
        resid2 = Inference.Anscombe_Poisson_residual(masked_model2, masked_data)
        resid3 = Inference.Anscombe_Poisson_residual(masked_model3, masked_data)
    elif residual == "linear":
        resid1 = Inference.linear_Poisson_residual(masked_model1, masked_data)
        resid2 = Inference.linear_Poisson_residual(masked_model2, masked_data)
        resid3 = Inference.linear_Poisson_residual(masked_model3, masked_data)
    else:
        raise ValueError("Unknown class of residual '%s'." % residual)
    pylab.plot(resid1, "-oc", label="allotet")
    pylab.plot(resid2, "-om", label="segtet")
    pylab.plot(resid3, "-ob", label="douglas")
    if plot_masked:
        pylab.plot(resid1.data, "--oc", mfc="w", zorder=-100)
        pylab.plot(resid2.data, "--om", mfc="w", zorder=-100)
        pylab.plot(resid3.data, "--ob", mfc="w", zorder=-100)

    ax.set_xlim(0, data.shape[0] - 1)
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.legend(loc="upper right")
    ax.set_ylabel("Allele count")
    ax2.legend(loc="lower right")
    ax2.set_xlabel("Sample allele frequency")
    ax2.set_ylabel("Anscombe residual")
    if show:
        pylab.show()


# Defining demographic models


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
    nu0, nuBot, T1, T2, misid = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["CbpA", "CbpB"]
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D, misid)
    fs_1D = fs_2D.combine_pops([1, 2])
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
    nu0, nuBot, T1, T2, dij, misid = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot, m12=dij, m21=dij)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["CbpA", "CbpB"]
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D, misid)
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


def douglas_et_al(params, ns, pts):
    """
    Implementation of Model C from Douglas et al. for the A subgenome of
    Capsella busra-pastoris and the B subgenome of Capsella bursa-pastoris.
    Includes exponential growth and is based on C. grandiflora as the reference
    population (even though it isn't in the model directly).

    Parameters
    ----------
    params : Tuple[Union[int, float]]
        A tuple of values for the following parameters:
          - nu_Co0: Population size for ancestral Co size before Cbp formation
          - nu_CbpA0: Initial population size for C. bursa-pastoris subgenome A
          - nu_CbpAF: Final population size for C. bursa-pastoris subgenome A
          - nu_CbpB0: Initial population size for C. bursa-pastoris subgenome B
          - nu_CbpBF: Final population size for C. bursa-pastoris subgenome B
          - T1: Time period of constant pop size for parental pops (Cg and Co)
          - T2: Time period of size change and exp. growth after Cbp formation

    ns : Tuple[int]
        The sample sizes for the modeled populations.
    pts : List[int]
        Grid points for numerical integration.

    Returns
    -------
    dadi.Spectrum
        The site frequency spectrum as a ``dadi.Spectrum`` object.
    """
    (nu_Co0, nu_CbpA0, nu_CbpAF, nu_CbpB0, nu_CbpBF, T1, T2) = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T2, 1, nu_Co0)

    nu_CbpA_func = lambda t: nu_CbpA0 * (nu_CbpAF / nu_CbpA0) ** (t / T1)
    nu_CbpB_func = lambda t: nu_CbpB0 * (nu_CbpBF / nu_CbpB0) ** (t / T1)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu_CbpA_func, nu_CbpB_func)

    fs_2D = dadi.Spectrum.from_phi(
        phi,
        new_ns,
        (xx, xx),
        pop_ids=["CbpA", "CbpB"],
    )
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


if __name__ == "__main__":
    data = dadi.Spectrum.from_file(
        "Capsella_intergene_4fold_corr_4pop_4_DSFS.fs"
    )
    L = data.data.sum()  # sequence length
    mu = 7.0e-9  # mutation rate
    g = 1  # generation time
    scalar = 4.0 * L * mu

    data.pop_ids = ["CbpA", "CbpB", "grand", "ori"]
    data = data.marginalize([2, 3]).combine_pops([1, 2])
    pts_l = [50, 60, 70]

    allotet_params = [100.000000, 0.2023164, 2.064304, 0.5392525, 0.01660832]
    segtet_params = [
        100.00000000,
        7.166250,
        100.00000,
        18.3723521,
        0.0009271943,
        0.031397990,
    ]

    # Parameters reported in the supplement of Douglas et al.
    ppaper = {
        "NeCgF": 840e3,
        "NeCoF": 4e3,
        "NeCbpAF": 37e3,
        "NeCbpBF": 75e3,
        "T1": 128e3,
        "T2": 931e3 - 128e3,
    }

    # Back-calculated starting pop sizes based on growth rates
    ppaper["NeCg0"] = ppaper["NeCgF"] * (1 - 4.1e-6) ** ppaper["T1"]
    ppaper["NeCbpA0"] = ppaper["NeCbpAF"] * (1 + 4.8e-7) ** ppaper["T1"]
    ppaper["NeCo0"] = ppaper["NeCoF"] * (1 + 2.6e-5) ** ppaper["T1"]
    ppaper["NeCbpB0"] = ppaper["NeCbpBF"] * (1 + 4.8e-7) ** ppaper["T1"]

    douglas_params = [
        ppaper["NeCo0"] / ppaper["NeCg0"],
        ppaper["NeCbpA0"] / ppaper["NeCg0"],
        ppaper["NeCbpAF"] / ppaper["NeCg0"],
        ppaper["NeCbpB0"] / ppaper["NeCg0"],
        ppaper["NeCbpBF"] / ppaper["NeCg0"],
        ppaper["T1"] / (2 * ppaper["NeCg0"]),
        ppaper["T2"] / (2 * ppaper["NeCg0"]),
    ]

    print(data.sample_sizes, data.pop_ids, L)

    func1_ex = dadi.Numerics.make_extrap_log_func(allotetraploid_bottleneck)
    model1 = func1_ex(allotet_params, data.sample_sizes, pts_l)
    func2_ex = dadi.Numerics.make_extrap_log_func(segtetraploid_bottleneck)
    model2 = func2_ex(segtet_params, data.sample_sizes, pts_l)
    func3_ex = dadi.Numerics.make_extrap_log_func(douglas_et_al)
    model3 = func3_ex(douglas_params, data.sample_sizes, pts_l)

    llik1 = -651.0186
    theta1 = 5393.922
    llik2 = -512.4066
    theta2 = 167.3206
    theta3 = scalar * ppaper["NeCg0"]
    print(theta3)

    plot_1d_comp_Poisson(
        model1 * theta1, model2 * theta2, model3 * theta3, data
    )
    sys.exit(0)

    print(scalar)
    Nref1 = theta1 / scalar
    Nref2 = theta2 / scalar
    print(Nref1, Nref2, sep="\t")
    X_sq = 2 * (llik2 - llik1)
    print(f"Likelihood ratio statistic = {X_sq}")
    p_val = dadi.Godambe.sum_chi2_ppf(X_sq, weights=(0.5, 0.5))
    print(f"p_val = {p_val}")

    allotet_uncerts = dadi.Godambe.FIM_uncert(
        func1_ex, pts_l, allotet_params, data, log=True
    )
    allotet_FIM = dadi.Godambe.get_godambe(
        func1_ex,
        pts_l,
        [],
        allotet_params,
        data,
        0.01,
        True,
        just_hess=True,
    )
    allotet_vcov = np.linalg.inv(allotet_FIM)
    allotet_uncerts2 = [
        np.sqrt(
            allotet_vcov[-1, -1] + allotet_vcov[0, 0] + 2 * allotet_vcov[0, -1]
        ),
        np.sqrt(
            allotet_vcov[-1, -1] + allotet_vcov[1, 1] + 2 * allotet_vcov[1, -1]
        ),
        np.sqrt(
            allotet_vcov[-1, -1] + allotet_vcov[2, 2] + 2 * allotet_vcov[2, -1]
        ),
        np.sqrt(
            allotet_vcov[-1, -1] + allotet_vcov[3, 3] + 2 * allotet_vcov[3, -1]
        ),
        allotet_uncerts[4],
    ]

    allotet_log_params = [
        np.log([theta1, allotet_params[0], 1 / scalar]).sum(),  # nu0
        np.log([theta1, allotet_params[1], 1 / scalar]).sum(),  # nuBot
        np.log([theta1, allotet_params[2], 2 / scalar]).sum(),  # T1
        np.log([theta1, allotet_params[3], 2 / scalar]).sum(),  # T2
    ]

    print(
        "\nParameter estimates and 95% confidence intervals for allotet model:"
    )
    print(
        "Nref = {} ({}--{})".format(
            Nref1,
            np.exp(np.log(theta1) - 1.96 * allotet_uncerts2[-1]) / scalar,
            np.exp(np.log(theta1) + 1.96 * allotet_uncerts2[-1]) / scalar,
        )
    )
    print(
        "nu0   = {} ({}--{})".format(
            allotet_params[0] * Nref1,
            np.exp(allotet_log_params[0] - 1.96 * allotet_uncerts2[0]),
            np.exp(allotet_log_params[0] + 1.96 * allotet_uncerts2[0]),
        )
    )
    print(
        "nuBot   = {} ({}--{})".format(
            allotet_params[1] * Nref1,
            np.exp(allotet_log_params[1] - 1.96 * allotet_uncerts2[1]),
            np.exp(allotet_log_params[1] + 1.96 * allotet_uncerts2[1]),
        )
    )
    print(
        "T1   = {} ({}--{})".format(
            allotet_params[2] * 2 * g * Nref1,
            np.exp(allotet_log_params[2] - 1.96 * allotet_uncerts2[2]),
            np.exp(allotet_log_params[2] + 1.96 * allotet_uncerts2[2]),
        )
    )
    print(
        "T2   = {} ({}--{})".format(
            allotet_params[3] * 2 * g * Nref1,
            np.exp(allotet_log_params[3] - 1.96 * allotet_uncerts2[3]),
            np.exp(allotet_log_params[3] + 1.96 * allotet_uncerts2[3]),
        )
    )
    print(
        "% misid    = {} ({}--{})".format(
            allotet_params[4] * 100,
            np.exp(np.log(allotet_params[4]) - 1.96 * allotet_uncerts2[4]) * 100,
            np.exp(np.log(allotet_params[4]) + 1.96 * allotet_uncerts2[4]) * 100,
        )
    )

    segtet_uncerts = dadi.Godambe.FIM_uncert(
        func2_ex, pts_l, segtet_params, data, log=True
    )
    segtet_FIM = dadi.Godambe.get_godambe(
        func2_ex,
        pts_l,
        [],
        segtet_params,
        data,
        0.01,
        True,
        just_hess=True,
    )
    segtet_vcov = np.linalg.inv(segtet_FIM)
    segtet_uncerts2 = [
        np.sqrt(
            segtet_vcov[-1, -1] + segtet_vcov[0, 0] + 2 * segtet_vcov[0, -1]
        ),
        np.sqrt(
            segtet_vcov[-1, -1] + segtet_vcov[1, 1] + 2 * segtet_vcov[1, -1]
        ),
        np.sqrt(
            segtet_vcov[-1, -1] + segtet_vcov[2, 2] + 2 * segtet_vcov[2, -1]
        ),
        np.sqrt(
            segtet_vcov[-1, -1] + segtet_vcov[3, 3] + 2 * segtet_vcov[3, -1]
        ),
        np.sqrt(
            segtet_vcov[-1, -1] + segtet_vcov[4, 4] + 2 * segtet_vcov[4, -1]
        ),
        segtet_uncerts[5],
    ]

    segtet_log_params = [
        np.log([theta2, segtet_params[0], 1 / scalar]).sum(),  # nu0
        np.log([theta2, segtet_params[1], 1 / scalar]).sum(),  # nuBot
        np.log([theta2, segtet_params[2], 2 / scalar]).sum(),  # T1
        np.log([theta2, segtet_params[3], 2 / scalar]).sum(),  # T2
        np.log([1 / theta2, segtet_params[4], scalar / 2]).sum(),  # e_ij
    ]

    print(
        "\nParameter estimates and 95% confidence intervals for segtet model:"
    )
    print(
        "Nref = {} ({}--{})".format(
            Nref2,
            np.exp(np.log(theta2) - 1.96 * segtet_uncerts2[-1]) / scalar,
            np.exp(np.log(theta2) + 1.96 * segtet_uncerts2[-1]) / scalar,
        )
    )
    print(
        "nu0   = {} ({}--{})".format(
            segtet_params[0] * Nref2,
            np.exp(segtet_log_params[0] - 1.96 * segtet_uncerts2[0]),
            np.exp(segtet_log_params[0] + 1.96 * segtet_uncerts2[0]),
        )
    )
    print(
        "nuBot   = {} ({}--{})".format(
            segtet_params[1] * Nref2,
            np.exp(segtet_log_params[1] - 1.96 * segtet_uncerts2[1]),
            np.exp(segtet_log_params[1] + 1.96 * segtet_uncerts2[1]),
        )
    )
    print(
        "T1   = {} ({}--{})".format(
            segtet_params[2] * 2 * g * Nref2,
            np.exp(segtet_log_params[2] - 1.96 * segtet_uncerts2[2]),
            np.exp(segtet_log_params[2] + 1.96 * segtet_uncerts2[2]),
        )
    )
    print(
        "T2   = {} ({}--{})".format(
            segtet_params[3] * 2 * g * Nref2,
            np.exp(segtet_log_params[3] - 1.96 * segtet_uncerts2[3]),
            np.exp(segtet_log_params[3] + 1.96 * segtet_uncerts2[3]),
        )
    )
    print(
        "e_ij   = {} ({}--{})".format(
            segtet_params[4] / (2 * g * Nref2),
            np.exp(segtet_log_params[4] - 1.96 * segtet_uncerts2[4]),
            np.exp(segtet_log_params[4] + 1.96 * segtet_uncerts2[4]),
        )
    )
    print(
        "% misid    = {} ({}--{})".format(
            segtet_params[5] * 100,
            np.exp(np.log(segtet_params[5]) - 1.96 * segtet_uncerts2[5]) * 100,
            np.exp(np.log(segtet_params[5]) + 1.96 * segtet_uncerts2[5]) * 100,
        )
    )
