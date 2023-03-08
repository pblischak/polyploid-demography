import matplotlib.pyplot as plt
import dadi


"""
DEMOGRAPHIC MODELS

Defining demographic models from our paper and from Douglas et al. to check that
we can reproduce their results and then make the appropriate comparisons.
"""


@dadi.Numerics.make_extrap_func
def Cg_CbpA(params, ns, pts):
    """
    Implementation of Model C from Douglas et al. for Capsella grandiflora
    and the A subgenome of Capsella bursa-pastoris. Includes exponential
    growth and C. grandiflora as the reference population.

    Parameters
    ----------
    params : Tuple[Union[int, float]]
        A tuple of values for the following parameters:
          - nu_CgF: Final population size for C. grandiflora
          - nu_CbpA0: Initial population size for C. bursa-pastoris A
          - nu_CbpAF: Final population size for C. bursa-pastoris A
          - T: Divergence time between Cg and CbpA

    ns : Tuple[int]
        The sample sizes for the modeled populations.
    pts : List[int]
        Grid points for numerical integration.

    Returns
    -------
    dadi.Spectrum
        The site frequency spectrum as a ``dadi.Spectrum`` object.
    """
    (nu_CgF, nu_CbpA0, nu_CbpAF, T) = params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu_Cg_func = lambda t: nu_CgF ** (t / T)
    nu_CbpA_func = lambda t: nu_CbpA0 * (nu_CbpAF / nu_CbpA0) ** (t / T)

    phi = dadi.Integration.two_pops(phi, xx, T, nu_Cg_func, nu_CbpA_func)

    return dadi.Spectrum.from_phi(
        phi, ns, (xx, xx), pop_ids=["C. grandiflora", "C. bursa-pastoris A"]
    )


@dadi.Numerics.make_extrap_func
def Cg_Co(params, ns, pts):
    """
    Implementation of Model C from Douglas et al. for Capsella grandiflora
    and Capsella orientalis. Includes exponential growth and C. grandiflora as
    the reference population.

    Parameters
    ----------
    params : Tuple[Union[int, float]]
        A tuple of values for the following parameters:
          - nu_CgF: Final population size for C. grandiflora
          - nu_Co0: Initial population size for C. orientalis
          - nu_CoF: Final population size for C. orientalis
          - T1: Time period of recent exponential growth after split
          - T2: Time period of initial constant pop size after split

    ns : Tuple[int]
        The sample sizes for the modeled populations.
    pts : List[int]
        Grid points for numerical integration.

    Returns
    -------
    dadi.Spectrum
        The site frequency spectrum as a ``dadi.Spectrum`` object.
    """
    (nu_CgF, nu_Co0, nu_CoF, T1, T2) = params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T2, 1, nu_Co0)

    nu_Cg_func = lambda t: nu_CgF ** (t / T1)
    nu_Co_func = lambda t: nu_Co0 * (nu_CoF / nu_Co0) ** (t / T1)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu_Cg_func, nu_Co_func)

    return dadi.Spectrum.from_phi(
        phi, ns, (xx, xx), pop_ids=["C. grandiflora", "C. orientalis"]
    )


@dadi.Numerics.make_extrap_func
def Co_CbpB(params, ns, pts):
    """
    Implementation of Model C from Douglas et al. for Capsella orientalis
    and the B subgenome of Capsella bursa-pastoris. Includes exponential growth
    and is based on C. grandiflora as the reference population (even though it
    isn't in the model directly).

    Parameters
    ----------
    params : Tuple[Union[int, float]]
        A tuple of values for the following parameters:
          - nu_Co0: Initial population size for C. orientalis
          - nu_CoF: Final population size for C. orientalis
          - nu_CbpB0: Initial population size for C. bursa-pastoris subgenome B
          - nu_CbpBF: Final population size for C. bursa-pastoris subgenome B
          - T1: Time period of exp. growth after split between Co and CbpB
          - T2: Time period of constant pop size for Co

    ns : Tuple[int]
        The sample sizes for the modeled populations.
    pts : List[int]
        Grid points for numerical integration.

    Returns
    -------
    dadi.Spectrum
        The site frequency spectrum as a ``dadi.Spectrum`` object.
    """
    (nu_Co0, nu_CoF, nu_CbpB0, nu_CbpBF, T1, T2) = params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T2, nu_Co0)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu_Co_func = lambda t: nu_Co0 * (nu_CoF / nu_Co0) ** (t / T1)
    nu_CbpB_func = lambda t: nu_CbpB0 * (nu_CbpBF / nu_CbpB0) ** (t / T1)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu_Co_func, nu_CbpB_func)

    return dadi.Spectrum.from_phi(
        phi, ns, (xx, xx), pop_ids=["C. orientalis", "C. bursa-pastoris B"]
    )


@dadi.Numerics.make_extrap_func
def CbpA_CbpB(params, ns, pts):
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
          - T1: Time period of size change and exp. growth after Cbp formation
          - T2: Time period of constant pop size for parental pops (Cg and Co)

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
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T2, 1, nu_Co0)

    nu_CbpA_func = lambda t: nu_CbpA0 * (nu_CbpAF / nu_CbpA0) ** (t / T1)
    nu_CbpB_func = lambda t: nu_CbpB0 * (nu_CbpBF / nu_CbpB0) ** (t / T1)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu_CbpA_func, nu_CbpB_func)

    return dadi.Spectrum.from_phi(
        phi,
        ns,
        (xx, xx),
        pop_ids=["C. bursa-pastoris A", "C. bursa-pastoris B"],
    )


@dadi.Numerics.make_extrap_func
def segtetraploid_bottleneck(params, ns, pts):
    """
    Parameters
    ----------
    params : Tuple[Union[int, float]]
        A tuple of values for the flollowing parameters:
          - nu0: population size after split
          - nuBot: population size after botttleneck
          - T1: Timer period between split and polyploid formation
          - T2: Time period since polyploid formation and bottleneck
          - eij: Effective homoeologous exchange rate
          - misid: probability of ancestral allele misidentification
    
    ns : Tuple[int]
        The sample sizes for the modeled populations.
    pts : List[int]
        Grid points for numerical integration.

    Returns
    -------
    dadi.Spectrum
        The site frequency spectrum as a ``dadi.Spectrum`` object.
    """
    nu0, nuBot, T1, T2, eij, misid = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot, m12=eij, m21=eij)

    fs_2D = dadi.Spectrum.from_phi(
        phi,
        new_ns,
        (xx, xx),
        pop_ids=["C. bursa-pastoris A", "C. bursa-pastoris B"],
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D, misid)
    return fs_2D


def compare_results():
    """
    """
    ns = (26, 10)  # Sample sizes
    pts_l = [30, 35, 40]  # Grid points

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

    douglas_data_orig = dadi.Spectrum.from_file(
        "Capsella_intergene_4fold_corr_4pop_4_DSFS.fs"
    )
    douglas_data_orig.pop_ids = [
        "C. bursa-pastoris A",
        "C. bursa-pastoris B",
        "C. grandiflora",
        "C. orientalis",
    ]
    theta = 4 * 7e-9 * ppaper["NeCg0"] * douglas_data_orig.data.sum()

    p1 = [
        ppaper["NeCgF"] / ppaper["NeCg0"],
        ppaper["NeCbpA0"] / ppaper["NeCg0"],
        ppaper["NeCbpAF"] / ppaper["NeCg0"],
        ppaper["T1"] / (2 * ppaper["NeCg0"]),
    ]
    model1 = Cg_CbpA(p1, ns, pts_l)

    p2 = [
        ppaper["NeCo0"] / ppaper["NeCg0"],
        ppaper["NeCoF"] / ppaper["NeCg0"],
        ppaper["NeCbpB0"] / ppaper["NeCg0"],
        ppaper["NeCbpBF"] / ppaper["NeCg0"],
        ppaper["T1"] / (2 * ppaper["NeCg0"]),
        ppaper["T2"] / (2 * ppaper["NeCg0"]),
    ]
    model2 = Co_CbpB(p2, (10, 10), pts_l)

    p3 = [
        ppaper["NeCgF"] / ppaper["NeCg0"],
        ppaper["NeCo0"] / ppaper["NeCg0"],
        ppaper["NeCoF"] / ppaper["NeCg0"],
        ppaper["T1"] / (2 * ppaper["NeCg0"]),
        ppaper["T2"] / (2 * ppaper["NeCg0"]),
    ]
    model3 = Cg_Co(p3, ns, pts_l)

    fig = plt.figure(1, figsize=(6, 4), dpi=250)
    fig.clear()
    ax = fig.add_subplot(2, 3, 1)
    dadi.Plotting.plot_single_2d_sfs(
        douglas_data_orig.filter_pops([1, 3]).reorder_pops([2, 1]),
        vmin=1,
        cmap="hsv",
    )
    ax = fig.add_subplot(2, 3, 2)
    dadi.Plotting.plot_single_2d_sfs(
        douglas_data_orig.filter_pops([2, 4]).reorder_pops([2, 1]),
        vmin=1,
        cmap="hsv",
    )
    ax = fig.add_subplot(2, 3, 3)
    dadi.Plotting.plot_single_2d_sfs(
        douglas_data_orig.filter_pops([3, 4]), vmin=1, cmap="hsv"
    )

    ax = fig.add_subplot(2, 3, 4)
    dadi.Plotting.plot_single_2d_sfs(theta * model1, cmap="hsv", vmin=1)
    ax = fig.add_subplot(2, 3, 5)
    dadi.Plotting.plot_single_2d_sfs(theta * model2, cmap="hsv", vmin=1)
    ax = fig.add_subplot(2, 3, 6)
    dadi.Plotting.plot_single_2d_sfs(theta * model3, cmap="hsv", vmin=1)
    fig.tight_layout()
    fig.savefig("reproduction.pdf")

    p4 = [
        ppaper["NeCo0"] / ppaper["NeCg0"],
        ppaper["NeCbpA0"] / ppaper["NeCg0"],
        ppaper["NeCbpAF"] / ppaper["NeCg0"],
        ppaper["NeCbpB0"] / ppaper["NeCg0"],
        ppaper["NeCbpBF"] / ppaper["NeCg0"],
        ppaper["T1"] / (2 * ppaper["NeCg0"]),
        ppaper["T2"] / (2 * ppaper["NeCg0"]),
    ]
    model4 = CbpA_CbpB(p4, (10, 10), pts_l)

    fig = plt.figure(34, figsize=(6, 4), dpi=250)
    fig.clear()
    ax = fig.add_subplot(2, 3, 1)
    dadi.Plotting.plot_single_2d_sfs(
        douglas_data_orig.filter_pops([1, 2]), vmin=1, cmap="hsv"
    )
    ax.set_title("data")
    ax = fig.add_subplot(2, 3, 2)
    ax.set_title("douglas model")
    dadi.Plotting.plot_single_2d_sfs(theta * model4, cmap="hsv", vmin=1)

    segtetraploid_params = [
        100.00000000,
        7.166250,
        100.00000,
        18.3723521,
        0.0009271943,
        0.031397990,
    ]
    our_theta = 167.3206

    our_model = segtetraploid_bottleneck(segtetraploid_params, [20], pts_l)
    ax = fig.add_subplot(2, 3, 3)
    ax.set_title("our model")
    dadi.Plotting.plot_single_2d_sfs(our_theta * our_model, cmap="hsv", vmin=1)

    douglas_data_collapsed = douglas_data_orig.filter_pops([1, 2]).combine_pops(
        [1, 2]
    )
    model4_collapsed = model4.combine_pops([1, 2])
    our_model_collapsed = our_model.combine_pops([1, 2])

    ax = fig.add_subplot(2, 1, 2)
    ax.semilogy(douglas_data_collapsed, "-o", label="data")
    ax.semilogy(theta * model4_collapsed, "-o", label="douglas model")
    ax.semilogy(our_theta * our_model_collapsed, "-o", label="our model")
    ax.legend(loc="upper right")
    ax.set_xticks([0, 5, 10, 15, 20])

    print("Likelihoods of 2D data")
    print(
        "douglas model: {0:.2f}".format(
            dadi.Inference.ll(
                theta * model4, douglas_data_orig.filter_pops([1, 2])
            )
        )
    )
    print(
        "our model: {0:.2f}".format(
            dadi.Inference.ll(
                our_theta * our_model, douglas_data_orig.filter_pops([1, 2])
            )
        )
    )

    print("Likelihoods of collapsed data")
    print(
        "douglas model: {0:.2f}".format(
            dadi.Inference.ll(theta * model4_collapsed, douglas_data_collapsed)
        )
    )
    print(
        "our model: {0:.2f}".format(
            dadi.Inference.ll(
                our_theta * our_model_collapsed, douglas_data_collapsed
            )
        )
    )

    fig.tight_layout()
    fig.savefig("comparison.pdf")


if __name__ == "__main__":
    compare_results()
