import dadi
import matplotlib.pyplot as plt
import numpy as np


# Define demographic models
@dadi.Numerics.make_extrap_func
def douglas_et_al(params, ns, pts) -> dadi.Spectrum:
    """
    params = (nu_Cg, nu_Co, nu_CbpA, nu_CbpB, T1, T2, eij, misid)
    """
    (
        nu_Cg,
        nu_Co,
        nu_CbpA0,
        nu_CbpAF,
        nu_CbpB0,
        nu_CbpBF,
        T1,
        T2,
        eij,
        misid,
    ) = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Here we simulate the parental populations: grandiflora and orientalis
    phi = dadi.Integration.two_pops(phi, xx, T1, nu_Cg, nu_Co)

    # Set up exponential growth lambdas
    def nu_CbpA_func(t):
        # return np.exp(np.log(nu_CbpA) * t / T2)
        return nu_CbpA0 * (nu_CbpAF / nu_CbpA0) ** (t / T2)

    def nu_CbpB_func(t):
        # return np.exp(np.log(nu_CbpB) * t / T2)
        return nu_CbpB0 * (nu_CbpBF / nu_CbpBF) ** (t / T2)

    phi = dadi.Integration.two_pops(
        phi, xx, T2, nu_CbpA_func, nu_CbpB_func, m12=eij, m21=eij
    )

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["CbpA", "CbpB"]
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D, misid)
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


@dadi.Numerics.make_extrap_func
def segtetraploid_bottleneck(params, ns, pts) -> dadi.Spectrum:
    """
    params = (nu0,nuBot,T1,T2,eij,misid)
    ns = (n1,n2)

    Split into two populations of specified size.

    nu0: Size of populations after split.
    nuBot: Size of populations after bottleneck
    T1: Time in the past between split and polyploid formation
        (in units of 2*Na generations)
    T2: Time since polyploid formation and bottleneck
    n1,n2: Sample sizes of resulting Spectrum
    eij: Effective homoeologous exchange rate (2Nm)
    misid: probability of ancestral allele misidentification
    pts: Number of grid points to use in integration
    """
    nu0, nuBot, T1, T2, eij, misid = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot, m12=eij, m21=eij)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["CbpA", "CbpB"]
    )
    if misid > 0.0:
        fs_2D = dadi.Numerics.apply_anc_state_misid(fs_2D, misid)
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


def main() -> None:
    pts_l = [50, 60, 70]
    douglas_data = dadi.Spectrum.from_file(
        "Capsella_intergene_4fold_corr_4pop_4_DSFS.fs"
    )
    douglas_data.pop_ids = ["CbpA", "CbpB", "grand", "ori"]
    douglas_data = douglas_data.marginalize([2, 3]).combine_pops([1, 2])

    # Constants used for parameter conversion
    L = np.array(douglas_data).sum()  # sequence length
    mu = 7.0e-9  # mutation rate
    g = 1  # generation time

    # These are values taken from Model C in the Douglas et al. paper
    # and are converted to dadi units
    douglas_real_params = [
        # Final size from model C, then multiplied by exponential growth rate
        # since the time of divergence (128k years)
        840000 * (1.0 - 4.4e-6) ** 128000,  # N_Cg
        4000 * (1.0 + 2.6e-5) ** 128000,  # N_Co
        37000 * (1.0 + 4.8e-7) ** 128000,  # N_CbpA0
        37000,  # N_CbpAF
        75000 * (1.0 + 1.0 + 4.8e-7),  # N_CbpB0
        75000,  # N_CbpBF
        (931 - 128) * 1000,  # T1
        128 * 1000,  # T2
    ]

    # Here we'll use the same theta and Nref between the models
    theta = 167.3206
    Nref = theta / (4 * g * mu * L)

    douglas_params = [
        douglas_real_params[0] / Nref,
        douglas_real_params[1] / Nref,
        douglas_real_params[2] / Nref,
        douglas_real_params[3] / Nref,
        douglas_real_params[4] / Nref,
        douglas_real_params[5] / Nref,
        douglas_real_params[6] / (2 * Nref),
        douglas_real_params[7] / (2 * Nref),
        0.0009271943,
        0.031397990,
    ]

    print(douglas_params)

    segtetraploid_params = [
        100.00000000,
        7.166250,
        100.00000,
        18.3723521,
        0.0009271943,
        0.031397990,
    ]

    print(segtetraploid_params)

    douglas_sfs = douglas_et_al(
        douglas_params, douglas_data.sample_sizes, pts_l
    )
    segtetraploid_sfs = segtetraploid_bottleneck(
        segtetraploid_params, douglas_data.sample_sizes, pts_l
    )

    # dadi.Plotting.plot_1d_comp_Poisson(douglas_sfs, segtetraploid_sfs)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.semilogy(theta * douglas_sfs, "-og")
    ax1.semilogy(theta * segtetraploid_sfs, "-ob")
    ax1.semilogy(douglas_data, "-or")
    ax2.plot(
        dadi.Inference.Anscombe_Poisson_residual(
            douglas_sfs, segtetraploid_sfs
        ),
        "-o",
    )
    plt.savefig("TestComparison.pdf")


if __name__ == "__main__":
    main()
