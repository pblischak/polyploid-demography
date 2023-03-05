#!/usr/bin/env python3

import dadi
import matplotlib.pyplot as plt


def allotetraploid_bottleneck(params, ns, pts):
    """
    params = (nu0,nuBot,T1,T2)
    ns = (n1,n2)
    Split into two populations of specifed size.
    nu0: Size of populations after split.
    nuBot: Size of populations after bottleneck.
    T1: Length of time between split and bottleneck.
    T2: Length of bottleneck.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu0, nuBot, T1, T2 = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["sub1", "sub2"]
    )
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


def segtetraploid_bottleneck(params, ns, pts):
    """
    params = (nu,nuBot,T1,T2,eij)
    ns = (n1,n2)
    Split into two populations of specifed size.
    nu: Size of populations after split.
    nuBot: Size of populations after bottleneck
    T1: Time in the past between split and polyploid
        formation(in units of 2*Na generations)
    T2: Time since polyploid formation and bottleneck
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu, nuBot, T1, T2, eij = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu, nu)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot, m12=eij, m21=eij)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["sub1", "sub2"]
    )
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


if __name__ == "__main__":
    """
    Make a single plot with two bottleneck sizes: 0.1 x nu0 and 0.5 x nu0

    These two bottleneck sizes will be in the rows and the three models
    will each have their own column.
    """
    pts_l = [60, 70, 80]
    theta = 5000.0

    # Allotetraploid Bottleneck
    allo_func_ex = dadi.Numerics.make_extrap_log_func(allotetraploid_bottleneck)

    allo_slim10 = dadi.Spectrum.from_file(
        "AveragedResults_Complete/AllotetraploidBottleneck/"
        + "avg_allotetraploid_bottleneck_0.1_1_0.25____.fs"
    )
    allo_dadi10 = theta * allo_func_ex([1.0, 0.1, 1.0, 0.25], [40], pts_l)

    allo_slim50 = dadi.Spectrum.from_file(
        "AveragedResults_Complete/AllotetraploidBottleneck/"
        + "avg_allotetraploid_bottleneck_0.5_1_0.25____.fs"
    )
    allo_dadi50 = theta * allo_func_ex([1.0, 0.5, 1.0, 0.25], [40], pts_l)

    # Segmental Allotetraploid Bottleneck
    seg_func_ex = dadi.Numerics.make_extrap_log_func(segtetraploid_bottleneck)

    seg_slim10 = dadi.Spectrum.from_file(
        "AveragedResults_Complete/SegtetraploidBottleneck/"
        + "avg_segtetraploid_0.1_1_0.25_5e-06.fs"
    )
    seg_dadi10 = theta * seg_func_ex([1.0, 0.1, 1.0, 0.25, 0.01], [40], pts_l)

    seg_slim50 = dadi.Spectrum.from_file(
        "AveragedResults_Complete/SegtetraploidBottleneck/"
        + "avg_segtetraploid_0.5_1_0.25_5e-06.fs"
    )
    seg_dadi50 = theta * seg_func_ex([1.0, 0.5, 1.0, 0.25, 0.01], [40], pts_l)

    # Autotetraploid Bottleneck
    auto_func_ex = dadi.Numerics.make_extrap_log_func(
        dadi.Demographics1D.two_epoch
    )

    auto_slim10 = dadi.Spectrum.from_file(
        "AveragedResults_Complete/AutotetraploidBottleneck/"
        + "autotetraploid_bottleneck_0.1_1.0_x.fs"
    )
    auto_dadi10 = theta * auto_func_ex([0.1, 1.0], [40], pts_l)

    auto_slim50 = dadi.Spectrum.from_file(
        "AveragedResults_Complete/AutotetraploidBottleneck/"
        + "autotetraploid_bottleneck_0.5_1.0_x.fs"
    )
    auto_dadi50 = theta * auto_func_ex([0.5, 1.0], [40], pts_l)

    ordered_slim_spectra = [
        auto_slim50,
        seg_slim50,
        allo_slim50,
        auto_slim10,
        seg_slim10,
        allo_slim10,
    ]

    ordered_dadi_spectra = [
        auto_dadi50,
        seg_dadi50,
        allo_dadi50,
        auto_dadi10,
        seg_dadi10,
        allo_dadi10,
    ]

    fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(7, 5))
    for i in range(2):
        for j in range(3):
            axes[i, j].semilogy(
                ordered_dadi_spectra[i * 3 + j],
                "-o",
                markersize=2,
                color="silver",
                linewidth=0.9,
            )
            axes[i, j].semilogy(
                ordered_slim_spectra[i * 3 + j],
                "-ok",
                markersize=2,
                linewidth=0.9,
            )

    plt.savefig("slim_dadi_comp.pdf")
