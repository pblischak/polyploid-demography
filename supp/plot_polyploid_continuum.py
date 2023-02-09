#!/usr/bin/env python3

import dadi
import dadi.NLopt_mod
import matplotlib
import matplotlib.pyplot as plt


# Defining demographic model
def segtetraploid_iso(params, ns, pts):
    """
    params = (nu,T1,T2,dij)
    ns = (n1,n2)

    Split into two populations of specifed size.

    nu: Size of populations after split.
    T: Time in the past of split (in units of 2*Na generations)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu, T1, T2, dij = params
    new_ns = [int(ns[0] / 2), int(ns[0] / 2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu, nu)
    phi = dadi.Integration.two_pops(phi, xx, T2, nu, nu, m12=dij, m21=dij)

    fs_2D = dadi.Spectrum.from_phi(
        phi, new_ns, (xx, xx), pop_ids=["sub1", "sub2"]
    )
    fs_1D = fs_2D.combine_pops([1, 2])
    return fs_1D


if __name__ == "__main__":
    matplotlib.rcParams.update({"font.size": 16})
    exchange_rates = [0.0, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2]
    N_ref = 1000.0
    T1 = 1.0
    T2 = 0.25

    pts_l = [60, 70, 80]
    func = segtetraploid_iso
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    fig, ax = plt.subplots(figsize=(7, 5))
    viridis = ["#fde725", "#7ad151", "#22a884", "#2a788e", "#414487", "#440154"]
    for i, e in enumerate(exchange_rates):
        print(e)
        sfs = func_ex([1.0, 1.0, 0.25, (2.0 * e * N_ref)], (40,), pts_l)
        ax.semilogy(sfs, "-o", color=viridis[i], label=str(e))
    ax.legend()
    plt.savefig("polyploid_continuum_sfs.svg")
