import matplotlib.pyplot as plt
import dadi


@dadi.Numerics.make_extrap_func
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


# For when I figure out SLiM scaling
# L = 1e6
# N,K = 1000,4
# theta = 5000
# mu = theta/(2*K*N*L)


theta = 5000

pbase = [1.0, 1.0, 0.25]
ns = (40,)
pts_l = [45, 50, 55]

fs_low = segtetraploid_iso(pbase + [0], ns, pts_l)
fs_high = segtetraploid_iso(pbase + [20], ns, pts_l)
fs_mix = 0.5 * fs_low + 0.5 * fs_high

fs_mid = segtetraploid_iso(pbase + [0.21], ns, pts_l)

data = fs_mix * theta

popt, llopt = dadi.Inference.opt(
    pbase + [0.3],
    data,
    segtetraploid_iso,
    pts_l,
    lower_bound=[1e-2, 0, 0, 0],
    upper_bound=[100, 5, 2, 20],
    verbose=1,
)
fs_fit = segtetraploid_iso(popt, ns, pts_l)

fig = plt.figure(2113, figsize=(6, 3), dpi=250)
fig.clear()
for ii in range(1, 3):
    ax = fig.add_subplot(1, 2, ii)
    ax.plot(fs_low, label=r"$E_{i \leftrightarrow j} = 0$")
    ax.plot(fs_high, label=r"$E_{i \leftrightarrow j} = 20$")
    ax.plot(
        fs_mix,
        label=(
            r"$\frac{1}{2} (E_{i \leftrightarrow j} = 0) + "
            + r"\frac{1}{2} (E_{i\leftrightarrow j} = 20)$"
        ),
    )
    ax.plot(fs_mid, label=r"$E_{i \leftrightarrow j} = 0.21$")
    # Plot fit result, accounting for multinomial scaling
    theta_scale = dadi.Inference.optimal_sfs_scaling(fs_fit, data / theta)
    ax.plot(theta_scale * fs_fit, label=r"Fit model")
    ax.set_xlabel("Sample allele frequency")

fig.axes[1].set_yscale("log")
fig.axes[0].legend(fontsize="x-small", loc="upper right")

fig.tight_layout(pad=1)

fig.savefig("mixture.pdf")
