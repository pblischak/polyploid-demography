#!/usr/bin/env python3

"""
<< run_segtetraploid_iso.py >>
"""

# Import system-level libraries
from sys import argv, exit
import argparse as ap

# Import installed libraries
import dadi
import dadi.NLopt_mod


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
    # Print out the script docstring if only the script name is given
    if len(argv) < 2:
        print(__doc__)
        exit(0)

    # Set up argument parsing
    parser = ap.ArgumentParser(
        description="Options for run_segtetraploid_iso.py", add_help=True
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-T1",
        "--div_time1",
        action="store",
        type=float,
        required=True,
        metavar="\b",
        help="Length of time after split but before polyploid formation",
    )
    required.add_argument(
        "-T2",
        "--div_time2",
        action="store",
        type=float,
        required=True,
        metavar="\b",
        help="Time since polyploid formation",
    )
    required.add_argument(
        "-dij",
        "--ex_rate",
        action="store",
        type=float,
        required=True,
        metavar="\b",
        help="homoeologous exchange rate (M=2Nm)",
    )
    required.add_argument(
        "-r",
        "--rep",
        action="store",
        type=int,
        required=True,
        metavar="\b",
        help="Replicate number",
    )
    additional = parser.add_argument_group("additional arguments")
    additional.add_argument(
        "--gbs", action="store_true", help="Run GBS-like simulation"
    )
    additional.add_argument(
        "--optimization_runs",
        action="store",
        type=int,
        default=50,
        metavar="\b",
        help="Desired number of independent optimizations",
    )
    additional.add_argument(
        "--max_failures",
        action="store",
        type=int,
        default=50,
        metavar="\b",
        help="Maximum number of failed optimization attempts",
    )

    # Get arguments and store
    args = parser.parse_args()
    T1 = args.div_time1
    T2 = args.div_time2
    dij = args.ex_rate
    rep = args.rep
    optimization_runs = args.optimization_runs
    max_failures = args.max_failures
    gbs = args.gbs

    if gbs:
        mode = "_gbs"
    else:
        mode = ""

    # Open output file to record optimization results
    f_out = open(
        f"segtetraploid_iso/segtetraploid_iso_{T1}_{T2}_{dij}_{rep}{mode}.csv",
        "w",
    )
    print(
        "rep",
        "loglik",
        "nu_true",
        "nu_est",
        "T1_true",
        "T1_est",
        "T2_true",
        "T2_est",
        "dij_true",
        "dij_est",
        "theta",
        sep=",",
        file=f_out,
    )

    # Get data, sample sizes, and extract T
    data = dadi.Spectrum.from_file(
        f"segtetraploid_iso/segtetraploid_iso_{T1}_{T2}_{dij}_{rep}{mode}.fs"
    )
    ns = data.sample_sizes

    # Setting up grid points for extrapolation
    pts_l = [60, 70, 80]
    func = segtetraploid_iso
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Set bounds
    upper_bound = [100, 100, 100, 10]
    lower_bound = [1e-4, 1e-4, 1e-4, 1e-4]

    # True parameter values
    p_true = [1.0, T1, T2, dij]

    # Optimization loop
    opt_successes = 0
    opt_failures = 0
    opt_rep = 0
    while opt_successes < optimization_runs and opt_failures < max_failures:
        opt_rep += 1
        print(f"Starting optimization rep {opt_rep}")

        # Perturb parameters to random starting point
        p0 = dadi.Misc.perturb_params(
            p_true, fold=2, upper_bound=upper_bound, lower_bound=lower_bound
        )

        # Attempt optimization and record success/failure
        try:
            popt, LLopt = dadi.Inference.opt(
                p0,
                data,
                func_ex,
                pts_l,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                verbose=10,
            )
            opt_successes += 1
        except Exception:
            print(f"WARNING: Optimization failed for round {opt_rep+1}...")
            opt_failures += 1
            continue

        # Compare data and model, record results
        model = func_ex(popt, ns, pts_l)
        theta = dadi.Inference.optimal_sfs_scaling(model, data)
        print(
            opt_rep,
            LLopt,
            "1.0",
            popt[0],
            T1,
            popt[1],
            T2,
            popt[2],
            dij,
            popt[3],
            theta,
            sep=",",
            file=f_out,
        )
        if opt_rep % 10 == 0:
            f_out.flush()

    # Finished; print optimization stats
    print(f"\n\nOptimization finished for segtetraploid_iso rep {rep}")
    print(f"  Total number of optimization:      {opt_rep}")
    print(
        "  Number of succesful optimizations: "
        + f"{opt_successes} ({round(opt_successes/opt_rep * 100, 3)}%)"
    )
    print(
        "  Number of failed optimizations:    "
        + f"{opt_failures}  ({round(opt_failures/opt_rep * 100, 3)}%)\n\n"
    )
