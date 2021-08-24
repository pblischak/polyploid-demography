#!/usr/bin/env python3

"""
<< run_allotetraploid_bottleneck.py >>


"""

# Import system-level libraries
from sys import argv,exit
import argparse as ap

# Import installed libraries
import dadi
import dadi.NLopt_mod
import nlopt
import numpy as np

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

if __name__ == "__main__":
    # Set up argument parsing
    # All command line flags are optional
    parser = ap.ArgumentParser(
        description = "Options for run_allotetraploid_iso.py",
        add_help = True
    )
    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        '--optimization_runs', action="store", type=int, default=100,
        metavar='\b', help="Desired number of independent optimizations"
    )
    optional.add_argument(
        '--max_failures', action="store", type=int, default=50,
        metavar='\b', help="Maximum number of failed optimization attempts"
    )

    # Get arguments and store
    args              = parser.parse_args()
    optimization_runs = args.optimization_runs
    max_failures      = args.max_failures

    # Open output file to record optimization results
    f_out = open("capsella_allotetraploid_bottleneck.csv", 'w')
    print(
        "loglik", "nu0", "nuBot", "T1", "T2", "misid", "theta",
        sep=",", file=f_out
    )

    # Get data, sample sizes, and extract T
    data = dadi.Spectrum.from_file(
        "Capsella_intergene_4fold_corr_4pop_4_DSFS.fs"
    )
    data.pop_ids = ["CbpA","CbpB","grand","ori"]
    data = data.marginalize([2,3]).combine_pops([1,2])
    ns = data.sample_sizes

    # Setting up grid points for extrapolation
    pts_l = [50,60,70]
    func = allotetraploid_bottleneck
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Set bounds
    # nu0,nuBot,T1,T2,misid
    upper_bound = [100,100,100,100,1.0]
    lower_bound = [1e-4,1e-4,1e-4,1e-4,0.0]

    # True parameter values
    p_true = [1.0,1.0,1.0,1.0,0.05]

    # Optimization loop
    opt_successes = 0
    opt_failures = 0
    opt_rep = 0
    while opt_successes < optimization_runs and opt_failures < max_failures:
        opt_rep += 1
        print(f"Starting optimization rep {opt_rep}")

        # Perturb parameters to random starting point
        p0 = dadi.Misc.perturb_params(
            p_true, fold=2, upper_bound=upper_bound,
            lower_bound=lower_bound
        )

        # Attempt optimization and record success/failure
        try:
            popt,LLopt = dadi.Inference.opt(
                p0, data, func_ex, pts_l,
                lower_bound=lower_bound, upper_bound=upper_bound,
                maxeval=int(1e3), maxtime=int(3600), verbose=10
            )
            opt_successes +=1
        except:
            print(f"WARNING: Optimization failed for round {opt_rep+1}...")
            opt_failures +=1
            continue

        # Compare data and model, record results
        model = func_ex(popt, ns, pts_l)
        theta = dadi.Inference.optimal_sfs_scaling(model, data)
        print(
            LLopt, popt[0], popt[1], popt[2], popt[3], popt[4],
            theta, sep=",", file=f_out
        )
        if opt_rep % 10 == 0:
            f_out.flush()

    # Finished; print optimization stats
    print(f"\n\nOptimization finished for Capsella allotetraploid_bottleneck run")
    print(f"  Total number of optimization:      {opt_rep}")
    print(f"  Number of succesful optimizations: {opt_successes} ({round(opt_successes/opt_rep * 100, 3)}%)")
    print(f"  Number of failed optimizations:    {opt_failures}  ({round(opt_failures/opt_rep * 100, 3)}%)\n\n")
