#!/usr/bin/env python3

"""
<< run_segtetraploid_bottleneck.py >>


"""

# Import system-level libraries
from sys import argv,exit
from os import system
import argparse as ap

# Import installed libraries
import dadi
import dadi.NLopt_mod
import nlopt
import numpy as np

# Function to match floats in SLiM output files
def match_slim_outfloat(flt):
    """
    SLiM likes to drop the decimal places from floats that could be integers
    (e.g., 1.0 becomes just 1). This causes issues which matching up file names
    so here we convert the floats appropriately here.
    """
    str_flt = str(round(flt,1))
    if str_flt[-1] == '0':
        return str_flt[0]
    else:
        return(str(flt))

# Setting up SLiM run
def run_slim(nuBot, T1, T2, eij, rep):
    """
    Takes arguments needed to run SLiM and makes a system call to generate
    a simulated SFS.
    """
    cmd = [
        "slim",
        f"-d nuBot={nuBot}",
        f'-d "T1={T1}"',
        f'-d "T2={T2}"',
        f'-d "eij={eij}"',
        f'-d "rep={rep}"',
        "./SLiM/segtetraploid_bottleneck.slim"
    ]
    print(" ".join(cmd))
    system(" ".join(cmd))

# Defining demographic model
def segtetraploid_bottleneck(params, ns, pts):
    """
    params = (nu,nuBot,T1,T2,Eij)
    ns = (n1,n2)

    Split into two populations of specifed size.

    nu: Size of populations after split.
    nuBot: Size of populations after bottleneck
    T1: Time in the past between split and polyploid formation(in units of 2*Na generations)
    T2: Time since polyploid formation and bottleneck
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,nuBot,T1,T2,eij = params
    new_ns = [int(ns[0]/2),int(ns[0]/2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu, nu)
    phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot, m12=eij, m21=eij)

    fs_2D = dadi.Spectrum.from_phi(phi, new_ns, (xx,xx), pop_ids=['sub1','sub2'])
    fs_1D = fs_2D.combine_pops([1,2])
    return fs_1D

if __name__ == "__main__":
    # Print out the script docstring if only the script name is given
    if len(argv) < 2:
        print(__doc__)
        exit(0)

    # Set up argument parsing
    parser = ap.ArgumentParser(
        description = "Options for run_segtetraploid_bottleneck.py",
        add_help = True
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        '-nu', '--nuBot', action="store", type=float, required=True,
        metavar='\b', help="Size of population after bottleneck"
    )
    required.add_argument(
        '-T1', '--div_time1', action="store", type=float, required=True,
        metavar='\b', help="Length of time after split but before polyploid formation"
    )
    required.add_argument(
        '-T2', '--div_time2', action="store", type=float, required=True,
        metavar='\b', help="Time since polyploid formation"
    )
    required.add_argument(
        '-eij', '--ex_rate', action="store", type=float, required=True,
        metavar='\b', help="homoeologous exchange rate (M=2Nm)"
    )
    required.add_argument(
        '-r', '--rep', action="store", type=int, required=True,
        metavar='\b', help="Replicate number"
    )
    additional = parser.add_argument_group("additional arguments")
    additional.add_argument(
        '--optimization_runs', action="store", type=int, default=50,
        metavar='\b', help="Desired number of independent optimizations"
    )
    additional.add_argument(
        '--max_failures', action="store", type=int, default=50,
        metavar='\b', help="Maximum number of failed optimization attempts"
    )
    additional.add_argument(
        '--skip_slim', action="store_true",
        help="Skip SLiM simulation (use if already complete)"
    )

    # Get arguments and store
    args              = parser.parse_args()
    nuBot             = args.nuBot
    T1                = args.div_time1
    T2                = args.div_time2
    eij               = args.ex_rate
    rep               = args.rep
    optimization_runs = args.optimization_runs
    max_failures      = args.max_failures
    skip_slim         = args.skip_slim

    # Run SLiM simulation to generate SFS
    if not skip_slim:
        run_slim(nuBot,T1,T2,eij,rep)

    # Deal with SLiM converting floats to integers
    T1_str = match_slim_outfloat(T1)
    T2_str = match_slim_outfloat(T2)

    # Open output file to record optimization results
    f_out = open(
        f'segtetraploid_bottleneck/segtetraploid_bottleneck_{nuBot}_'+T1_str+"_"+T2_str+f'_{eij}_{rep}.csv', 'w'
    )
    print(
        "rep","loglik","nu_true","nu_est","nuBot_true","nuBot_est","T1_true",
        "T1_est","T2_true","T2_est","eij_true","eij_est","theta",
        sep=",", file=f_out
    )

    # Get data, sample sizes, and extract T
    data = dadi.Spectrum.from_file(
        f'segtetraploid_iso/segtetraploid_nottlenecl_{nuBot}_'+T1_str+"_"+T2_str+f'_{eij}_{rep}.fs'
    )
    ns = data.sample_sizes

    # Setting up grid points for extrapolation
    pts_l = [60,70,80]
    func = segtetraploid_iso
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Set bounds
    upper_bound = [100, 100, 100, 100, 10]
    lower_bound = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4]

    # True parameter values
    p_true = [1.0,nuBot,T1,T2,eij]

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
                verbose=10
            )
            opt_successes +=1
        except:
            print(f"WARNING: Optimization failed for round {opt_rep+1}...")
            opt_failures +=1
            continue

        # Compare data and model, record results
        model = func_ex(popt, ns, pts_l)
        theta = dadi.Inference.optimal_sfs_scaling(model, data)
        print(opt_rep,LLopt,"1.0",popt[0],nuBot,popt[1],T1,popt[2],T2,popt[3],eij,popt[4],theta,sep=",",file=f_out)
        if opt_rep % 10 == 0:
            f_out.flush()

    # Finished; print optimization stats
    print(f"\n\nOptimization finished for segtetraploid_iso rep {rep}")
    print(f"  Total number of optimization:      {opt_rep}")
    print(f"  Number of succesful optimizations: {opt_successes} ({round(opt_successes/opt_rep * 100, 3)}%)")
    print(f"  Number of failed optimizations:    {opt_failures}  ({round(opt_failures/opt_rep * 100, 3)}%)\n\n")
