#!/usr/bin/env python3

"""
<< sim_allotetraploid_iso.py >>


"""

# Import system-level libraries
from sys import argv,exit
import subprocess
import argparse as ap

# Import installed libraries
import dadi
import numpy as np

# Setting up SLiM run
def run_slim(T,rep,mode=""):
    """
    Takes arguments needed to run SLiM and makes a system call to generate
    a simulated SFS.
    """
    cmd = [
        "slim",
        f'-d "T1={T}"',
        f'-d "rep={rep}"',
        "./SLiM/allotetraploid_iso"+mode+".slim"
    ]
    print(" ".join(cmd))
    output = subprocess.Popen(" ".join(cmd), shell=True, stdout=subprocess.PIPE)
    fs = dadi.Spectrum([
        float(i) for i in output.stdout.read().decode().split("\n")[-2].split()
    ])
    return fs

if __name__ == "__main__":
    # Print out the script docstring if only the script name is given
    if len(argv) < 2:
        print(__doc__)
        exit(0)

    # Set up argument parsing
    parser = ap.ArgumentParser(
        description = "Options for run_allotetraploid_iso.py",
        add_help = True
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        '-T', '--div_time', action="store", type=float, required=True,
        metavar='\b', help="Divergence time between subgenomes"
    )
    required.add_argument(
        '-r', '--rep', action="store", type=int, required=True,
        metavar='\b', help="Replicate number"
    )
    additional = parser.add_argument_group("additional arguments")
    additional.add_argument(
        '-l', '--nloci', type=int, default=5000,
        metavar='\b', help="Number of GBS loci. Not used for WGS simulation"
    )
    additional.add_argument(
        '--gbs', action="store_true",
        help="Run GBS-like simulation"
    )

    # Get arguments and store
    args              = parser.parse_args()
    T                 = args.div_time
    rep               = args.rep
    nloci             = args.nloci
    gbs               = args.gbs

    # Run SLiM simulation to generate SFS
    fs = dadi.Spectrum(np.zeros((41,)))
    if gbs:
        mode = "_gbs"
    else:
        mode = ""
        nloci = 10

    for i in range(nloci):
        fs += run_slim(T,rep,mode=mode)

    fs.to_file(
        f'allotetraploid_iso/allotetraploid_iso_{T}_{rep}{mode}.fs'
    )
