#!/usr/bin/env python3

"""
<< segtetraploid_bottleneck.py >>


"""

# Import system-level libraries
from sys import argv,exit
from os import system
import dadi
import numpy as np

# Define demographic model
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
    with open('segtetraploid_bottleneck.txt', 'w') as f_out:
        # Setting up grid points for extrapolation
        ns = [40]
        pts_l = [60,70,80]
        func = segtetraploid_bottleneck
        func_ex = dadi.Numerics.make_extrap_log_func(func)

        # True parameter values
        for nuBot in [0.1, 0.25, 0.5]:
            for T1 in [0.5, 1.0, 1.5, 2.0]:
                for T2 in [0.25, 0.5, 1.0]:
                    for eij in [0.001, 0.01, 0.1]:
                        p_true = [1.0,nuBot,T1,T2,eij]
                        print(p_true)
                        H = dadi.Godambe.get_godambe(
                            func_ex, pts_l, [], p_true,
                            func_ex(p_true, ns, pts_l),
                            eps=0.01, log=True, just_hess=True
                        )
                        print(f"## Eigendecomposition of Hessian for params nuBot,T1,T2,eij = {p_true[1:]}:\n", file=f_out)
                        eval,evec = np.linalg.eig(H)
                        print("\t".join([str(x) for x in eval]), "\n", sep = "", file=f_out)
                        for row in range(evec.shape[0]):
                            print("\t".join([str(x) for x in evec[row,:]]), file = f_out)
                        print("\n\n", file = f_out)
