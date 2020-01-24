#!/usr/bin/env python3

import dadi
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit

def combine_pops(fs, idx=[0,1]):
    """
    Combine the frequency spectra of two populations.
        fs:  Spectrum object (2D or 3D).
        idx: Indices for populations being collapsed. (defaul=[0,1])
    
    The function will always return the combined populations along the
    first axis of the sfs.
    """
    ns = fs.sample_sizes
    if len(ns) == 3:
        fs_tmp = np.array(fs)
        if idx == [0,1]:
            fs2 = np.zeros((ns[0]+ns[1]+1,ns[2]+1))
            for ii in range(ns[0]+1):
                for jj in range(ns[1]+1):
                    for kk in range(ns[2]+1):
                        fs2[ii+jj,kk] += fs_tmp[ii,jj,kk]
        elif idx == [0,2]:
            fs2 = np.zeros((ns[0]+ns[2]+1,ns[1]+1))
            for ii in range(ns[0]+1):
                for jj in range(ns[2]+1):
                    for kk in range(ns[1]+1):
                        fs2[ii+kk,jj] += fs_tmp[ii,jj,kk]
        elif idx == [1,2]:
            fs2 = np.zeros((ns[1]+ns[2]+1,ns[0]+1))
            for ii in range(ns[1]+1):
                for jj in range(ns[2]+1):
                    for kk in range(ns[0]+1):
                        fs2[jj+kk,ii] += fs_tmp[ii,jj,kk]
        else:
            print("Error: did not recognize population indices: {}".format(idx))
            exit(-1)
    elif len(ns) == 2:
        fs_tmp = np.array(fs)
        fs2    = np.zeros((ns[0]+ns[1]+1,))
        for ii in range(ns[0]+1):
            for jj in range(ns[1]+1):
                fs2[ii+jj] += fs_tmp[ii,jj]
    else:
        print("Error: could not combine populations.")
        exit(-1)
    return dadi.Spectrum(fs2)

def subgenome_model_2D_1(params, ns, pts):
    """
    A model for subgenomes
    """
    delta, T = params
    """
    if delta < 0.0:
        delta = 0.0
    elif delta > 1.0:
        delta = 1.0
    else:
        pass
    """
    
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    if delta == 0.0:
        phi = dadi.Integration.two_pops(phi, xx, T, 1.0, 1.0)
    else:
        phi = dadi.Integration.two_pops(phi, xx, T, 1.0, 1.0, m12=delta, m21=delta)
    
    
    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, [xx,xx], [0.0001,0.0001], [2,2])
    fs2 = combine_pops(fs)
    """
    fs_tmp = np.array(fs)
    fs2 = np.zeros((ns[0]+ns[1]+1,))
    for i in range(ns[0]+1):
        for j in range(ns[1]+1):
            fs2[i+j] += fs_tmp[i,j]
    """
    return fs2

def subgenome_model_2D_2(params, ns, pts):
    """
    A model for subgenomes
    """
    delta, T = params
    """
    if delta < 0.0:
        delta = 0.0
    elif delta > 1.0:
        delta = 1.0
    else:
        pass
    """
    
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    if delta == 0.0:
        phi = dadi.Integration.two_pops(phi, xx, T, 2.0, 1.0)
    else:
        phi = dadi.Integration.two_pops(phi, xx, T, 2.0, 1.0, m12=delta, m21=delta)
    
    
    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, [xx,xx], [0.0001,0.0001], [4,2])
    fs2 = combine_pops(fs)
    """
    fs_tmp = np.array(fs)
    fs2 = np.zeros((ns[0]+ns[1]+1,))
    for i in range(ns[0]+1):
        for j in range(ns[1]+1):
            fs2[i+j] += fs_tmp[i,j]
    """
    return fs2

def subgenome_model_3D(params, ns, pts):
    """
    A 3D model for 3 subgenomes.
    """
    delta, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)
    if delta == 0.0:
        delta = 0.0001
    phi = dadi.Integration.three_pops(phi, xx, T, 1.0, 1.0, 1.0, m12=delta, m13=delta,
                                      m23=delta, m21=delta, m31=delta, m32=delta)
    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, [xx,xx,xx], [0.0001, 0.0001, 0.0001], [2,2,2])
    fs_tmp = np.array(fs)
    fs2 = np.zeros((ns[0]+ns[1]+ns[2]+1,))
    for i in range(ns[0]+1):
        for j in range(ns[1]+1):
            for k in range(ns[2]+1):
                fs2[i+j+k] += fs_tmp[i,j,k]
    #fs2 = combine_pops(fs)
    #fs3 = combine_pops(fs2)
    return dadi.Spectrum(fs2)

if __name__ == "__main__":
    """
    Let's go.
    """
    D = float(argv[1]) # get delta from the command line
    t = float(argv[2]) # get integration time from the command line
    
    pts_l = [80,90,100]
    N = 20
    #func = subgenome_model
    #func_ex = dadi.Numerics.make_extrap_log_func(func)
    sfs  = subgenome_model_2D_1([D,t], [N,N], 80)
    sfs2 = subgenome_model_2D_2([D,t], [2*N,N], 80)
    sfs3 = subgenome_model_3D([D,t], [N,N,N], 80)
    #sfs = func_ex([D,t], [N,N], pts_l)
    
    # Increase the font size
    plt.rcParams.update({'font.size': 20})
    
    dadi.Plotting.plot_1d_fs(sfs2)
    plt.savefig("pgrp2/1D_2and2_{}_{}_2.svg".format(D,t))
    plt.close()
    
    dadi.Plotting.plot_1d_fs(sfs2)
    plt.savefig("pgrp2/1D_4and2_{}_{}.svg".format(D,t))
    plt.close()
    
    dadi.Plotting.plot_1d_fs(sfs3)
    plt.savefig("pgrp2/1D_2and2and2_{}_{}.svg".format(D,t))
    plt.close()
