# %%
import dadi
import numpy as np

# Define demographic model (Note the use of the Python decorator here)
@dadi.Numerics.make_extrap_func
def allotetraploid_iso(params, ns, pts):
    nu,T = params
    # Fix theta, to keep analysis simpler (for now)
    theta0 = 1e3
    new_ns = [int(ns[0]/2),int(ns[0]/2)]

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx, theta0=theta0)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu, nu, theta0=theta0)

    fs_2D = dadi.Spectrum.from_phi(phi, new_ns, (xx,xx), pop_ids=['sub1','sub2'])
    fs_1D = fs_2D.combine_pops([1,2])
    return fs_1D

# %%

# How do I explore for more parameters?
# For 3 parameters, can make a 3D plot.
# For more, project sloppy axes downward. Think I need to formulate this in terms of marginalizing 
#     the ND Gaussian.

# Godambe/hessian works well for shorter T. Not so well for long T.
# It's getting directions correct, but magnitude of sloppy axes too small.
# Could it be an issue with pts?
p_true = [2.0,4]

# Running with small ns to keep testing quick
ns = [10]
pts_l = [12,14,16]
func_ex = allotetraploid_iso
fs = func_ex(p_true, ns, pts_l)

# %%
# Generate bootstrap data sets. Fixing seed to make analysis repeatable.
Nboot = 50
np.random.seed(2131)
boots = [fs.sample() for ii in range(Nboot)]

# %%
# Calculate Godambe uncertainties and Hessian matrix
uncerts, _, H = dadi.Godambe.GIM_uncert(func_ex, pts_l, boots, p_true, fs, log=True,
                      multinom=False, eps=1e-3, return_GIM=True)

# Take Singular Value Decomposition of Hessian
u, s, vh = np.linalg.svd(H)

# %%
# Fit bootstrap data sets
boot_fits = []
for boot_ii, data in enumerate(boots):
    p0 = dadi.Misc.perturb_params(p_true, fold=0.5)
    # Only doing one fit per data set, for a quick rough analysis
    popt, llopt = dadi.Inference.opt(p0, data, func_ex, pts_l, 
                                     lower_bound=[1e-1, 0], upper_bound=[10, 10])
    print(boot_ii, popt)
    boot_fits.append(popt)
boot_fits = np.array(boot_fits)


# %%
print('Godambe uncertainties:', uncerts)
print('Bootstrap uncertainties:', np.std(np.log(boot_fits), axis=0))

print('Hessian eigvenalues:', s)
print('Stiffest parameter combination:', u[:,0])
print('Sloppiest parameter combination:', u[:,-1])

import matplotlib.pyplot as plt
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1, aspect='equal')

# Plot bootstraps (in log parameter space)
ax.plot(np.log(boot_fits[:,0]), np.log(boot_fits[:,1]), 'o')

# Plot scaled principal axes of ellipses (out to 3 sigma)
# For two-parameter model...
x0, y0 = np.log(p_true)
ax.plot([x0 - 3/np.sqrt(s[0])*u[0,0], x0 + 3/np.sqrt(s[0])*u[0,0]], [y0 - 3/np.sqrt(s[0])*u[1,0], y0+3/np.sqrt(s[0])*u[1,0]], '-g')
ax.plot([x0 - 3/np.sqrt(s[1])*u[0,1], x0 + 3/np.sqrt(s[1])*u[0,1]], [y0 - 3/np.sqrt(s[1])*u[1,1], y0+3/np.sqrt(s[1])*u[1,1]], '-c')

ax.set_xlabel('log(nu)')
ax.set_ylabel('log(T)')

plt.show()

# %%
