# Parameter Inference and Optimization

The code in this folder contains code for simulating SFS data using SLiM, running
maximum likelihood parameter inference on the simulated data using dadi, and
plotting the results in R. The four models tested are:

 - allotetraploid_iso
 - allotetraploid_bottleneck
 - segtetraploid_iso
 - segtetraploid_bottleneck

The ``SLiM/`` directory contains code to implement the four models used for
testing inference in dadi. The ``sim_*.py`` scripts run the simulations,
the ``run_*.py`` scripts run the inference, and the ``plot_*.R`` scripts find
the maximum likelihood parameters across replicates and plot the results.

Required Python packages:
 - dadi
 - numpy

Required R packagesL
 - tidyverse
 - patchwork