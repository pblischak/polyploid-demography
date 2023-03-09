.. _Simulations:

Simulation-Based Model Validation
=================================

We used `SLiM 3 <https://messerlab.org/slim/>`__ to perform forward simulations
as a means to validate the diffusion approximation for collapsed polyploid spectra
in dadi. Below we give details on how to run each of the validation scripts.
The scripts can be found in the ``validation/`` folder in the GitHub repo. To
run multiple replicates, each call to SLiM can be wrapped with a Bash ``for``
loop with the loop index variable being passed to the ``rep`` parameter,
followed by any other parameters needed by the model and then the name of the
SLiM script:

.. code-block:: bash
   :caption: Example run for SLiM simulation.
   
   for i in {1..100}
   do
       slim -d "rep=${i}" [other parameters] <script-name>
   done

When run, each script generates a single realization of the demographic model
and outputs a site frequency spectrum that is the format required by dadi.
The corresponding Python code for comparing results from the diffusion approximation
in dadi can be found in the ``supp/SupplementalSfsData/plot_slim_dadi_comp.py``
script.

SLiM scripts
------------

``autotetraploid_snm.slim``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The standard neutral model for an autotetraploid.

.. code-block:: bash

   # Parameters: rep
   slim -d "rep=1" autotetraploid_snm.slim

``autotetraploid_bottleneck.slim``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates a bottleneck to population size ``nuBot`` times the original
population size.

.. code-block:: bash

   # Parameters: rep, nuBot
   slim -d "rep=1" -d "nuBot=0.5" autotetraploid_bottleneck.slim

``allotetraploid_iso.slim``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates an equilibrium population that splits at time ``T1`` in the past
before being sampled.

.. code-block:: bash
   
   # Parameters: rep, T1
   slim -d "rep=1" -d "T1=0.5" allotetraploid_iso.slim

``allotetraploid_bottleneck.slim``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates an equilibrium population that splits at time ``T1 + T2`` in the past
before being sampled, with a bottleneck of size ``nuBot`` occurring at time
``T2`` in the past.

.. code-block:: bash
   
   # Parameters: rep, T1, T2, nuBot
   slim -d "rep=1" -d "T1=0.5" -d "T2=0.25" -d "nuBot=0.5" allotetraploid_iso.slim

``segtetraploid_iso.slim``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates an equilibrium population that splits at time ``T1 + T2`` in the past
before being sampled, with homoeologous exchange occurring at a rate of ``dij``
starting at time ``T2`` in the past. Here ``dij`` is equivalent to the exchange
parameter :math:`e_{i \leftrightarrow j}` described in the manuscript. The use
of *iso* to describe this model is a bit of a misnomer since the subgenomes
are not isolated, but we kept it for naming consistency with the allotetraploid
models.

.. code-block:: bash
   
   # Parameters: rep, T1, T2, dij
   slim -d "rep=1" -d "T1=0.5" -d "T2=0.25" -d "dij=0.001" segtetraploid_iso.slim

``segtetraploid_bottleneck.slim``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates an equilibrium population that splits at time ``T1 + T2`` in the past
before being sampled, with a bottleneck of size ``nuBot`` and homoeologous
exchange occurring at a rate of ``dij`` starting at time ``T2`` in the past.
Here ``dij`` is equivalent to the exchange parameter
:math:`e_{i \leftrightarrow j}` described in the manuscript.

.. code-block:: bash
   
   # Parameters: rep, T1, T2, dij, nuBot
   slim -d "rep=1" -d "T1=0.5" -d "T2=0.25" -d "dij=0.001" -d "nuBot=0.5" segtetraploid_iso.slim

Inference and plotting scripts
------------------------------

In the ``inference/optimization/`` folder, there are Python and R scripts to
simulate data under the allotetraploid and segmental allotetraploid models, run
inference on the simulated data with dadi, and then collate and plot the maximum
likelihood parameter estimates across all replicates for each combination of
parameters. The subfolders for each model contain the simulated frequency spectra
and optimization results (in CSV format) for all combinations of parameters tested.

**References**

Haller, B. C. and P. W. Messer. 2019. SLiM 3: Forward genetic simulations beyond
the Wright–Fisher model. *Molecular Biology and Evolution* 36:632–-637.
