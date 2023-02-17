.. _Capsella:

Analysis of *Capsella bursa-pastoris*
=====================================

Code for all analyses of *Capsella* can be found in the ``capsella/`` folder
in the GitHub repo.

Main Scripts
------------

- ``run_*.py``: Scripts used to run demographic inference under the
  ``allotetraploid_bottleneck`` and ``segtetraploid_bottleneck`` models.
- ``capsella_*.csv``: Results of the independent optimization runs for each model.
- ``analyze_capsella_results.py``: Script to analyze all of the results for the
  maximum likelihood parameter estimates under each model, as well as code for
  estimating confidence intervals using the Fisher information Matrix and 
  propogation of uncertainty.

The SFS data file is named ``Capsella_intergene_4fold_corr_4pop_4_DSFS.fs``. It
can be read into Python as a ``Spectrum`` object using the code below:

.. code-block:: python
   :caption: Read a SFS file into Python as a ``dadi.Spectrum`` object.
   
   import dadi
   fs = dadi.Spectrum.from_file("Capsella_intergene_4fold_corr_4pop_4_DSFS.fs")