.. _Tutorial:

Tutorial
========

In this tutorial we'll walk through building demographic models in both SLiM
and dadi, focusing on how to capture the appropriate types of chromosomal
dynamics for the different types of polyploids.
We'll also cover how to think about converting polyploid parameter values from 
"real" units to the scaled, genetic units used by dadi. Finally, we'll go through
a brief overview of how to analyze real polyploid data from a VCF file.

Simulating Polyploid Data in SLiM
---------------------------------

.. note::

   The scripts for simulating frequency spectra for various types of polyploids
   with SLiM can be found in the ``validation/`` folder in the main GitHub repo.

Here, we will give an overview of how SLiM models are constructed to accommodate
different types of inheritance patterns (e.g., polysomic vs. disomic) that
are typically associated with the various categories of polyploids. The main
differences between the models revolve around how the initial populations
reach approximate equilibrium through the burn-in phase before undergoing
any changes associated with demographic events or polyploid formation.

Tetrasomic Autopolyploid
~~~~~~~~~~~~~~~~~~~~~~~~

For the tetrasomic autopolyploid example used in the paper, we simulated
a population at approximate equilibrium, meaning that the ancestral population
before any demographic events was a tetrasomic autopolyploid as well. Although
we don't consider it here, it would also be possible to do a more direct
simulation of autotetraploid formation by assuming a diploid ancestral population
that doubles at some :math:`T` generations in the past, and that the doubled
population undergoes tetrasomic inheritance for the remaining :math:`T`
generations of the simulation.

Building Polyploid Demographic Models in Dadi
---------------------------------------------



Converting Parameter Values
---------------------------

To best understand the relationship between the units used in forward simulators
like SLiM, where we simulate the actual number of individuals, generations, etc.,
the original description of dadi's genetic units in the manual is a good place
to start: `link <https://dadi.readthedocs.io/en/latest/user-guide/specifying-a-model/>`_
(see the subsection titled **Units**). 

Working with Real Data
----------------------

.. note::

   There is a ``scripts/`` folder in the main GitHub repo with a script named
   ``sfs_from_vcf.py``. It contains more information and comments for each step
   described below.

To use real data stored in a VCF file, dadi has several built-in functions that
allow users to read in and format their data into a site frequency spectrum for
use in downstream analyses. These steps include reading in the VCF file using
population information (individual assignment to populations)

.. code-block:: python
   :linenos:
   :caption: Generating a frequency spectrum from a VCF file.

   import dadi

   # Put name of VCF file here. Can be gzipped or not.
   vcf_file = "variants.vcf.gz"

   # Put the name of your population info file
   population_info = "popinfo.txt"

   # Read the VCF data into a data dictionary using the population info
   data_dict = dadi.Misc.data_dict_from_vcf(
       vcf_file,
       population_info

   )

   # Convert the data dictionary into a Spectrum object. A Spectrum object is
   # a wrapper around a NumPy array and is what dadi uses to work with site
   # frequency spectra
   sfs = dadi.Spectrum.from_data_dict(
       data_dict=data_dict,
       pop_ids=["pop1", "pop2"] # names for populations in info file
       projections=(20,40) # number of chromosomes sampled from each population
       polarized=False # Can we determine ancestral vs. derived allelic states
   )

   # Write the SFS to file
   sfs.to_file("spectrum.fs")
