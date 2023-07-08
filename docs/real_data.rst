.. _RealData:

Working with Real Data
======================

.. note::

   There is a ``scripts/`` folder in the main GitHub repo with a script named
   ``sfs_from_vcf.py``. It contains more information and comments for each step
   described below.

To use data stored in a VCF file, dadi has several built-in functions that
allow users to read in and format their data into a site frequency spectrum for
use in downstream analyses. These steps include reading in the VCF file using
population information (individual assignment to populations), which can handle
individuals of any ploidy level based on the number of alleles in their
genotype (GT) field.

.. code-block:: python
   :linenos:
   :caption: Generating a frequency spectrum from a VCF file.

   import dadi

   # Put name of VCF file here. Can be gzipped or not.
   vcf_file = "example.vcf.gz"

   # Put the name of your population info file
   popinfo_file = "popinfo.txt"

   # Read the VCF data into a data dictionary using the population info
   data_dict = dadi.Misc.make_data_dict_vcf(
       vcf_file,
       popinfo_file
   )

   # Convert the data dictionary into a Spectrum object. A Spectrum object is
   # a wrapper around a NumPy array and is what dadi uses to work with site
   # frequency spectra
   sfs = dadi.Spectrum.from_data_dict(
       data_dict=data_dict,
       pop_ids=["pop1", "pop2"]  # names for populations in info file
       projections=(20, 40)  # number of chromosomes sampled from each pop
       polarized=False  # Can we determine ancestral vs. derived allelic states
   )

   # Write the SFS to file
   sfs.to_file("spectrum.fs")
