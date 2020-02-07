.. _ArenosaDemog:

Demographic Analyses in *Arabidopsis arenosa*
=============================================

.. code-block:: python
  :linenos:

  #!/usr/bin/env python3

  """ Code for file run_arenosa_demog.py """
  import dadi
  import numpy as np
  from glob import glob

  if __name__ == "__main__":
    vcf_dicts = [dadi.Misc.make_data_dict_vcf(f) for f in glob("*.vcf")]
    spectra   = [dadi.Spectrum.from_data_dict(dd) for dd in vcf_dicts]
    fs        = spectra[0]
    for s in range(1,len(spectra)):
        fs += spectra[s]

**References**

Monnahan, P. *et al*. 2019. Pervasive population genomic consequences of genome duplication in *Arabidopsis arenosa*.
*Nature Ecology & Evolution* 3:457--468.
