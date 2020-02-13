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
    vcf_dicts = [dadi.Misc.make_data_dict_vcf(f, "arenosa_popinfo.txt") for f in glob("*.vcf")]
    spectra   = [dadi.Spectrum.from_data_dict(dd, ['arenosa'], projections=[20]) for dd in vcf_dicts]
    fs        = specta[0]
    for s in range(1,len(spectra)):
        fs += spectra[s]
