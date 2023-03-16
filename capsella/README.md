# Analysis of *Capsella bursa-pastoris*

**Main Scripts**

 - `run_*.py`: Scripts used to run demographic inference under the
   `allotetraploid_bottleneck` and `segtetraploid_bottleneck` models.
 - `capsella_*.csv`: Results of the independent optimization runs for each model.
 - `analyze_capsella_results.py`: Script to analyze all of the results for the
   maximum likelihood parameter estimates under each model, as well as code for
   estimating confidence intervals using the Fisher information Matrix and 
   propogation of uncertainty.
 - `douglas_et_al_comparison.py`: Code for running a comparison of the
   `segtetraploid_bottleneck` model with the model used in the original
   Douglas et al. (2015) paper. The comparison is run with both 2D frequency
   spectra and the collapsed version of the SFS we introduce in our paper.

**SFS Data**

The SFS data file is also in the repo and is named
`Capsella_intergene_4fold_corr_4pop_4_DSFS.fs`. It can be read into Python as
a `dadi.Spectrum` object using the code below:

```python
# Read the *Capsella* SFS file into Python as a `dadi.Spectrum` object.

import dadi
fs = dadi.Spectrum.from_file(
    "Capsella_intergene_4fold_corr_4pop_4_DSFS.fs"
)
```
