[![Documentation Status](https://readthedocs.org/projects/polyploid-demography/badge/?version=latest)](https://polyploid-demography.readthedocs.io/en/latest/?badge=latest)

# Demographic Inference in Polyploids

## [Read the Docs](http://polyploid-demography.rtfd.io/)

This repository hosts all the code for testing models for demographic inference
in polyploid species. The models are implemented in the open-source Python package
[dadi](https://bitbucket.org/gutenkunstlab/dadi).

Each folder corresponds to a different analysis from the accompanying manuscript
(<a href="" target="_blank">link</a>). More detailed documentation can be found
on our **Read the Docs** website (link above), as well as within the comments of
the scripts in this repository.

> **Note**
> 
> While working on the manuscript we changed the homoeologous exchange parameter
> from $d_{i \leftrightarrow j}$ to $e_{i \leftrightarrow j}$. Because, of this
> some of the code may still have `dij` as a parameter but it represents the
> same thing.

**Preprint information**

Blischak, PD, M Sajan, MS Barker, RN Gutenkunst. 2022. Demographic history
inference and the polyploid continuum. *bioRxiv* doi:
https://doi.org/10.1101/2022.09.15.508148.
