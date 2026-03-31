# TrioBEASTIE

TrioBEASTIE is a Bayesian graphical model for detecting allele-specific activity in a familial trio.
Please cite the publication: https://www.biorxiv.org/content/10.64898/2026.03.28.714974v1

## Description

TrioBEASTIE is a Stan model. Code is provided to run the model using Rstan through python with rpy2. 
TrioBEASTIE can be applied to read counts from RNA-seq to detect allele-specific gene expression (ASE) or
to ATAC-seq to detect allele-specific chromatin accessibility (ASA). The model takes in data in the essex format. 
In addition to read count data, genotypes are also required. These do not need to be phased, phasing is done using trio information with code developed in the lab.

## Getting Started

### Dependencies

* Originally run with Python version 3.13.2
* requires python packages numpy, scipy, and rpy2 (which requires htslib, available from anaconda)
* requires lab developed python libraries Rex and EssexParser available from https://github.com/bmajoros/python
* to phase essex files, use lab developed phaser: https://github.com/bmajoros/TrioBEAST/blob/main/phase-trio.C

### Installing

* the only necessary installation is of dependencies; stan and python files for running the model just need to be downloaded.
* after installing python dependencies, use python with rpy2 to install R packages:
```
python
> import rpy2
> from rpy2.robjects.packages import importr
> utils = importr('utils')
> utils.install_packages('rstan')
> utils.install_packages('codetools')
```

### Executing program

* First phase the essex file
```
./phase-trio input.essex input.phased.essex
```
* Then run the model
```
AFFECTED=0.04, NUM_MCMC=5000, NUM_GENES=4999; ./refactored_11_mode_model.py Refactored_constant_singleprior_no_triplehets input.phased.essex $NUM_MCMC 0-$NUM_GENES $AFFECTED trio_beastie.out
```

## Help

Contact stephanie.hoyt@duke.edu or bmajoros@duke.edu with any issues.

## Authors

Stephanie H. Hoyt: stephanie.hoyt@duke.edu \
William H. Majoros: bmajoros@duke.edu

## Version History

* 0.1
    * Initial Release March 27, 2026

## License

This is OPEN SOURCE SOFTWARE governed by the Gnu General Public License (GPL) version 3, as described at www.opensource.org.
Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu) and Stephanie Hoyt (stephanie.hoyt@duke.edu)

## Acknowledgments

Additional project advisors and funding:
* Andrew S. Allen
* Raluca Gordan
* Tim E. Reddy
* Research reported in this publication was supported in part by the National Institute of General Medical Sciences of the National Institutes of Health under award number 1R35-GM150404 to W.H.M., and by NIH under award number RM1-HG011123 to T.E.R. and A.S.A. and R.G. Content is solely the responsibility of the authors.
