# TrioBEASTIE

TrioBEASTIE is a Bayesian graphical model for detecting allele-specific activity in a familial trio.
Please cite the publication: https://www.biorxiv.org/content/10.64898/2026.03.28.714974v1

## Description

TrioBEASTIE is a Stan model (https://mc-stan.org/). Code is provided to run the model using Rstan through python with rpy2. 
TrioBEASTIE can be applied to read counts from RNA-seq to detect allele-specific gene expression (ASE) or
to ATAC-seq to detect allele-specific chromatin accessibility (ASA). The model takes in data in the essex format. 
In addition to read count data, genotypes are also required. These do not need to be phased, phasing is done using trio information with code developed in the lab.

## Getting Started

### Dependencies

* Originally run with Python version 3.13.2
* requires python packages numpy, scipy, and rpy2 (which requires htslib, available from anaconda)
    * ```conda create --name rpy2 conda-forge::python numpy scipy rpy2 htslib```
* requires lab developed python libraries Rex and EssexParser available from https://github.com/bmajoros/python
* to phase essex files, use lab developed phaser: https://github.com/bmajoros/TrioBEAST/blob/main/phase-trio.C
   * Examples of essex files (original and phased) are provided in this repo
   * In this example file, GENE0 is simulated to show no ASE in any individual, GENE1 is simulated to show the father and child affected, and GENE2 is simulated to show ASE in the father only, however this is not visible without making use of the triple heterozygous site (with the folded model). 

### Installing

* the only necessary installation is of dependencies; stan and python files for running the model just need to be downloaded.
    * ```conda create --name rpy2 conda-forge::python numpy scipy rpy2 htslib```
    * Set up other python library files from Bill Majoros:
        * ```git clone https://github.com/bmajoros/python.git```
        * ```export PYTHONPATH="path/to/bmajoros/python/"```
* after installing python dependencies, use python with rpy2 to install R packages:
```
python
> import rpy2
> from rpy2.robjects.packages import importr
> utils = importr('utils')
> utils.install_packages('pak')
> pak = importr('pak')
> pak.pak('StanHeaders@2.32.10')
> pak.pak('rstan@2.32.7')
> utils.install_packages('codetools')
```
* to check that all dependencies have been installed correctly: 
    * the expected output of ```./refactored_11_mode_model.py``` with no inputs is: ``` refactored_11_mode_model.py[-c continue] <model> <input.essex> <#MCMC-samples> <firstGene-lastGene> <P(affected)> <P(recomb)> <P(denovo)> <outFile>
  gene range is zero-based and inclusive ```
    * the expected version of stan via ```rstan.stan_version()``` is 2.32.2. Later versions may give errors when compiling stan files

### Executing program

* First copy over the provided essex file so you can edit freely and compare output, and phase the essex file
```
cp example.essex input.essex
./phase-trio input.essex input.phased.essex
```
* Then run the model
```
PROBAFFECTED=0.04; PROBRECOMB=0.01; PROBDENOVO=0.001; NUM_MCMC=5000; NUM_GENES=4999; ./refactored_11_mode_model.py stan_files/TrioBEASTIE input.phased.essex $NUM_MCMC 0-$NUM_GENES $PROBAFFECTED $PROBRECOMB $PROBDENOVO trio_beastie.out
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

This is OPEN SOURCE SOFTWARE governed by the MIT License.
Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu) and Stephanie H. Hoyt (stephanie.hoyt@duke.edu)

## Acknowledgments

Additional project advisors and funding:
* Andrew S. Allen
* Raluca Gordan
* Tim E. Reddy
* Research reported in this publication was supported in part by the National Institute of General Medical Sciences of the National Institutes of Health under award number 1R35-GM150404 to W.H.M., and by NIH under award number RM1-HG011123 to T.E.R. and A.S.A. and R.G. Content is solely the responsibility of the authors.
