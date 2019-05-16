
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--05--16-brightgreen.svg)](https://github.com/pat-s/pathogen-modeling/commits/master)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)

## Reproducible analysis for paper:

##### tailfindr: Alignment-free poly(A) length measurement for Oxford Nanopore RNA and DNA sequencing

Maximilian Krause, Adnan M. Niazi, Kornel Labun, Yamila N. Torres
Cleuren, Florian S. Müller, Eivind Valen

## About the repo

This repository is organized as an R package, providing functions and
the raw data to reproduce and extend the analysis reported in the
publication. By raw data, we mean the output of tools such as tailfindr
and Nanopolish etc.

This project is setup with a drake workflow, ensuring reproducibility.
Intermediate targets/objects will be stored in a hidden .drake
directory.

The R library of this project is managed by packrat. This makes sure
that the exact same package versions are used when recreating the
project. When calling packrat::restore(), all required packages will be
installed with their specific version.

Please note that this project was built with R version 3.6.0 on a MAC
OSx Mojave operating system. The packrat packages from this project are
not compatible with R versions prior version 3.6.0 (In general, it
should be possible to reproduce the analysis on any other operating
system.)

## Getting started

To clone the project, a working installation of git is required. Open a
terminal in the directory of your choice and execute:

``` sh
git clone https://github.com/adnaniazi/krauseNiazi2019Analyses.git
```

Then go into the `krauseNiazi2019Analyses` directory using:

``` sh
cd krauseNiazi2019Analyses
```

Now start R in this directory and run:

``` r
packrat::restore() # restores all R packages with their specific version
drake::r_make()    # recreates the analysis
```
