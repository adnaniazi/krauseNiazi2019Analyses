# krauseNiazi2019Analyses

[![Build Status](https://travis-ci.org/adnaniazi/krauseNiazi2019Analyses.png?branch=master)](https://travis-ci.org/adnaniazi/krauseNiazi2019Analyses)  [![codecov](https://codecov.io/gh/adnaniazi/krauseNiazi2019Analyses/branch/master/graph/badge.svg)](https://codecov.io/gh/adnaniazi/krauseNiazi2019Analyses)
[![cran](https://www.r-pkg.org/badges/version-last-release/krauseNiazi2019Analyses)](https://cran.r-project.org/package=krauseNiazi2019Analyses)

## How to finish setting up your new package

Now that you've got a working package skeleton, there are a few steps to finish setting up all the integrations:

### 1. Git(Hub)

Go to https://github.com/adnaniazi and create a new repository. Then, in the directory where this package is, create your git repository from the command line, add the files, and push it to GitHub:

    git init
    git add --all
    git commit -m "Initial commit of package skeleton"
    git remote add origin git@github.com:adnaniazi/krauseNiazi2019Analyses.git
    git push -u origin master

### 2. Travis

Now you can go to [Travis](https://travis-ci.org/profile/adnaniazi) and turn on continuous integration for your new package. You may need to click the "Sync account" button to get your new package to show up in the list.

If you have a codecov.io account, running your tests on Travis will trigger the code coverage job. No additional configuration is necessary

### 3. Appveyor

Go to [Appveyor's new project page](https://ci.appveyor.com/projects/new) and select your new repository from the list. Then you can go to the [badges](https://ci.appveyor.com/project/adnaniazi/krauseNiazi2019Analyses/settings/badges) page, copy the markdown code it provides, and paste it up with the other badges above. (Their badge API has a random token in it, so `skeletor` can't include it in the template for you.)

### 4. Delete this "How to finish setting up your new package" section from your README.md

## Installing

<!-- If you're putting `krauseNiazi2019Analyses` on CRAN, it can be installed with

    install.packages("krauseNiazi2019Analyses") -->

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/r-lib/devtools) package:

    # install.packages("devtools")
    devtools::install_github("adnaniazi/krauseNiazi2019Analyses", build_vignettes=TRUE)

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](http://testthat.r-lib.org/) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `testthat::test_package()` will do a regular-expression pattern match within the file names (ignoring the `test-` prefix and the `.R` file extension). 

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
