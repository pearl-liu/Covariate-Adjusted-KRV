# Covariate-adjusted kernel RV test

## Overview
The covariate-adjusted kernel RV (KRV) test evaluates the generalized association between two multivariate (potentially high-dimensional) variables, while adjusting for covariates. It is implemented as part of the `KRV()` function in the [MiRKAT R package](https://CRAN.R-project.org/package=MiRKAT), which is a comprehensive R package for association analysis of microbiome data. We introduce the usage of the `KRV()` function and demonstrate its application in genetic association analysis of microbiome composition.

## System Requirements

### Operating systems
The MiRKAT package (v1.1.3) works in R with version 3.5.0 or higher.
The package is successfully tested on the following systems:

Mac OSX: Mojave version 10.14.6 (R version 4.0.2)

Windows 10 version 1809 (R version 3.6.1) 

### Package dependencies
The following packages are required for functions and examples in the MiRKAT package: MASS, GUniFrac, CompQuadForm, quantreg, PearsonDS, propr, lme4, Matrix, permute, survival, and stats. All of these packages are available on [CRAN](https://cran.r-project.org/).

## Installation guide
The MiRKAT package can be downloaded and installed from CRAN.

```{r message=FALSE, warning=FALSE, eval=FALSE}
install.packages("MiRKAT")
```
The install time is approximately 2 seconds on a laptop (2.7 GHz Intel Core i7).

We can then load the package:
```{r message=FALSE, warning=FALSE}
library(MiRKAT)
```

##  Demo and instructions for use
See the [vignette](https://pearl-liu.github.io/Covariate-Adjusted-KRV/vignette.html).
