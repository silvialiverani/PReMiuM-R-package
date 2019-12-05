# Welcome to the documentation for PReMiuM!

This repository is home to _.md_ files rendered to html at [PReMiuM](https://premium-r-package.readthedocs.io/en/latest/).

See [Updating the Docs](/updating/) for instructions to modify content in this repository.


### EnviroTyping

Dirichlet Process Bayesian Clustering, Profile Regression

### Core Contributors

- Silvia Liverani : Queen Mary University of London, School of Mathematical Sciences

### Purpose

Bayesian clustering using a Dirichlet process mixture model. This model is an alternative to regression models, non-parametrically linking a response vector to covariate data through cluster membership. The package allows Bernoulli, Binomial, Poisson, Normal, survival and categorical response, as well as Normal and discrete covariates. It also allows for fixed effects in the response model, where a spatial CAR (conditional autoregressive) term can be also included. Additionally, predictions may be made for the response, and missing values for the covariates are handled. Several samplers and label switching moves are implemented along with diagnostic tools to assess convergence. A number of R functions for post-processing of the output are also provided. In addition to fitting mixtures, it may additionally be of interest to determine which covariates actively drive the mixture components. This is implemented in the package as variable selection. 

### Reference

The main reference for the package is Liverani, Hastie, Azizi, Papathomas and Richardson (2015) <doi:10.18637/jss.v064.i07>.

