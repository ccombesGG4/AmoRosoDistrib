# AmoRosoDistrib
The R package "AmoRosoDistrib" provides functions for fitting four-parameter Generalized Gamma univariate distribution to different types of data (continuous non-censored  non-truncated data and discrete data) from nine different approaches based on the maximum likelihood estimation method and the minimum distance estimation (MDE) method.
Package: AmoRoDistrib

R package for fitting Amoroso probability distribution family

-----------------
See the article: 
-----------------
"On Parameter Estimation for Amoroso Family of Distributions". Catherine Combes and Hon Keung Tony Ng.
Volume 191, January 2022, Pages 309-327, in Mathematics and Computers in Simulation - Journal - Elsevier. https://doi.org/10.1016/j.matcom.2021.07.004

-------------------
Citation to BibTeX: 
-------------------
@article{COMBES2022309,
title = {On parameter estimation for Amoroso family of distributions},
journal = {Mathematics and Computers in Simulation},
volume = {191},
pages = {309-327},
year = {2022},
issn = {0378-4754},
doi = {https://doi.org/10.1016/j.matcom.2021.07.004},
url = {https://www.sciencedirect.com/science/article/pii/S0378475421002548},
author = {Catherine Combes and Hon Keung Tony Ng},
keywords = {Maximum likelihood estimation, Minimum distance estimation, Optimization, Generalized gamma distribution},
abstract = {The four-parameter generalized gamma (GΓ) distribution, also known as the Amoroso family of distributions, is a flexible and versatile statistical distribution that encapsulates many well-known lifetime distributions, including the exponential, Weibull, lognormal, and gamma distributions as special instances. The four-parameter GΓ distribution is shown to be appropriate for fitting skewed and heavy-tailed data sets. However, even though the GΓ distribution is very useful and flexible, it remains less studied than its counterparts, probably due to the difficulty in estimating the parameters of the distribution. In this paper, we explore several novel iterative parameter estimation approaches for the four-parameter GΓ distribution, which includes the maximum likelihood estimation and minimum distance estimation approaches. Standard error and confidence interval of a function of the parameter estimates based on bootstrap method are also discussed. An R package is developed based on the proposed estimation methods. Numerical examples and Monte Carlo simulations are used to illustrate the usefulness of the proposed approaches for fitting the four-parameter GΓ distribution.}
}

------------
Installation
------------

On RStudio, install package devtools and curl:

    install.packages("devtools") 
    install.packages("curl")

Afterwards, run these commands:

    library(devtools)
    install_github("ccombesGG4/AmoRosoDistrib")

------------
How to use?
------------
See 'vignettes' directory
