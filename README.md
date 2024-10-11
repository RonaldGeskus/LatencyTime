# Repository for Latency Time paper [![DOI](https://zenodo.org/badge/865307690.svg)](https://doi.org/10.5281/zenodo.13912402)

Code and data to obtain the estimates as reported in the paper "The latency time of SARS-CoV-2 Delta
variant in infection- and vaccine-naive individuals from Vietnam". 

The file Descriptives.R contains the code to obtain the growth rate of the 2021 Ho Chi Minh city
SARS-CoV-2 outbreak. It uses the data set incidence.csv

The file Estimation.R contains the code to obtain the estimate of the latency time distribution. It
uses the data set intervals.csv. The JAGS code is provided for both the uniform infection risk and
the exponential increase in infection risk. The code assumes a generalized gamma distribution, but
that can be changed to Weibull or gamma (similar code can be found in the function Estimate_doublIn
in the doublIn package on CRAN).

