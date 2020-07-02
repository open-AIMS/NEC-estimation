<!-- README.md is generated from README.Rmd. Please edit that file -->

jagsNEC
=======

‘jagsNEC’ is an R package to fit concentration(dose) - response curves
to toxicity data, and derive No-Effect-Concentration (NEC),
No-Significant-Effect-Concentration (NSEC), and Effect-Concentration (of
specified percentage ‘x’, ECx) thresholds from non-linear models fitted
using Bayesian MCMC fitting methods via the R2jags package and jags.

Background
==========

Bayesian model fitting can be difficult to automate across a broad range
of usage cases, particularly with respect to specifying valid initial
values and appropriately vague priors. This is one reason the use of
Bayesian statistics for NEC estimation is not currently widely adopted
across the broader ecotoxicological community, who rarely have access to
specialist statistical expertise. The jagsNEC package attempts to
provide an accessible interface to the R2jags package specifically for
fitting NEC models, with a range of models specified based on the known
distribution of the “concentration” or “dose” variable (the predictor,
x) as well as the “response” (y) variable. The model formula, including
priors and the required init function required to call a jags model are
automatically generated based on information contained in the supplied
data. While the distribution of the x and y variables can be specified
directly, jagsNEC will automatically ‘guess’ the correct distribution to
use, based on the characteristics of the provided data. See
?write.jags.NECmod for the currently available x and y data types.

This project started with an implementation of the NEC model based on
that described in \[@Fox2010\]. That package has been further
generalised to allow a large range of response variables to be modelled
using the appropriate statistical distribution, and the current
implementation supports gaussian, poisson, binomial, gamma, negbin and
beta response data. We have further added a range of alternative NEC
model types, as well as a range of typically used concentration-response
models (such as 4-parameter logistic and weibull models) that have no
NEC ‘step’ function but simply model response as a smooth function of
concentration.

Models can be fit directly using fit.jagsNEC, or alternatively using the
function fit.jagsMANEC it is possible to fit a specific set or all of
the available models. The fit.jagsMANEC function returns a model
weighted estimate of predicted posterior values, based on DIC model
weights. It is also possible to obtain all individual model fits from
the fitted jagsMANECfit model object if required.

An additional endpoint has also been derived using Bayesian posterior
predicted values to estimate the “No-Statistical-Effect-Concentration”
as the concentration at which predicted values for each MCMC chain fall
below a lower percentile bound (defined as sig.val) of the control
(assumed to be the lowest treatment (x.var) concentration in the data.
NSEC estimates are currently used to approximate NEC for models without
a specific NEC step parameter (in jagsNEC these have the prefix ECx in
their model name).

Important information on the current package is contained in the jagsNEC
helpfile (see ?jagsNEC).

This package is currently under development. We are keen on any feedback
regarding usage, and especially bug reporting that includes an easy self
contained reproducible example of either unexpected behaviour or example
model fits that fail to converge (have poor chain mixing) or yield other
errors. Such information will hopefully help us towards building a more
robust package. We cannot help troublshoot issues if an easy to run
reproducible example is not supplied.

Installation
============

To install the latest version from github
(<a href="https://github.com/AIMS/NEC-estimation" class="uri">https://github.com/AIMS/NEC-estimation</a>)

    install.packages("remotes")
    remotes::install_github("AIMS/NEC-estimation")
