% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predict_fitted_vals.R
\name{predict_WB2mod}
\alias{predict_WB2mod}
\title{predict_WB2mod}
\usage{
predict_WB2mod(x.vec, EC50, top, beta, bot = 0)
}
\arguments{
\item{x.vec}{the x vector over which to calculate}

\item{EC50}{the 50 percent effect concentration}

\item{top}{the upper plateau}

\item{beta}{the exponential decay rate (hillslope)}

\item{bot}{the lower plateau}
}
\value{
A list containing x and fitted y, with up and lw values
}
\description{
Calculates predicted y (response) values for a supplied vector of x (concentration) values
}
