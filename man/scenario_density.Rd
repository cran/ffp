% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scenario_density.R, R/scenario_histogram.R
\name{scenario_density}
\alias{scenario_density}
\alias{scenario_histogram}
\title{Plot Scenarios}
\usage{
scenario_density(x, p, n = 10000)

scenario_histogram(x, p, n = 10000)
}
\arguments{
\item{x}{An univariate marginal distribution.}

\item{p}{A probability from the \code{ffp} class.}

\item{n}{An \code{integer} scalar with the number of scenarios to be generated.}
}
\value{
A \code{ggplot2} object.
}
\description{
This functions are designed to make it easier to visualize the impact of a
\emph{View} in the P&L distribution.
}
\details{
To generate a scenario-distribution the margins are bootstrapped using
\code{\link{bootstrap_scenarios}}. The number of resamples can be controlled
with the \code{n} argument (default is \code{n = 10000}).
}
\examples{
x <- diff(log(EuStockMarkets))[, 1]
p <- exp_decay(x, 0.005)

scenario_density(x, p, 500)

scenario_histogram(x, p, 500)
}
