\name{VariableSelectionUnknowCompFMR}
\alias{VariableSelectionUnknowCompFMR}
\title{Variable selection in the mixture model with unknown number of components.}
\description{Using RJMCMC procedure to generate posterior samples.}
\usage{VariableSelectionUnknowCompFMR(Num.of.iterations, Num.of.iterations.Inside, Max.num.of.groups, list.Data, list.InitialValues, list.GivenValues, bool.Robust, bool.RJMCMC)
}
\arguments{
\item{Num.of.iterations}{Number of iterations.}
\item{Num.of.iterations.Inside}{Num of iterations to run when the demension changes.}
\item{Max.num.of.groups}{The maximum number of groups.}
\item{list.Data}{\emph{List} of data.}
\item{list.InitialValues}{\emph{List} of initial values for parameters.}
\item{list.GivenValues}{\emph{List} of given hyperparameters in priors.}
\item{bool.Robust}{\emph{List} of given hyperparameters in priors.}
\item{bool.RJMCMC}{Determine which parameter will be updated.}
}
\value{Lists of posterior samples.}
\examples{
\dontrun{
library(MASS)
library(UnknownCompFMR)
rm(list=ls(all=TRUE))
set.seed(1)
}
}