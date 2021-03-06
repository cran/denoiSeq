% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Class_denoiSeq.R
\docType{methods}
\name{setOutput}
\alias{setOutput}
\alias{setOutput,readsData-method}
\title{Generic  for the altering output slot.}
\usage{
setOutput(object, outval)

\S4method{setOutput}{readsData}(object, outval)
}
\arguments{
\item{object}{a readsData object.}

\item{outval}{A list of the output from Bayesian inference.}
}
\value{
The same readsData object with the output slot updated.
}
\description{
Updates the value of the output slot for the readsData object supplied.
}
\section{Methods (by class)}{
\itemize{
\item \code{readsData}: Alters the value of the output slot of
the readsData object.
}}

\keyword{internal}
