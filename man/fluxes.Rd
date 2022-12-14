% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fluxes}
\alias{fluxes}
\title{Calculated fluxes}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 1180 rows and 11 columns.
}
\usage{
fluxes
}
\description{
A data set containing metabolite fluxes.
}
\details{
\describe{
\item{metabolite}{name of measurement}
\item{abbreviation}{abbreviated metabolite name}
\item{cell_type}{
\code{lf} = lung fibroblasts \cr
\code{pasmc} = pulmonary artery smooth muscle cells}
\item{experiment}{
\code{02} = 0.2\% oxygen for hypoxia \cr
\code{05} = 0.5\% oxygen for hypoxia \cr
\code{bay} = molidustat treatment \cr
\code{05-bay} = 0.5\% oxygen plus molidustat}
\item{batch}{a group of biological replicates analyzed similarly}
\item{date}{start date of an experiment}
\item{oxygen}{ambient oxygen level}
\item{treatment}{
\code{None} = no treatment \cr
\code{DMSO} = 0.1\% DMSO \cr
\code{BAY} = 10 μM molidustat}
\item{flux}{fmol / cell / h, positive fluxes indicate secretion, negative fluxes
indicate uptake}
}
}
\keyword{datasets}
