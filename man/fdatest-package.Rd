\name{fdatest-package}
\alias{fdatest-package}
\alias{fdatest}
\docType{package}
\title{
Interval Testing Procedure for functional data
}
\description{
The package implements the Interval Testing Procedure for functional data in different frameworks (i.e., one or two-population frameworks) by means of different basis expansions (i.e., B-spline, Fourier, and phase-amplitude Fourier). The current version of the package requires functional data evaluated on a uniform grid; it automatically projects each function on a chosen functional basis; it performs the entire family of multivariate tests; and, finally, it provides the matrix of the p-values of the previous tests and the vector of the corrected p-values. The functional basis, the coupled or uncoupled scenario, and the kind of test can be chosen by the user. The package provides also a plotting function creating a graphical output of the procedure: the p-value heat-map, the plot of the corrected p-values, and the plot of the functional data.
}
\details{
\tabular{ll}{
Package: \tab fdatest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-02-27\cr
License: \tab GPL-2\cr
}

}
\author{
Alessia Pini, Simone Vantini

Maintainer: Alessia Pini <alessia.pini@mail.polimi.it>
}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

}
\keyword{ package }
\seealso{
See also \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}}, and \code{\link{ITPimage}}.
}
\examples{
data(NASAtemp)
# Performing the ITP for one population with the Fourier basis
ITP.result <- ITP1fourier(NASAtemp$milan,maxfrequency=15,B=1000)
# Plotting the results of the ITP
\dontrun{

ITPimage(ITP.result)

}
# Selecting the significant coefficients
which(ITP.result$corrected.pval < 0.05)

# Performing the ITP for two populations with the B-spline basis
ITP.result <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=30,B=1000)
# Plotting the results of the ITP
\dontrun{

ITPimage(ITP.result,abscissa.range=c(0,12))

}
# Selecting the significant components for the radius at 5% level
which(ITP.result$corrected.pval < 0.05)

}
