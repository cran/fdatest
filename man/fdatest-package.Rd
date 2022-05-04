\name{fdatest-package}
\alias{fdatest-package}
\alias{fdatest}
\docType{package}
\title{
Interval Testing Procedure for Functional Data
}
\description{
Implementation of the Interval Testing Procedure for functional data in different frameworks (i.e., one or two-population frameworks, functional linear models) by means of different basis expansions (i.e., B-spline, Fourier, and phase-amplitude Fourier). The current version of the package requires functional data evaluated on a uniform grid; it automatically projects each function on a chosen functional basis; it performs the entire family of multivariate tests; and, finally, it provides the matrix of the p-values of the previous tests and the vector of the corrected p-values. The functional basis, the coupled or uncoupled scenario, and the kind of test can be chosen by the user. The package provides also a plotting function creating a graphical output of the procedure: the p-value heat-map, the plot of the corrected p-values, and the plot of the functional data.
}
\details{
\tabular{ll}{
Package: \tab fdatest\cr
Type: \tab Package\cr
Version: \tab 2.1.1\cr
Date: \tab 2022-05-04\cr
License: \tab GPL-2\cr
}

}
\author{
Alessia Pini, Simone Vantini

Maintainer: Alessia Pini <alessia.pini@polimi.it>
}
\encoding{utf-8}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

K. Abramowicz, S. De Luna, C. Häger, A. Pini, L. Schelin, and S. Vantini (2015). Distribution-Free Interval-Wise Inference for Functional-on-Scalar Linear Models. \emph{MOX-report 3/2015}, Politecnico di Milano.
}
\keyword{ package }
\seealso{
See also \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}}, \code{\link{ITPlmbspline}}, \code{\link{ITPaovbspline}}, and \code{\link{ITPimage}}.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

# Example 1: 
# Performing the ITP for one population with the Fourier basis
ITP.result <- ITP1fourier(NASAtemp$milan,maxfrequency=10,B=1000)
# Plotting the results of the ITP
plot(ITP.result)

# Plotting the p-value heatmap
ITPimage(ITP.result)

# Selecting the significant coefficients
which(ITP.result$corrected.pval < 0.05)

# Example 2: 
# Performing the ITP for two populations with the B-spline basis
ITP.result <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=20,B=1000)
# Plotting the results of the ITP
plot(ITP.result)

# Plotting the p-values heatmap
ITPimage(ITP.result)

# Selecting the significant components for the radius at 5% level
which(ITP.result$corrected.pval < 0.05)


# Example 3:
# Fitting and testing a functional-on-scalar linear model  
# Defining data and covariates 
temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the ITP
ITP.result <- ITPlmbspline(temperature ~ groups,B=1000,nknots=20,order=3)
# Summary of the ITP results
summary(ITP.result)

# Plot of the ITP results
layout(1)
plot(ITP.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))

# All graphics on the same device
layout(matrix(1:6,nrow=3,byrow=FALSE))
plot(ITP.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))


}
