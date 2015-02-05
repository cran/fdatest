\name{summary.ITPlm}
\alias{summary.ITPlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing Functional-on-Scalar Linear Model Fits
}
\description{
\code{summary} method for class "\code{ITPlm}".
}
\usage{
\method{summary}{ITPlm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "\code{ITPlm}", usually, a result of a call to \code{\link{ITPlmbspline}}.
}
  \item{\dots}{
Further arguments passed to or from other methods.
}
}
\value{
The function \code{summary.ITPlm} computes and returns a list of summary statistics of the fitted functional-on-scalar linear model given in \code{object}, using the component "\code{call}" from its arguments, plus:

\item{ttest}{A \code{L+1 x 1} matrix with columns for the functional regression coefficients, and corresponding (two-sided) ITP-adjusted minimum p-values of t-tests (i.e., the minimum p-value over all \code{p} basis components used to describe functional data). }
\item{R2}{Range of the functional R-squared.}
\item{ftest}{ITP-adjusted minimum p-value of functional F-test.}
}
\encoding{utf-8}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

K. Abramowicz, S. De Luna, C. Häger, A. Pini, L. Schelin, and S. Vantini (2015). Distribution-Free Interval-Wise Inference for Functional-on-Scalar Linear Models. \emph{MOX-report 3/2015}, Politecnico di Milano.

}
\author{
Alessia Pini, Simone Vantini
}
\seealso{
See \code{\link{ITPlmbspline}} for fitting and testing the functional linear model and \code{\link{plot.ITPlm}} for plots.
See also  \code{\link{ITPaovbspline}}, \code{\link{ITP1bspline}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}}.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the ITP
ITP.result <- ITPlmbspline(temperature ~ groups, B=1000,nknots=20)

# Summary of the ITP results
summary(ITP.result)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
