\name{plot.ITPlm}
\alias{plot.ITPlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ITP results for functional-on-scalar linear model testing
}
\description{
\code{plot} method for class "\code{ITPlm}".
Plotting function creating a graphical output of the ITP for the test on a functional-on-scalar linear model: functional data, functional coefficients and ITP-adjusted p-values for the F-test and t-tests are plotted.
}
\usage{
\method{plot}{ITPlm}(x, xrange = c(0, 1), alpha1 = 0.05, alpha2 = 0.01, 
     plot.adjpval = FALSE, col = c(1, rainbow(dim(x$corrected.pval.t)[1])), 
     ylim = range(x$data.eval), ylab = "Functional Data", 
     main = NULL, lwd = 1, pch = 16, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The object to be plotted.
An object of class "\code{ITPlm}", usually, a result of a call to \code{\link{ITPlmbspline}}.
}
  \item{xrange}{
Range of the \code{x} axis.
}
  \item{alpha1}{
First level of significance used to select and display significant effects. Default is \code{alpha1 = 0.05}.
}
  \item{alpha2}{
Second level of significance used to select and display significant effects. Default is \code{alpha1 = 0.01}. \code{alpha1} and \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are switched.
}
  \item{plot.adjpval}{
A logical indicating wether the plots of adjusted p-values have to be done. Default is \code{plot.adjpval = FALSE}.
}
  \item{col}{
Vector of colors for the plot of functional data (first element), and functional coefficients (following elements). 

Default is \code{col = c(1, rainbow(dim(x$corrected.pval.t)[1]))}.
}
  \item{ylim}{
Range of the \code{y} axis. Default is \code{ylim = range(x$data.eval)}.
}
  \item{ylab}{
Label of \code{y} axis of the plot of functional data. Default is "\code{Functional Data}".
}
  \item{main}{
An overall title for the plots (it will be pasted to "\code{Functional Data and F-test}" for the first plot and "\code{t-test}" for the other plots). 
}
  \item{lwd}{
Line width for the plot of functional data. Default is \code{lwd=16}.
}
  \item{pch}{
Point character for the plot of adjusted p-values. Default is \code{pch=16}.
}
  \item{\dots}{
Additional plotting arguments that can be used with function \code{plot}, such as \code{\link{graphical parameters}} (see \code{\link{par}}).
}
}
\value{
No value returned. 
The function produces a graphical output of the ITP results:  the plot of the functional data, functional regression coefficients, and ITP-adjusted p-values for the F-test and t-tests. 
The basis components selected as significant by the tests at level \code{alpha1} and \code{alpha2} are highlighted in the plot of the corrected p-values and in the one of functional data by gray areas (light and dark gray, respectively). 
The plot of functional data reports the gray areas corresponding to a significant F-test. The plots of functional regression coefficients report the gray areas corresponding to significant t-tests for the corresponding covariate.
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
See also \code{\link{ITPlmbspline}} to fit and test a functional-on-scalar linear model applying the ITP, and \code{\link{summary.ITPlm}} for summaries. 
See \code{\link{plot.ITPaov}}, \code{\link{plot.ITP1}}, and \code{\link{plot.ITP2}} for the plot method applied to the ITP results of functional analysis of variance, one-population and two-population, respectively.

}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

data <- rbind(NASAtemp$milan,NASAtemp$paris)
lab <- c(rep(0,22),rep(1,22))

# Performing the ITP
ITP.result <- ITPlmbspline(data ~ lab,B=1000,nknots=20,order=3)
# Summary of the ITP results
summary(ITP.result)

# Plot of the ITP results
layout(1)
plot(ITP.result,main='NASA data',xlab='Day',xrange=c(1,365))

# Plots of the adjusted p-values
plot(ITP.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))

# To have all plots in one device
layout(matrix(1:6,nrow=3,byrow=FALSE))
plot(ITP.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
