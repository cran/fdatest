\name{plot.ITP2}
\alias{plot.ITP2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ITP results for two-population tests
}
\description{
\code{plot} method for class "\code{ITP2}".
Plotting function creating a graphical output of the ITP for the test of comparison between two populations: functional data and ITP-adjusted p-values are plotted.
}
\usage{
\method{plot}{ITP2}(x, xrange = c(0, 1), alpha1 = 0.05, alpha2 = 0.01, 
     ylab = "Functional Data", main = NULL, lwd = 1, 
     col = c(1, 2), pch = 16, ylim = range(object$data.eval), ...)
}
\arguments{
  \item{x}{
  The object to be plotted.
An object of class "\code{ITP2}", that is, a result of an ITP for comparison between two populations. Usually a call to \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}} or \code{\link{ITP2pafourier}}.
}
  \item{xrange}{
Range of the \code{x} axis.
}
  \item{alpha1}{
First level of significance used to select and display significant differences. Default is \code{alpha1 = 0.05}.
}
  \item{alpha2}{
Second level of significance used to select and display significant differences. Default is \code{alpha1 = 0.01}. \code{alpha1} and \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are switched.
}
  \item{ylab}{
Label of \code{y} axis of the plot of functional data. Default is "\code{Functional Data}".
}
  \item{main}{
An overall title for the plots (it will be pasted to "\code{Functional Data}" for the first plot and "\code{adjusted p-values}" for the second plot). 
}
  \item{lwd}{
Line width for the plot of functional data. 
}
  \item{col}{
Color used to plot the functional data.
}
  \item{pch}{
Point character for the plot of adjusted p-values.
}
  \item{ylim}{
Range of the \code{y} axis.
}
  \item{\dots}{
Additional plotting arguments that can be used with function \code{plot}, such as \code{\link{graphical parameters}} (see \code{\link{par}}).
}
}
\value{
No value returned. 
The function produces a graphical output of the ITP results:  the plot of the functional data and the one of the adjusted p-values. 
The basis components selected as significant by the test at level \code{alpha1} and \code{alpha2} are highlighted in the plot of the corrected p-values and in the one of functional data (in case the test is based on a local basis, such as B-splines) by gray areas (light and dark gray, respectively). 
In the case of a Fourier basis with amplitude and phase decomposition, two plots of adjusted p-values are done, one for phase and one for amplitude.
}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

}
\author{
Alessia Pini, Simone Vantini
}
\seealso{
\code{\link{ITPimage}} for the plot of p-values heatmaps. 

See also \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}} to perform the ITP to test for differences between two populations. 
See \code{\link{plot.ITP1}} and \code{\link{plot.ITPlm}} for the plot method applied to the ITP results of one-population tests and a linear models, respectively.
}
\examples{

# Importing the NASA temperatures data set
data(NASAtemp)

# Performing the ITP for two populations with the B-spline basis
ITP.result.bspline <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=30,B=1000)
# Plotting the results of the ITP
plot(ITP.result.bspline,xlab='Day',xrange=c(1,365),main='NASA data')

# Selecting the significant components for the radius at 5% level
which(ITP.result.bspline$corrected.pval < 0.05)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
