\name{plot.ITPaov}
\alias{plot.ITPaov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ITP results for functional analysis of variance testing
}
\description{
\code{plot} method for class "\code{ITPaov}".
Plotting function creating a graphical output of the ITP for the test on a functional analysis of variance: functional data, and ITP-adjusted p-values of the F-tests on the whole model and on each factor are plotted.
}
\usage{
\method{plot}{ITPaov}(x,xrange=c(0,1), alpha1=0.05, alpha2=0.01,
     plot.adjpval=FALSE,ylim=range(x$data.eval),col=1,
     ylab='Functional Data',main=NULL,lwd=1,pch=16,...)
}
\arguments{
  \item{x}{
  The object to be plotted.
An object of class "\code{ITPaov}", usually, a result of a call to \code{\link{ITPaovbspline}}.
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
Colors for the plot of functional data. Default is \code{col = 1}.
}
  \item{ylim}{
Range of the \code{y} axis. Default is \code{ylim = range(x$data.eval)}.
}
  \item{ylab}{
Label of \code{y} axis of the plot of functional data. Default is "\code{Functional Data}".
}
  \item{main}{
An overall title for the plots (it will be pasted to "\code{Functional Data and F-test}" for the first plot and "\code{factor}" for the other plots). 
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
The function produces a graphical output of the ITP results:  the plot of the functional data, functional regression coefficients, and ITP-adjusted p-values of the F-tests on the whole model and on each factor. 
The basis components selected as significant by the tests at level \code{alpha1} and \code{alpha2} are highlighted in the plot of the corrected p-values and in the one of functional data by gray areas (light and dark gray, respectively). 
The first plot reports the gray areas corresponding to a significant F-test on the whole model. The remaining plots report the gray areas corresponding to significant F-tests on each factor (with colors corresponding to the levels of the factor).
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
See also \code{\link{ITPaovbspline}} to fit and test a functional analysis of variance applying the ITP, and \code{\link{summary.ITPaov}} for summaries. 
See \code{\link{plot.ITPlm}}, \code{\link{plot.ITP1}}, and \code{\link{plot.ITP2}} for the plot method applied to the ITP results of functional-on-scalar linear models, one-population and two-population, respectively.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the ITP
ITP.result <- ITPaovbspline(temperature ~ groups,B=1000,nknots=20,order=3)

# Summary of the ITP results
summary(ITP.result)

# Plot of the ITP results
layout(1)
plot(ITP.result)

# All graphics on the same device
layout(matrix(1:4,nrow=2,byrow=FALSE))
plot(ITP.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
