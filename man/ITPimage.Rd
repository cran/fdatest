\name{ITPimage}
\alias{ITPimage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the Interval Testing Procedure results
}
\description{
Plotting function creating a graphical output of the ITP: the p-value heat-map, the plot of the corrected p-values, and the plot of the functional data.
}
\usage{
ITPimage(ITP.result, alpha = 0.05, abscissa.range = c(0, 1), nlevel = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ITP.result}{
Results of the ITP, as created by \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, and \code{\link{ITP2pafourier}}.
}
  \item{alpha}{
Level of the hypothesis test. The default is \code{alpha}=0.05.
}
  \item{abscissa.range}{
Range of the plot abscissa. The default is \code{c(0,1)}.
}
  \item{nlevel}{
Number of desired color levels for the p-value heatmap. The default is \code{nlevel=20}.
}
}

\value{
No value returned. 
The function produces a graphical output of the ITP results: the p-value heatmap, a plot of the corrected p-values and the plot of the functional data. 
The basis components selected as significant by the test at level \code{alpha} are highlighted in the plot of the corrected p-values by a gray area. 
}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

}
\author{
Alessia pini, Simone Vantini
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{plot.ITP1}}, \code{\link{plot.ITP2}}, \code{\link{plot.ITPlm}}, and \code{\link{plot.ITPaov}} for the plot method applied to the ITP results of one- and two-population tests, linear models, and ANOVA, respectively.

See also \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, and \code{\link{ITP2pafourier}} for applying the ITP.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

# Performing the ITP for two populations with the B-spline basis
ITP.result <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=20,B=1000)

# Plotting the results of the ITP
ITPimage(ITP.result,abscissa.range=c(0,12))

# Selecting the significant components for the radius at 5% level
which(ITP.result$corrected.pval < 0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
