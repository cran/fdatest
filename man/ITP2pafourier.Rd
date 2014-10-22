\name{ITP2pafourier}
\alias{ITP2pafourier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two populations Interval Testing Procedure with Fourier basis (phase-amplitude decomposition)
}
\description{
The function implements the Interval Testing Procedure for testing the difference between two functional populations evaluated on a uniform grid. Data are represented by means of the Fourier basis expansion with the phase-amplitude decomposition and the significance of the amplitude and phase of each frequency is tested with an interval-wise control of the Family Wise Error Rate.
}
\usage{
ITP2pafourier(data1, data2, 
              maxfrequency=floor(dim(data1)[2]/2), B = 10000, paired = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
Pointwise evaluations of the first population's functional data set on a uniform grid. \code{data1} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
}
  \item{data2}{
Pointwise evaluations of the second population's functional data set on a uniform grid. \code{data2} is a matrix of dimensions \code{c(n2,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
}
  \item{maxfrequency}{
The maximum frequency to be used in the Fourier basis expansion of data. The default is \code{floor(dim(data1)[2]/2)}, leading to an interpolating expansion.
}
  \item{B}{
The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=10000}.
}
  \item{paired}{
A logical indicating whether the test is paired. The default is \code{FALSE}.
}
}

\value{
\code{ITP2pafourier} returns an object of \code{\link{class}} "\code{ITP2}".

An object of class "\code{ITP2}" is a list containing at least the following components:

\item{basis}{String vector indicating the basis used for the first phase of the algorithm. Equal to \code{"paFourier"}.}
\item{test}{String vector indicating the type of test performed. Equal to \code{"2pop"}.}
\item{paired}{Logical indicating whether the test is paired (as entered by the user).}
\item{coeff_phase}{Matrix of dimensions \code{c(n,p)} of the \code{p} phases of the Fourier basis expansion. Rows are associated to units and columns to frequencies: the first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients of the second population units.}
\item{coeff_amplitude}{Matrix of dimensions \code{c(n,p)} of the \code{p} amplitudes of the Fourier basis expansion. Rows are associated to units and columns to frequencies: the first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients of the second population units.}
\item{pval_phase}{Uncorrected p-values of the phase tests for each frequency.}
\item{pval_amplitude}{Uncorrected p-values of the amplitude tests for each frequency.}
\item{pval.matrix_phase}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests on phase. The element \code{(i,j)} of matrix \code{pval.matrix_phase} contains the p-value of the joint NPC test of the frequencies \code{(j,j+1,...,j+(p-i))}.}
\item{pval.matrix_amplitude}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests on amplitude. The element \code{(i,j)} of matrix \code{pval.matrix_amplitude} contains the p-value of the joint NPC test of the frequencies \code{(j,j+1,...,j+(p-i))}.}
\item{corrected.pval_phase}{Corrected p-values of the phase tests for each frequency.}
\item{corrected.pval_amplitude}{Corrected p-values of the amplitude tests for each frequency.}
\item{labels}{Labels indicating the population membership of each data.}
\item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
\item{heatmap.matrix_phase}{Heatmap matrix of p-values for phase (used only for plots).}
\item{heatmap.matrix_amplitude}{Heatmap matrix of p-values for amplitude (used only for plots).}
}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

}
\author{
Alessia Pini, Simone Vantini
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{ITP2fourier}},  \code{\link{ITP2bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP1bspline}}, and \code{\link{ITPimage}}.
}
\examples{

# Importing the NASA temperatures data set
data(NASAtemp)
# Performing the ITP
ITP.result <- ITP2pafourier(NASAtemp$milan,NASAtemp$paris,maxfrequency=20,B=1000,paired=TRUE)
# Plotting the results of the ITP
#\dontrun{

plot(ITP.result,main='NASA data',xrange=c(1,365),xlab='Day')

#}
# Selecting the significant coefficients
which(ITP.result$corrected.pval < 0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
