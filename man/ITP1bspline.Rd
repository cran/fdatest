\name{ITP1bspline}
\alias{ITP1bspline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One population Interval Testing Procedure with B-spline basis
}
\description{
The function implements the Interval Testing Procedure for testing the center of symmetry of a functional population evaluated on a uniform grid. Data are represented by means of the B-spline expansion and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate. The default parameters of the basis expansion lead to the piece-wise interpolating function.
}
\usage{
ITP1bspline(data, mu = 0, order = 2, nknots = dim(data)[2], B = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Pointwise evaluations of the functional data set on a uniform grid. \code{data} is a matrix of dimensions \code{c(n,J)}, with \code{J} evaluations on columns and \code{n} units on rows.
}
  \item{mu}{
The center of symmetry under the null hypothesis: either a constant (in this case, a constant function is used) or a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated. The default is \code{mu=0}.
}
  \item{order}{
Order of the B-spline basis expansion. The default is \code{order=2}.
}
  \item{nknots}{
Number of knots of the B-spline basis expansion. The default is \code{nknots=dim(data)[2]}.
}
  \item{B}{
The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=10000}.
}
}

\value{
\code{ITP1bspline} returns an object of \code{\link{class}} "\code{ITP1}".

An object of class "\code{ITP1}" is a list containing at least the following components:

\item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"B-spline"}.}
\item{test}{String vector indicating the type of test performed. In this case equal to \code{"1pop"}.}
\item{mu}{Center of symmetry under the null hypothesis (as entered by the user).}
\item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion. Rows are associated to units and columns to the basis index.}
\item{pval}{Uncorrected p-values for each basis coefficient.}
\item{pval.matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
\item{corrected.pval}{Corrected p-values for each basis coefficient.}
\item{labels}{Labels indicating the population membership of each data (in this case always equal to \code{1}).}
\item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
\item{heatmap.matrix}{Heatmap matrix of p-values (used only for plots).}
}
\references{
A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

}
\author{
Alessia Pini, Simone Vantini
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}}, and \code{\link{ITPimage}}.
}
\examples{

# Importing the NASA temperatures data set
data(NASAtemp)
# Performing the ITP for two populations with the B-spline basis
ITP.result <- ITP1bspline(NASAtemp$paris,mu=4,nknots=50,B=1000)
# Plotting the results of the ITP
plot(ITP.result,xrange=c(0,12),main='Paris temperatures')

# Plotting the p-value heatmap
ITPimage(ITP.result,abscissa.range=c(0,12))


# Selecting the significant components for the radius at 5% level
which(ITP.result$corrected.pval < 0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ functional data analysis }
\keyword{ local inference }
\keyword{ permutation test }
\keyword{ one-population test }
