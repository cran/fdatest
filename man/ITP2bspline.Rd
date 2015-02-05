\name{ITP2bspline}
\alias{ITP2bspline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two populations Interval Testing Procedure with B-spline basis
}
\description{
The function implements the Interval Testing Procedure for testing the difference between two functional populations evaluated on a uniform grid. Data are represented by means of the B-spline basis and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate. The default parameters of the basis expansion lead to the piece-wise interpolating function.
}
\usage{
ITP2bspline(data1, data2, mu = 0, 
            order = 2, nknots = dim(data1)[2], B = 10000, paired = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
Pointwise evaluations of the first population's functional data set on a uniform grid. \code{data1} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
}
  \item{data2}{
Pointwise evaluations of the second population's functional data set on a uniform grid. \code{data2} is a matrix of dimensions \code{c(n2,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
}
  \item{mu}{
The difference between the first functional population and the second functional population under the null hypothesis. Either a constant (in this case, a constant function is used) or a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated. The default is \code{mu=0}.
}
  \item{order}{
Order of the B-spline basis expansion. The default is \code{order=2}.
}
  \item{nknots}{
Number of knots of the B-spline basis expansion. The default is \code{nknots=dim(data1)[2]}.
}
  \item{B}{
The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=10000}.
}
  \item{paired}{
A logical indicating whether the test is paired. The default is \code{FALSE}.
}
}

\value{
\code{ITP2bspline} returns an object of \code{\link{class}} "\code{ITP2}".

An object of class "\code{ITP2}" is a list containing at least the following components:

\item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"B-spline"}.}
\item{test}{String vector indicating the type of test performed. In this case equal to \code{"2pop"}.}
\item{mu}{Difference between the first functional population and the second functional population under the null hypothesis (as entered by the user).}
\item{paired}{Logical indicating whether the test is paired (as entered by the user).}
\item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion, with \code{n=n1+n2}. Rows are associated to units and columns to the basis index. The first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients of the second population units}
\item{pval}{Uncorrected p-values for each basis coefficient.}
\item{pval.matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
\item{corrected.pval}{Corrected p-values for each basis coefficient.}
\item{labels}{Labels indicating the population membership of each data.}
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
For tests of comparison between two populations, see \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}}.
For different types of ITP-based tests, see \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITPlmbspline}}, \code{\link{ITPaovbspline}} and \code{\link{ITPimage}}.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)
# Performing the ITP
ITP.result <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=50,B=1000)

# Plotting the results of the ITP
plot(ITP.result,main='NASA data',xrange=c(1,365),xlab='Day')

# Plotting the p-values heatmap
ITPimage(ITP.result,abscissa.range=c(0,12))

# Selecting the significant components at 5% level
which(ITP.result$corrected.pval < 0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
