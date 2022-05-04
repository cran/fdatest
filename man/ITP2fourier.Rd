\name{ITP2fourier}
\alias{ITP2fourier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two populations Interval Testing Procedure with Fourier basis
}
\description{
The function implements the Interval Testing Procedure for testing the difference between two functional populations evaluated on a uniform grid. Data are represented by means of the Fourier basis and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate. 
}
\usage{
ITP2fourier(data1, data2, mu = 0, 
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
  \item{mu}{
The difference between the first functional population and the second functional population under the null hypothesis. Either a constant (in this case, a constant function is used) or a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated. The default is \code{mu=0}.
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
\code{ITP2fourier} returns an object of \code{\link{class}} "\code{ITP2}".

An object of class "\code{ITP2}" is a list containing at least the following components:

\item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"Fourier"}.}
\item{test}{String vector indicating the type of test performed. in this case equal to \code{"2pop"}.}
\item{mu}{Difference between the first functional population and the second functional population under the null hypothesis (as entered by the user).}
\item{paired}{Logical indicating whether the test is paired (as entered by the user).}
\item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the Fourier basis expansion. Rows are associated to units and columns to the basis index: the first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients of the second population units; the first column is a0, the following \code{(p-1)/2} columns are the ak coefficients (sine coefficients) and the last \code{(p-1)/2} columns the bk coefficients (cosine coefficients).}
\item{pval}{Uncorrected p-values for each frequency.}
\item{pval.matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the frequencies \code{(j,j+1,...,j+(p-i))}.}
\item{corrected.pval}{Corrected p-values for each frequency.}
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
See also \code{\link{ITP2pafourier}},  \code{\link{ITP2bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP1bspline}}, and \code{\link{ITPimage}}.
}
\examples{

# Importing the NASA temperatures data set
data(NASAtemp)
# Performing the ITP
ITP.result <- ITP2fourier(NASAtemp$milan,NASAtemp$paris,maxfrequency=20,B=1000,paired=TRUE)

# Plotting the results of the ITP
plot(ITP.result,main='NASA data',xrange=c(1,365),xlab='Day')

# Plotting the p-value heatmap
ITPimage(ITP.result,abscissa.range=c(1,365))

# Selecting the significant coefficients
which(ITP.result$corrected.pval < 0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ functional data analysis }
\keyword{ local inference }
\keyword{ permutation test }
\keyword{ two-population test }