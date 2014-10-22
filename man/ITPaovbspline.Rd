\name{ITPaovbspline}
\alias{ITPaovbspline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interval Testing Procedure for testing Functional analysis of variance with B-spline basis
}
\description{
\code{ITPaovbspline} is used to fit and test functional analysis of variance.
The function implements the Interval Testing Procedure for testing for significant differences between several functional population evaluated on a uniform grid. Data are represented by means of the B-spline basis and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate. The default parameters of the basis expansion lead to the piece-wise interpolating function.
}
\usage{
ITPaovbspline(formula, order = 2, 
              nknots = dim(model.response(model.frame(formula)))[2], 
              B = 10000, method = "residuals")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{order}{
Order of the B-spline basis expansion. The default is \code{order=2}.
}
  \item{nknots}{
Number of knots of the B-spline basis expansion. 

The default is \code{dim(model.response(model.frame(formula)))[2]}.
}
  \item{B}{
The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=10000}.
}
  \item{method}{
Permutation method used to calculate the p-value of permutation tests. Choose "\code{residuals}" for the permutations of residuals under the reduced model, according to the Freedman and Lane scheme, and "\code{responses}" for the permutation of the responses, according to the Manly scheme. 
}
}
\value{
\code{ITPaovbspline} returns an object of \code{\link{class}} "\code{ITPaov}".

The function \code{summary} is used to obtain and print a summary of the results. 

An object of class "\code{ITPlm}" is a list containing at least the following components:

\item{call}{The matched call.}
\item{design.matrix}{The design matrix of the functional-on-scalar linear model.}
\item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"B-spline"}.}
\item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion. Rows are associated to units and columns to the basis index.}
\item{coeff.regr}{Matrix of dimensions \code{c(L+1,p)} of the \code{p} coefficients of the B-spline basis expansion of the intercept (first row) and the \code{L} effects of the covariates specified in \code{formula}. Columns are associated to the basis index.}
\item{pval.F}{Uncorrected p-values of the functional F-test for each basis coefficient.}
\item{pval.matrix.F}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate F-tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
\item{corrected.pval.F}{Corrected p-values of the functional F-test for each basis coefficient.}
\item{pval.factors}{Uncorrected p-values of the functional F-tests on each factor of the analysis of variance, separately (rows) and each basis coefficient (columns).}
\item{pval.matrix.factors}{Array of dimensions \code{c(L+1,p,p)} of the p-values of the multivariate F-tests on factors. The element \code{(l,i,j)} of array \code{pval.matrix} contains the p-value of the joint NPC test on factor \code{l} of the components \code{(j,j+1,...,j+(p-i))}.}
\item{corrected.pval.factors}{Corrected p-values of the functional F-tests on each factor of the analysis of variance (rows) and each basis coefficient (columns).}
\item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
\item{coeff.regr.eval}{Evaluation on a fine uniform grid of the functional regression coefficients.}
\item{fitted.eval}{Evaluation on a fine uniform grid of the fitted values of the functional regression.}
\item{residuals.eval}{Evaluation on a fine uniform grid of the residuals of the functional regression.}
\item{R2.eval}{Evaluation on a fine uniform grid of the functional R-squared of the regression.}
\item{heatmap.matrix.F}{Heatmap matrix of p-values of functional F-test (used only for plots).}
\item{heatmap.matrix.factors}{Heatmap matrix of p-values of functional F-tests on each factor of the analysis of variance (used only for plots).}
}
\references{
D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported Significance Levels. \emph{Journal of Business & Economic Statistics} 1.4, 292-298.

B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods in Biology}. Vol. 70. CRC Press.

A. Pini and S. Vantini (2013). The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. \emph{MOX-report 13/2013}, Politecnico di Milano.

}
\author{
Alessia Pini, Simone Vantini
}


\seealso{
See \code{\link{summary.ITPaov}} for summaries and \code{\link{plot.ITPaov}} for plotting the results. 

See also \code{\link{ITPlmbspline}} to fit and test a functional-on-scalar linear model applying the ITP, and \code{\link{ITP1bspline}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}} for one-population and two-population tests.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the ITP
\dontrun{

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
