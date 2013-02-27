\name{NASAtemp}
\alias{NASAtemp}
\docType{data}
\title{
NASA daily temperatures data set
}
\description{
It contains the daily mean temperatures registered from July 1983 to June 2005 and stored in the NASA database Earth Surface Meteorology for Solar Energy of two different geographical locations: the region (45-46 North, 9-10 East), including the city of Milan (Italy), and the region (48-49 North, 2-3 East), including the city of Paris (France).
}
\usage{data(NASAtemp)}
\format{
List of 2 elements:
 \itemize{
	\item{\code{milan}}{ Matrix of dimensions \code{c(22,365)} containing the daily mean temperatures of the region (45-46 North, 9-10 East), including the city of Milan (Italy) registered from July 1983 to June 2005 (22 years).}
	\item{\code{paris}}{ Matrix of dimensions \code{c(22,365)} containing the daily mean temperatures of the region (48-49 North, 2-3 East), including the city of Paris (France) registered from July 1983 to June 2005 (22 years).}
 }
}

\source{
These data were obtained from the NASA Langley Research Center Atmospheric Science Data Center Surface meteorological and Solar Energy (SSE) web portal supported by the NASA LaRC POWER Project. Data are freely available at: NASA Surface Meteorology and Solar Energy, A Renewable Energy Resource web site (release 6.0): http://eosweb.larc.nasa.gov
}

\examples{
data(NASAtemp)
\dontrun{

matplot(t(NASAtemp$milan),type='l')
matplot(t(NASAtemp$paris),type='l')

}
}
\keyword{datasets}
