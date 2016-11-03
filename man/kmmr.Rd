\name{kmmr}
\alias{kmmr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial mark correlation function
%%  ~~function to do ... ~~
}
\description{Hola a todos
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kmmr(xyt,s.region,s.lambda,ds,ks="epanech",hs,correction="none",approach="simplified")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyt}{Spatial coordinates and times \eqn{(x,y,t)} of the point pattern.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If s.region is missing, the Ripley-Rasson estimate of the spatial domain is considered.}
  \item{s.lambda}{Vector of values of the spatial intensity function evaluated at the points \eqn{(x,y)} in \eqn{W}. If s.lambda is missing, the estimate of the spatial mark correlation function is computed as for the homogeneous case, i.e. considering \eqn{n/|W|} as an estimate of the spatial intensity.}
  \item{ds}{A vector of distances \code{u} at which \eqn{k_[mm](u)} is computed.}
  \item{ks}{A kernel function for the spatial distances. The default is the \code{"epanech"} kernel. It can also be \code{"box"} for the uniform kernel, or \code{"biweight"}.}
  \item{hs}{A bandwidth of the kernel function \code{ks}.}
  \item{correction}{A character vector specifying the edge-correction(s) to be applied among "isotropic","border", "modified.border", "translate", "setcovf" and "none". The default is "none".}
  \item{approach}{A character vector specifying the approach to use for the estimation to be applied among "simplified" or "standardised".}
%%     ~~Describe \code{x} here~~
}
\details{By default, this command calculates an estimate of Stoyan's mark correlation k_[mm](r) for the point pattern.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{kmmr}{A vector containing the values of \eqn{k_[mm](u)} estimated.}
\item{ds}{If \code{ds} is missing, a vector of distances \code{u} at which \eqn{k_{mm}(u)} is computed from 0 to until quarter of the maximum distance between the points in the pattern.}
\item{kernel}{A vector of names and bandwidth of the spatial kernel.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Baddeley, A., Rubak, E., Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.

Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J. and Wilfried, G. (2016). Mark variograms for spatio-temporal
point processes, Submitted .
}
\author{Francisco J. Rodriguez Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
Rodriguez-Cortes, F. J.(2015). \code{mvstpp}: Mark Variogram for Spatio-Temporal Point Processes. GitHub repository. URL \url{https://github.com/frajaroco/mvstpp}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")%
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
