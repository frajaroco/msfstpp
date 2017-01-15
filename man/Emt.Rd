\name{Emt}
\alias{Emt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Temaporal mean mark function
%%  ~~function to do ... ~~
}
\description{Computes an estimator of the temporal mean mark function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{Emt(xyt,t.region,t.lambda,dt,kt="epanech",ht,correction="none",approach="simplified")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyt}{Spatial coordinates and times \eqn{(x,y,t)} of the point pattern.}
  \item{t.region}{Vector containing the minimum and maximum values of the time interval.}
  \item{t.lambda}{Vector of values of the temporal intensity function evaluated at the points \eqn{t} in \eqn{T}.}
  \item{dt}{A vector of times \code{v} at which \code{Emt(v)} is computed.}
  \item{kt}{A kernel function for the temporal distances. The default is the \code{"epanech"} kernel. It can also be \code{"box"} for the uniform kernel, or \code{"biweight"}.}
  \item{ht}{A bandwidth of the kernel function \code{kt}.}
  \item{correction}{A character vector specifying the edge-correction(s) to be applied among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"}, \code{"setcovf"} and \code{"none"}. The default is \code{"none"}.}
  \item{approach}{A character vector specifying the approach to use for the estimation to be applied among "simplified" or \code{"standardised"}. If approach is missing, \code{"simplified"} is considered by default.}
%%     ~~Describe \code{x} here~~
}
\details{By default, this command calculates an estimate of the temporal mean mark function \eqn{E_[te](t)} for a spatio-temporal point pattern.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{eEmt}{A vector containing the values of \eqn{E_[te](v)} estimated.}
\item{dt}{Parameter passed in argument. If \code{dt} is missing, a vector of temporal distances \code{v} at which \code{Emt(v)} is computed from 0 to until quarter of the maximum distance between the times in the temporal pattern.}
\item{kernel}{Parameters passed in argument. A vector of names and bandwidth of the spatial kernel.}
\item{t.region}{If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
\item{t.lambda}{If \code{t.lambda} is missing, the estimate of the temporal mark correlation function is computed as for the homogeneous case, i.e. considering \eqn{n/|T|} as an estimate of the temporal intensity under the parameter \code{approach="standardised"}.}
}
%% ...
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

\examples{
## Not run:
#################

# A realisation of spatio-temporal homogeneous Poisson point processes
hpp <- rpp(lambda = 100, replace = FALSE)$xyt

# R plot
plot(hpp)

# This function provides an kernel estimator of the temporal mean mark function
out <- Emt(hpp,ht=0.12/sqrt(100))

# Temporal mean mark constant under homogeneous Poisson point processes
a <- 1
b <- 1

ete <- ((a^5+b^5-sqrt(a^2+b^2)*(a^4-3*(b*a)^2+b^4))/(15*(a*b)^2))+
((b^2*log((a+sqrt(a^2+b^2))/b))/(6*a))+((a^2*log((b+sqrt(a^2+b^2))/a))/(6*b))

# R plot - Temporal mean mark function
par(mfrow=c(1,1))
yl <- c(min(out$eEmt,ete),max(out$eEmt,ete))
plot(out$dt,out$eEmt,type="l",xlab="t = time",ylab=expression(E[te](t)),
     ylim=yl,col=1,cex.lab=1.5,cex.axis=1.5)
lines(out$dt,rep(ete,length(out$dt)),col=11)

## End(Not run)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")%
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
