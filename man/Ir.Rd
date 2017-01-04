\name{Ir}
\alias{Ir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial Moran's \eqn{I}-statistics
%%  ~~function to do ... ~~
}
\description{Computes an estimator of the spatial Moran's \eqn{I}-statistics.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Ir(xyt, s.region, s.lambda, ds, ks = "epanech", hs, correction = "none", approach = "simplified")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyt}{Spatial coordinates and times \eqn{(x,y,t)} of the point pattern.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. }
  \item{s.lambda}{Vector of values of the spatial intensity function evaluated at the points \eqn{(x,y)} in \eqn{W}.}
  \item{ds}{A vector of distances \code{u} at which \eqn{Ir(u)} is computed.}
  \item{ks}{A kernel function for the spatial distances. The default is the \code{"epanech"} kernel. It can also be \code{"box"} for the uniform kernel, or \code{"biweight"}.}
  \item{hs}{A bandwidth of the kernel function \code{ks}.}
  \item{correction}{A character vector specifying the edge-correction(s) to be applied among "isotropic","border", "modified.border", "translate", "setcovf" and "none". The default is "none".}
  \item{approach}{A character vector specifying the approach to use for the estimation to be applied among "simplified" or "standardised". If approach is missing, "simplified" is considered by default.}
%%     ~~Describe \code{x} here~~
}
\details{By default, this command calculates an estimate of the spatial Moran's \eqn{I}-statistics \eqn{Ir(r)} for a spatio-temporal point pattern.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{eImr}{A vector containing the values of \eqn{Ir(u)} estimated.}
\item{Imr0}{A single value of \eqn{Ir(0)} estimated at \eqn{u=0}.}
\item{ds}{If \code{ds} is missing, a vector of distances \code{u} at which \eqn{Ir(u)} is computed from 0 to until quarter of the maximum distance between the points in the pattern.}
\item{kernel}{A vector of names and bandwidth of the spatial kernel.}
\item{s.region}{If \code{s.region} is missing, the Ripley-Rasson estimate of the spatial domain is considered under the parameter \code{approach="standardised"}. The bounding box of \code{xyt[,1:2]} is considered under the parameter \code{approach="simplified"}.}
\item{s.lambda}{If \code{s.lambda} is missing, the estimate of the spatial mark correlation function is computed as for the homogeneous case, i.e. considering \eqn{n/|W|} as an estimate of the spatial intensity under the parameter \code{approach="standardised"}.}
}
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
## Not run:
#################

# Realisations of the homogeneous Poisson processes
hpp <- rpp(lambda = 100, replace = FALSE)$xyt

# R plot
plot(hpp)

# This function provides an kernel estimator of the spatial mark correlation function
out <- Ir(hpp,hs=0.12/sqrt(100))

# R plot - Spatial mark correlation function
par(mfrow=c(1,1))
xl <- c(0,0.25)
yl <- c(min(1,out$eImr,out$Imr0),max(1,out$eIr,out$Imr0))
plot(out$ds,out$eImr,type="l",xlab="r = distance",ylab=expression(I(r)),
                 xlim=xl,ylim=yl,col=1,cex.lab=1.5,cex.axis=1.5)
lines(out$ds,rep(0,length(out$ds)),col=11)
points(0,out$Imr0,col=3,pch=13,cex=1.5)

## End(Not run)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")%
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
