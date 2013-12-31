\name{cens.d}
\alias{cens.d}

\title{Censored probability density function of a gamlss.family distribution}
\description{
 Creates a probability density function from a current \code{\link[gamlss.dist]{gamlss.family}}  
 distribution to be used for fitting a censored or interval response variable.
}
\usage{
cens.d(family = "NO", type = c("right", "left", "interval"), ...)
}
\arguments{
  \item{family}{ a \code{\link[gamlss.dist]{gamlss.family}} object, which is used to 
          define the distribution and the link functions of the various parameters. 
          The distribution families supported by \code{gamlss()} 
        can be found in \code{\link[gamlss.dist]{gamlss.family}} and in the package \code{gamlss.dist}.}
  \item{type}{ whether \code{right}, \code{left}  or in \code{interval} censoring is required, (right is the default) }
  \item{\dots}{for extra arguments }
}

\details{
This function is used to calculate the likelihood function for censored data. 
This function is not supposed to be used on its own but it is used in function \code{gen.cens}. 
}
\value{
Returns a modified d family function.  
The argument of the original function \code{d} function are the same.  
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
 }
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}}
\note{For an example see \code{gen.cens()}   
}
\seealso{  \code{\link{cens.p}},  \code{\link{gen.cens}}}
\examples{#see the help for function cens for an example}

\keyword{distribution}
\keyword{regression}
