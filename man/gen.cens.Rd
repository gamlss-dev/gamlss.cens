\name{gen.cens}
\alias{gen.cens}

\title{A Function to Generate Appropriate Functions to Be Used to Fit a Censored Response variable in GAMLSS}
\description{
The \code{gen.cens()} function allows the user to generate a
\code{d}, \code{p}, (dummy) \code{q} and fitting  \code{gamlss} functions for  
censor and interval response variables.  The function can take any  \code{\link[gamlss.dist]{gamlss.family}} distribution. 
}
\usage{
gen.cens(family = "NO", name = "cens", 
          type = c("right", "left", "interval"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family}{a \code{\link[gamlss.dist]{gamlss.family}} object, which is 
        used to define the distribution and the link functions of the various parameters.
        The distribution families supported by \code{gamlss()} 
        can be found in \code{\link[gamlss.dist]{gamlss.family}} and in the package \code{gamlss.dist}.}
  \item{name}{the characters you want to add to the name of new functions, by default is the first letter of \code{type} and \code{c}
               i.e WEIic for WEI (weibull) interval response variable}
   \item{type}{whether \code{right}, \code{left}  or in \code{interval} censoring is required, 
              (right is the default) }
  \item{\dots}{for extra arguments}
}

\value{
  Returns  the \code{d}, \code{p},  (dummy) \code{q} and the fitting used in the fitting gamlss algorithm (The one 
         used in the fitting gamlss algorithm)  of a \code{gamlss.family} distribution. 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
 }
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby}

\seealso{ \code{\link{cens.d}}, \code{\link{cens.p}}, \code{\link{cens}}  }
\examples{
library(gamlss.dist)
data(lip)
gen.cens(WEI,type="interval") 
WEIic
gg1<- gamlss(y ~ poly(Tem,2)+poly(pH,2)+poly(aw,2), data=lip, 
     family=WEIic, c.crit=0.00001, n.cyc=200, trace=FALSE)
}
\keyword{distribution}
\keyword{regression}
