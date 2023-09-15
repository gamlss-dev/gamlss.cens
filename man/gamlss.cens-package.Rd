\name{gamlss.cens-package}
\alias{gamlss.cens-package}
\alias{gamlss.cens}
\docType{package}
\title{
\packageTitle{gamlss.cens}
}
\description{
\packageDescription{gamlss.cens}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{gamlss.cens}
\packageIndices{gamlss.cens}
}
\author{
\packageAuthor{gamlss.cens}

Maintainer: \packageMaintainer{gamlss.cens}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{https://www.gamlss.com/}).

}

\keyword{ package }
\seealso{\code{\link[gamlss]{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}}  
}
\examples{
library(survival)
library(gamlss)
library(gamlss.dist)
# comparing results with package survival
# fitting the exponential distribution
ms1<-survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, 
             dist='exponential')
mg1<-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
             family=cens(EXP),c.crit=0.00001)
if(abs(-2*ms1$loglik[2]-deviance(mg1))>0.001) stop(paste("descrepancies in exp")) 
if(sum(coef(ms1)-coef(mg1))>0.001) warning(paste("descrepancies in coef in exp")) 
summary(ms1)
summary(mg1)
# fitting the Weibull distribution
ms2 <-survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull')
mg2 <-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
           family=cens(WEI, delta=c(0.001,0.001)), c.crit=0.00001)
if(abs(-2*ms2$loglik[2]-deviance(mg2))>0.005) 
     stop(paste("descrepancies in deviance in WEI")) 
summary(ms2);summary(mg2)
# compare the scale parameter
 1/exp(coef(mg2,"sigma"))
# now fit the Weibull in different parameterrazions  
mg21<-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
             family=cens(WEI2), method=mixed(2,30)) 
mg21<-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
             family=cens(WEI3)) 
}
