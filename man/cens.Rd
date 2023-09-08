\name{cens}
\alias{cens}

\title{Function to Fit Censored Data Using a gamlss.family Distribution}

\description{ This function can be used to fit censored or interval response variables. 
 It takes as an argument an existing \code{\link[gamlss.dist]{gamlss.family}} distribution  
 and  generates 
 a new \code{\link[gamlss.dist]{gamlss.family}} object which then can be used to fit 
 right, left or interval censored data. 
}
\usage{
cens(family = "NO", type = c("right", "left", "interval"), name = "cens", 
       local = TRUE, delta = NULL, ...)
}

\arguments{
   \item{family}{ a \code{\link[gamlss.dist]{gamlss.family}} object, which is used to 
          define the distribution and the link functions of the various parameters. 
          The distribution families supported by \code{gamlss()} 
        can be found in \code{\link[gamlss.dist]{gamlss.family}} and in the package \code{gamlss.dist}.}
  \item{name}{the characters you want to add to the name of new functions, by default is \code{cens}}
  \item{type}{what type of censoring is required, \code{right}, \code{left} or \code{interval}.}
  \item{local}{if TRUE the function will try to find the environment of \code{gamlss} to generate 
               the \code{d} and \code{p} functions required for the fitting, 
               if FALSE the functions will be generated in the global environment}
  \item{delta}{the delta increment used in the numerical derivatives} 
  \item{\dots}{for extra arguments }
}
\details{
This function is created to help users to fit censored data using an existing 
\code{gamlss.family} distribution.
It does this by taking an existing \code{\link[gamlss.dist]{gamlss.family}} and changing 
some of the components of the distribution to help the fitting process. 
It particular it (i) creates a (\code{d}) function (for calculating the censored 
likelihood) and a (\code{p}) function (for generating the quantile residuals) 
within \code{gamlss}, 
(ii) changes  the global deviance function \code{G.dev.incr}, 
the first derivative functions (see note below) 
   and other quantities from the original distribution.   
}
\value{
 It returns a \code{\link[gamlss.dist]{gamlss.family}} object which has all the components needed for fitting a distribution in \code{gamlss}.
}
\references{ Rigby, R. A. and Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), Appl. Statist., 54, part 3, 1-38.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}}

\note{
This function is experimental and could be changed in the future. 
The function \code{cens} changes the first derivatives of the original gamlss family 
\code{d} function to numerical derivatives for the new censored \code{d} function. 
The default increment \code{delta}, for this numerical derivatives function, 
is \code{eps * pmax(abs(x), 1)} where  \code{eps<-sqrt(.Machine$double.eps)}. 
The default \code{delta} could be inappropriate for 
specific applications and can be overwritten by using the argument \code{delta}.

Note that in order to get the correct standard errors you have to generate the "d" function by using
\code{gen.cens()}.
}

\seealso{\code{\link{cens.d}}, \code{\link{cens.p}},  \code{\link{gen.cens}} 
}
\examples{
# comparing output with the survreg() of package survival
library(gamlss.dist)
library(survival)
#--------------------------------------------------------------------
# right censoring example 
# example from survreg() 
# fitting the exponential distribution
mexp<-survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='exponential')
gexp<-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
             family=cens(EXP), c.crit=0.00001)
if(abs(-2*mexp$loglik[2]-deviance(gexp))>0.001) 
         stop(paste("descrepancies in exponential models")) 
if(sum(coef(mexp)-coef(gexp))>0.001) 
        warning(paste("descrepancies in coef in exponential models")) 
summary(mexp)
gen.cens(EXP)
summary(gexp)
# fitting different distributions
# weibull 
mwei <-survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull')
gwei<-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
             family=cens(WEI, delta=c(0.0001,0.0001)), c.crit=0.00001)
if(abs(-2*mwei$loglik[2]-deviance(gwei))>0.005) 
        stop(paste("descrepancies in deviance in WEI")) 
scoef <- sum(coef(mwei)-coef(gwei))
if(abs(scoef)>0.005) 
         warning(cat("descrepancies in coef in WEI of ", scoef, "\n")) 
# WEI3 is weibull parametrised with mu as the mean
gwei3 <- gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
                 family=cens(WEI3)) 
# log normal
mlogno <-survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, 
                  dist='lognormal')
glogno<-gamlss(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
                family=cens(LOGNO, delta=c(0.001,0.001)), c.cyc=0.00001)
if(abs(-2*mlogno$loglik[2]-deviance(glogno))>0.005) 
          stop(paste("descrepancies in deviance in LOGNO")) 
coef(mlogno);coef(glogno) 
#-------------------------------------------------------------------- 
# now interval response variable 
data(lip)
with(lip, y)
mg1<-survreg(y ~ poly(Tem,2)+poly(pH,2)+poly(aw,2), data=lip, dist="weibull")
gg1<- gamlss(y ~ poly(Tem,2)+poly(pH,2)+poly(aw,2), data=lip, 
      family=cens(WEI,type="interval"), c.crit=0.00001, n.cyc=200, trace=FALSE)
summary(mg1)
gen.cens(WEI,type="interval")
summary(gg1)
#--------------------------------------------------------------------
# now fitting discretised continuous distribution to count data
# fitting discretised Gamma
data(species)
# first generate the distributions
gen.cens(GA, type="interval")
gen.cens(IG, type="interval")
 mGA<-gamlss(Surv(fish,fish+1,type= "interval2")~log(lake)+I(log(lake)^2), 
       sigma.fo=~log(lake), data=species, family=GAic)
# fitting discretised inverse Gaussian
 mIG<-gamlss(Surv(fish,fish+1,type= "interval2")~log(lake)+I(log(lake)^2), 
      sigma.fo=~log(lake), data=species, family=IGic)
AIC(mGA,mIG)
plot(fish~log(lake), data=species)
with(species, lines(log(lake)[order(lake)], fitted(mIG)[order(lake)]))             
#--------------------------------------------------------------------
}
\keyword{distribution}
\keyword{regression}
