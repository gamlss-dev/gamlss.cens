################################################################################
################################################################################
################################################################################
################################################################################
cens <- function (
              family = "NO", 
                type = c("right", "left",  "interval"),
                name = "cens", 
               local = TRUE,
               delta = NULL, 
                ...)
{
#------------------------------------------
## dummy name
     TEST <- "TEST" 
     type <- match.arg(type)   
    fname <- if (is.name(family)) as.character(family)
             else if (is.character(family)) family
             else if (is.call(family)) as.character(family[[1]])
             else if (is.function(family)) deparse(substitute(family))
             else if (is(family, "gamlss.family"))  family$family[1]
             else stop("the family must be a character or a gamlss.family name")
     fam1 <- eval(parse(text=fname)) # the family to output
      fam <- as.gamlss.family(family) # this is created so I can get things
   family <- c("None", "None") 
   dorfun <- paste("d",fname,sep="") # say dNO
   porfun <- paste("p",fname,sep="") # say pNO
     dfun <- paste(paste("d",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say dNOrc
     pfun <- paste(paste("p",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say pNOrc
    nopar <- fam$nopar # or fam1$nopar
#-------------------------------------------------------------------------------
if (local)
 {
#--trying to get gamlss sys.frame--  
     rexpr<-regexpr("gamlss",sys.calls())
for (i in 1:length(rexpr)){ 
    position <- i 
    if (rexpr[i]==1) break}
gamlss.environment <- sys.frame(position)      
#--end here------------------------
 }
 else gamlss.environment <- sys.frame(0)
#   generate d within gamlss
    eval(dummy <- cens.d(family = fname, type = type, ...))
    eval(call("<-",as.name(dfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# generate p within gamlss
    eval(dummy <- cens.p(family = fname, type = type, ...))
    eval(call("<-",as.name(pfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
#-------------------------------------------------------------------------------
# rename the family 
   family[[1]] <- paste(fname, substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="")
   family[[2]] <- paste(type, "censored",fam$family[[2]])
    fam$family <- family
body(fam1)[[nopar+2]][[2]]$family <- family # and in fam1
# Global deviance increment  
           sGD <- gsub(dorfun, dfun, deparse(body(fam$G.dev.incr)))
  body(fam$G.dev.incr) <- parse(text=sGD)
body(fam1)[[nopar+2]][[2]]$G.dev.incr <- fam$G.dev.incr
# check for the delta
 if (length(delta)==0) delta <- rep(NA,nopar) 
 if (length(delta)==1) delta <- rep(delta,nopar)
 if (length(delta)!=nopar)  stop("delta should be the same length the parameters in the family ") 
#-------------------------------------------------------------------------------
# now change the first derivatives
  switch(nopar,  
          {
          # 1 parameter 
      fam$dldm <- function(y,mu) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, log=TRUE), "mu", delta=NULL), "gradient")) 
           # mu
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)])
body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
          # residuals
      cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
     sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres[[1]])) 
          #sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)
body(fam1)[[nopar+2]][[2]]$rqres <- fam$rqres
         # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial[[1]]))
         #  inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
body(fam1)[[nopar+2]][[2]]$mu.initial <- fam$mu.initial
          #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$y.valid) <- parse(text=yval)[[1]] 
body(fam1)[[nopar+2]][[2]]$y.valid <- fam$y.valid
          }, 
          { 
          # 2 parameters  
      fam$dldm <- function(y,mu,sigma) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, log=TRUE), "sigma", delta=NULL), "gradient"))
          # mu 
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
         body(fam$dldm) <- parse(text=sMU[length(sMU)]) 
         body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
         # sigma  
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA)
         body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])
         body(fam1)[[nopar+2]][[2]]$dldd  <- fam$dldd
          #d2ldmdd 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdd)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldmdd) <- parse(text=yval)[[1]] 
          body(fam1)[[nopar+2]][[2]]$d2ldmdd  <- fam$d2ldmdd
           #d2ldm2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldm2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldm2) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldm2  <- fam$d2ldm2
          #d2ldd2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldd2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldd2) <- parse(text=yval)[[1]] 
          body(fam1)[[nopar+2]][[2]]$d2ldd2  <- fam$d2ldd2
          #residuals
         #cen <- paste("censored=", "\"", type, "censored" , "\"")
         # sres <- gsub(porfun, paste(pfun,cen),  deparse(fam$rqres)) 
           cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres[[1]])) 
         #sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)
          body(fam1)[[nopar+2]][[2]]$rqres <- fam$rqres
          # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial[[1]]))
           #inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
           body(fam1)[[nopar+2]][[2]]$mu.initial <- fam$mu.initial
         # initial sigma
           inisigma <- gsub("y", "y[,1]",  deparse(fam$sigma.initial[[1]]))
         # inisigma <- gsub("expression", "",  inisigma)
           fam$sigma.initial <- parse(text=inisigma)
           body(fam1)[[nopar+2]][[2]]$sigma.initial <- fam$sigma.initial
         #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           body(fam$y.valid) <- parse(text=yval)[[1]]
           body(fam1)[[nopar+2]][[2]]$y.valid <- fam$y.valid
           }, 
           {
           # 3 parameters   
      fam$dldm <- function(y,mu,sigma,nu) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma,nu) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function(y,mu,sigma,nu) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "nu", delta=NULL), "gradient"))
           # mu
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1]))sMU <- sub("NULL",  as.character(delta[1]), sMU)          
           body(fam$dldm) <- parse(text=sMU[length(sMU)]) 
           body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
           # sigma 
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
           body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])
           body(fam1)[[nopar+2]][[2]]$dldd  <- fam$dldd
           # nu   
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)
           body(fam$dldv) <- parse(text=sNU[length(sNU)])
           body(fam1)[[nopar+2]][[2]]$dldv  <- fam$dldv
           # d2ldmdd 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdd)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldmdd) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldmdd  <- fam$d2ldmdd
           #d2ldmdv 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldmdv) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldmdv  <- fam$d2ldmdv
           #d2ldddv 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldddv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldddv) <- parse(text=yval)[[1]] 
          body(fam1)[[nopar+2]][[2]]$d2ldddv  <- fam$d2ldddv
            #d2ldm2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldm2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldm2) <- parse(text=yval)[[1]] 
          body(fam1)[[nopar+2]][[2]]$d2ldm2  <- fam$d2ldm2
          #d2ldd2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldd2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldd2) <- parse(text=yval)[[1]] 
          body(fam1)[[nopar+2]][[2]]$d2ldd2  <- fam$d2ldd2
           #d2ldv2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldv2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldv2) <- parse(text=yval)[[1]] 
          body(fam1)[[nopar+2]][[2]]$d2ldv2  <- fam$d2ldv2
           #residuals
      cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres[[1]])) 
         #sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres) 
          body(fam1)[[nopar+2]][[2]]$rqres  <- fam$rqres
           # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial[[1]]))
           # inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
           body(fam1)[[nopar+2]][[2]]$mu.initial  <- fam$mu.initial
          # initial sigma
           inisigma <- gsub("y", "y[,1]",  deparse(fam$sigma.initial[[1]]))
          # inisigma <- gsub("expression", "",  inisigma)
           fam$sigma.initial <- parse(text=inisigma)
           body(fam1)[[nopar+2]][[2]]$sigma.initial  <- fam$sigma.initial
           # initial nu
           ininu <- gsub("y", "y[,1]",  deparse(fam$nu.initial[[1]]))
           #ininu <- gsub("expression", "",  ininu)
           fam$nu.initial <- parse(text=ininu)
           body(fam1)[[nopar+2]][[2]]$nu.initial  <- fam$nu.initial
          #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           body(fam$y.valid) <- parse(text=yval)[[1]]
           body(fam1)[[nopar+2]][[2]]$y.valid  <- fam$y.valid
           },
           {
           # 4 parameters
      fam$dldm <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=NULL), "gradient"))
      fam$dldt <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=NULL), "gradient"))
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU)      
         body(fam$dldm) <- parse(text=sMU[length(sMU)])
         body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
         body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
         body(fam1)[[nopar+2]][[2]]$dldd  <- fam$dldd
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)           
         body(fam$dldv) <- parse(text=sNU[length(sNU)])
         body(fam1)[[nopar+2]][[2]]$dldv  <- fam$dldv
          sTAU <- sub("TEST", dfun, body(fam$dldt))
if (!is.na(delta[4])) sTAU <- sub("NULL",  as.character(delta[4]), sTAU)
         body(fam$dldt) <- parse(text=sTAU[length(sTAU)]) 
         body(fam1)[[nopar+2]][[2]]$dldt  <- fam$dldt
           # d2ldmdd 1
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdd)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldmdd) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldmdd  <- fam$d2ldmdd
           # d2ldmdv 2
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldmdv) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldmdv  <- fam$d2ldmdv
           # d2ldmdt 3 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdt)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldmdt) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldmdt  <- fam$d2ldmdt
           # d2ldddv 4
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldddv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldddv) <- parse(text=yval)[[1]] 
         body(fam1)[[nopar+2]][[2]]$d2ldddv  <- fam$d2ldddv
           # d2ldddt 5 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldddt))) 
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldmdt) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldmdt  <- fam$d2ldmdt
           #d2ldvdt 6 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldvdt)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldvdt) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldvdt  <- fam$d2ldvdt
           #d2ldm2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldm2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldm2) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldm2  <- fam$d2ldm2
           #d2ldd2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldd2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
          body(fam$d2ldd2) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldd2  <- fam$d2ldd2
           #d2ldv2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldv2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldv2) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldv2  <- fam$d2ldv2
           #d2ldt2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldt2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
          body(fam$d2ldt2) <- parse(text=yval)[[1]]
          body(fam1)[[nopar+2]][[2]]$d2ldt2  <- fam$d2ldt2
          # residuals 
      cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres[[1]])) 
        # sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)
          body(fam1)[[nopar+2]][[2]]$rqres  <- fam$rqres
           # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial[[1]]))
           #inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
           body(fam1)[[nopar+2]][[2]]$mu.initial  <- fam$mu.initial
          # initial sigma
           inisigma <- gsub("y", "y[,1]",  deparse(fam$sigma.initial[[1]]))
           #inisigma <- gsub("expression", "",  inisigma)
           fam$sigma.initial <- parse(text=inisigma) 
           body(fam1)[[nopar+2]][[2]]$sigma.initial  <- fam$sigma.initial
           # initial nu
           ininu <- gsub("y", "y[,1]",  deparse(fam$nu.initial[[1]]))
           #ininu <- gsub("expression", "",  ininu)
           fam$nu.initial <- parse(text=ininu)
           body(fam1)[[nopar+2]][[2]]$nu.initial  <- fam$nu.initial
          # initial tau
           initau <- gsub("y", "y[,1]",  deparse(fam$tau.initial[[1]]))
           #initau <- gsub("expression", "",  initau)
           fam$tau.initial <- parse(text=initau)
           body(fam1)[[nopar+2]][[2]]$tau.initial  <- fam$tau.initial
          #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           body(fam$y.valid) <- parse(text=yval)[[1]]
           body(fam1)[[nopar+2]][[2]]$y.valid  <- fam$y.valid
           })
##          nfam <- function() fam
## formals(nfam) <- formals(fam1) 
 fam1
}
################################################################################
################################################################################
################################################################################
################################################################################