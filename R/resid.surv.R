################################################################################
################################################################################
################################################################################
################################################################################
# resid.surv() function 
# to take residualas from a fitted GAMLSS survival model
# and plot its Kaplan Meir surival function with Normal dustribution 
# survival added  
################################################################################
################################################################################
################################################################################
################################################################################
resid.surv <- function(model)
{
gamlss.bi.list <- .binom
if (!is.null(model)) family <- model[["family"]][1]
    fname <- if (is.name(family)) as.character(family)
        else if (is.character(family)) family
        else if (is.call(family)) as.character(family[[1]])
        else if (is.function(family)) deparse(substitute(family))
        else if (is(family, "gamlss.family"))  family$family[1]
        else stop("the family must be a character or a gamlss.family name")
        fam1 <- eval(parse(text=fname)) # the family to output
         fam <- as.gamlss.family(family) # this is created so I can get things
      dorfun <- paste("p",fname,sep="") # say dNO
       nopar <- fam$nopar # or fam1$nopar
        type <- fam$type  
#-------------------------------------------------------------------------------
# get the information of the survival response   

if (!is(model$y, "Surv")) stop("the response is not Surv object")
       TYPE <- attr(model$y, "type") 
if (TYPE=="right")
       {
       time <-  model$y[,"time"]
     status <-  model$y[,"status"]
       pfun <- substr(dorfun, 1, nchar(dorfun)-2)
       } 
if (TYPE=="left")
       {
      time <-  model$y[,1]
    status <-  model$y[,"status"]# need checking
      pfun <- substr(dorfun, 1, nchar(dorfun)-2)
       }  
if (TYPE=="interval")
       {
     time <- model$y[,1]
    time2 <- model$y[,2]
   status <- model$y[,"status"]
     pfun <- substr(dorfun, 1, nchar(dorfun)-2)
       }   
################################################################################
# get the cdf of the model  to get the PIT
     pit<- switch(nopar, 
          eval(call(pfun, q = time, mu = fitted(model,"mu"))),# 1
          eval(call(pfun, q = time, mu = fitted(model,"mu"), 
                    sigma = fitted(model,"sigma"))),        # 2
          eval(call(pfun, q = time,  mu = fitted(model,"mu"), 
                      sigma = fitted(model,"sigma"),  
                       nu = fitted(model,"nu"))),           # 3                   
          eval(call(pfun, q = time,  mu = fitted(model,"mu"), 
                    sigma = fitted(model,"sigma"),
                       nu = fitted(model,"nu"), 
                      tau = fitted(model,"tau"))))
#  pit <- pGA(colcancer$followup/12, mu=fitted(gamfit), 
#             sigma=fitted(gamfit, "sigma"))
pit <- ifelse(pit==1, 0.999, pit)  
# get the normalised residuals
     z.scores <- qNO(pit)
if (TYPE=="interval")
{
  pit2 <- switch(nopar, 
               eval(call(pfun, q= time2, mu=fitted(model,"mu"))),# 1
               eval(call(pfun, q= time2,  mu=fitted(model,"mu"), 
                         sigma=fitted(model,"sigma"))),        # 2
               eval(call(pfun, q= time2,  mu=fitted(model,"mu"), 
                         sigma=fitted(model,"sigma"),  
                         nu=fitted(model,"nu"))),           # 3                   
               eval(call(pfun, q= time2,  mu=fitted(model,"mu"), 
                         sigma=fitted(model,"sigma"),
                         nu=fitted(model,"nu"), 
                         tau=fitted(model,"tau"))))
  z.scores2 <- qNO(pit2)
}  
  
################################################################################  
#  get  the Kaplan Meier survival from the normalised residuals
if (TYPE=="right")
     {
  kmpit <- survival::survfit(Surv(z.scores, status, type="right") ~ 1)       
     } 
if (TYPE=="left")
     {
       kmpit <- survival::survfit(Surv(z.scores, status) ~ 1)  
     } 
if (TYPE=="interval")
     {
  second <- ifelse(status==0, max(z.scores[z.scores!=Inf]), 0)+ 
            ifelse(status==1, z.scores, 0)+
            ifelse(status==2, min(z.scores[z.scores!=-Inf]), 0)+
            ifelse(status==3,z.scores2, 0)
   kmpit <- survival::survfit(Surv(z.scores, second, status, type="interval")                           ~ 1)  
     } 
          
plot(kmpit)
# add the normal distribution survival
  lines(z.scores[order(z.scores)], pNO(z.scores[order(z.scores)], lower.tail = FALSE), col="blue")
# if the normal curve in the middle 
# it does not contradict the original distribution assumption 
# the residual could have come from a censored normal distribution
# this plot can be generalised  to any distribution and in ggplot 
} 