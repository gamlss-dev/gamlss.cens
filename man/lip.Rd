\name{lip}
\alias{lip}
\docType{data}
\title{Data for lip}
\description{
  The data set used in this package are collected by Dr Peggy Braun (University of Leipzig) and 
  passed on to use by professor Jane Sutherland of London Metropolitan University. 

  It consists of experimental  enzymology  results from a research project 
  which attempted to develop a generic food spoilage model.
  
  The data set contains a column called NAMES, which shows the experiment name, 
  three columns with values of the environmental conditions: temperature (\code{Tem}),
   \code{pH} and water activity (\code{aw}),  and the rest of the columns 
   contains the activity of the cocktails, observed at certain days. 

  The researchers recorded the activity of proteases and lipases in each cocktail and 
  were interested in predicting the time when the activity started given the  
  environmental conditions. 
  The activity is a positive integer and enzymes are considered inactive 
  when activity=0. 

}
\usage{data(lip)}
\format{
  A data frame with 120 observations on the following 14 variables.
  \describe{
    \item{\code{name}}{a factor with levels the different experiment} 
    \item{\code{Tem}}{a numeric vector showing the temperature }
    \item{\code{pH}}{a numeric vector PH}
    \item{\code{aw}}{a numeric vector water activity}
    \item{\code{X0.d}}{a numeric vector if enzyme reacted at day 0}
    \item{\code{X1.d}}{a numeric vector if enzyme reacted at day 1}
    \item{\code{X2.d}}{a numeric vector if enzyme reacted at day 2}
    \item{\code{X4.d}}{a numeric vector if enzyme reacted at days 3 or 4}
    \item{\code{X11.d}}{a numeric vector if enzyme reacted at days 5 to 11}
    \item{\code{X18d.}}{a numeric vector if enzyme reacted at days 12 to q18}
    \item{\code{X25.d}}{a numeric vector if enzyme reacted at days 19 to 25}
    \item{\code{X32.d}}{a numeric vector if enzyme reacted at days 26 to 32}
    \item{\code{X39.d}}{a numeric vector if enzyme reacted at days 33 to 39}
    \item{\code{y}}{a matrix with 3 columns: this is a \code{Surv()} 
    object indicating the start the finish and censored  indicator as defined in 
    function \code{Surv()} of \pkg{survival}.}
  }
}

\source{
Prof.  Jane Sutherland, London Metropolitan University  
}
\examples{
data(lip)
with(lip, y)
}
\keyword{datasets}
