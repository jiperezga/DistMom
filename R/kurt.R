#' @title Kurtosis function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical kurtosis or excess kurtosis of any continuous or discrete distribution.
#' @param dist density or mass name for the distribution. The created density or mass functions must have a name of the form \code{dxxx}. To understand its use see details and examples.
#' @param param are the parameters of the distribution. The name of each parameter must be specified. To understand its use see examples.
#' @param domain defines the domain of the distribution function. The type of domain of distribution to be tried see details.
#' @param excess logical; if TRUE, the excess kurtosis are calculated. FALSE is the default value.
#' @details The \code{kurt} function supports probability distribution functions of a large number of libraries.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @details In the \code{dist} argument, you must enter the name of the distribution of interest, for example, you can enter \code{"gamma"} or \code{"dgamma"}, both will produce the same result.
#' @details If \eqn{f(x)} has no parameters, then do \code{param = NULL}.
#' @details The following are the different \code{domain} argument: \itemize{
#' \item \code{binom}: for discrete distributions of binomial type.
#' \item \code{counts}: for discrete distributions of counting type.
#' \item \code{realline}: for continuous distributions defined between -\eqn{\infty} and \eqn{\infty}. 
#' \item \code{realplus}: for continuous distributions defined between 0 and \eqn{\infty}.
#' \item \code{real0to1}: for continuous distributions defined between 0 and 1.
#' \item \code{real-1to1}: for continuous distributions defined between -1 and 1.
#' \item \code{c(lower = a, upper = b)}: for continuous distributions defined between \code{a} and \code{b}.
#' }
#' @return \code{kurt} gives the theorical kurtosis or excess kurtosis of any continuous or discrete probability distribution function. 
#' @note Many continuous distributions support \code{domain = "realline"} even though they are not defined from -\eqn{\infty} to \eqn{\infty} because of their programming. 
#' @note In the same way, many discrete distributions support \code{domain = "counts"} even though they are not defined from \eqn{0} to \eqn{\infty} or \eqn{1} to \eqn{\infty} because of their programming.
#' @note It is recommended to try initially with this argument.
#' @note Discrete distributions require the existence of the quantile function, of the form qxxx.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @seealso \code{\link{moments}}, \code{\link{cumulants}}, \code{\link{kurt}}.
#' @examples
#' # Let's try first with distributions of the library stats
#' kurt(dist = "dchisq", param = c(df = 3), domain = "realplus")
#' # or
#' kurt(dist = "chisq", param = c(df = 4), domain = "realplus")
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # The name of the created density functions must have a name
#' # of the form dxxx. Also, how does it not have parameters
#' # then \code{param = NULL}
#' dmyfunction <- function(x) x^3/4 
#' # so that it integrates to 1, x must be between 0 to 2.
#' kurt(dist = "dmyfunction", param = NULL, domain = c(0, 2))
#' kurt(dist = "dmyfunction", param = NULL, domain = c(0, 2), excess = TRUE)
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try distributions from other libraries
#' if(!require("extraDistr")) install.packages("extraDistr") # to install the package
#' # The same result is obtained with the diferent domain (see 'Note')
#' 
#' kurt(dist = "dpareto", param = c(a = 5, b = 10),
#'      domain = "realline")
#'      
#' kurt(dist = "dpareto", param = c(a = 5, b = 10),
#'      domain = "realline", excess = TRUE)    
#'  
#'         
#' # In this case, no moments are calculated for k > 4, because the
#' # parameter of the pareto distribution is \code {a = 5}, and
#' # therefore, the moments are defined for \eqn{E (X ^ k) < a}.
#' # Read about pareto distribution for more information.       
#' 
#' kurt(dist = "dbhatt", param = c(mu = 3, sigma = 7),
#'      domain = "realline")
#'  
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try distributions from other libraries
#' if(!require("gamlss.dist")) install.packages("gamlss.dist") # to install the package
#' kurt(dist = "PE", param = c(mu = -25, sigma = 7, nu = 4),
#'      domain = "realline") 
#' kurt(dist = "BEOI", param = c(mu = 0.3, sigma = 7, nu = 0.3),
#'      domain = "real0to1") 
#' kurt(dist = "BCT", param = c(mu = 12, sigma = 0.2, nu = 3,
#'      tau = 5), domain = "realplus")
#' kurt(dist = "SEP2", param = c(mu = 0.5, sigma = 3,
#'      nu = 0, tau = 5), domain = "realline")
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try with a discrete counting distribution
#' if(!require("gamlss.dist")) install.packages("gamlss.dist") # to install the package
#' kurt(dist = "DEL", param = c(mu = 2, sigma = 3, nu = 0.5),
#'      domain = "counts")
#' kurt(dist = "NBF", param = c(mu = 4, sigma = 3, nu = 2),
#'      domain = "counts")
#'         
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try with a discrete binomial type distribution
#' kurt(dist = "binom", param = c(size = 15, prob = 0.3),
#'         domain = "binom")
#' kurt(dist = "dhyper", param = c(m = 10, n = 7, k = 8),
#'         domain = "binom")
#' @export

kurt <- function(excess = FALSE, dist, param, domain){
  CM4 <- tryCatch(expr = moments(k = 4, dist = dist, param = param, domain = domain,
                                 central = TRUE), error = function(e) "error")
  if(CM4 == "error") stop("The asymptotic method does not converge, the value of the fourth moment is very large or the Fourth moment of the distribution does not exist. \n Read about the conditions for the existence of the moments of the distribution of interest")
  
  CM2 <- moments(k = 2, dist = dist, param = param, domain = domain, central = TRUE)
  if(excess) {
    kurt <- CM4/CM2^(2) - 3
    names(kurt) <- "Ex.Kurtosis"
    } else {
      kurt <- CM4/CM2^(2)
      names(kurt) <- "Kurtosis"
    }
  
  return(kurt)
}

