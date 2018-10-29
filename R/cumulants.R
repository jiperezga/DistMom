#' @title Cumulant function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical cumulants of continuous or discrete probability distribution function.
#' @param k order of the moment of interest
#' @param dist density or mass name for the distribution. The created density or mass functions must have a name of the form \code{dxxx}. To understand its use see details and examples.
#' @param param are the parameters of the distribution. The name of each parameter must be specified. To understand its use see examples.
#' @param domain defines the domain of the distribution function. The type of domain of distribution to be tried see details.
#' @details The \code{cumulants} function supports probability distribution functions of a large number of libraries.
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
#' @return \code{cumulants} gives the theorical k-th cumulant of any continuous or discrete probability distribution function. 
#' @note Many continuous distributions support \code{domain = "realline"} even though they are not defined from -\eqn{\infty} to \eqn{\infty} because of their programming. 
#' @note In the same way, many discrete distributions support \code{domain = "counts"} even though they are not defined from \eqn{0} to \eqn{\infty} or \eqn{1} to \eqn{\infty} because of their programming.
#' @note It is recommended to try initially with this argument.
#' @note Discrete distributions require the existence of the quantile function, of the form qxxx.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @seealso \code{\link{moments}}, \code{\link{skew}}, \code{\link{kurt}}.
#' @examples
#' # Let's try first with distributions of the library stats
#' cumulants(k = 1:4, dist = "dchisq", param = c(df = 3), domain = "realplus")
#' # or
#' cumulants(k = 1:4, dist = "chisq", param = c(df = 4), domain = "realplus")
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # The name of the created density functions must have a name
#' # of the form dxxx. Also, how does it not have parameters
#' # then \codeparam = NULL
#' dmyfunction <- function(x) x^3/4 
#' # so that it integrates to 1, x must be between 0 to 2.
#' cumulants(k = 1:4, dist = "dmyfunction", param = NULL, domain = c(0, 2))
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try distributions from other libraries
#' if(!require("extraDistr")) install.packages("extraDistr") # to install the package
#' # The same result is obtained with the diferent domain (see 'Note')
#' cumulants(k = 1:2, dist = "dpareto", param = c(a = 3, b = 7),
#'         domain = "realline")
#' # or
#' cumulants(k = 1:2, dist = "dpareto", param = c(a = 3, b = 7),
#'         domain = c(7, Inf))
#' # In this case, no moments are calculated for k> 2, because the
#' # parameter of the pareto distribution is a = 3, and
#' # therefore, the moments are defined for E (X ^ k) < a.
#' # Read about pareto distribution for more information.
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try distributions from other libraries
#' if(!require("gamlss.dist")) install.packages("gamlss.dist") # to install the package
#' cumulants(k = 3, dist = "PE", param = c(mu = -25, sigma = 7, nu = 4),
#'         domain = "realline") 
#' cumulants(k = 1:3, dist = "BEOI", param = c(mu = 0.3, sigma = 7, nu = 0.3),
#'         domain = "real0to1") 
#' cumulants(k = 1:6, dist = "BCT", param = c(mu = 12, sigma = 0.2, nu = 3, tau = 5),
#'         domain = "realplus")
#' cumulants(k = 2, dist = "PE", param = c(mu = 24, sigma = 18, nu = 2),
#'         domain = "realline")
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try with a discrete counting distribution
#' if(!require("gamlss.dist")) install.packages("gamlss.dist") # to install the package
#' cumulants(k = 1:4, dist = "DEL", param = c(mu = 2, sigma = 3, nu = 0.5),
#'         domain = "counts")
#' cumulants(k = 1:4, dist = "NBF", param = c(mu = 4, sigma = 3, nu = 2),
#'         domain = "counts")
#'         
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try with a discrete binomial type distribution
#' cumulants(k = 1:4, dist = "binom", param = c(size = 15, prob = 0.3),
#'         domain = "binom")
#' cumulants(k = 1:4, dist = "dhyper", param = c(m = 10, n = 7, k = 8),
#'         domain = "binom")
#' @export
cumulants <- function(k, dist, param, domain){
  if(k < 1) stop("the value k must be greater than or equal to 1")
  if(k != as.integer(k)) stop("k must be an integer")
  
  mom <- moments(k = 1:k, dist = dist, param = param, domain = domain)
  cum <- mom[1]
  names(cum) <- paste0("cum_", 1)
  if(k == 1) return(cum)
  
  i <- 2
  while(i <= k){
    comb <- choose(n = i - 1, k = 0:(i-2))
    j <- 1
    aux <- numeric()
    while(j <= length(comb)){
      aux[j] <- comb[j] * cum[j] * mom[i - j]
      j <- j + 1
    }
    cum[i] <- mom[i] - sum(aux)
    names(cum)[i] <- paste0("cum_", i)
    i = i + 1
  }
  
  return(cum[k])
}
cumulants <- Vectorize(cumulants, vectorize.args = "k")