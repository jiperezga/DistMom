#' @title Moment function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical raw, central, absolute or absolute central moments of any continuous probability distribution function.
#' @param k order of the moment of interest
#' @param dist density name for the distribution. The created density functions must have a name of the form \code{dxxx}. To understand its use see details and examples.
#' @param param are the parameters of the distribution. The name of each parameter must be specified. To understand its use see examples.
#' @param domain defines the domain of the distribution function. The type of domain of distribution to be tried see details.
#' @param central logical; if TRUE, the k-th central moments are given as \eqn{E[(X-\mu)^k]}.
#' @param absolute logical; if TRUE, the k-th absolute moments are given as \eqn{E|X|^k]}.
#' @param rel.tol relative accuracy requested. By default it's value is \eqn{1e-6.}.
#' @details The \code{moments} function supports probability distribution functions of a large number of libraries.
#' @details In the \code{dist} argument, you must enter the name of the distribution of interest, for example, you can enter \code{"gamma"} or \code{"dgamma"}, both will produce the same result.
#' @details If \eqn{f(x)} has no parameters, then do \code{param = NULL}.
#' @details The following are the different \code{domain} argument: \itemize{
#' \item \code{realline}: the function is defined between -\eqn{\infty} and \eqn{\infty}. 
#' \item \code{realplus}: the function is defined between 0 and \eqn{\infty}.
#' \item \code{real0to1}: the function is defined between 0 and 1.
#' \item \code{real-1to1}: the function is defined between -1 and 1.
#' \item \code{c(lower = a, upper = b)}: the function is defined between \code{a} and \code{b}.
#' }
#' @details If \code{central = TRUE} and \code{absolute = TRUE} are selected, the k-th central absolute moments is calculated and given as \eqn{E|(X-\mu)^k|}.
#' @details the value of \code{rel.tol} can be modified if there are convergence problems or if it is believed that the obtained precision is not enough.
#' @details \code{rel.tol} cannot be less than max(50*.Machine$double.eps, 0.5e-28)
#' @return \code{moments} gives the theorical k-th raw, central, absolute or absolute central moments of any continuous probability distribution function.
#' @note many distributions support \code{domain = "realline"} even though they are not defined from -\eqn{\infty} to \eqn{\infty} because of their programming. It is recommended to try initially with this argument.
#' @importFrom extraDistr dpareto
#' @importFrom gamlss.dist dPE
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @seealso \code{\link{cumulants}}, \code{\link{skew}}, \code{\link{kurt}}.
#' @examples
#' # Let's try first with distributions of the library stats
#' moments(k = 1:4, dist = "dchisq", param = c(df = 3), domain = "realplus")
#' # or
#' moments(k = 1:4, dist = "chisq", param = c(df = 4), domain = "realplus")
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # the name of the created density functions must have a name
#' # of the form dxxx. Also, how does it not have parameters
#' # then \code{param = NULL}
#' dmyfunction <- function(x) x^3/4 
#' # so that it integrates to 1, x must be between 0 to 2.
#' moments(k = 1:4, dist = "dmyfunction", param = NULL, domain = c(0, 2))
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try distributions from other libraries
#' if(!require("extraDistr")) install.packages("extraDistr") # to install the package
#' # The same result is obtained with the diferent domain (see 'Note')
#' moments(k = 1:2, dist = "dpareto", param = c(a = 3, b = 7),
#'         domain = "realline")
#' # or
#' moments(k = 1:2, dist = "dpareto", param = c(a = 3, b = 7),
#'         domain = c(7, Inf))
#' # In this case, no moments are calculated for k> 2, because the
#' # parameter of the pareto distribution is \code {a = 3}, and
#' # therefore, the moments are defined for \eqn{E (X ^ k) < a}.
#' # Read about pareto distribution for more information.
#' 
#' #---------------------------------------------------------------------------------------
#' 
#' # Let's try distributions from other libraries to calculated central
#' # and absolute moments
#' if(!require("gamlss.dist")) install.packages("gamlss.dist") # to install the package
#' moments(k = 3, dist = "PE", param = c(mu = -25, sigma = 7, nu = 4),
#'         central = TRUE) 
#' moments(k = 3, dist = "PE", param = c(mu = -25, sigma = 7, nu = 4),
#'         absolute = TRUE)
#' moments(k = 3, dist = "PE", param = c(mu = -25, sigma = 7, nu = 4),
#'         central = TRUE, absolute = TRUE)
#' 
#' @export
moments <- function(k, dist, param, domain = "realline", central = FALSE, absolute = FALSE, rel.tol = 1e-6){
  if(k < 1) stop("the value k must be greater than or equal to 1")
  if(k != as.integer(k)) stop("k must be an integer")
  dist <- ifelse(exists(paste0("d", dist)), paste0("d", dist),
                 ifelse(exists(dist), dist, stop("The probability distribution entered does not exist or the library to which it belongs is not loaded")))
  if(is.character(domain)){
    if(domain == "realline") {lowerDomain = -Inf; upperDomain = Inf}
    if(domain == "realplus") {lowerDomain = 0; upperDomain = Inf}
    if(domain == "real0to1") {lowerDomain = 0; upperDomain = 1}
    if(domain == "real-1to1") {lowerDomain = -1; upperDomain = 1}
  } else {
    if(length(domain) == 2) {lowerDomain = min(domain); upperDomain = max(domain)}
  }
  
  if(!central & !absolute){
    mom <- if(!is.null(names(param))){
      if(all(names(param) %in% names(formals(dist)))){
        integrate(function(x) x^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
      }  else {
        dif1 <- setdiff(names(param), names(formals(dist)))
        dif2 <- setdiff(names(formals(dist)), names(param))
        dif2 <- dif2[!dif2 == "x" & !dif2 == "log"]
        stop(paste0("The name of the parameters entered ", '"',  paste(as.character(dif1), collapse = '" or "'), '"', " does not match the name of the parameters ", '"', paste(as.character(dif2), collapse = '" or "'), '"', " of the probability distribution. \n See ", '"', dist, '"', " help for more info."))
      }
    } else {
      if(is.null(param)){
        integrate(function(x) x^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
      } else {
        dif3 <- setdiff(names(formals(dist)), names(param))
        dif3 <- dif3[!dif3 == "x" & !dif3 == "log"]
        stop(paste0("The parameters of the distribution ", '"', dist, '"', " are ", '"', paste(as.character(dif3), collapse = '" or "'), '"', ". \n Please name what value belongs to which parameter. \n See the help of the ", '"', dist, '"', " function in ", getAnywhere(dist)$where[1], " for more information."))
      }
    }
    res <- mom$value
    names(res) <- paste0("E[X^", k, "]")
  }
  
  if(central & !absolute){
    if(k == 1){
      mom <- if(!is.null(names(param))){
        if(all(names(param) %in% names(formals(dist)))){
          integrate(function(x) x^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        }  else {
          dif1 <- setdiff(names(param), names(formals(dist)))
          dif2 <- setdiff(names(formals(dist)), names(param))
          dif2 <- dif2[!dif2 == "x" & !dif2 == "log"]
          stop(paste0("The name of the parameters entered ", '"',  paste(as.character(dif1), collapse = '" or "'), '"', " does not match the name of the parameters ", '"', paste(as.character(dif2), collapse = '" or "'), '"', " of the probability distribution. \n See ", '"', dist, '"', " help for more info."))
        }
      } else {
        if(is.null(param)){
          integrate(function(x) x^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        } else {
          dif3 <- setdiff(names(formals(dist)), names(param))
          dif3 <- dif3[!dif3 == "x" & !dif3 == "log"]
          stop(paste0("The parameters of the distribution ", '"', dist, '"', " are ", '"', paste(as.character(dif3), collapse = '" or "'), '"', ". \n Please name what value belongs to which parameter. \n See the help of the ", '"', dist, '"', " function in ", getAnywhere(dist)$where[1], " for more information."))
        }
      }
    } else {
      mom <- if(!is.null(names(param))){
        if(all(names(param) %in% names(formals(dist)))){
          med <- integrate(function(x) x * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
          integrate(function(x) (x - med$value)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        }  else {
          dif1 <- setdiff(names(param), names(formals(dist)))
          dif2 <- setdiff(names(formals(dist)), names(param))
          dif2 <- dif2[!dif2 == "x" & !dif2 == "log"]
          stop(paste0("The name of the parameters entered ", '"',  paste(as.character(dif1), collapse = '" or "'), '"', " does not match the name of the parameters ", '"', paste(as.character(dif2), collapse = '" or "'), '"', " of the probability distribution. \n See ", '"', dist, '"', " help for more info."))
        }
      } else {
        if(is.null(param)){
          med <- integrate(function(x) x * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
          integrate(function(x) (x - med$value)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        } else {
          dif3 <- setdiff(names(formals(dist)), names(param))
          dif3 <- dif3[!dif3 == "x" & !dif3 == "log"]
          stop(paste0("The parameters of the distribution ", '"', dist, '"', " are ", '"', paste(as.character(dif3), collapse = '" or "'), '"', ". \n Please name what value belongs to which parameter. \n See the help of the ", '"', dist, '"', " function in ", getAnywhere(dist)$where[1], " for more information."))
        }
      }
    }
    res <- mom$value
    names(res) <- paste0("E[(X-EX)^", k, "]")
  }
  
  if(!central & absolute){
    mom <- if(!is.null(names(param))){
      if(all(names(param) %in% names(formals(dist)))){
        integrate(function(x) abs(x)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
      }  else {
        dif1 <- setdiff(names(param), names(formals(dist)))
        dif2 <- setdiff(names(formals(dist)), names(param))
        dif2 <- dif2[!dif2 == "x" & !dif2 == "log"]
        stop(paste0("The name of the parameters entered ", '"',  paste(as.character(dif1), collapse = '" or "'), '"', " does not match the name of the parameters ", '"', paste(as.character(dif2), collapse = '" or "'), '"', " of the probability distribution. \n See ", '"', dist, '"', " help for more info."))
      }
    } else {
      if(is.null(param)){
        integrate(function(x) abs(x)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
      } else {
        dif3 <- setdiff(names(formals(dist)), names(param))
        dif3 <- dif3[!dif3 == "x" & !dif3 == "log"]
        stop(paste0("The parameters of the distribution ", '"', dist, '"', " are ", '"', paste(as.character(dif3), collapse = '" or "'), '"', ". \n Please name what value belongs to which parameter. \n See the help of the ", '"', dist, '"', " function in ", getAnywhere(dist)$where[1], " for more information."))
      }
    }
    res <- mom$value
    names(res) <- paste0("E|X^", k, "|")
  }
  
  if(central & absolute){
    if(k == 1){
      mom <- if(!is.null(names(param))){
        if(all(names(param) %in% names(formals(dist)))){
          integrate(function(x) abs(x)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        }  else {
          dif1 <- setdiff(names(param), names(formals(dist)))
          dif2 <- setdiff(names(formals(dist)), names(param))
          dif2 <- dif2[!dif2 == "x" & !dif2 == "log"]
          stop(paste0("The name of the parameters entered ", '"',  paste(as.character(dif1), collapse = '" or "'), '"', " does not match the name of the parameters ", '"', paste(as.character(dif2), collapse = '" or "'), '"', " of the probability distribution. \n See ", '"', dist, '"', " help for more info."))
        }
      } else {
        if(is.null(param)){
          integrate(function(x) abs(x)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        } else {
          dif3 <- setdiff(names(formals(dist)), names(param))
          dif3 <- dif3[!dif3 == "x" & !dif3 == "log"]
          stop(paste0("The parameters of the distribution ", '"', dist, '"', " are ", '"', paste(as.character(dif3), collapse = '" or "'), '"', ". \n Please name what value belongs to which parameter. \n See the help of the ", '"', dist, '"', " function in ", getAnywhere(dist)$where[1], " for more information."))
        }
      }
    } else {
      mom <- if(!is.null(names(param))){
        if(all(names(param) %in% names(formals(dist)))){
          med <- integrate(function(x) x * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
          integrate(function(x) abs(x - med$value)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ", ", paste0(sapply(X = 1:length(param), FUN = function(i) paste(names(param)[i], " = ", param[i])), collapse = ', '), ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        }  else {
          dif1 <- setdiff(names(param), names(formals(dist)))
          dif2 <- setdiff(names(formals(dist)), names(param))
          dif2 <- dif2[!dif2 == "x" & !dif2 == "log"]
          stop(paste0("The name of the parameters entered ", '"',  paste(as.character(dif1), collapse = '" or "'), '"', " does not match the name of the parameters ", '"', paste(as.character(dif2), collapse = '" or "'), '"', " of the probability distribution. \n See ", '"', dist, '"', " help for more info."))
        }
      } else {
        if(is.null(param)){
          med <- integrate(function(x) x * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
          integrate(function(x) abs(x - med$value)^k * eval(parse(text = paste0(dist ,"(", names(formals(paste0(dist)))[1], " = ", "x", ")"))), lower = lowerDomain, upper = upperDomain, rel.tol = rel.tol)
        } else {
          dif3 <- setdiff(names(formals(dist)), names(param))
          dif3 <- dif3[!dif3 == "x" & !dif3 == "log"]
          stop(paste0("The parameters of the distribution ", '"', dist, '"', " are ", '"', paste(as.character(dif3), collapse = '" or "'), '"', ". \n Please name what value belongs to which parameter. \n See the help of the ", '"', dist, '"', " function in ", getAnywhere(dist)$where[1], " for more information."))
        }
      }
    }
    res <- mom$value
    names(res) <- paste0("E|(X-EX)^", k, "|")
  }
  
  return(res)
}
moments <- Vectorize(moments, vectorize.args = "k")
