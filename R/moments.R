#' @title Moment function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical raw, central, absolute or absolute central moments of any continuous probability distribution function.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' Enter example D: D: D:
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
