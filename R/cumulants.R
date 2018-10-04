#' @title Moment function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical raw moments, central moments, absolute moments or absolute central moments of any continuous distribution.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' # Enter example D: D: D:
#' @rdname cumulants
#' @export

cumulants <- function(k, dist, param, domain = "realline", rel.tol = 1e-6){
  if(k < 1) stop("the value k must be greater than or equal to 1")
  if(k != as.integer(k)) stop("k must be an integer")
  
  mom <- moments(k = 1:k, dist = dist, param = param, domain = domain, rel.tol = rel.tol)
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