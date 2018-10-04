#' @title Skewness function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical skewness of any continuous distribution.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' # Enter example D: D: D:
#' @rdname skew
#' @export

skew <- function(dist, param, domain = "realline", rel.tol = 1e-6){
  CM3 <- moments(k = 3, dist = dist, param = param, domain = domain,
                central = TRUE, rel.tol = rel.tol)
  
  CM2 <- moments(k = 2, dist = dist, param = param, domain = domain,
                 central = TRUE, rel.tol = rel.tol)
  
  skew <- CM3/CM2^(3/2)
  
  return(skew)
}

