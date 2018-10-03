#' @title Skewness function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical skewness of any continuous distribution.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' Enter example D: D: D:
#' @rdname moment
#' @export

skew <- function(dist, param, lowerDomain = -Inf, upperDomain = Inf, rel.tol = 1e-6){
  CM3 <- moment(k = 3, dist = dist, param = param, lowerDomain = lowerDomain,
                        upperDomain = upperDomain, central = TRUE, rel.tol = rel.tol)
  
  CM2 <- moment(k = 2, dist = dist, param = param, lowerDomain = lowerDomain,
                upperDomain = upperDomain, central = TRUE, rel.tol = rel.tol)
  
  skew <- CM3/CM2^(3/2)
  
  return(skew)
}

