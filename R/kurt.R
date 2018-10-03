#' @title Kurtosis function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical kurtosis or excess kurtosis of any continuous distribution.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' Enter example D: D: D:
#' @rdname moment
#' @export

kurt <- function(excess = TRUE, dist, param, lowerDomain = -Inf, upperDomain = Inf, rel.tol = 1e-6){
  CM4 <- moment(k = 4, dist = dist, param = param, lowerDomain = lowerDomain,
                        upperDomain = upperDomain, central = TRUE, rel.tol = rel.tol)
  
  CM2 <- moment(k = 2, dist = dist, param = param, lowerDomain = lowerDomain,
                upperDomain = upperDomain, central = TRUE, rel.tol = rel.tol)
  if(excess) kurt <- CM3/CM2^(2) - 3 else kurt <- CM3/CM2^(2)
  
  return(kurt)
}

