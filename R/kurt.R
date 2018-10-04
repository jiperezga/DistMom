#' @title Kurtosis function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule theoretical kurtosis or excess kurtosis of any continuous distribution.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' # Enter example D: D: D:
#' @rdname kurt
#' @export

kurt <- function(excess = TRUE, dist, param, domain = "realline", rel.tol = 1e-6){
  CM4 <- momentss(k = 4, dist = dist, param = param, domain = domain, central = TRUE, rel.tol = rel.tol)
  
  CM2 <- moments(k = 2, dist = dist, param = param, domain = domain, central = TRUE, rel.tol = rel.tol)
  if(excess) kurt <- CM3/CM2^(2) - 3 else kurt <- CM3/CM2^(2)
  
  return(kurt)
}

