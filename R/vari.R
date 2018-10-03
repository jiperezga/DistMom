#' @title Variance function
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @description Calcule variance of any continuous distribution.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' Enter example D: D: D:
#' @rdname moment
#' @export

vari <- function(dist, param, lowerDomain = -Inf, upperDomain = Inf, rel.tol = 1e-6){
  vari <- moment(k = 2, dist = dist, param = param, lowerDomain = lowerDomain,
                upperDomain = upperDomain, central = TRUE, rel.tol = rel.tol)
  return(vari)
}

