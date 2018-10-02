#' @title The McCullagh Distribution
#' @author Jorge Iván Pérez, \email{jivan.perez@udea.edu.co}
#' @author Juan Carlos Correa, \email{jccorrea@unal.edu.co}
#' @description Density, distribution function, quantile function and random generation for the McCullagh distribution with parameters \code{theta} and \code{nu}.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param theta parameter. Must be between -1 and 1, inclusive.
#' @param nu parameter. Must be greater than -0.5.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#' @details If \code{theta} or \code{nu} are not specified they assume the default values of 0 and 4, respectively.
#' @details The McCullagh distribution with parameters \eqn{\theta} (\code{theta}) and \eqn{\nu} (\code{nu}) has density
#' @details \deqn{f(x) = (1 - x^2)^(\nu - 0.5) / [ (1 - 2x\theta + \theta^2)^v * B(\nu + 0.5, 0.5) ]}
#' @details for -1 ≤ \eqn{x} ≤ 1, -1 < \eqn{\theta} < 1, \eqn{\nu} > - 0.5. The mean, variance, skewness and kurtosis are
#' @details \eqn{E(X) = \nu\theta / (\nu + 1)}
#' @details \eqn{Var(X) = (1 - [\nu * (\nu-1) * \theta^2] / [ (\nu + 1) * (\nu + 2) ] ) / [2 * (\nu + 1) ]}
#' @details \eqn{\kappa_3(X) = (\nu / [2 * (\nu + 1)^2 * (\nu + 2) ] ) * [ ( [ (\nu - 1) * (3\nu - 1) * \theta^3] / [ (\nu - 1) * (\nu + 3) ] ) - 3\theta ]}
#' @details \eqn{\kappa_4(X) = (-3 / [4 * (\nu + 1)^2 * (\nu + 2) ] ) * (1 - [4\nu * (3\nu - 1) * \theta^2] / [ (\nu + 1) * (\nu + 3) ] + ( [\nu (\nu - 1) * (11\nu^3 + 16\nu^2 - 17\nu + 2) * \theta^4] / [ (\nu + 1)^2 * (\nu + 2) * (\nu + 3) * (\nu + 4) ] ) )}
#' @note Because the cumulative distribution function Mccullagh does not have a closed form, the values of the functions \code{pmcullagh}, \code{qmcullagh} and \code{rmcullagh} are obtained through the solution of systems of nonlinear equations, which could delay the algorithm.
#' @note We will try to improve the times, especially in the function \code{rmccullagh}.
#' @return \code{dnorm} gives the density, \code{pnorm} gives the distribution function, \code{qnorm} gives the quantile function, and \code{rnorm} generates random deviates.
#' @return The length of the result is determined by \code{n} for \code{rnorm}, and is the maximum of the lengths of the numerical arguments for the other functions.
#' @return The numerical arguments other than \code{n} are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' @import GoFKernel
#' @references McCullagh, P. (1989). Some statistical properties of a family of continuous univariate distributions. Journal of the American Statistical Association, 84(405), 125-129.
#' @references Seshadri, V. (1991). A family of distributions related to the McCullagh family. Statistics & probability letters, 12(5), 373-378.
#' @seealso \code{\link{Distributions}} for other standard distributions.
#' @examples
#' dmccullagh(x = 0.3, theta = -0.8, nu = 2.3)
#'
#' system.time(sim <- rmccullagh(n = 1000, theta = 0.2, nu = -0.15)) # it could take several seconds
#' plot(ecdf(sim))
#' curve(pmccullagh(q = x, theta = 0.2, nu = -0.15), add = TRUE, col = "red", lwd = 2)
#'
#' qmccullagh(p = 0.2, theta = 0, nu = 0.41, lower.tail = FALSE)
#' @rdname mccullagh
#' @export
dmccullagh <- function(x, theta = 0, nu = 4, log = FALSE){
  if(x < -1 | x > 1) stop(return(0))
  if(theta <= -1 | theta >= 1) stop("theta must be between -1 and 1")
  if(nu <= -0.5) stop("nu must be greater than or equal to -0.5")
  d <- ((1 - (x^2))^(nu - 0.5)) / (((1 - (2 * theta * x) + theta^2)^nu) * beta(a = nu + 0.5, b = 0.5))
  if(log) d <- log(d)
  return(d)
}
dmccullagh <- Vectorize(dmccullagh)

#' @rdname mccullagh
#' @export
pmccullagh <- function(q, theta = 0, nu = 4, lower.tail = TRUE, log.p = FALSE){
  if(q <= -1) stop(return(0))
  if(q >= 1) stop(return(1))
  p <- integrate(f = dmccullagh, lower = -1, upper = q, theta = theta, nu = nu, rel.tol = 1e-4)$value
  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}
pmccullagh <- Vectorize(pmccullagh)

#' @rdname mccullagh
#' @export
qmccullagh <- function(p, theta = 0, nu = 4, lower.tail = TRUE, log.p = FALSE){
  if(!log.p & (p < 0 | p > 1)) {print(NaN); stop("NaNs produced")}
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  f <- function(x) pmccullagh(x, theta = theta, nu = nu)
  q.inv <- inverse(f = f, lower = -1, upper = 1)
  q <- q.inv(p)
  return(q)
}
qmccullagh <- Vectorize(qmccullagh)

#' @rdname mccullagh
#' @export
rmccullagh <- function(n, theta = 0, nu = 4){
  u <- runif(n)
  r <- qmccullagh(p = u, theta = theta, nu = nu)
  return(r)
}
rmccullagh <- Vectorize(rmccullagh)
