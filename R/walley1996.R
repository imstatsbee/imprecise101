#'
#' @rdname idm
#' @title imprecise Dirichlet model
#' @description lower and upper posterior probabilities
#' @param nj number of observations in the j-th category
#' @param s learning parameter
#' @param N a total number of drawings
#' @param tj proability
#' @examples
#' idm(nj=1, N=6, s=2)
#' idm(nj=1, N=6, s=1)
#' @export
idm <- function(nj, s=1, N, tj=NA){
  stopifnot(s >= 1)

  # section 2.3
  ## point probability p.10, 0 <= tj <= 1
  p <- tj.star <- (nj+s*tj)/(N+s)

  ## upper bound (limiting as tj -> 1)
  p.u <- (nj+s)/(N+s)

  ## lower bound (limiting as tj -> 0)
  p.l <- (nj)/(N+s)

  p.delta <- s/(N+s)

  # uppder bound of Variance (p.17)
  n.u <- nj*(N-nj) + 1/4*(N+s)*(s+1)^2
  d.u <- (N+s)*(N+s+1)^2
  v.u <- n.u/d.u

  # lower bound of Variance (p.17)
  n.l <- nj*(N-nj) + s*min(c(nj, N-nj))
  d.l <- (N+s)^2*(N+s+1)
  v.l <- n.l/d.l

  robj <- list(p.lower=p.l, p.upper=p.u, v.lower=v.l, v.upper=v.u, s.lower=sqrt(v.l), s.upper=sqrt(v.u), p=p, p.delta=p.delta)
  return(robj)
}


#' @rdname idm
#' @title Find HPD
#' @description Find the highest posterior density interval
#' @param maxiter maximum number of iterations
#' @example
#' x <- hpd(alpha=3, beta=5, p=0.95) # c(0.0031, 0.6587) when s=2
#' # round(x,4); x*(1-x)^5
#' x <- hpd(alpha=2, beta=5, p=0.95) # c(0.0076, 0.5834) when s=1
#' # round(x,4); x*(1-x)^5
#' x <- hpd(alpha=3, beta=5, p=0.9) # c(0.0066, 0.5962) when s=2
#' # round(x,4); x*(1-x)^5
#' x <- hpd(alpha=2, beta=5, p=0.9) # c(0.0150, 0.5141) when s=1
#' # round(x,4); x*(1-x)^5
#' x <- hpd(alpha=3, beta=5, p=0.5) # c(0.0481, 0.3656) when s=2 (strange/error)
#' # round(x,4); x*(1-x)^5
#' x <- hpd(alpha=2, beta=5, p=0.5) # c(0.0761, 0.2958) when s=1
#' # round(x,4); x*(1-x)^5
#' @export
hpd <- function(alpha=3, beta=5, p=0.95, runs=1e1, tolerance=1e-4, maxiter=1e2){

  # objective function 1
  fn1 <- function(a, b){
    y1 <- a*(1-a)^5
    y2 <- b*(1-b)^5
    dif <- abs(y1-y2)
    return(dif)
  }

  # objective function 2
  fn2 <- function(b, a){ # let fix a
    # fn0 <- function(theta) 105*theta^2*(1-theta)^4
    fn0 <- function(theta) dbeta(theta, alpha, beta) # 105*theta^2*(1-theta)^4
    dif <- abs(integrate(f=fn0, lower=a, upper=b)$value - p)
    return(dif)
  }

  # initialization
  a <- 1e-8
  b <- 1-a
  dif <- 1e-2
  niter <- 0

  while( dif > tolerance){

    # for fixed b, searching for a
    op1 <- optimize(f=fn1, b=b, lower=0, upper=qbeta(1-p, alpha, beta))
    a <- op1$minimum

    # for fixed a, searching for b
    op2 <- optimize(f=fn2, a=a, lower=qbeta(p, alpha, beta), upper=1)
    b <- op2$minimum

    dif <- fn1(a=a, b=b)
    print(c(a=a,b=b, dif=dif))
    niter <- niter + 1
    if(niter == maxiter) break
  }

  robj <- c(a=a, b=b)
  return(robj)
}




#' @rdname idm
#
# Eqn 2. P(D(l)|n) = \int_a^b 105 t^2 (1-t)^4 dt = 0.95
# computing the 95% credible interval (highest posterior density interval)
#
# library(pscl)
# x <- betaHPD(alpha=3, beta=5, p=0.95, plot=TRUE)

alpha <- 3
beta <- 5
p <- 0.95

fn2 <- function(x0, alpha, beta, p) {
  y0 <- dbeta(x0, alpha, beta)
  p0 <- pbeta(x0, alpha, beta)
  x1 <- qbeta(p0 + p, alpha, beta)
  y1 <- dbeta(x1, alpha, beta)
  v <- abs(y0 - y1)
  return(v)
}
op2 <- optimize(f = fn2, alpha=alpha, beta=beta, p=p, interval = c(0, qbeta(1 - p, alpha, beta)))
c(a=op2$minimum, b=qbeta(pbeta(op2$minimum, a, b) + p, a, b))
