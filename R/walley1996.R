#'
#'
#' @title Means and Standard Deviations of $\theta_R$
#' @description Formulae in Section 3
#' @export
iddm.inf <- function(nj, s=1, N, tj=NA){
  stopifnot(s >= 1)

  # point probability p.10, 0 <= tj <= 1
  p <- (nj+s*tj)/(N+s)
  # upper bound of expectation (probability) p.10
  p.u <- (nj+s)/(N+s) # limit as tj -> 1
  # lower bound of expectation (probability) p.10
  p.l <- (nj)/(N+s) # limit as tj -> 0
  # p.delta <- p.u - p.l
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

# 4.2. Means and standard deviations of theta_R
iddm.inf(nj=1, N=6, s=2) # pass
iddm.inf(nj=1, N=6, s=1) # pass

f <- function(theta){
  v <- 105*theta^2*(1-theta)^4
  return(v)
}

a1 <- a*(1-a)^5
b1 <- b*(1-b)^5
a1 == b1

library(pscl)
betaHPD(alpha=2, beta=6, p=0.95, plot=TRUE)

