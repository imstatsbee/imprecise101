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

# library(pscl)
# x <- betaHPD(alpha=3, beta=5, p=0.95, plot=TRUE)

alpha <- 3
beta <- 5

a0 <- .Machine$double.eps
b0 <- 1 - a0

for(i in 1:10){

  b <- b0
  fn1 <- function(a){ # let fix b
    y1 <- a*(1-a)^5
    y2 <- b*(1-b)^5
    v <- abs(y1-y2)
    return(v)
  }
  op1 <- optimize(f=fn1, lower=a0, upper=qbeta(1-0.95, 3, 5))

  a <- op1$minimum
  fn2 <- function(b){ # let fix a
    fn0 <- function(theta) 105*theta^2*(1-theta)^4
    v <- integrate(f=fn0, lower=a, upper=b)$value - 0.95
    return(abs(v))
  }
  op2 <- optimize(f=fn2, lower=a, upper=qbeta(0.95, 3, 5))

  b <- op2$minimum

  b0 <- b
  a0 <- a
}

c(a,b)


# Eqn 2. P(D(l)|n) = \int_a^b 105 t^2 (1-t)^4 dt = 0.95
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


