#'
#'
#' @title Means and Standard Deviations of $\theta_R$
#' @description Formulae in Section 3
#' @export
iddm.inf <- function(nj, s=1, N){
  stopifnot(s >= 1)

  # uppder bound of Variance (p.17)
  n.u <- nj*(N-nj) + 1/4*(N+s)*(s+1)^2
  d.u <- (N+s)*(N+s+1)^2
  v.u <- n.u/d.u

  # lower bound of Variance (p.17)
  n.l <- nj*(N-nj) + s*min(c(nj, N-nj))
  d.l <- (N+s)^2*(N+s+1)
  v.l <- n.l/d.l

  robj <- list(v.lower=v.l, v.upper=v.u)
  return(robj)
}

# 4.2. Means and standard deviations of theta_R
lapply(iddm.inf(nj=1, N=6, s=2), sqrt)
