#'
#' @rdname idm
#' @title Imprecise Dirichlet Model
#' @description This function computes lower and upper posterior probabilities under an imprecise Dirichlet model when prior information is not available.
#' @param nj number of observations in the j th category
#' @param s learning parameter
#' @param N total number of drawings
#' @param k number of elements in the sample space
#' @param tj mean probability associated with the j th category
#' @param cA the number of elements in the event A
#' @examples
#' idm(nj=1, N=6, s=2, k=4)
#' @references
#' Walley, P. (1996), Inferences from Multinomial Data: Learning About a Bag of Marbles. Journal of the Royal Statistical Society: Series B (Methodological), 58: 3-34. https://doi.org/10.1111/j.2517-6161.1996.tb02065.x
#' @export
idm <- function(nj, s=1, N, tj=NA_real_, k, cA=1){
  stopifnot(s >= 0)
  stopifnot(nj >= 0)
  stopifnot(N >= 0)

  if(is.na(tj)) tj <- 0.5

  # sec2.3/
  # p <- tj.star <- (nj+s*tj)/(N+s)
  # sec2.4/
  p <- (nj + s/k*cA)/(N+s)

  ## sec2.3/upper bound (limiting as tj -> 1)
  p.u <- (nj+s)/(N+s)

  ## sec2.3/lower bound (limiting as tj -> 0)
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

#' @rdname betabinom
#' @title Beta-Binomial Distribution
#' @description This function computes the predictive posterior density of the outcome of interest under the imprecise Dirichlet prior distribution. It follows a beta-binomial distribution.
#' @param x number of occurrence of event A in the N previous trials
#' @param N total number of previous trials
#' @param M number of future trials
#' @param i number of occurrences of event A in the M future trials
#' @param y maximum number of occurrences of event A in the M future trials
#' @param tA prior probability of event A under the Dirichlet prior
#' @param s learning parameter
#' @examples
#' pbetabinom(M=6, x=1, s=1, N=6, y=0)
#' @export
dbetabinom <- function(i, M, x, s, N, tA){
  # stopifnot(0 >= i & i <= M)
  # stopifnot(0 >= tA & tA <= 1)
  alpha <- s*tA
  beta <- s-s*tA
  pi <- choose(n=M,k=i)*beta(alpha+x+i, beta+N-x+M-i)/beta(alpha+x, beta+N-x)
  return(pi)
}

#' @rdname betabinom
#' @export
pbetabinom <- function(M, x, s, N, y){
  p.l <- 0
  for(i in 0:y){
    p0 <- dbetabinom(i=i, M=M, x=x, s=s, N=N, tA=1) # minimized by limiting tA to 1
    p.l <- p.l + p0
  }
  p.u <- 0
  for(i in 0:y){
    p0 <- dbetabinom(i=i, M=M, x=x, s=s, N=N, tA=0) # maximized by limiting tA to 0
    p.u <- p.u + p0
  }
  robj <- list(p.l=p.l, p.u=p.u)
  return(robj)
}



#' @rdname idm
#' @title Highest Posterior Density Interval
#' @description This function searches for the lower and upper bounds of a given level of the highest posterior density interval under the imprecise Dirichlet prior.
#' @param maxiter maximum number of iterations
#' @param tolerance level of error allowed
#' @param alpha shape1 parameter of beta distribution
#' @param beta shape2 parameter of beta distribution
#' @param p level of credible interval
#' @param verbose logical option suppressing messages
#' @examples
#' x <- hpd(alpha=3, beta=5, p=0.95) # c(0.0031, 0.6587) when s=2
#' # round(x,4); x*(1-x)^5
#' @export
hpd <- function(alpha=3, beta=5, p=0.95, tolerance=1e-4, maxiter=1e2, verbose=FALSE){

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
    fn0 <- function(theta) stats::dbeta(theta, alpha, beta) # 105*theta^2*(1-theta)^4
    dif <- abs(stats::integrate(f=fn0, lower=a, upper=b)$value - p)
    return(dif)
  }

  # initialization
  a <- 1e-8
  b <- 1-a
  dif <- 1e-2
  niter <- 0

  while( dif > tolerance){

    # for fixed b, searching for a
    op1 <- stats::optimize(f=fn1, b=b, lower=0, upper=stats::qbeta(1-p, alpha, beta))
    a <- op1$minimum

    # for fixed a, searching for b
    op2 <- stats::optimize(f=fn2, a=a, lower=stats::qbeta(p, alpha, beta), upper=1)
    b <- op2$minimum

    dif <- fn1(a=a, b=b)
    if(verbose) print(c(a=a,b=b, dif=dif))
    niter <- niter + 1
    if(niter == maxiter) break
  }

  robj <- c(a=a, b=b)
  return(robj)
}




#' @rdname ibm
#' @title Impreicse Beta Model
#' @description This function computes lower and upper posterior probabilities under an imprecise Beta model when prior information is not available.
#' @param n total of trials
#' @param m number of observations realized
#' @param s0 learning parameter
#' @param xlab1 x axis text
#' @param main1 main title text
#' @references
#' Walley, P. (1996), Inferences from Multinomial Data: Learning About a Bag of Marbles. Journal of the Royal Statistical Society: Series B (Methodological), 58: 3-34. https://doi.org/10.1111/j.2517-6161.1996.tb02065.x
#' @examples
#' tc <- seq(0,1,0.1)
#' s <- 2
#' ibm(n=10, m=6)
#' @export
ibm <- function(n=10, m=6, s0=2, xlab1=NA, main1=NA){

  # a set of priors
  t0 <- seq(from=0.0, to=1, by=0.1)

  # translation of mean parameters
  sn <- s0 + n
  tn <- (s0*t0 + m)/(s0+n)

  # parameter space for theta
  theta <- seq(from=0, to=1, by=0.02)

  # Updating the parameters of posterior distributions
  alpha.n <- sn*tn
  beta.n <- sn*(1-tn)

  # Calculating the CDF of posteriors
  fn <- function(x, y, ...) stats::pbeta(q=theta, shape1=x, shape2=y)

  posteriors <- as.data.frame(mapply(fn, x=alpha.n, y=beta.n))

  plot(x=theta, y=rep(0, length(theta)), type='n', ylim=c(0,1), ylab="F(z)", xlab=xlab1, main=main1)
  lapply(posteriors, graphics::lines, x=theta, type='l', lty=3)
  graphics::text(x=0.1, y=0.9, bquote(paste("n=", .(n), ", m=", .(m))))

  graphics::lines(x=theta, y=posteriors[[1]], col='blue') # upper probabilities
  graphics::lines(x=theta, y=posteriors[[length(t0)]], col='red') # lower probabilities
  graphics::abline(v=m/n)
}


#' @rdname betadif
#' @title Distribution of Difference of Two Proportions
#' @param x difference of two beta distributions
#' @param a1 shape 1 parameter of Beta distribution with control
#' @param b1 shape 2 parameter of Beta distribution with control
#' @param a2 shape 1 parameter of Beta distribution with treatment
#' @param b2 shape 2 parameter of Beta distribution with treatment
#' @references
#' Chen, Y., & Luo, S. (2011). A few remarks on 'Statistical distribution of the difference of two proportions' by Nadarajah and Kotz, Statistics in Medicine 2007; 26 (18): 3518-3523. Statistics in Medicine, 30(15), 1913-1915.
#' @export
dbetadif <- function(x, a1, b1, a2, b2){

  stopifnot(-1<= x & x <= 1)
  stopifnot(a1 > 0 & b1 > 0 & a2 > 0 & b2 > 0)

  if(x <= 0 & x > -1){
    v <- ((-x)^(b1+b2-1) * (1+x)^(a1+b2-1)) / (gamma(b1)*gamma(a2)*gamma(a1+b2)) * tolerance::F1(a=b2, b=a1+a2+b1+b2-2, b.prime=1-a2, c=a1+b2, x=1+x, y=1-x^2)
  } else { # x > 0 & x <= 1
    v <- x^(b1+b2-1)*(1-x)^(a2+b1-1) / (gamma(b2)*gamma(a1)*gamma(a2+b1)) * tolerance::F1(a=b1, b=a1+a2+b1+b2-2, b.prime=1-a1, c=a2+b1, x=1-x, y=1-x^2)
  }
  v <- v*gamma(a1+b1)*gamma(a2+b2)
  return(v)
}
