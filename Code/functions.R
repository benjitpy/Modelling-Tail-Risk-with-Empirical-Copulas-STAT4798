## Functions ###

## First Function: s.C.n, pbetab, pbeta2, kbetab - UNUSED!, kbeta2 - UNUSED!
##' @title A class of smooth, data-adaptive nonparametric copula estimators derived from EBC (Kojadinovic and Yi), and some modified distributions
##' @param u matrix of d-dimensional evaluation points (points correspond to rows)
##' @param r matrix of d-dimensional ranks from data
##' @param marg smoothing margins: scaled binomial, scaled beta-binomial or beta
##' @param rho smoothing parameter for the scaled beta-binomial and beta margins
##' @param cop the smoothing d-dimensional survival copula
##' @return values of smooth empcop
##' @author Ivan Kojadinovic and Bingqing Yi, Modified by Benjamin Tong

pbetab <- function(x, m, u, rho, lower.tail) {
  if (rho <= 1) {
    warning("rho is <= 1; using binomial instead")
    return(pbinom(x, size = m, prob = u, lower.tail = lower.tail))
  }
  ifelse((u == 0) | (u == 1),
         pbinom(x, size = m, prob = u, lower.tail = lower.tail), # degenerate cases
         pbbinom(x, size = m, alpha = (rho - m) * u / (1 - rho),
                 beta = (rho - m) / (1 - rho) * (1 - u),
                 lower.tail = lower.tail)
  )
}

pbeta2 <- function(x, m, u, rho, ...)
  pbeta(x, shape1 = (m - rho) / rho * u, shape2 = (m - rho) / rho * (1 - u), ...)

kbetab.inv <- function(u, r, m, rho = 4) {
  f <- function(v)
    pbetab(x = r - 1, m = m, u = v, rho = rho, lower.tail = FALSE) - u
  uniroot(f, interval = 0:1)$root
}

kbeta2.inv <- function(u, r, m, rho = 4) {
  f <- function(v)
    pbeta2(x = (r - 0.5)/m, m = m, u = v, rho = rho) - u
  uniroot(f, interval = 0:1)$root
} 

s.C.n <- function(u, r, marg = c("binomial", "betabinomial", "beta"), rho ,
                  cop = indepCopula(dim = ncol(r)))

{
if(any(u < 0, 1 < u))
  stop("'u' must be in [0,1].")

  stopifnot((d <- ncol(u)) == ncol(r) && d == dim(cop))
  #depends on marg chosen
  
  marg <- match.arg(marg) #depends on marg chosen
  n <- nrow(r) #number of row of rank matrix
  m <- nrow(u) #number of evaluation points
  v <- matrix(NA, nrow = m, ncol = d) #create a m*d matrix to construct copula
  ec <- numeric(m) #construct a empty vector as ec (as the smooth copula is weighted avg of ec)
  for (i in 1:n) {
    for (j in 1:d) {
      v[,j] <- switch(marg,
                      binomial = pbinom(r[i,j] - 1, size = n, prob = u[,j], lower.tail = FALSE),
                      betabinomial = pbetab(r[i,j] - 1, m = n, u = u[,j], rho = rho, lower.tail = FALSE),
                      beta = pbeta2((r[i,j] - 0.5)/n, m = n, u = u[,j], rho = rho, lower.tail = FALSE))
    }
    ec <- ec + pCopula(v, copula = cop)
  }
  ec / n
}  

## Second Function: rBetaBCopula, rBeta2Copula - UNUSED!
##' @title Samping from the "BetaB4" and "Beta2" smooth empirical copula
##' @param n sample size
##' @param x matrix of d-dimensional observations (or multivariate ranks)
##'          from which the "BetaB4" smooth empirical copula will be computed
##' @param ranks logical indicating whether x contains observations
##'              or corresponding multivariate ranks
##' @param rho parameter of the margins of the smoothing distributions;
##'            set to 4 by default
##' @return a random sample of size n from the "BetaB4" smooth empirical copula
##' @author Ivan Kojadinovic and Bingqing Yi, modified by Benjamin Tong

rBetaBCopula <- function(n, x, ranks = FALSE, rho = 4, ebcop) {
  
  ## Checks
  if(!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  stopifnot(n >= 1L)
  
  ## Data dimensions
  m <- nrow(x)
  d <- ncol(x)
  
  ## Multivariate ranks from x
  if (!ranks)
    x <- apply(x, 2, rank)
  
  ## Generate n realizations from the empirical
  ## beta copula first
  
  u <- rCopula(n, copula = ebcop)
  
  ## Transform these realizations marginally
  for (i in 1:n) {
    I <- sample(1:m, size = 1) # select a row
    #u[i,] = vapply(seq_len(d), function(k) # iterate over rows k of u
    #  kbetab.inv(u[i,k], r = x[I, k], m = m, rho = rho), NA_real_)
    for (j in 1:d) {
      u[i,j] <- kbetab.inv(u[i,j], r = x[I, j], m = m, rho = rho)
    }
  }
  u
}

rBeta2Copula <- function(n, x, ranks = FALSE, rho = 4, ebcop) {
  
  ## Checks
  if(!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  stopifnot(n >= 1L)
  
  ## Data dimensions
  m <- nrow(x)
  d <- ncol(x)
  
  ## Multivariate ranks from x
  if (!ranks)
    x <- apply(x, 2, rank)
  
  ## Generate n realizations from the empirical
  ## beta copula first
  u <- rCopula(n, copula = ebcop)
  
  ## Transform these realizations marginally
  for (i in 1:n) {
    I <- sample(1:m, size = 1) # select a row
    for (j in 1:d) {
      r <- x[I, j]
      u[i,j] <- kbeta2.inv(u[i,j], r = r, m = m, rho = rho)
    }
  }
  u
}

## Third Function: empirical_scopula
##' @title compute survival empirical copula (for common empirical classes) without using numerical approximation
##' @author Marius Hofert, modified by Benjamin Tong
empirical_scopula <- function(u, U, smoothing = c("none", "pbeta", "dbeta", "checkerboard"),
                             offset = 0, log = FALSE, ...)
{
  stopifnot(0 <= u, u <= 1, 0 <= U, U <= 1)
  if(!is.matrix(u)) u <- rbind(u)
  if(!is.matrix(U)) U <- rbind(U)
  m <- nrow(u) # number of evaluation points
  n <- nrow(U) # number of points based on which the empirical copula is computed
  d <- ncol(U) # dimension
  stopifnot(ncol(u) == d)
  R <- apply(U, 2, rank) # (n, d)-matrix of ranks
  switch(match.arg(smoothing),
         "none" = {
           R. <- t(R) / (n + 1) # (d, n)-matrix
           vapply(seq_len(m), function(k) # iterate over rows k of u
             sum(colSums(R. >= u[k,]) == d) / (n + offset), NA_real_)
         },
         "pbeta" = {
           vapply(seq_len(m), function(k) { # iterate over rows k of u
             sum( # sum() over i
               vapply(seq_len(n), function(i)
                 prod( pbeta(u[k,], shape1 = R[i,], shape2 = n + 1 - R[i,], lower.tail = FALSE) ), # prod() over j
                 NA_real_)) / (n + offset)
           }, NA_real_)
         },
         "dbeta" = {
           if(log) {
             vapply(seq_len(m), function(k) { # iterate over rows k of u
               lsum( # lsum() over i
                 vapply(seq_len(n), function(i) {
                   ## k and i are fixed now
                   lx.k.i <- sum( dbeta(u[k,], shape1 = R[i,], shape2 = n + 1 - R[i,], log = TRUE) ) # log(prod()) = sum(log()) over j for fixed k and i
                 },
                 NA_real_)) - log(n + offset)
             }, NA_real_)
           } else { # as for 'pbeta', just with dbeta()
             vapply(seq_len(m), function(k) { # iterate over rows k of u
               sum( # sum() over i
                 vapply(seq_len(n), function(i)
                   prod( dbeta(u[k,], shape1 = R[i,], shape2 = n + 1 - R[i,], lower.tail = FALSE) ), # prod() over j
                   NA_real_)) / (n + offset)
             }, NA_real_)
           }
         },
         "checkerboard" = {
           R. <- t(R) # (d, n)-matrix
           vapply(seq_len(m), function(k) # iterate over rows k of u
             sum(apply(pmin(pmax(n * u[k,] - R. + 1, 0), 1), 2, prod)) / (n + offset),
             NA_real_) # pmin(...) = (d, n)-matrix
         },
         stop("Wrong 'smoothing'"))
}

## Fourth Function: power_set_fn, sempcop_analytical
##' @title Calculate survival copula by first principles for empirical copula, power set function to generate J
##' @param d dimension
##' @param copula desired true copula
##' @author Benjamin Tong

power_set_fn <- function(set) {
  if (length(set) == 0) {
    return(list(list()))
  }
  
  first <- set[1]
  remaining <- set[-1]
  
  subsets <- power_set_fn(remaining)
  
  new_subsets <- lapply(subsets, function(subset) c(first, subset))
  
  return(c(subsets, new_subsets))
}


sempcop_analytical <- function(u, U, copula, smoothing) {
  d = ncol(u)
  scop = 0
  if (smoothing == "none") {
    U = pobs(U)
    scop = empirical_scopula(u, U, smoothing = "none")
    return(scop)
  }
  if (smoothing == "beta") {
    U = pobs(U)
    ### The below is the algorithm used to verify that the modified definition is equivalent. 
    #  u = matrix(rep(1 - rev(u[,1]), d), ncol = d)
    #   power_set <- power_set_fn(1:d)
    #   for (i in 1:length(power_set)) {
    #     if (length(power_set[[i]]) == 0) {
    #       scop = scop + 1
    #       next
    #     }
    #     if(length(power_set[[i]]) == 1) {
    #       scop = scop - (1 - u[,unlist(power_set[[i]])])
    #       next
    #     }
    #     else {
    #       column_target = u[,unlist(power_set[[i]])]
    #       #target_cop = empCopula(pobs(U[,unlist(power_set[[i]])]), smoothing = smoothing)
    #       #target_cop = gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = length(power_set[[i]]))
    #       scop = scop + (-1)^(length(power_set[[i]]))*C.n(1 - column_target, U[,unlist(power_set[[i]])], smoothing = "beta")
    #       #scop = scop + (-1)^(length(power_set[[i]]))*pCopula(1 - column_target, target_cop)
    # 
    #
    #    }
    #  }
    scop = empirical_scopula(u, U, smoothing = "pbeta")
    return(scop)
    #return(rev(scop))
  }
  if (smoothing == "checkerboard") {
    U = pobs(U)
    u = matrix(rep(1 - rev(u[,1]), d), ncol = d)
    power_set <- power_set_fn(1:d)
    for (i in 1:length(power_set)) {
      if (length(power_set[[i]]) == 0) {
        scop = scop + 1
        next
      }
      if(length(power_set[[i]]) == 1) {
        scop = scop - (1 - u[,unlist(power_set[[i]])])
        next
      }
      else {
        column_target = u[,unlist(power_set[[i]])]
        scop = scop + (-1)^(length(power_set[[i]]))*C.n(1 - column_target, pobs(U[,unlist(power_set[[i]])]), smoothing = "checkerboard")
      }
    }
    return(rev(scop))
  }
  if (smoothing == "betabinomial4") {       
    
    u = matrix(rep(1 - rev(u[,1]), d), ncol = d)
    power_set <- power_set_fn(1:d)
    for (i in 1:length(power_set)) {
      if (length(power_set[[i]]) == 0) {
        scop = scop + 1
        next
      }
      if(length(power_set[[i]]) == 1) {
        scop = scop - (1 - u[,unlist(power_set[[i]])])
        next
      }
      else {
        column_target = u[,unlist(power_set[[i]])]
        target_cop = empCopula(pobs(U[,unlist(power_set[[i]])]), smoothing = "beta")
        #target_cop = gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = length(power_set[[i]]))
        scop = scop + (-1)^(length(power_set[[i]]))*s.C.n(1 - column_target, apply(U[,unlist(power_set[[i]])],2,rank), marg = "betabinomial", rho = 4, cop = target_cop)
      }
    }
    return(rev(scop))
  }
  if (smoothing == "binomial4") {       
    
    u = matrix(rep(1 - rev(u[,1]), d), ncol = d)
    power_set <- power_set_fn(1:d)
    for (i in 1:length(power_set)) {
      if (length(power_set[[i]]) == 0) {
        scop = scop + 1
        next
      }
      if(length(power_set[[i]]) == 1) {
        scop = scop - (1 - u[,unlist(power_set[[i]])])
        next
      }
      else {
        column_target = u[,unlist(power_set[[i]])]
        target_cop = empCopula(pobs(U[,unlist(power_set[[i]])]), smoothing = "beta")
        #target_cop = gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = length(power_set[[i]]))
        scop = scop + (-1)^(length(power_set[[i]]))*s.C.n(1 - column_target, apply(U[,unlist(power_set[[i]])],2,rank), marg = "binomial", rho = 4, cop = target_cop)
      }
    }
    return(rev(scop))
  }
  if (smoothing == "beta4") {       
    
    u = matrix(rep(1 - rev(u[,1]), d), ncol = d)
    power_set <- power_set_fn(1:d)
    for (i in 1:length(power_set)) {
      if (length(power_set[[i]]) == 0) {
        scop = scop + 1
        next
      }
      if(length(power_set[[i]]) == 1) {
        scop = scop - (1 - u[,unlist(power_set[[i]])])
        next
      }
      else {
        column_target = u[,unlist(power_set[[i]])]
        target_cop = empCopula(pobs(U[,unlist(power_set[[i]])]), smoothing = "beta")
        #target_cop = gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = length(power_set[[i]]))
        scop = scop + (-1)^(length(power_set[[i]]))*s.C.n(1 - column_target, apply(U[,unlist(power_set[[i]])],2,rank), marg = "beta", rho = 4, cop = target_cop)
      }
    }
    return(rev(scop))
  }
}

## Fifth Function: CvM
##' @title Calculate Discretized Cramer von Mises Statistic
##' @param res_mat matrix as returned by emp_prob()
##' @param n number of observations that created res_mat
##' @return ...
##' @author Benjamin Tong

CvM <- function(res_mat, n) {
  statistic = numeric(length(cop.names))
  realscop_stat <- res_mat[,1]
  for (i in 1:length(cop.names)) {
    empscop_stat <- res_mat[,-1 + 3*i]
    statistic[i] <- sum((empscop_stat[1:n] - realscop_stat[1:n])^2)
  }
  statistic
}

## Sixth Function: emp_survivalprob
##' @title Compute Simulated Mean Estimates and CIs of exceedance probability
##' @param copula true underlying copula that is sampled from
##' @return (length(u), 1+3*m)-matrix
##'         column 1: true copula evaluated at (u,...,u)
##'         column 1+3k-2, 1+3k-1, 1+3k: simulated mean estimates and CIs for kth empirical copula
##' @author Benjamin Tong and Marius Hofert

emp_survivalprob <- function(copula) 
{

  d <- dim(copula) #dimension
  u_mat = matrix(rep(u_uptail, d), ncol = d) #(length(u), d) evaluation matrix
  
  length.u_mat <- length(u_uptail) #length of evaluation matrix
  length.cop.names <- length(cop.names) #length of ecop.nms
  res_mat = matrix(NA, nrow = length.u_mat, ncol = 1 + 3 * length.cop.names) #result matrix
  colnames(res_mat) = c("True", sapply(cop.names, function(name) c(name, "Low", "Up")))
  arr <- array(NA, dim = c(length.cop.names, B, length.u_mat), dimnames = list(cop.names, 1:B, u_uptail))
  #array, for each ecop, do it for B replications and for all u evaluation points
  
  #MC simulation
  
  res_mat[,"True"] = pCopula(1 - u_mat, copula = rotCopula(copula = copula))

  for (b in 1:B) {
    
    U <- rCopula(n, copula = copula) #sample pseudo-observations from the true copula
     
    arr["Empirical copula", b,] <- sempcop_analytical(u_mat, U, copula, smoothing = "none")
    arr["Empirical beta copula", b,] <- sempcop_analytical(u_mat, U, copula, smoothing = "beta")
    arr["Empirical checkerboard copula", b,] <- sempcop_analytical(u_mat, U, copula, smoothing = "checkerboard")
    arr["Empirical binomial 4", b, ] <- sempcop_analytical(u_mat, U, copula, smoothing = "binomial4")
    arr["Empirical betabinomial 4", b, ] <- sempcop_analytical(u_mat, U, copula, smoothing = "betabinomial4")
    arr["Empirical beta2 4", b, ] <- sempcop_analytical(u_mat, U, copula, smoothing = "beta4")
    
    ## Progress report
    if(b %% ceiling(B/20) == 0)
      cat(sprintf("%2d%% ", ceiling(b/B * 100)))
    if(b == B) cat("\n")
  }
  
  #compute mean estimate and CIs
  mean_estimate <- apply(arr, c(1,3), mean) # (m, n.u)-matrix #take mean across all B replications
  CIs <- apply(arr, c(1,3), quantile, probs = c(0.025, 0.975), names = FALSE) #find quantile across all B replications
  
  #copy into res matrix
  for (k in 1:length.cop.names) {
    res_mat[, 1+3*k-2] <- mean_estimate[k,]
    res_mat[, 1+3*k-1] <- CIs[1,k,]
    res_mat[, 1+3*k] <- CIs[2,k,]
  }
  
  #return result_matrix
  res_mat
}

## Seventh Function: emp_cprob
##' @title Compute Simulated Mean Estimates and CIs of lower-tail probability
##' @param copula true underlying copula that is sampled from
##' @return (length(u), 1+3*m)-matrix
##'         column 1: true copula evaluated at (u,...,u)
##'         column 1+3k-2, 1+3k-1, 1+3k: simulated mean estimates and CIs for kth empirical copula
##' @author Benjamin Tong and Marius Hofert

emp_cprob <- function(copula) 
{
  
  d <- dim(copula) #dimension
  u_mat = matrix(rep(u_lowtail, d), ncol = d) #(length(u), d) evaluation matrix
  length.u_mat <- length(u_lowtail) #length of evaluation matrix
  length.cop.names <- length(cop.names) #length of ecop.nms
  res_mat = matrix(NA, nrow = length.u_mat, ncol = 1 + 3 * length.cop.names) #result matrix
  colnames(res_mat) = c("True", sapply(cop.names, function(name) c(name, "Low", "Up")))
  arr <- array(NA, dim = c(length.cop.names, B, length.u_mat), dimnames = list(cop.names, 1:B, u_lowtail))
  #array, for each ecop, do it for B replications and for all u evaluation points
  
  ## Evaluate true copula
  res_mat[,"True"] <- pCopula(u_mat, copula = copula) # n.u-vector of true probabilities under 'copula'
  
  #MC simulation
  
  for (b in 1:B) {
    
    U <- rCopula(n, copula = copula) #sample pseudo-observations from the true copula
    
    ecop <- empCopula(pobs(U)) #create empirical copula using the pseudo-observations
    ebcop <- empCopula(pobs(U), smoothing = "beta") #create empirical beta copula using the pseudo-observations
    eccop <- empCopula(pobs(U), smoothing = "checkerboard") #create empirical checkerboard copula using the pseudo-observations
    
    arr["Empirical copula", b,] <- pCopula(u_mat, copula = ecop)
    arr["Empirical beta copula", b,] <- pCopula(u_mat, copula = ebcop)
    arr["Empirical checkerboard copula", b, ] <- pCopula(u_mat, copula = eccop)    
    arr["Empirical binomial 4", b, ] <- s.C.n(u_mat, apply(U, 2, rank), marg = "binomial", rho = 4, cop = ebcop)
    arr["Empirical betabinomial 4", b, ] <- s.C.n(u_mat, apply(U, 2, rank), marg = "betabinomial", rho = 4, cop = ebcop)
    arr["Empirical beta2 4", b, ] <- s.C.n(u_mat, apply(U, 2, rank), marg = "beta", rho = 4, cop = ebcop)
    
    ## Progress report
    if(b %% ceiling(B/100) == 0)
      cat(sprintf("%2d%% ", ceiling(b/B * 100)))
    if(b == B) cat("\n")
    
  }
  
  #compute mean estimate and CIs
  mean_estimate <- apply(arr, c(1,3), mean) # (m, n.u)-matrix #take mean across all B replications
  CIs <- apply(arr, c(1,3), quantile, probs = c(0.025, 0.975), names = FALSE) #find quantile across all B replications
  
  #copy into res matrix
  for (k in 1:length.cop.names) {
    res_mat[, 1+3*k-2] <- mean_estimate[k,]
    res_mat[, 1+3*k-1] <- CIs[1,k,]
    res_mat[, 1+3*k] <- CIs[2,k,]
  }
  
  #return result_matrix
  res_mat
}

## Eighth Function: plot_res_mat
##' @title Plot an Object as Returned by emp_prob()
##' @param res_mat matrix as returned by emp_prob()
##' @param u low-tailed or high-tailed probability for x-axis

plot_res_mat <- function(res_mat, u, tau)
{
  plot(u, res_mat[,1], type = "l",
       #xlim = c(0.9,1), 
       ylim = c(min(res_mat),max(res_mat)),
       xlab = "u of copula evaluation point (u,...,u)", ylab = "")
  title(paste("tau = ", tau))
  m <- (ncol(res_mat)-1)/3 # number of empirical copulas
  for(k in 1:m) {
    lines(u, res_mat[,1+3*k-2], col = k+1) # estimate
    lines(u, res_mat[,1+3*k-1], lty = 3, col = k+1) # lower CI
    lines(u, res_mat[,1+3*k],   lty = 3, col = k+1) # upper CI
  }
  #lgd <- c("True copula value",
  #         sapply(1:m, function(k) c(cop.names[k], "Lower/upper 95%-CI")))
  #legend("topleft", bty = "n", lty = c(1, rep(c(1,3), m)),
  #       col = c(1, rep(2:(m+1), each = 2)), legend = lgd)
}







