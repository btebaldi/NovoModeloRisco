
dthanssen <- function(x, lambda = 0, eta = 3, log.p = FALSE){
  # dnorm(x, mean = 0, sd = 1, log = FALSE)
  # pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # rnorm(n, mean = 0, sd = 1
  
  if((lambda > 1) | (lambda < -1)){
    stop("erro: -1 < lambda < 1")
  }
  
  if(eta <= 2){
    stop("erro: 2 < eta < infty")
  }
  
  # Calculating hanssen variables
  c = gamma((eta+1)/2)/(((pi*(eta-2))^0.5) * gamma(eta/2))
  # cat(sprintf("c:%f\n", c))
  
  a = 4*lambda*c*((eta-2)/(eta-1))
  # cat(sprintf("a:%f\n", a))
  
  b = (1+ 3*lambda^2 - a^2)^0.5
  # cat(sprintf("b:%f\n", b))
  
  f = as.numeric(NA)
  
  L = ifelse(test = (x < -a/b), yes = -lambda, no = +lambda)
  f = b*c*(1 + (1/(eta-2)) * ((b*x+a)/(1+L))^2 )^(-(eta+1)/2)
  
  # 
  # if(x < -a/b){
  #   f = b*c*(1 + (1/(eta-2)) * ((b*x+a)/(1-lambda))^2 )^(-(eta+1)/2)
  # } else {
  #   f = b*c*(1 + (1/(eta-2)) * ((b*x+a)/(1+lambda))^2 )^(-(eta+1)/2)
  # }
  
  if(log.p){f= log(f)}
  
  return(f)
}

pthanssen <- function(q, lambda = 0, eta = 3, lower.tail = TRUE, log.p = FALSE){
  # dnorm(x, mean = 0, sd = 1, log = FALSE)
  # pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # rnorm(n, mean = 0, sd = 1
  base = 0.0001
  x = seq(from = -100, to = q, by = base)
  
  y = dthanssen(x, lambda = lambda, eta = eta)
  
  LT = sum(y*base  + c(y[1], base::diff(y)) * base /2)
  
  if(lower.tail){
    f = LT
  } else { 
    f = 1-LT
  }
  
  if(log.p){f = log(f)}
  
  return(f)    
}

qthanssen <- function(p, lambda = 0, eta = 3, lower.tail = TRUE){
  # dnorm(x, mean = 0, sd = 1, log = FALSE)
  # pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # rnorm(n, mean = 0, sd = 1
  
  base = 0.0001
  x = seq(from = -100, to = 100, by = base)
  
  y = dthanssen(x, lambda = lambda, eta = eta)
  
  LT = cumsum(y*base  + c(y[1], base::diff(y)) * base /2)
  
  if(lower.tail) { 
    id = max(which(LT <= p))
  }
  else {
    id = min(which(LT >= 1-p))
  }
  
  return(x[id])
}

rthanssen <- function(n, lambda = 0, eta = 3){
  # dnorm(x, mean = 0, sd = 1, log = FALSE)
  # pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # rnorm(n, mean = 0, sd = 1
  
  sampleCDF <- runif(n)
  
  base = 0.0001
  x = seq(from = -100, to = 100, by = base)
  
  y = dthanssen(x, lambda = lambda, eta = eta)
  
  LT = cumsum(y*base  + c(y[1], base::diff(y)) * base /2)

  return(x[sapply(sampleCDF, function(x){max(which(LT <= x))})])
}

mle.thanssen <- function(data, lambda = 0, eta = 2.1){
  
  param <- c("lambda" = lambda, "eta" = eta)
  # param <- c("lambda" = 1, "eta" = 1.1)
  
  # Restriction coefficients 
  ui = matrix(data = c(1,0,-1,0,0,1), nrow = 3, ncol = 2, byrow = TRUE)
  
  # Restriction matrix
  ci = matrix(data = c(-0.9999999999,-0.9999999999,2.0000000001), nrow = 3, ncol = 1, byrow = TRUE)
  
  # theta = matrix(data = param, nrow = 2, ncol = 1, byrow = TRUE)
  # ui %*% param - ci >= 0.
  
  mOpt <- constrOptim(theta = param,
                      f = nll.thanssen,
                      grad=NULL, ui=ui, ci=ci,
                      control = list(),
                      outer.iterations = 100,
                      outer.eps = 1e-05, 
                      data=data)
  return(mOpt)    
}


nll.thanssen <- function(param, data){
  lambda <- param[1]
  eta <- param[2]
  ll <- sum(dthanssen(x = data, lambda, eta, log.p = TRUE))
  return(-ll)
}