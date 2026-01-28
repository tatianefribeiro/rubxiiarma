# cumulative distribution function
cdf_RUBXII <- function(y,q,c, tau = tau)   #RUBXII cdf (useful for qqplot)
{
  1-(1+(log(1/(1-y)))^c)^(log(1-tau)/log(1+(log(1/(1-q)))^c))
}


# inversion method for randon generation
r_RUBXII <- function(n,q_t,c,tau = tau)   #It generates occurences of Y_i ~ RUBXII (q_t, c)
{
  u = runif(n)
  y_RUBXII =1-exp(-((1-u)^(log(1+(log(1/(1-q_t)))^c)/log(1-tau))-1)^(1/c)) #RUBXII qf
  return(y_RUBXII)
}

