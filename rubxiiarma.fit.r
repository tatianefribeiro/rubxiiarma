# Algumas informacoes:
# diag = 0 : nao plota graficos
# diag = 1 : plota graficos na tela
# diga = 2 : gera graficos em pdf e plota graficos na tela

rubxiiarma.fit<-function (y, 
                          ar = NA,
                          ma = NA,
                          tau = 0.5,
                          link = "logit",
                          h = 0,
                          diag = 0,
                          X = NA,
                          X_hat = NA,
                          check_hessian = FALSE)
{
  source("rubxii-funcs.r")
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if(is.ts(y)==T)  freq<-frequency(y) else stop("data can be a time-series object")
  
  z<-c()
  maxit1<-50
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  #n1 <- n-h
  m <- max(p,q,na.rm=T)
  y1 <- y[(m+1):n]
  p1 <- length(ar)  
  q1 <- length(ma)
  error <- rep(0,n) 
  eta <- rep(NA,n)
  
  y_prev <- c(rep(NA,(n+h))) 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog"))){
    stats <- make.link(linktemp)
  }  else {
    stop(paste(linktemp, "link not available, available links are \"logit\", ","\"probit\" and \"cloglog\""))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ynew = linkfun(y) 
  ynew_ar <- suppressWarnings(matrix(ynew,(n-1),max(p,1,na.rm=T)))
  
  ###########################################################3  
  if(any(is.na(ar)) == F) {
    names_phi <- c(paste("phi", ar, sep = ""))
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1),])} else {
      ar = p1<-0; Z <- NA  
    } 
  
  if(any(is.na(ma)) == F) {
    names_theta <- c(paste("theta", ma, sep = ""))
  } else ma = q1 <- 0 
  
  if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n, ]       
    k = ncol(X)
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
  }
  
  ###FB  recorrences   #*
  q_1 <- max(q1, 1)
  R <- matrix(rep(NA, (n-m)*q_1), ncol = q_1)  # pq precisamos repetir NA (nlinhas*ncol), nesse caso
  k_i <- q1/q_1                                 # útil para o vetor escore (termo R); vai ser sempre 0 (pois q1 pode ser) ou 1?
  #                                             # não queremos multip a dim pelo num par med moveis, mas sim por 0 ou 1, ie, ou tem ou não tem
  
  deta.dalpha <- rep(0, n)
  deta.dbeta <- matrix(0, ncol=max(k,1), nrow=n)
  deta.dphi <- matrix(0, ncol=p1, nrow=n)
  deta.dtheta <- matrix(0, ncol=q_1, nrow=n)
  
  Xstart <- (cbind(rep(1, (n-m)), Xm, Z))         #Conferir se essa Xm está certa mesmo, pois pode ser isso o erro no BETA-ARMA, USAMOS ELA NO CHUTE....
  Xstart <- matrix(apply(Xstart, 1, na.omit),nrow = (n-m),byrow = T)
  ols <- lm.fit(Xstart, ynew[(m+1) : n])$coef
  initial <- c(rep(0, k+p1+q1+1),1)
  initial[1 : (k+p1+1)] <- ols
  
  loglik <- function(z) 
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    c_par <- z[length(z)]
    
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    q_t <- linkinv(eta[(m+1):n])
    
    
    ll <- log(
      c_par*log(1/(1-tau))*(log(1/(1-y1)))^(c_par-1)*
        (1+(log(1/(1-y1)))^c_par)^(log(1-tau)/log(1+(log(1/(1-q_t)))^c_par)-1)/  
        ((1-y1)*log(1+(log(1/(1-q_t)))^c_par))    
    )
    sum(ll)
  } 
  
  #######################################################################
  escore.RUBXIIarma <- function(z)
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    c_par <- z[length(z)]
    
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1),ar] - Xbeta_ar[(i-1),ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
      #mu[i]   <- linkinv(eta[i])
      #error[i]<- qnorm(cdf_RUBXII(y[i],mu[i],c_par, tau = 0.1)) # qr
    }
    
    q_t <- linkinv(eta[(m+1):n])
    
    Xbeta <- X%*%beta
    for(i in 1:(n-m)){
      R[i,] <- error[i+m-ma]*k_i}
    
    for(i in (m+1):n)
    {
      deta.dalpha[i] <- 1 - deta.dalpha[i-ma]%*%theta
      deta.dbeta[i,] <- X[i,] - t(phi)%*%X[i-ar,] - t(theta)%*%deta.dbeta[i-ma,]
      deta.dphi[i,] <- ynew_ar[i-ar]- Xbeta[i-ar] - t(theta)%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- R[(i-m),] - t(theta)%*%deta.dtheta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rM <- deta.dbeta[(m+1):n,]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT <- diag(mu.eta(eta[(m+1):n]))
    
    r <- function(x){log(1+(log(1/(1-x)))^c_par)}
    s <- function(x){log(log(1/(1-x)))}
    
    a_t <- -(c_par*(log(1/(1-q_t)))^(c_par-1))/((1-q_t)*r(q_t)*exp(r(q_t)))+
      (-c_par*log(1-tau))*(log(1/(1-q_t)))^(c_par-1)/((1-q_t)*(r(q_t))^2*exp(r(q_t)))*r(y1)
    
    b_t <-  1/c_par+s(y1)+(s(y1)*(log(1/(1-y1)))^c_par*(log(1-tau)/r(q_t)-1))/exp(r(y1))-
      ((log(1/(1-q_t)))^c_par*s(q_t))/(r(q_t)*exp(r(q_t)))-
      (log(1-tau)*s(q_t)*(log(1/(1-q_t)))^c_par*r(y1))/((r(q_t))^2*exp(r(q_t)))
    
    ones <- rep(1,(n-m))
    
    Ualpha <- t(v) %*% mT %*% a_t
    Ubeta <- t(rM) %*% mT %*% a_t
    Uphi <-   t(rP) %*% mT %*% a_t
    Utheta <- t(rR) %*% mT %*% a_t
    Uc <- t(b_t) %*% ones
    
    rval <- c(Ualpha,Ubeta,Uphi,Utheta,Uc)
    return(rval[rval!=0])
  }
  #############################################
  
  opt<-optim(initial, loglik, 
             escore.RUBXIIarma, 
             method = "BFGS", #hessian = TRUE,
             control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
  
  
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
    opt<-optim(initial, loglik, 
               #  escore.RUBXIIarma, 
               method = "BFGS", hessian = TRUE,
               control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
    }else{
      warning("IT WORKS WITH NUMERICAL GRADIENT!")
    }
  }
  
  z$conv <- opt$conv
  coef <- (opt$par)[1:(p1+q1+k+2)]
  alpha <- coef[1]
  if(k==0) beta=names_beta=NULL else z$beta <- coef[2:(k+1)]
  if(p1==0) phi=names_phi=NULL else z$phi <- coef[(k+2):(k+p1+1)]
  if(q1==0) theta=names_theta=NULL else z$theta <- coef[(k+p1+2):(k+p1+q1+1)]
  z$c_par <- coef[length(coef)]
  
  # names_par <- c("alpha",names_beta,names_phi,names_theta,"c_par")
  # names(coef)<-names_par
  # z$coeff <- coef
  # J_inv <- solve(-(opt$hessian))
  # z$stderror<-sqrt(diag(J_inv))
  # z$zstat <- abs(z$coef/z$stderror)      # ??? pq módulo aqui? Não seria melhor colocar depois no p-valor?
  # z$pvalues <- 2*(1 - pnorm(z$zstat) )  #2*(1-pt(abs(t_value),df_mod))
  # # z$zstat <- z$coef/z$stderror  
  # # z$pvalues <- 2*(1 - pnorm(abs(z$zstat)) ) 
  # z$loglik <- opt$value
  # z$counts <- as.numeric(opt$counts[1])
  # 
  # if(any(is.na(X)==F))
  # {
  #   z$k<- (p1+q1+k+2)
  #   z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
  #   z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
  #   z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  # }else{
  #   z$k<- (p1+q1+2)
  #   z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
  #   z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
  #   z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  # }
  # 
  # model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  # colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  # z$model <- model_presentation
  
  
  #***************************************************************************************************
  #Fitted values   (NO REGRESSORS)
  if(k==0){                                     ##########(NO REGRESSORS)
    alpha <- as.numeric(coef[1])
    phi <- as.numeric(coef[2:(p1+1)])
    theta <- as.numeric(coef[(p1+2):(p1+q1+1)])
    c_par <- as.numeric(coef[p1+q1+2]) # precision parameter
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$c_par <- c_par
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    if(p1==0) {phi = as.matrix(0);ar=1}
    if(q1==0) {theta = as.matrix(0);ma=1}
    for(i in (m+1):n)
    {
      # etahat[i]<-alpha + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      etahat[i]<-alpha + ynew_ar[(i-1),ar]%*%as.matrix(phi) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    
    q_hat <- linkinv(etahat[(m+1):n])   ### fitted values 
    
    z$fitted <- ts(c(rep(NA,m),q_hat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
    
    z$forecast <- y_prev[(n+1):(n+h)]
  } else{                                        ########### ##########(with REGRESSORS)
    X_hat <- as.matrix(X_hat)
    
    alpha <- as.numeric(coef[1])
    beta <- as.numeric(coef[2:(k+1)])
    phi <- as.numeric(coef[(k+2):(k+p1+1)])
    theta <- as.numeric(coef[(k+p1+2):(k+p1+q1+1)])
    c_par <- as.numeric(coef[length(coef)]) # precision parameter
    
    z$alpha <- alpha
    z$beta <- beta
    z$phi <- phi
    z$theta <- theta
    z$c_par <- c_par
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    if(p1==0) {phi = as.matrix(0);ar=1}
    if(q1==0) {theta = as.matrix(0);ma=1}
    
    Xbeta <- X%*%as.matrix(beta)
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      etahat[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar])%*%as.matrix(phi) +
        t(as.matrix(theta))%*%error[i-ma]
      errorhat[i] <- ynew[i] - etahat[i] 
    }
    q_hat <- linkinv(etahat[(m+1):n])   ### fitted values 
    
    z$fitted <- ts(c(rep(NA,m),q_hat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h))
    y_prev[1:n] <- z$fitted
    
    X_prev <- rbind(X,X_hat)
    
    for(i in 1:h)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) +
        (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta))) +
        (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
    
    
    z$forecast <- y_prev[(n+1):(n+h)]
    
  }
  
  ########################################################
  #      J matrix: Observed information matrix 
  ########################################################
  
  if(q1 == 0){
    R <- matrix(rep(0,(n-m)*1), ncol = 1)
  }else{
    R <- matrix(rep(NA,(n-m)*q1), ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
  }
  
  p_1 = max(p1,1)
  if(p1 == 0){
    P <- matrix(rep(0,(n-m)*p_1), ncol = p_1)
  }else{
    P <- matrix(rep(NA,(n-m)*p1), ncol = p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }
  }
  
  
  if(k == 0){ # No regressors
    k1 = 1
    M = matrix(rep(0,(n-m)*k1), ncol = k1)
  }else{     # with regressors
    k1 <- length(beta)
    M  <- matrix(rep(NA,(n-m)*k1), ncol = k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
    }
  }
  
  # First derivatives of eta
  deta.dalpha = rep(0, n)
  deta.dphi   = matrix(0,  ncol = p_1, nrow = n)
  deta.dtheta = matrix(0,  ncol = q_1, nrow = n)
  deta.dbeta  =  matrix(0, ncol = k1, nrow = n)
  
  for(i in (m+1):n)
  {
    deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
    deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
    deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
  }
  
  v <- deta.dalpha[(m+1):n]
  rP <- deta.dphi[(m+1):n,]
  rR <- deta.dtheta[(m+1):n,]
  rM <- deta.dbeta[(m+1):n,]
  
  mT <- diag(mu.eta(etahat[(m+1):n]))
  
  # useful quantities
  q_t = q_hat
  
  # dell.dq
  r <- function(x){log(1+(log(1/(1-x)))^c_par)}
  s <- function(x){log(log(1/(1-x)))}
  
  a_t <- -(c_par*(log(1/(1-q_t)))^(c_par-1))/((1-q_t)*r(q_t)*exp(r(q_t)))+
    (-c_par*log(1-tau))*(log(1/(1-q_t)))^(c_par-1)/((1-q_t)*(r(q_t))^2*exp(r(q_t)))*r(y1)
  
  # Second-order derivatives of eta
  
  
  if(p1 == 0){
    mI.phi.beta     = matrix(0, ncol = p_1, nrow = k1)
  }else{
    mI.phi.beta     = matrix(0, ncol = p_1, nrow = k1)
    for(t in 1:k1){
      deta.dbeta.dphi <- matrix(0,ncol=p1,nrow=n)
      for(j in 1:p1)
      {
        for(i in (m+1):n)
        {
          deta.dbeta.dphi[i,j] <- -(X[i-j,t]) - theta%*%deta.dbeta.dphi[i-ma,j]
        }
      }   
      mI.phi.beta[t,] <- t(a_t)%*%mT%*%deta.dbeta.dphi[(m+1):n,]
    }
  }
  
  
  if(q1 == 0){  #AR
    deta.dalpha.dtheta = matrix(0, ncol = q_1, nrow = n)
    mA = deta.dalpha.dtheta[(m+1):n,]
    
    mJ.theta.beta   = matrix(0, ncol = q_1, nrow = k1)
    mK.theta.phi    = matrix(0, ncol = q_1, nrow = p_1)
    mL.theta.theta  = matrix(0, ncol = q_1, nrow = q_1)
  }else{
    
    deta.dalpha.dtheta = matrix(0, ncol = q_1, nrow = n)
    mJ.theta.beta   = matrix(0, ncol = q_1, nrow = k1)
    mK.theta.phi    = matrix(0, ncol = q_1, nrow = p_1)
    mL.theta.theta  = matrix(0, ncol = q_1, nrow = q_1)
    
    for(j in 1:q1)
    {
      for(i in (m+1):n)
      {
        deta.dalpha.dtheta[i,j] <- -deta.dalpha[i-j] - 
          theta%*%deta.dalpha.dtheta[i-ma,j] 
      }
    }
    mA = deta.dalpha.dtheta[(m+1):n,]
    
    for(t in 1:k1){
      deta.dtheta.dbeta <- matrix(0,ncol=q_1,nrow=n)
      for(j in 1:q_1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dbeta[i,j] <- -deta.dbeta[i-j, t] - theta%*%deta.dtheta.dbeta[i-ma,j]  
        }
      }   
      mJ.theta.beta[t,] <- t(a_t)%*%mT%*%deta.dtheta.dbeta[(m+1):n,]
    }
    
    for(t in 1:p_1){  # MA model
      deta.dphi.dtheta <- matrix(0,ncol=q_1,nrow=n)
      for(j in 1:q_1)
      {
        for(i in (m+1):n)
        {
          deta.dphi.dtheta[i,j] <- -deta.dphi[i-j, t] - theta%*%deta.dphi.dtheta[i-ma,j]  
        }
      }   
      mK.theta.phi[t,] <- t(a_t)%*%mT%*%deta.dphi.dtheta[(m+1):n,]
    }
    
    for(t in 1:q_1){
      deta.dtheta.dtheta <- matrix(0,ncol = q_1, nrow = n)
      for(j in 1:q_1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dtheta[i,j] <- -deta.dtheta[i-t,j] - deta.dtheta[i-j, t] - theta%*%deta.dtheta.dtheta[i-ma,j]  
        }
      }   
      mL.theta.theta[t,] <- t(a_t)%*%mT%*%deta.dtheta.dtheta[(m+1):n,]
    }
  }
  
  mI = mI.phi.beta    
  mJ = mJ.theta.beta  
  mK = mK.theta.phi   
  mL = mL.theta.theta 
  
  # Second-order and cross derivatives 
  # dell2.dq2
  h_func <- function(x,z1){
    1/x*(x-1)*log(1-z1)
  }
  l <- function(x){ exp(r(x))*1/r(x)}
  w <- function(x){log(1+(log(1/(1-x)))^c_par)}
  
  delta_t_qq =  ((c_par^2 *(exp(s(q_t)))^(2*c_par-2))/
                   ((1-q_t)^2*(exp(w(q_t)))^2*w(q_t)))*
    (
      1+(1/w(q_t))-(1/c_par)*(c_par-1)*(exp(s(q_t)))^(-c_par)*exp(w(q_t))-
        (1/c_par)*(exp(s(q_t)))^(1-c_par)*exp(w(q_t))+
        2*log(1-tau)*(w(q_t))^(-2)*w(y1)+
        log(1-tau)*(1/w(q_t))*w(y1)-
        h_func(c_par,tau)*(exp(s(q_t)))^(-c_par)*l(q_t)*w(y1)-
        h_func(c_par,tau)*(1/(c_par-1))*
        (exp(s(q_t)))^(1-c_par)*l(q_t)*w(y1)
    )
  
  
  # dell2.dqc
  h_star <- function(x,z1){ log(1-z1)/w(x)}
  
  delta_t_qc =   -((c_par*s(q_t)*(exp(s(q_t)))^(c_par - 1))/
                     ((1 - q_t)*exp(w(q_t))* w(q_t)) )*(
                       1+(1/c_par)*(s(q_t))^(-1)-(exp(s(q_t)))^c_par*(exp(w(q_t)))^(-1)-
                         (exp(s(q_t)))^c_par*(exp(w(q_t)))^(-1)*(w(q_t))^(-1)+
                         log(1-tau)*(
                           (s(q_t))^(-1)*s(y1)*(exp(s(y1)))^c_par*(w(q_t))^(-1)*(exp(w(y1)))^(-1)+
                             (1/c_par)*(s(q_t))^(-1)*w(y1)*(w(q_t))^(-1)+
                             w(y1)*(w(q_t))^(-1)-
                             2*(exp(s(q_t)))^c_par*w(y1)*(exp(w(q_t)))^(-1)*(w(q_t))^(-2)-
                             (exp(s(q_t)))^c_par*w(y1)*(exp(w(q_t)))^(-1)*(w(q_t))^(-1)
                         ) 
                     ) 
  
  
  # dell2.dcc
  delta_t_cc =  -((s(q_t)*log(1 - tau)*s(y1)* (exp(s(q_t)))^c_par *(exp(s(y1)))^c_par)/
                    (exp(r(q_t))*(r(q_t))^2* exp(r(y1))))*(
                      2-s(q_t)*(log(1-tau))^(-1)*(s(y1))^(-1)*(exp(s(y1)))^(-c_par)*
                        (  (exp(s(q_t)))^c_par-w(q_t)   )*(exp(w(q_t)))^(-1)*exp(w(y1))+
                        (s(y1))^-1*(exp(s(y1)))^(-c_par)*s(q_t)*w(y1)*exp(w(y1))-
                        (s(y1))^(-1)*(exp(s(y1)))^(-c_par)*s(q_t)*(exp(s(q_t)))^c_par*w(y1)*exp(w(y1))*(exp(w(q_t)))^(-1)-
                        2*(s(y1))^(-1)*(exp(s(y1)))^(-c_par)*s(q_t)*(exp(s(q_t)))^c_par*
                        (exp(w(q_t)))^(-1)*w(q_t)^(-1)*exp(w(y1))*w(y1)
                    )-
    c_par^(-2)+h_star(q_t,tau)*(s(y1))^2*(exp(s(y1)))^c_par*(exp(w(y1)))^(-1)*
    (
      1-(exp(s(y1)))^c_par*(exp(w(y1)))^(-1)
    )
  
  #mathcal A = mA
  g2lin  <- function(x){return( (2*x-1)/((x^2)*(1-x)^2) )}
  kappa =   (-g2lin(q_t))*(mu.eta(etahat[(m+1):n]))^3
  w_t = delta_t_qq * (diag(mT))^2 + a_t * kappa
  W = diag(w_t)
  
  # J matrix elements
  # Line 1 - alpha
  J.alpha.alpha   = -t(v) %*% W %*% v
  J.alpha.beta    = -t(v) %*% W %*% rM
  J.beta.alpha    = t(J.alpha.beta)
  J.alpha.phi     = -t(v) %*% W %*% rP
  J.phi.alpha     = t(J.alpha.phi)
  J.alpha.theta   = -t(v) %*% W %*% rR- t(a_t) %*% mT %*% mA
  J.theta.alpha   = t(J.alpha.theta)
  J.alpha.c       = -t(v) %*% mT %*% ( diag(delta_t_qc) %*% rep(1,length(delta_t_qc))  ) #para deixar D1 #ver porque é a cruzada aqui
  J.c.alpha       = t(J.alpha.c)
  
  # Line 2  - beta
  J.beta.beta     = -t(rM) %*% W %*% rM
  J.beta.phi      = -t(rM) %*% W %*% rP - mI
  J.phi.beta      = t(J.beta.phi) 
  J.beta.theta    = -t(rM) %*% W %*% rR - mJ
  J.theta.beta    = t(J.beta.theta)
  J.beta.c        = -t(rM) %*% mT %*% delta_t_qc
  J.c.beta        = t(J.beta.c)
  
  # Line 3  
  J.phi.phi       = -t(rP) %*% W %*% rP
  J.phi.theta     = -t(rP) %*% W %*% rR- mK
  J.theta.phi     = t(J.phi.theta)
  J.phi.c         = -t(rP) %*% mT %*% delta_t_qc
  J.c.phi         = t(J.phi.c)
  
  # Line 4
  J.theta.theta   = -t(rR) %*% W %*% rR- mL
  J.theta.c       = -t(rR) %*% mT %*% delta_t_qc
  J.c.theta       = t(J.theta.c)
  
  # Line 5
  J.c.c = -sum(delta_t_cc)
  
  # alpha
  line.alpha = cbind(J.alpha.alpha ,
                     J.alpha.beta  ,
                     J.alpha.phi   ,
                     J.alpha.theta ,
                     J.alpha.c 
  )
  line.alpha = line.alpha[, colSums(line.alpha)!= 0]
  
  # beta
  line.beta = cbind(J.beta.alpha,
                    J.beta.beta ,
                    J.beta.phi  ,
                    J.beta.theta,
                    J.beta.c   
  )
  line.beta = line.beta[, colSums(line.beta)!= 0]
  
  # phi
  line.phi = cbind(J.phi.alpha,
                   J.phi.beta,
                   J.phi.phi  ,
                   J.phi.theta,
                   J.phi.c    )
  line.phi = line.phi[, colSums(line.phi)!= 0]
  
  # theta
  line.theta = cbind(J.theta.alpha,
                     J.theta.beta,
                     J.theta.phi,
                     J.theta.theta,
                     J.theta.c    )
  line.theta = line.theta[, colSums(line.theta)!= 0]
  
  # c_par
  line.c = cbind(J.c.alpha,
                 J.c.beta,
                 J.c.phi,
                 J.c.theta,
                 J.c.c)
  line.c = line.c[, colSums(line.c)!= 0]
  
  ##############################
  # Observed information matrix
  ##############################
  J = rbind(line.alpha,
            line.beta,
            line.phi,
            line.theta,
            line.c)
  
  # round(J,2)
  # library(numDeriv)
  # round(-jacobian(escore.UDARMA,opt$par),2)
  # 
  # round(J-(-jacobian(escore.UDARMA,opt$par)),2)
  ########################################################
  
  names_par   = c("alpha",names_beta,names_phi,names_theta,"c_par")
  names(coef) = names_par
  z$coeff     =   coef
  J_inv       = solve(J) #solve(-(opt$hessian))
  z$stderror  = sqrt(diag(J_inv))
  # z$zstat     = abs(z$coef/z$stderror)      
  # z$pvalues   = 2*(1 - pnorm(z$zstat) )   
  z$zstat = z$coef/z$stderror
  z$pvalues = 2*(1 - pnorm(abs(z$zstat)) )
  z$loglik    = opt$value
  z$counts    = as.numeric(opt$counts[1])
  
  if(any(is.na(X) == F))
  {
    z$k   = (p1+q1+k+2)
    z$aic = -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic = -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq  = -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k   = (p1+q1+2)
    z$aic = -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic = -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq  = -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }
  
  # if(any(is.na(X)==F))
  # {
  #   z$k= (p1+q1+k+2)
  #   z$aic = -2*(z$loglik)+2*(z$k)
  #   z$bic = -2*(z$loglik)+log(n)*(z$k)
  #   z$hq = -2*(z$loglik)+log(log(n))*(z$k)
  # }else{
  #   z$k= (p1+q1+2)
  #   z$aic = -2*(z$loglik)+2*(z$k)
  #   z$bic = -2*(z$loglik)+log(n)*(z$k)
  #   z$hq = -2*(z$loglik)+log(log(n))*(z$k)
  # }
  # 
  model_presentation = cbind(round(z$coef,4),
                             round(z$stderror,4),
                             round(z$zstat,4),
                             round(z$pvalues,4))
  
  colnames(model_presentation) = c("Estimate",
                                   "Std. Error",
                                   "z value",
                                   "Pr(>|z|)")
  z$model = model_presentation
  
  z$J <- J
  
  # if (isTRUE(check_hessian)) {
  #   H_num <- optimHess(opt$par, loglik)
  #   J_num <- -H_num
  #   hess_diff <- J - J_num
  #   max_abs <- max(abs(hess_diff))
  #   max_rel <- max(abs(hess_diff)/(abs(J_num) + 1e-8))
  #   z$hessian_check <- list(
  #     J_numeric = J_num,
  #     diff = hess_diff,
  #     max_abs = max_abs,
  #     max_rel = max_rel
  #   )
  #   if (diag > 0) {
  #     message("Hessian check (analytical vs numeric): max_abs = ",
  #             signif(max_abs, 3), ", max_rel = ", signif(max_rel, 3))
  #   }
  # }
  # 
  
  # randomized quantile residuals 
  z$residuals <- as.vector(qnorm(cdf_RUBXII(y[(m+1):n],z$fitted[(m+1):n],
                                            z$c_par, tau =tau)))
  residc <- z$residuals 
  
  
  #*********************************************************************************************
  ###################################################
  ######### GRAPHICS ################################
  
  if(diag>0)
  {
    
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," SIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    t<-seq(-5,n+h,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",
         ylim=c(-4,4))
    lines(t,rep(-3,n+h+6),lty=2,col=1)
    lines(t,rep(3,n+h+6),lty=2,col=1)
    lines(t,rep(-2,n+h+6),lty=3,col=1)
    lines(t,rep(2,n+h+6),lty=3,col=1)
    
    
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="Fitted values",ylab="Observed data",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
         xlab="Fitted values",ylab="Residuals")
    #lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    
    densidade<-density(residc)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(residc,ylab="ACF",xlab="Lag") 
    
    pacf(residc,ylab="PACF",xlab="Lag") 
    
    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)
    qqnorm(residc, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="Normal quantiles",ylab="Empirical quantiles")
    #qqnorm(residc,pch=1,frame=T, main = "QQ-plot",
    #      make.plot = T, lwd=1)
    #lines(c(-10,10),c(-10,10),lty=2)
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="Serie",xlab="Time",ylim=c(min(y),max(y)))
    lines(z$fitted,col="blue",lty=2)
    legend("topright",c("Observed data","Predicted median"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
    
    # fim<-end(y)[1]+end(y)[2]/12                   
    # 
    # y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    # par(mfrow=c(1,1))
    # plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
    # abline(v=fim,lty=2)
    # lines(y)
    
    w1<-5
    h1<-4
    
    if(diag>1)
    {
      postscript(file = "histogram.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        hist(y,ylab="Frequency",xlab="y",main="",col=8) 
      }
      dev.off()
      
      postscript(file = "sazonality.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        monthplot(y,ylab="Stocked hydroelectric energy",xlab="Months") 
      }
      dev.off()
      
      postscript(file = "y_FAC.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        Acf(residc,ylab="ACF",xlab="Lag") 
      }
      dev.off()
      
      postscript(file = "y_FACP.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        Pacf(y,ylab="PACF",xlab="Lag")
      }
      dev.off()
      
      postscript(file = "qq_plot.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(residc, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="Normal quantiles",ylab="Empirical quantiles")
        #abline(0,1,col=2)
      }
      dev.off()
      
      postscript(file = "resid_v_ind.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",
             ylim=c(-4,4))
        lines(t,rep(-3,n+h+6),lty=2,col=1)
        lines(t,rep(3,n+h+6),lty=2,col=1)
        lines(t,rep(-2,n+h+6),lty=3,col=1)
        lines(t,rep(2,n+h+6),lty=3,col=1)
      }
      dev.off()
      
      
      postscript(file = "resid_FAC.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        acf(residc,ylab="ACF",xlab="Lag") 
      }
      dev.off()
      
      postscript(file = "resid_FACP.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        pacf(residc,ylab="PACF",xlab="Lag")
      }
      dev.off()
      
      
      postscript(file = "adjusted.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="Serie",xlab="Time")
        lines(z$fitted,col=2,lty=2)
        legend("topright",c("Observed data","Predicted median"),
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,2), cex=.8)
        
      }
      dev.off()
      
    }    
  }  
  #*****************************************************************************************************
  return(z)
}
