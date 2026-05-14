

CP_MTS = function(Y, xi = NULL, Rank = NULL, lag.k = 20, lag.ktilde = 10,
                  method = c("CP.Direct","CP.Refined","CP.Unified"),
                  thresh1 = FALSE, thresh2 = FALSE, thresh3 = FALSE,
                  delta1 = 2 * sqrt(log(dim(Y)[2] * dim(Y)[3]) / dim(Y)[1]),
                  delta2 = delta1, delta3 = delta1,c1=0,c2=0,c3=0){
  n <- dim(Y)[1]; p <- dim(Y)[2]; q <- dim(Y)[3];
  if(is.null(xi)){
    xi <- est.xi(Y)$xi
  }
  method <- match.arg(method)
  if(method == "CP.Direct"){
    S_yxi_1 <- Autocov_xi_Y(Y, xi, k = 1, thresh = thresh1, delta = delta1)
    S_yxi_2 <- Autocov_xi_Y(Y, xi, k = 2, thresh = thresh1, delta = delta1)
    if(p > q){
      ##(1) estimation of d
      K1 <- t(S_yxi_1) %*% S_yxi_1
      eg1 <- eigen(K1)
      w <- eg1$values
      ww <- w[-1] / w[-length(w)]
      d <- which(ww == min(ww[1:floor(0.75 * q)]))
      if(!is.null(Rank)){
        if(!is.null(Rank$d)){
          d <- Rank$d
        }
        else{stop("List Rank without d, use Rank=list(d=?)")}
      }
      if (d > 1){
        K1 <- eg1$vectors[, 1:d]%*%diag(eg1$values[1:d])%*%t(eg1$vectors[, 1:d]);
      }else{
        K1 <- eg1$vectors[, 1]%*%diag(eg1$values[1], 1)%*%t(eg1$vectors[, 1]);
      }
      K2 <- t(S_yxi_1) %*% S_yxi_2;
      
      ##(2) estimation of A and B
      Geg <- geigen::geigen(K2, K1);
      evalues <- Geg$values[which(Mod(Geg$values) <= 10^5 & Geg$values != 0)]
      
      Bl <- Geg$vectors[, which(Geg$values %in% evalues), drop = FALSE]
      A <- apply(S_yxi_1 %*% Bl, 2, l2s)
      Al <- t(MASS::ginv(A))
      B <- apply(t(S_yxi_1) %*% Al, 2, l2s)
    }else{
      
      ##(1) estimation of d
      K1 <- S_yxi_1 %*% t(S_yxi_1)
      eg1 <- eigen(K1)
      w <- eg1$values
      ww <- w[-1] / w[-length(w)]
      d <- which(ww == min(ww[1:floor(0.75 * p)]))
      if(!is.null(Rank)){
        if(!is.null(Rank$d)){
          d <- Rank$d
        }
        else{stop("List Rank without d, use Rank=list(d=?)")}
      }
      
      if (d > 1){
        K1 <- eg1$vectors[, 1:d] %*% diag(eg1$values[1:d]) %*% t(eg1$vectors[, 1:d]);
      }else{
        K1 <- eg1$vectors[, 1] %*% diag(eg1$values[1], 1) %*% t(eg1$vectors[, 1]);
      }
      K2 <- S_yxi_1 %*% t(S_yxi_2);
      ##(2) estimation of A and B
      Geg <- geigen::geigen(K2, K1);
      evalues <- Geg$values[which(Mod(Geg$values) <= 10^5 & Geg$values!=0)]
      Al <- Geg$vectors[, which(Geg$values %in% evalues), drop = FALSE]
      B <- apply(t(S_yxi_1) %*% Al, 2, l2s)
      Bl <- t(MASS::ginv(B))
      A <- apply(S_yxi_1 %*% Bl, 2, l2s)
    }
    ##(3) estimation of Xt
    H <- matrix(NA, p * q, d)
    for (ii in 1:d) {
      H[, ii] <- B[, ii] %x% A[, ii]
    }
    f <- Vec.tensor(Y) %*% H %*% MASS::ginv(t(H) %*% H)
    
    if(is.complex(A) == T || is.complex(B) == T ){
      A <- Complex2Real(A)
      B <- Complex2Real(B)
      f <- Complex2Real(f)
    }
    # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:", method))
    con <- structure(list(A = A,B = B,f = f,Rank = list(d = d), method = method),
                     class = "mtscp")
    return(con)
  }
  if(method == "CP.Refined"){
    ##(1) estimation of P,Q and d
    M1 = M2 <- 0
    dmax <- round(min(p, q) * 0.75)
    
    for (kk in 1:lag.k){
      S_yxi_k <- Autocov_xi_Y(Y, xi, k = kk, thresh = thresh1, delta = delta1)
      M1 <- M1 + S_yxi_k %*% t(S_yxi_k)
      M2 <- M2 + t(S_yxi_k) %*% S_yxi_k
    }
    ev_M1 <- eigen(M1)
    ev_M2 <- eigen(M2)
    
    d1 <-  which.max(ev_M1$values[1:dmax] / ev_M1$values[2:(dmax + 1)])
    d2 <-  which.max(ev_M2$values[1:dmax] / ev_M2$values[2:(dmax + 1)])
    
    d <- ifelse(p > q, d1, d2)
    
    if(!is.null(Rank)){
      if(!is.null(Rank$d)){
        d <- Rank$d
      }
      else{stop("List Rank without d, use Rank=list(d=?)")}
    }
    P <- ev_M1$vectors[, 1:d, drop = FALSE]
    Q <- ev_M2$vectors[, 1:d, drop = FALSE]
    
    if(d == 1){
      A <- as.matrix(P)
      B <- as.matrix(Q)
      
      f <- vector()
      for (tt in 1:n) {
        f[tt] <- t(A) %*% Y[tt, , ] %*% B
      }
      f <- as.matrix(f)
      
    }else{
      ##(2) estimation of U and V
      Z <- array(NA, dim = c(n, d, d))
      for (tt in 1:n) {
        Z[tt, , ] <- t(P) %*% Y[tt, , ] %*% Q
      }
      xi <- est.xi(Z)
      if(thresh2){
        w_hat <- xi$w_hat
        Xi <- diag(1, p) %x% ((Q %x% P) %*% as.matrix(w_hat))
        sigma_ycheck_1 <- Sigma_Ycheck(Y, 1)
        sigma_ycheck_1 <- thresh_C(sigma_ycheck_1, delta2)
        sigma_ycheck_2 <- Sigma_Ycheck(Y, 2)
        sigma_ycheck_2 <- thresh_C(sigma_ycheck_2, delta2)
        S_Zxi_1 <- t(P) %*% t(Xi) %*% sigma_ycheck_1 %*% Q
        S_Zxi_2 <- t(P) %*% t(Xi) %*% sigma_ycheck_2 %*% Q
      }
      else{
        S_Zxi_1 <- Autocov_xi_Y(Z, xi$xi, k = 1)
        S_Zxi_2 <- Autocov_xi_Y(Z, xi$xi, k = 2)
      }
      
      vl <- eigen(MASS::ginv(t(S_Zxi_1) %*% S_Zxi_1) %*% t(S_Zxi_1) %*% S_Zxi_2)$vectors ##MASS
      ul <- eigen(MASS::ginv(S_Zxi_1 %*% t(S_Zxi_1)) %*% S_Zxi_1 %*% t(S_Zxi_2))$vectors
      
      U <- apply(S_Zxi_1 %*% vl, 2, l2s)
      V <- apply(t(S_Zxi_1) %*% ul, 2, l2s)
      
      ##(3) estimation of A and B
      A <- P %*% U
      B <- Q %*% V
      
      ##(4) estimation of Xt
      W <- matrix(NA, d^2, d)
      
      for (ii in 1:d) {
        W[, ii] <- V[, ii] %x% U[, ii]
      }
      
      f <- Vec.tensor(Z) %*% W %*% solve(t(W) %*% W)
      
      if(is.complex(A) == T || is.complex(B) == T ){
        A <- Complex2Real(A)
        B <- Complex2Real(B)
        f <- Complex2Real(f)
      }
    }
    # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:", method))
    con <- structure(list(A = A,B = B,f = f,Rank = list(d = d), method = method),
                     class = "mtscp")
    return(con)
  }
  if(method == "CP.Unified"){
    ##(1) estimation of P,Q and d1,d2
    if(is.null(Rank)){
      
      PQ_hat_tol <- est.d1d2.PQ(Y, xi, K = lag.k, thresh = thresh1, delta = delta1,c1 = c1,c2 = c2)
      
      d1 <- PQ_hat_tol$d1_hat
      d2 <- PQ_hat_tol$d2_hat
      d  <- NULL
      P  <- PQ_hat_tol$P_hat
      Q  <- PQ_hat_tol$Q_hat
      
      if(d1 == 1 || d2 == 1){d <- d1 * d2}
      
    }else{
      if(all(!is.null(Rank$d1), !is.null(Rank$d1), !is.null(Rank$d2))){
        d  <- Rank$d
        d1 <- Rank$d1
        d2 <- Rank$d2
      }
      else{stop("List Rank without d, d1 and d2, use Rank=list(d=?, d1=?, d2=?)")}
      
      PQ_hat_tol <- est.PQ(Y, xi, d1, d2, K = lag.k, thresh = thresh1, delta = delta1)
      
      P   <- PQ_hat_tol$P_hat
      Q   <- PQ_hat_tol$Q_hat
      
    }
    ##(2) estimation of W* = (v1*u1,v2*u2,...,vd*ud)H = WH
    if(d1 == 1 & d2 == 1){
      d <- 1
      f <- vector()
      for (tt in 1:n) {
        f[tt] = t(P) %*% Y[tt, , ] %*% Q
      }
      A <- P
      B <- Q
      
      # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      rank <- list(d = d, d1 = d1, d2 = d2)
      con <- structure(list(A = as.matrix(A), B = as.matrix(B), 
                            f = as.matrix(f), Rank = rank, method = method),
                       class = "mtscp")
      
      return(con)
    }else{
      if(is.null(d)){
        W_hat_tol  <-  est.d.Wf(Y, P, Q, Ktilde = lag.ktilde, thresh = thresh3, delta = delta3,c3 = c3)
        d          <-  W_hat_tol$d_hat
        W          <-  W_hat_tol$W_hat
        f          <-  W_hat_tol$f_hat
      }else{
        d          <-  d
        W_hat_tol  <-  est.Wf(Y, P, Q, d, Ktilde = lag.ktilde, thresh = thresh3, delta = delta3)
        W          <-  W_hat_tol$W_hat
        f          <-  W_hat_tol$f_hat
      }
      
      ##(3) estimation of U and V
      if(d1 == 1 || d2 == 1){
        Theta <- NULL
        if(d1 == 1){
          U <- 1;
          V <- W;
          warning("d1 equal to 1, V cannot be identified uniquely!")
        }
        if(d2 == 1){
          U <- W;
          V <- 1;
          warning("d2 equal to 1, U cannot be identified uniquely!")
        }
        if(d1 == 1 & d2 == 1){
          U <- 1;
          V <- 1;
        }
        U <- as.matrix(U)
        V <- as.matrix(V)
      }else{
        
        UV_hat_tol <- est.UV.JAD.BSel(W, d1, d2, d)
        
        U          <- UV_hat_tol$U
        V          <- UV_hat_tol$V
        Theta      <- UV_hat_tol$Theta
      }
      
      
      ##(4) estimation of A and B
      A <- P %*% U
      B <- Q %*% V
      # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      rank <- list(d = d, d1 = d1, d2 = d2)
      con <- structure(list(A = as.matrix(A), B = as.matrix(B),
                            f = as.matrix(f), Rank = rank, method = method),
                       class = "mtscp")
      
      return(con)
    }
    
  }
}

rho2.loss = function(A_hat,A){
  max(apply(1-(t(A_hat)%*%A)^2,2,min))
}

l2s = function(x){x/sqrt(sum(x^2))}

fnorm = function(x){sqrt(sum(x^2))}

MAE = function(x){sum(abs(x))/length(x)}


Mat.tensor = function(Y,j){ # fold on j-th mode
  t(apply(Y,j,c))
}


Complex2Real = function(A){
  REA = round(Re(A),8)
  IMA = round(Im(A),8)
  
  real.index  =  which( colSums(IMA)   == 0)
  
  if(length(real.index) == 0){
    complex_real  =  REA
    complex_image =  IMA
    
    complex_take  = which(duplicated(complex_real[1,]) == T)
    
    real = as.matrix(complex_real[,complex_take])
    
    img  = as.matrix(complex_image[,complex_take])
    
    new_A = cbind(real,img)
  }else{
    real.vector   = REA[,real.index]
    
    complex_real  =  REA[,-real.index]
    complex_image =  IMA[,-real.index]
    
    complex_take  = which(duplicated(complex_real[1,]) == T)
    
    real = as.matrix(complex_real[,complex_take])
    
    img  = as.matrix(complex_image[,complex_take])
    
    
    new_A = cbind(real.vector,real,img)
  }
  
  return(new_A)
}


Vec.tensor = function(Y){
  p = dim(Y)[2];q = dim(Y)[3];
  
  if(p == q & q == 1){
    Y_tilde = apply(Y,1,FUN = as.vector)
  }else{
    Y_tilde = t(apply(Y,1,FUN = as.vector))
  }
  return(Y_tilde)
  
}


DGP.CP = function(n,p,q,d,d1,d2){
  
  par_A = c(-3,3)
  par_B = c(-3,3)
  par_X = c(0.6,0.95)
  par_E = 1
  
  Input = list(n = n,p = p,q = q,d1 = d1,d2 = d2,d = d)
  
  A_inl = matrix(runif(p*d,par_A[1],par_A[2]),p,d)
  B_inl = matrix(runif(q*d,par_B[1],par_B[2]),q,d)
  
  svd_A = svd(A_inl)
  P = svd_A$u[,1:d1]
  
  A = (svd_A$u[,1:d1]) %*% diag(svd_A$d[1:d1],nrow = d1,ncol = d1) %*%t(svd_A$v[,1:d1])
  
  U = t(P)%*%apply(A,2,l2s)
  
  A_s = P%*%U
  
  
  svd_B = svd(B_inl)
  Q = svd_B$u[,1:d2]
  B = (svd_B$u[,1:d2]) %*% diag(svd_B$d[1:d2],nrow = d2,ncol = d2) %*% (t(svd_B$v[,1:d2]))
  V = t(Q)%*%apply(B,2,l2s)
  
  B_s = Q%*%V
  
  
  W = matrix(NA,d1*d2,d)
  
  for (ii in 1:d) {
    W[,ii] = V[,ii]%x%U[,ii]
  }
  
  X = array(0,c(n,d,d))
  X_m = matrix(NA,n,d)
  
  signal = runif(d,-1,1)
  signal[signal>=0] <-  1
  signal[signal< 0] <- -1
  par_ar = runif(d,par_X[1],par_X[2])
  for (ii in 1:d) {
    xd = arima.sim(model = list(ar = par_ar[ii]*signal[ii]),n = n)
    X_m[,ii] =  xd*(apply(A, 2, fnorm)*apply(B, 2, fnorm))[ii]
    X[,ii,ii] <- xd*(apply(A, 2, fnorm)*apply(B, 2, fnorm))[ii]
    
  }
  
  S_m = X_m%*%t(W)
  W_star = svd(S_m)$v[,1:d]
  
  Y = S = array(NA,dim = c(n,p,q))
  for (tt in 1:n) {
    S[tt,,] <- A_s%*%X[tt,,]%*%t(B_s)
    Y[tt,,] <- S[tt,,] + matrix(rnorm(p*q,0,par_E),p,q)
  }
  
  return(list(Y      =  Y,
              A      =  A_s,
              B      =  B_s,
              X      =  X
  ))
  
}


Autocov_xi_Y = function(Y, xi, k, thresh = FALSE, delta = NULL){
  
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  q <- dim(Y)[3]
  # k <- lag.k
  
  Y_mean <- 0
  xi_mean <- 0
  for (ii in 1:n) {
    Y_mean <- Y_mean + Y[ii,,]
    xi_mean <- xi_mean + xi[ii]
  }
  Y_mean  <- Y_mean/n
  xi_mean <- xi_mean/n
  
  Sigma_Y_xi_k <- 0
  for (ii in (k+1):n) {
    Sigma_Y_xi_k <- Sigma_Y_xi_k + (Y[ii,,] - Y_mean)*(xi[ii-k] - xi_mean)
  }
  Sigma_Y_xi_k <- Sigma_Y_xi_k/(n-k)
  if(thresh){
    Sigma_Y_xi_k <- thresh_C(Sigma_Y_xi_k, delta)
  }
  return(Sigma_Y_xi_k)
}


est.eta  = function(Y,thresh_per = 0.99){
  
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  eta.mat = Vec.tensor(Y)
  
  if(n > p*q){
    eig_eta.mat = eigen(MatMult(t(eta.mat),eta.mat))
    cfr =  cumsum(eig_eta.mat$values)/sum(eig_eta.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = d_hat
    w_hat = eig_eta.mat$vectors[,1:d_fin]
    
    eta.f = eta.mat%*%w_hat
    
    eta   = rowMeans(eta.f)
    
  }else{
    eig_eta.mat = eigen(MatMult(eta.mat,t(eta.mat)))
    cfr = cumsum(eig_eta.mat$values)/sum(eig_eta.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = d_hat
    eta.f1 = as.matrix(eig_eta.mat$vectors[,1:d_fin])
    
    weight = sqrt(eig_eta.mat$values[1:d_fin])
    
    eta.f1 = eta.f1%*%diag(weight)
    
    eta = rowMeans(eta.f1)
    
  }
  return(eta)
}
est.xi  = function(Y, thresh_per = 0.99, d_max = 20){
  
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  xi.mat = Vec.tensor(Y)
  
  if(n > p*q){
    eig_xi.mat = eigen(MatMult(t(xi.mat),xi.mat))
    cfr =  cumsum(eig_xi.mat$values)/sum(eig_xi.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    w_hat = eig_xi.mat$vectors[,1:d_fin, drop=FALSE]
    if(d_fin == 1){
      sign_value = adjust_sign(w_hat[, 1])
      w_hat = w_hat * sign_value
    }else{
      column_signs = apply(w_hat, 2, adjust_sign)
      w_hat = w_hat %*% diag(column_signs)
    }
    xi.f = xi.mat%*%w_hat
    xi   = rowMeans(xi.f)
    w_hat = rowMeans(w_hat)
  }else{
    eig_xi.mat = eigen(MatMult(xi.mat,t(xi.mat)))
    cfr = cumsum(eig_xi.mat$values)/sum(eig_xi.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    xi.f1 = as.matrix(eig_xi.mat$vectors[,1:d_fin, drop=FALSE])
    if(d_fin == 1){
      sign_value = adjust_sign(xi.f1[, 1])
      xi.f1 = xi.f1 * sign_value
    }else{
      column_signs = apply(xi.f1, 2, adjust_sign)
      xi.f1 = xi.f1 %*% diag(column_signs)
    }
    weight = sqrt(eig_xi.mat$values[1:d_fin])
    xi.f1 = xi.f1%*%diag(weight)
    xi = rowMeans(xi.f1)
    w_hat = as.matrix(t(t(xi.f1)%*%xi.mat))
    w_hat = rowMeans(w_hat)
  }
  return(list(xi=xi, w_hat = w_hat))
}


tensor.est.xi = function(Y,d_max = 10,thresh_per = 0.99, random = F){
  
  n = dim(Y)[1];D = dim(Y)[-1]
  
  Y.mat = Mat.tensor(Y,1)
  
  if(n > prod(D)){
    eig_Y.mat = eigen(HDTSA:::MatMult(t(Y.mat),Y.mat))
    cfr =  cumsum(eig_Y.mat$values)/sum(eig_Y.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    
    if(random == T){
      w_hat = (eig_Y.mat$vectors[,1:d_fin])%*%pracma:::randortho(d_fin)
    }else{
      w_hat = (eig_Y.mat$vectors[,1:d_fin])
    }
    
    xi.f = scale(Y.mat%*%w_hat)
    xi   = rowMeans(xi.f)
    
  }else{
    eig_Y.mat = eigen(HDTSA:::MatMult(Y.mat,t(Y.mat)))
    cfr = cumsum(eig_Y.mat$values)/sum(eig_Y.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    
    if(random == T){
      xi.f = (as.matrix(eig_Y.mat$vectors[,1:d_fin]))%*%pracma:::randortho(d_fin)
    }else{
      xi.f = as.matrix(eig_Y.mat$vectors[,1:d_fin])
    }
    
    weight = sqrt(eig_Y.mat$values[1:d_fin])
    
    xi.f = xi.f   #%*%diag(weight)
    
    xi = rowMeans(xi.f)
    
  }
  return(xi)
}



est.d1d2.PQ = function(Y,xi,K = 10, thresh = FALSE, delta = NULL,c1 = 0,c2 = 0){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  d2_list = d1_list =vector()
  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  P_list = Q_list = list()
  for (kk in 1:K){
    
    S_yxi_k = Autocov_xi_Y(Y,xi, k = kk, thresh = thresh, delta = delta)
    
    M1 = M1 + S_yxi_k%*%t(S_yxi_k)
    M2 = M2 + t(S_yxi_k)%*%S_yxi_k
    
    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)
    
    d1_list[kk] =  which.max((ev_M1$values[1:dmax]+c1)/(ev_M1$values[2:(dmax+1)]+c1))
    d2_list[kk] =  which.max((ev_M2$values[1:dmax]+c2)/(ev_M2$values[2:(dmax+1)]+c2))
    
    P_list[[kk]] = ev_M1$vectors[,1:(d1_list[kk]), drop = FALSE]
    Q_list[[kk]] = ev_M2$vectors[,1:(d2_list[kk]), drop = FALSE]
    
  }
  
  d1_list[1] = 0
  d2_list[1] = 0
  
  
  d1_hat =  d1_list[K]
  d2_hat =  d2_list[K]
  P_hat  =  P_list[[K]]
  Q_hat  =  Q_list[[K]]
  
  return(list(d1_hat  = d1_hat,
              d2_hat  = d2_hat,
              P_hat   = P_hat,
              Q_hat   = Q_hat,
              d1_list = d1_list,
              d2_list = d2_list,
              P_list  = P_list,
              Q_list  = Q_list))
}

est.PQ = function(Y,xi,d1,d2,K = 20, thresh = FALSE, delta = NULL){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  P_list = Q_list = list()
  for (kk in 1:K){
    
    S_yxi_k = Autocov_xi_Y(Y,xi, k = kk, thresh = thresh, delta = delta)
    
    M1 = M1 + S_yxi_k%*%t(S_yxi_k)
    M2 = M2 + t(S_yxi_k)%*%S_yxi_k
    
    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)
    
    P_list[[kk]] = ev_M1$vectors[ ,1:(d1), drop = FALSE]
    Q_list[[kk]] = ev_M2$vectors[ ,1:(d2), drop = FALSE]
    
  }
  
  P_hat  =  P_list[[K]]
  Q_hat  =  Q_list[[K]]
  
  return(list(P_hat   = P_hat,
              Q_hat   = Q_hat,
              P_list  = P_list,
              Q_list  = Q_list))
}


est.d.Wf = function(Y,P,Q, Ktilde = 10, thresh = FALSE, delta = NULL,c3 = 0){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  d1 = NCOL(P);d2 = NCOL(Q);
  
  Z = array(NA,dim = c(n,d1,d2))
  for (tt in 1:n) {
    Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
  }
  
  Z_tilde = Vec.tensor(Z)
  
  M = 0
  W_list = f_list = list()
  d_list = vector()
  dmax = d1*d2
  dstar = max(d1,d2)
  
  for (kk in 1:Ktilde) {
    
    if(thresh){
      Y_2d <- Vec.tensor(Y)
      S_ytilde_k <- sigmak(t(Y_2d), as.matrix(colMeans(Y_2d)), n = n, k = kk)
      S_ytilde_k <- thresh_C(S_ytilde_k, delta)
      S_ztilde_k <- MatMult(MatMult(t(Q) %x% t(P), S_ytilde_k), Q %x% P)
    }
    else{
      S_ztilde_k = sigmak(t(Z_tilde), as.matrix(colMeans(Z_tilde)),n = n, k= kk)
    }
    
    M = M + S_ztilde_k%*%t(S_ztilde_k)
    
    ev_M = eigen(M)
    
    evalues = ev_M$values
    
    
    d_list[kk] =  max(which.max((evalues[1:(dmax-1)]+c3)/(evalues[2:(dmax)]+c3)),dstar)
    
    W_list[[kk]] = ev_M$vectors[,1:(d_list[kk]), drop = FALSE]
    
    f_list[[kk]] = Z_tilde%*%(W_list[[kk]])
  }
  
  d_list[1]  = 0
  
  d_hat  =  d_list[Ktilde]
  W_hat  =  W_list[[Ktilde]]
  f_hat  =  f_list[[Ktilde]]
  
  return(list(d_hat   = d_hat,
              W_hat   = W_hat,
              f_hat   = f_hat,
              d_list  = d_list,
              W_list  = W_list,
              f_list  = f_list))
}


est.d.Wf.nPQ = function(Z, Ktilde = 10,c3 = 0){
  n = dim(Z)[1];d1 = dim(Z)[2];d2 = dim(Z)[3];
  
  Z_tilde = Vec.tensor(Z)
  
  M = 0
  W_list = f_list = list()
  d_list = vector()
  dmax = d1*d2
  dstar = max(d1,d2)
  for (kk in 1:Ktilde){
    
    S_ztilde_k = sigmak(t(Z_tilde),as.matrix(colMeans(Z_tilde)),n = n, k= kk)
    
    M = M + S_ztilde_k%*%t(S_ztilde_k)
    
    ev_M = eigen(M)
    
    evalues = ev_M$values
    
    d_list[kk] =  max(which.max((evalues[1:(dmax-1)]+c3)/(evalues[2:(dmax)]+c3)),dstar)
    
    W_list[[kk]] = ev_M$vectors[,1:d_hat, drop = FALSE]
    
    f_list[[kk]] = Z_tilde%*%(W_list[[kk]])
  }
  
  d_list[1]  = 0
  
  
  d_hat  =  d_list[Ktilde]
  W_hat  =  W_list[[Ktilde]]
  f_hat  =  f_list[[Ktilde]]
  
  
  return(list(d_hat   = d_hat,
              W_hat   = W_hat,
              f_hat   = f_hat,
              d_list  = d_list,
              W_list  = W_list,
              f_list  = f_list))
}


est.Wf = function(Y,P,Q,d,Ktilde = 10, thresh = FALSE, delta = NULL,c3 = 0){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  d1 = NCOL(P);d2 = NCOL(Q);
  
  Z = array(NA,dim = c(n,d1,d2))
  for (tt in 1:n) {
    Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
  }
  
  Z_tilde = matrix(NA,n,d1*d2)
  for (tt in 1:n) {
    Z_tilde[tt,] = as.vector(Z[tt,,])
  }
  
  M = 0
  W_list = f_list = list()
  
  for (kk in 1:Ktilde){
    if(thresh){
      Y_2d <- Vec.tensor(Y)
      S_ytilde_k <- sigmak(t(Y_2d), as.matrix(colMeans(Y_2d)), n = n, k = kk)
      S_ytilde_k <- thresh_C(S_ytilde_k, delta)
      S_ztilde_k <- MatMult(MatMult(t(Q) %x% t(P), S_ytilde_k), Q %x% P)
    }
    else{S_ztilde_k = sigmak(t(Z_tilde),as.matrix(colMeans(Z_tilde)),n = n, k= kk)}
    
    M = M + S_ztilde_k%*%t(S_ztilde_k)
    
    ev_M = eigen(M)
    
    W_list[[kk]] = ev_M$vectors[,1:d, drop = FALSE]
    f_list[[kk]] = Z_tilde%*%(W_list[[kk]])
  }
  
  W_hat  =  W_list[[Ktilde]]
  f_hat  =  f_list[[Ktilde]]
  return(list(W_hat   = W_hat,
              f_hat   = f_hat,
              W_list  = W_list,
              f_list  = f_list))
}


est.UV.JAD = function(W,d1,d2,d){
  
  W_tilde_tol = array(NA,dim= c(d,d1,d2))
  
  for (jj in 1:d){
    W_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        W_tilde_i[mm,pp] <- W[mm + (pp - 1)*d1,jj]
      }
    }
    W_tilde_tol[jj,,] = W_tilde_i
  }
  
  P_tol = vector()
  for (ss in 1:d) {
    for(rr in ss:d){
      if(ss == rr){
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }else{
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }
      P_tol = cbind(P_tol,P_rs)
    }
  }
  
  M_tol    = svd(P_tol)$v
  dt       = NCOL(M_tol)
  M        = M_tol[,(dt-d+1):dt]
  M_tensor = array(NA,dim = c(d,d,d))
  eigen_gap = vector()
  for (ii in 1:d) {
    M_tensor[,,ii] = Vech2Mat_new(M[,ii], d)
    eigen_gap[ii] = min(abs(eigen(M_tensor[,,ii])$values))
  }
  
  
  ##construct H_star
  Ms = M_tensor[,,which.max(eigen_gap)]
  
  HH0 = vector()
  HH1 = vector()
  HH2 = vector()
  
  for (ii in 1:d) {
    HH0 = cbind(HH0,c(M_tensor[,,ii]))
    HH1 = cbind(HH1,c(solve(Ms)%*%M_tensor[,,ii]))
    HH2 = cbind(HH2,c(M_tensor[,,ii]%*%solve(Ms)))
  }
  
  PP = (t(HH0)%*%HH2)%*%solve(t(HH1)%*%HH2 + t(HH2)%*%HH1)%*%(t(HH2)%*%HH0)
  
  EVD = eigen(PP)
  if( min(EVD$values) < 0){
    R_NOJD = EVD$vectors
  }else{
    R_NOJD = sqrt(2)/2*(EVD$vectors)%*%diag((EVD$values)^{-1/2})%*%t(EVD$vectors)
  }
  
  M_1 = M%*%R_NOJD
  
  M_tensor_1 = array(NA,dim = c(d,d,d))
  for (ii in 1:d) {
    M_tensor_1[,,ii] = Vech2Mat_new(M_1[,ii],d)
    
  }
  
  H = jointDiag::ffdiag(M_tensor_1)$B ## jointDiag::ffdiag
  
  Theta = apply(MASS::ginv(H), 2, l2s)
  
  Wt = W%*%Theta
  
  U = matrix(NA,d1,d)
  V = matrix(NA,d2,d)
  
  for (jj in 1:d){
    Wt_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        Wt_tilde_i[mm,pp] <- Wt[mm + (pp - 1)*d1,jj]
      }
    }
    svdi = svd(Wt_tilde_i)
    U[,jj] = svdi$u[,1]
    V[,jj] = svdi$v[,1]
  }
  
  return(list(U = U,V = V,Theta = Theta))
  
}

est.UV.EVD = function(W,d1,d2,d){
  
  W_tilde_tol = array(NA,dim= c(d,d1,d2))
  
  for (jj in 1:d){
    W_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        W_tilde_i[mm,pp] <- W[mm + (pp - 1)*d1,jj]
      }
    }
    W_tilde_tol[jj,,] = W_tilde_i
  }
  
  P_tol = vector()
  for (ss in 1:d) {
    for(rr in ss:d){
      if(ss == rr){
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }else{
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }
      P_tol = cbind(P_tol,P_rs)
    }
  }
  Theta = 1 + 1i
  M_tol    = svd(P_tol)$v
  dt       = NCOL(M_tol)
  M_org        = M_tol[,(dt-d+1):dt]
  
  
  M  = M_org
  
  M_tensor = array(NA,dim = c(d,d,d))
  
  for (ii in 1:d) {
    M_tensor[,,ii] = Vech2Mat_new(M[,ii],d)
  }
  
  L0 = L1 = 0
  for (tt in 1:(d-1)) {
    L0 = L0 + M_tensor[,,tt]
    L1 = L1 + M_tensor[,,tt + 1]
  }
  
  L0 = L0/(d-1)
  L1 = L1/(d-1)
  
  theta.l = eigen(solve(L1)%*%L0)$vectors
  
  Theta = apply(L0%*%theta.l,2,l2s)
  
  if(is.complex(Theta)){
    Theta = Complex2Real(Theta)
  }
  
  Wt = W%*%Theta
  
  U = matrix(NA,d1,d)
  V = matrix(NA,d2,d)
  
  for (jj in 1:d){
    Wt_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        Wt_tilde_i[mm,pp] <- Wt[mm + (pp - 1)*d1,jj]
      }
    }
    svdi = svd(Wt_tilde_i)
    U[,jj] = svdi$u[,1]
    V[,jj] = svdi$v[,1]
  }
  
  return(list(U = U,V = V,Theta = Theta))
  
}

est.UV.JAD.BSel = function(W,d1,d2,d){
  
  W_tilde_tol = array(NA,dim= c(d,d1,d2))
  
  for (jj in 1:d){
    W_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        W_tilde_i[mm,pp] <- W[mm + (pp - 1)*d1,jj]
      }
    }
    W_tilde_tol[jj,,] = W_tilde_i
  }
  
  if(d1 == 1 || d2 == 1){
    P_tol = vector()
    for (ss in 1:d) {
      for(rr in ss:d){
        if(ss == rr){
          P_rs = minor_P_vector(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
        }else{
          P_rs = minor_P_vector(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
        }
        P_tol = cbind(P_tol,P_rs)
      }
    } 
  }else{
    P_tol = vector()
    for (ss in 1:d) {
      for(rr in ss:d){
        if(ss == rr){
          P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
        }else{
          P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
        }
        P_tol = cbind(P_tol,P_rs)
      }
    }
  }
  
  
  M_tol    = svd(P_tol)$v
  dt       = NCOL(M_tol)
  M        = M_tol[,(dt-d+1):dt]  
  M_tensor = array(NA,dim = c(d,d,d))
  eigen_gap = vector()
  for (ii in 1:d) {
    M_tensor[,,ii] = Vech2Mat_new(M[,ii],d)
    eigen_gap[ii] = min(abs(eigen(M_tensor[,,ii])$values))
  }
  
  
  ##construct H_star
  
  Ms = M_tensor[,,which.max(eigen_gap)] # + M_tensor[,,2]
  
  if(0 %in% eigen(Ms)$values){
    Ms = 0
    for (jj in 1:d) {
      Ms = Ms + M_tensor[,,jj]
    }
  }
  
  HH0 = vector()
  HH1 = vector()
  HH2 = vector()
  
  for (ii in 1:d) {
    HH0 = cbind(HH0,c(M_tensor[,,ii])) 
    HH1 = cbind(HH1,c(solve(Ms)%*%M_tensor[,,ii])) 
    HH2 = cbind(HH2,c(M_tensor[,,ii]%*%solve(Ms))) 
  }
  
  PP = (t(HH0)%*%HH2)%*%solve(t(HH1)%*%HH2 + t(HH2)%*%HH1)%*%(t(HH2)%*%HH0)
  
  EVD = eigen(PP)
  if( min(EVD$values) < 0){
    R_NOJD = EVD$vectors 
  }else{
    R_NOJD = sqrt(2)/2*(EVD$vectors)%*%diag((EVD$values)^{-1/2})%*%t(EVD$vectors)
  }
  
  M_1 = M%*%R_NOJD 
  
  M_tensor_1 = array(NA,dim = c(d,d,d))
  for (ii in 1:d) {
    M_tensor_1[,,ii] = Vech2Mat_new(M_1[,ii],d)
    
  }
  
  H = jointDiag::ffdiag(M_tensor_1)$B
  
  Theta = apply(MASS::ginv(H), 2, l2s)
  
  Wt = W%*%Theta 
  
  U = matrix(NA,d1,d)
  V = matrix(NA,d2,d)
  
  for (jj in 1:d){
    Wt_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        Wt_tilde_i[mm,pp] <- Wt[mm + (pp - 1)*d1,jj]
      }
    }
    svdi = svd(Wt_tilde_i)
    U[,jj] = svdi$u[,1]
    V[,jj] = svdi$v[,1]
  }
  
  return(list(U = U,V = V,Theta = Theta))
  
}

Sigma_Ycheck <- function(Y, k){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  q <- dim(Y)[3]
  sigmaYk <- matrix(0, nrow = p^2 * q, ncol = q)
  Y_mean <- apply(Y, c(2, 3), mean)
  for (t in (k + 1):n) {
    A <- Y[t, , ] - Y_mean # Y_t - Y_bar (p x q)
    B <- Y[t - k, , ] - Y_mean # Y_{t-k} - Y_bar (p x q)
    
    # vec(B): 将 B 转换为列向量
    B_vec <- as.vector(B)
    
    # Kronecker product
    kron_prod <- A %x% B_vec #  (p*q) x (p*q)
    sigmaYk <- sigmaYk + kron_prod
  }
  return(sigmaYk/(n-k))
}

adjust_sign <- function(column) {
  first_nonzero_idx <- which(column != 0)[1]
  if (!is.na(first_nonzero_idx)) {
    return(sign(column[first_nonzero_idx]))
  } else {
    return(1)
  }
}


FAC.estimation = function(Y,h0  = 1,rank = NULL){
  
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3]; 
 
  M1 = 0
  for (h in 1:h0) {
  for (ii in 1:p) {
    for (jj in 1:p) {
        Omega_ij_h =  cov(Y[1:(n-h),ii,], Y[(h+1):(n),jj,])
        M1 = M1 + Omega_ij_h%*%t(Omega_ij_h)
    }
  }
}
  
  M2 = 0
  for (h in 1:h0) {
  for (ii in 1:q) {
    for (jj in 1:q) {
      Omega_ij_h = cov(Y[1:(n-h),,ii], Y[(h+1):(n),,jj] )
      M2 = M2 + Omega_ij_h%*%t(Omega_ij_h)
    }
  }
}
  
  if(is.null(rank)){
    k1 = which.max(eigen(M1)$values[1:(p/2)]/eigen(M1)$values[2:(p/2+1)])
    k2 = which.max(eigen(M2)$values[1:(p/2)]/eigen(M2)$values[2:(p/2+1)])
  }else{
    k1 = rank[1]
    k2 = rank[2] 
  }
  
  P = eigen(M1)$vectors[,1:k1]
  Q = eigen(M2)$vectors[,1:k2]
 
  Ft = array(NA,c(n,k1,k2))
  
  for (tt in 1:n) {
    Ft[tt,,] = t(P)%*%Y[tt,,]%*%Q
  }
  
  
  return(list(P = P ,Q = Q,Ft = Ft, rank = list(k1,k2), M1 = M1 , M2 = M2))
}


forecast.factor = function(f,forecast.step = 1){
  
  if(NCOL(f) > 1){
    colnames(f) <- paste0("y",1:NCOL(f))
    var_f = vars::VAR(f, type = "const", lag.max = 6, ic = "HQ")
    pre_f = predict(var_f)
    f_tenstep = vector()
    for (jj in 1:NCOL(f)) {
      f_tenstep = cbind(f_tenstep,pre_f$fcst[[jj]][,1]) 
    }
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step")  
  }else{
    ar_f  = ar(f)
    pre_f = predict(ar_f,n.ahead = 10)
    f_tenstep = as.matrix(pre_f$pred)
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step") 
  }
  
  f_forecast = f_tenstep[forecast.step,]
  
  return(list(f.forecast = f_forecast, f_tenstep = f_tenstep))
  
}  


Chang2022.forecast = function(Y,eta,Kmax = 10,Rank = NULL, forecast.step = 1,...){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3]; 
  
  
  ##estimation of P,Q and d1,d2
  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  
  for (kk in 1:Kmax){
    
    S_yeta_k = Autocov_xi_Y(Y,eta,k = kk)
    
    M1 = M1 + S_yeta_k%*%t(S_yeta_k)
    M2 = M2 + t(S_yeta_k)%*%S_yeta_k
    
  }
  
  ev_M1 = eigen(M1)
  ev_M2 = eigen(M2)
  
  if(is.null(Rank)){
    d1 =  which.max(ev_M1$values[1:dmax]/ev_M1$values[2:(dmax+1)])
    d2 =  which.max(ev_M2$values[1:dmax]/ev_M2$values[2:(dmax+1)])
    
    d   = ifelse(p>q,d1,d2)  
  }else{
    d   = Rank
  }
  
  
  P = ev_M1$vectors[,1:d]
  Q = ev_M2$vectors[,1:d]
  
  if(d == 1){
    
    A = as.matrix(P)
    B = as.matrix(Q)
    U = V = W = as.matrix(1)
    f = vector()
    for (tt in 1:n) {
      f[tt] = t(A)%*%Y[tt,,]%*%B
    }
    f  = as.matrix(f)
    
    Z = array(NA,dim = c(n,d,d))
    for (tt in 1:n) {
      Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
    }
    
  }else{
    Z = array(NA,dim = c(n,d,d))
    for (tt in 1:n) {
      Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
    }
    
    xi =  est.xi(Z)$xi
    
    S_Zxi_1 = Autocov_xi_Y(Z,xi,k = 1)
    S_Zxi_2 = Autocov_xi_Y(Z,xi,k = 2)  
    
    vl = eigen(MASS::ginv(t(S_Zxi_1)%*%S_Zxi_1)%*%t(S_Zxi_1)%*%S_Zxi_2)$vectors
    ul = eigen(MASS::ginv(S_Zxi_1%*%t(S_Zxi_1))%*%S_Zxi_1%*%t(S_Zxi_2))$vectors
    
    U = apply(S_Zxi_1%*%vl,2,l2s)
    V = apply(t(S_Zxi_1)%*%ul,2,l2s)
    
    if(is.complex(U)){
      U = Complex2Real(U)
      V = Complex2Real(V)
    }
    
    A = as.matrix(P%*%U)
    B = as.matrix(Q%*%V)
    
  }
  
  W = matrix(NA,d^2,d)
  
  for (ii in 1:d) {
    W[,ii] = V[,ii]%x%U[,ii]
  }
  f = Vec.tensor(Z)%*%W%*%MASS::ginv(t(W)%*%W)
  
  if(NCOL(f) > 1){ #when d >1, we use VAR model 
    colnames(f) <- paste0("y",1:NCOL(f))
    var_f = vars::VAR(f, type = "const", lag.max = 6, ic = "HQ")
    pre_f = predict(var_f)
    f_tenstep = vector()
    for (jj in 1:NCOL(f)) {
      f_tenstep = cbind(f_tenstep,pre_f$fcst[[jj]][,1]) 
    }
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step")  
  }else{ #when d = 1, we use AR model 
    ar_f  = ar(f)
    pre_f = predict(ar_f,n.ahead = 10)
    f_tenstep = as.matrix(pre_f$pred)
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step") 
  }
  
  f_forecast=Z_forecast=Y_forecast=list()
  for (gg in forecast.step) {
    f_forecast[[gg]] = f_tenstep[gg,]
    
    Z_tilde_forecast = W%*%(f_forecast[[gg]])
    
    Z_forecast_i = matrix(NA,d,d)
    for (pp in 1:d){
      for (mm in 1:d) {
        Z_forecast_i[mm,pp] <- Z_tilde_forecast[mm + (pp - 1)*d,]
      }
    }
    Z_forecast[[gg]] = Z_forecast_i
    Y_forecast[[gg]] = P%*%Z_forecast_i%*%t(Q)
  }
  
  con = list(Y.forecast = Y_forecast,Z.forecast = Z_forecast, f.forecast = f_forecast, A = A,B = B, W = W, f = f, P = P, Q = Q, Rank = d)  
  
  return(con) 
  
}


Wang.forecast = function(Y,r,h0,forecast.step = c(1,2)){
   
  est_mat_fac = FAC.estimation(Y,h0,rank = r)
  
  r = est_mat_fac$rank
  
  if(max(r) == 1){
    f = est_mat_fac$Ft[,1,1]
    ar_f  = ar(f)
    pre_f = predict(ar_f,n.ahead = 10)
    f_tenstep = as.matrix(pre_f$pred)
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step")
    Ql = est_mat_fac$P
    Qr = est_mat_fac$Q
    Y.forecast = list()
    for (gg in forecast.step) {
      Y.forecast[[gg]] <- Ql%*%as.matrix(f_tenstep[gg])%*%t(Qr)
    }
    
  }else{
    
    f = est_mat_fac$Ft
    Ql = est_mat_fac$P
    Qr = est_mat_fac$Q
    k1 = dim(f)[2]
    k2 = dim(f)[3]
    
    
    fore_mar = UniARMA.forecast(f,forecast.step = 1:10)
    
    Y.forecast = list()
    for (gg in forecast.step) {
      Y.forecast[[gg]] <- Ql%*%as.matrix(fore_mar$Y.forecast[[gg]])%*%t(Qr)
    }
    
  }
  
  return(list(Y.forecast = Y.forecast, Ql = Ql, Qr = Qr,f = f))
  
  
}



HDMTS.CP.Unified.forecast = function(Y,xi,Rank = NULL,K = 20, Ktilde = 10,
                                     solve.UV = NULL, thresh1 = FALSE, thresh2 = FALSE, thresh3 = FALSE,
                                     delta1 = 0,delta2 = delta1, delta3 = delta1,
                                     c1=0,c2=0,c3=0,forecast.step = c(1,2),All_out = F,...){
  
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3]; 
  
  ##estimation of P,Q and d1,d2
  if(is.null(Rank)){ # if (d,d1,d2) is not given, we need determine d1 and d2 then estimate P and Q.
    PQ_hat_tol = est.d1d2.PQ(Y, xi, K = K, thresh = thresh1, delta = delta1,c1 = c1,c2 = c2) 
    
    d1  = PQ_hat_tol$d1_hat
    d2  = PQ_hat_tol$d2_hat
    d   = NULL
    
    if(d1 == 1 || d2 == 1){
      d = d1*d2
    }
    
    P   = PQ_hat_tol$P_hat
    Q   = PQ_hat_tol$Q_hat
    
  }else{# (d,d1,d2) is given, we estimate P and Q directly
    
    d   = Rank$d
    d1  = Rank$d1
    d2  = Rank$d2
    
    PQ_hat_tol = est.PQ(Y, xi, d1, d2, K = K, thresh = thresh1, delta = delta1)
    
    P   = PQ_hat_tol$P_hat
    Q   = PQ_hat_tol$Q_hat
    
  }
  
  ## estimation of W* = (v1*u1,v2*u2,...,vd*ud)H = WH
  if(d1 == 1 & d2 == 1){# when d1 = 1 and d2 = 1, P = A and Q = B, we can forecast based on ft = t(P)*Yt*Q
    d = 1
    
    f = vector()
    for (tt in 1:n) {
      f[tt] = t(P)%*%Y[tt,,]%*%Q
    }
    
    f = as.matrix(f)
    W = as.matrix(1)
    
    ar_f  = ar(f)
    pre_f = predict(ar_f,n.ahead = 10) # default to 10 ahead forecast
    f_tenstep = as.matrix(pre_f$pred)
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step") 
    
    f_forecast = Y_forecast = Z_forecast = list()
    for (gg in forecast.step) {
      f_forecast[[gg]] = f_tenstep[gg,]
      
      Z_forecast[[gg]] = as.matrix(f_forecast[[gg]])
      
      Y_forecast[[gg]] = P%*%(Z_forecast[[gg]])%*%t(Q)
      
    }
    
    if(All_out == F){ # tuning results are not reported when All_out = FALSE
      
      con = list(Y.forecast = Y_forecast,Z.forecast = Z_forecast, f.forecast = f_forecast, W = 1, f = f, P = P, Q = Q, Rank = list(d=d,d1=d1,d2=d2))  
      
    }else{
      
      con = list(Y.forecast = Y_forecast,Z.forecast = Z_forecast, f.forecast = f_forecast, W = 1, f = f, P = P, Q = Q, Rank = list(d=d,d1=d1,d2=d2),
                 P_list  = PQ_hat_tol$P_list ,
                 Q_list  = PQ_hat_tol$Q_list ,
                 W_list  = NULL,
                 d1_list = PQ_hat_tol$d1_list,
                 d2_list = PQ_hat_tol$d2_list,
                 d_list  = NULL,
                 f_tenstep = f_tenstep
      )  
      
      
    }
    
    return(con) 
    
  }else{# when d1 > 1 or d2 > 1, we need estimate W* and ft
    
    if(is.null(d)){
      W_hat_tol  =  est.d.Wf(Y, P, Q, Ktilde = Ktilde, thresh = thresh3, delta = delta3,c3 = c3)
      d          =  W_hat_tol$d_hat
      W          =  W_hat_tol$W_hat
      f          =  W_hat_tol$f_hat
    }else{
      d          =  d
      W_hat_tol  =  est.Wf(Y, P, Q, d, Ktilde = Ktilde, thresh = thresh3, delta = delta3)
      W          =  W_hat_tol$W_hat
      f          =  W_hat_tol$f_hat
    }
    
    if(!is.null(solve.UV)){ #solve U and V when they can be identified uniquely
      
      if(d > 1){
        if(solve.UV == "OPT"){
          UV_hat_tol = est.UV.OPT(W,d1,d2,d)
        }
        if(solve.UV == "EVD"){
          UV_hat_tol = est.UV.EVD(W,d1,d2,d)
        }
        if(solve.UV == "JAD"){
          UV_hat_tol = est.UV.JAD(W,d1,d2,d)
        }
        if(solve.UV == "JAD.BSel"){
          UV_hat_tol = est.UV.JAD.BSel(W,d1,d2,d)
        }
        
        Theta      = UV_hat_tol$Theta
        U          = UV_hat_tol$U
        V          = UV_hat_tol$V
        
        Z = array(NA,dim = c(n,d1,d2))
        for (tt in 1:n) {
          Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
        }
        
        W = matrix(NA,d1*d2,d)
        
        for (ii in 1:d) {
          W[,ii] = V[,ii]%x%U[,ii]
        }
        
        f = Vec.tensor(Z)%*%W%*%MASS::ginv(t(W)%*%W) 
        
      }
      
    }
    
    # modeling ft with VAR
    colnames(f) <- paste0("y",1:NCOL(f))
    var_f = try(vars::VAR(f, type = "const", lag.max = 6, ic = "HQ"),silent = T)
    
    if(class(var_f) == "try-error"){ #sometimes when ft are almostly perfect correlated, VAR report error, we need do AR for each ft
      f_tenstep = vector()
      colnames(f) <- paste0("y",1:NCOL(f))
      
      for (jj in 1:d) {
        ar_f  = ar(f[,jj])
        pre_f = predict(ar_f,n.ahead = 10)
        f_tenstep_i = as.matrix(pre_f$pred)
        f_tenstep = cbind(f_tenstep,f_tenstep_i)  
      }
      
      colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
      row.names(f_tenstep) <- paste(1:10,"step")   
      
    }else{ #if VAR is successful
      
      pre_f = predict(var_f)
      f_tenstep = vector()
      for (jj in 1:NCOL(f)) {
        f_tenstep = cbind(f_tenstep,pre_f$fcst[[jj]][,1]) 
      }
      colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
      row.names(f_tenstep) <- paste(1:10,"step")  
      
      if(c(is.na(f_tenstep))[1]){ #VAR may generate NA, if NA exist, we do AR instead
        f_tenstep = vector()
        colnames(f) <- paste0("y",1:NCOL(f))
        
        for (jj in 1:d) {
          ar_f  = ar(f[,jj])
          pre_f = predict(ar_f,n.ahead = 10)
          f_tenstep_i = as.matrix(pre_f$pred)
          f_tenstep = cbind(f_tenstep,f_tenstep_i)  
        }
        
        colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
        row.names(f_tenstep) <- paste(1:10,"step")   
      }    
    }
    
    f_forecast = Y_forecast = Z_forecast = list()  
    
    for (gg in forecast.step) {
      f_forecast[[gg]] = f_tenstep[gg,]
      
      Z_tilde_forecast = W%*%(f_forecast[[gg]])
      
      Z_forecast_i = matrix(NA,d1,d2)
      for (pp in 1:d2){
        for (mm in 1:d1) {
          Z_forecast_i[mm,pp] <- Z_tilde_forecast[mm + (pp - 1)*d1,]
        }
      }
      Z_forecast[[gg]] = Z_forecast_i
      Y_forecast[[gg]] = P%*%Z_forecast_i%*%t(Q)  
    }
    
    if(All_out == F){
      
      con = list(Y.forecast = Y_forecast,Z.forecast = Z_forecast, f.forecast = f_forecast, W = W, f = f, P = P, Q = Q, Rank = list(d=d,d1=d1,d2=d2))  
      
    }else{
      
      con = list(Y.forecast = Y_forecast,Z.forecast = Z_forecast, f.forecast = f_forecast, W = W, f = f, P = P, Q = Q, Rank = list(d=d,d1=d1,d2=d2),
                 P_list     = PQ_hat_tol$P_list ,
                 Q_list     = PQ_hat_tol$Q_list ,
                 W_list     = W_hat_tol$W_list,
                 d1_list    = PQ_hat_tol$d1_list,
                 d2_list    = PQ_hat_tol$d2_list,
                 d_list     = W_hat_tol$d_list,
                 f_tenstep  = f_tenstep
      )  
      
      
    }
    
    return(con) 
    
  } 
  
}

 
PCATS.forecast = function(Y,permatation = "max",forecast.step = c(1,2,6),...){
  p = dim(Y)[2];q = dim(Y)[3];
  Y_m =  Vec.tensor(Y)
  
  pcats = HDTSA::PCA_TS(Y_m,permutation = permatation,beta = 0.1)
  
  Z     = pcats$X
  B     = pcats$B
  group = pcats$Groups
  
  Z_fore = matrix(0,max(forecast.step),NCOL(Z))
  for (ii in 1:NCOL(group)) {
    
    Z_i = Z[,group[,ii]]
    
    Z_fore_i = forecast.factor(Z_i,forecast.step = forecast.step)$f.forecast
    
    
    Z_fore[forecast.step,group[,ii]] = Z_fore_i
  }
  
  Y_fore = solve(B)%*%t(Z_fore)
  
  Y.forecast = list()
  
  for (jj in forecast.step) {
    Y.forecast[[jj]] <- matrix(Y_fore[,jj],p,q)
  }
  
  return(list(Y.forecast = Y.forecast, Z = Z, B = B))
  
}

UniARMA.forecast = function(Y,forecast.step = c(1,2)){
  
  p = dim(Y)[2]
  q = dim(Y)[3]
  
  Y_vect = Vec.tensor(Y)
  
  Y_tenstep = apply(Y_vect, 2, function(x){
    res = ar(x)
    arma_fore = as.numeric(predict(res,n.ahead = 10)$pred)
    return(arma_fore)
  })
  
  
  Y.forecast = list()
  for (gg in forecast.step) {
    Y.forecast[[gg]] <- matrix(Y_tenstep[gg,],p,q)
  }
  
  return(list(Y.forecast = Y.forecast))
}


MAR.forecast = function(Y,method = "RRLSE",k1, k2,tol = 10^-4,forecast.step = c(1,2)){
  
  est_mat_ar1  = tensorTS::matAR.RR.est(Y,method = method, k1 = k1, k2 = k2, tol = tol)
  
  fore_mat_ar1 = tensorTS:::predict.tenAR(est_mat_ar1,n.ahead = 10)
  
  Y.forecast = list()
  for (gg in forecast.step) {
    Y.forecast[[gg]] <- fore_mat_ar1[gg,,]
  }
  
  return(list(Y.forecast = Y.forecast))
  
}

FAC.forecast = function(Y,r = NULL,h0,method,forecast.step = c(1,2)){
  require(tensorTS)
  
  if(is.null(r)){
    r =  tensorTS::tenFM.rank(Y, rank = "IC",method = "TIPUP",penalty = 3)$factor.num
  }
  
  est_mat_fac = tensorTS::tenFM.est(Y,r = r,h0 = h0,method = method,tol = 10^-6)
  
  if(max(r) == 1){
    f = est_mat_fac$Ft[,1,1]
    ar_f  = ar(f)
    pre_f = predict(ar_f,n.ahead = 10)
    f_tenstep = as.matrix(pre_f$pred)
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step")
    Ql = est_mat_fac$Q[[1]]
    Qr = est_mat_fac$Q[[2]]
    Y.forecast = list()
    for (gg in forecast.step) {
      Y.forecast[[gg]] <- Ql%*%as.matrix(f_tenstep[gg])%*%t(Qr)
    }
    
  }else{
    
    f = est_mat_fac$Ft
    Ql = est_mat_fac$Q[[1]]
    Qr = est_mat_fac$Q[[2]]
    k1 = dim(f)[2]
    k2 = dim(f)[3]
    
    
    fore_mar = UniARMA.forecast(f,forecast.step = 1:10)
    
    Y.forecast = list()
    for (gg in forecast.step) {
      Y.forecast[[gg]] <- Ql%*%as.matrix(fore_mar$Y.forecast[[gg]])%*%t(Qr)
    }
    
  }
  
  return(list(Y.forecast = Y.forecast, Ql = Ql, Qr = Qr,f = f))
  
  
}


