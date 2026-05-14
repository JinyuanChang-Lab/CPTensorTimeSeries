l2s = function(x){x/sqrt(sum(x^2))}
fnorm = function(x){sqrt(sum(x^2))}

cp_residuals_general <- function(Y, f, A_list) {
  # Y: array of dimension n x d1 x d2 x ... x dm
  # f: matrix of dimension n x r
  # A_list: list(A1, ..., Am), where Aj is dj x r
  
  if (!is.array(Y)) {
    stop("Y must be an array with dimensions n x d1 x ... x dm.")
  }
  if (!is.matrix(f)) {
    stop("f must be an n x r matrix.")
  }
  if (!is.list(A_list) || length(A_list) < 1) {
    stop("A_list must be a non-empty list: list(A1, ..., Am).")
  }
  
  dims <- dim(Y)
  if (length(dims) < 2) {
    stop("Y must have at least 2 dimensions: n x d1.")
  }
  
  n <- dims[1]
  mode_dims <- dims[-1]
  m <- length(mode_dims)
  
  if (length(A_list) != m) {
    stop("length(A_list) must equal length(dim(Y)) - 1.")
  }
  
  r <- ncol(f)
  if (nrow(f) != n) {
    stop("nrow(f) must equal dim(Y)[1].")
  }
  
  # Check each loading matrix
  for (j in seq_len(m)) {
    Aj <- A_list[[j]]
    if (!is.matrix(Aj)) {
      stop(sprintf("A_list[[%d]] must be a matrix.", j))
    }
    if (nrow(Aj) != mode_dims[j]) {
      stop(sprintf("nrow(A_list[[%d]]) must equal dim(Y)[%d].", j, j + 1))
    }
    if (ncol(Aj) != r) {
      stop(sprintf("ncol(A_list[[%d]]) must equal ncol(f).", j))
    }
  }
  
  # Precompute rank-1 basis tensors (vectorized)
  total_dim <- prod(mode_dims)
  basis_mat <- matrix(0, nrow = total_dim, ncol = r)
  
  for (i in seq_len(r)) {
    # Start from the first mode loading vector
    comp_vec <- A_list[[1]][, i]
    
    # Sequentially build the tensor product
    if (m >= 2) {
      for (j in 2:m) {
        comp_vec <- as.vector(outer(comp_vec, A_list[[j]][, i]))
      }
    }
    
    basis_mat[, i] <- comp_vec
  }
  
  # Fitted values in matricized form: n x (d1*...*dm)
  Y_hat_mat <- f %*% t(basis_mat)
  
  # Convert back to array
  Y_hat <- array(Y_hat_mat, dim = dims)
  E <- Y - Y_hat
  
  list(
    residual = E,
    fitted = Y_hat,
    basis = basis_mat
  )
}
svd_inverse <- function(mat, threshold = 1e-6) {
  # Input validation
  if (!is.matrix(mat)) {
    stop("Input must be a matrix!")
  }
  
  if (!is.numeric(mat)) {
    stop("Matrix must be a numeric matrix!")
  }
  
  if (nrow(mat) != ncol(mat)) {
    stop("Only square matrices are supported for inversion!")
  }
  
  if (!is.numeric(threshold) || threshold <= 0) {
    stop("Threshold must be a positive number!")
  }
  
  # Perform singular value decomposition (SVD)
  svd_result <- svd(mat)
  
  # Extract singular value vector
  d <- svd_result$d
  
  # Replace singular values smaller than threshold with the threshold
  d_corrected <- pmax(d, threshold)
  
  # Construct the inverse of the diagonal matrix with corrected singular values
  d_inv <- diag(1 / d_corrected, nrow = length(d_corrected))
  
  # Calculate the corrected inverse matrix (V * D^{-1} * U^T)
  inv_mat <- svd_result$v %*% d_inv %*% t(svd_result$u)
  
  return(inv_mat)
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

# Function: ensure the first element of each column is positive
make_first_row_positive <- function(mat) {
  stopifnot(is.matrix(mat))  # ensure input is a matrix
  
  for (j in seq_len(ncol(mat))) {
    if (mat[1, j] < 0) {
      mat[, j] <- -mat[, j]
    }
  }
  return(mat)
}

generate_tensor_ar1_indep <- function(n, dims, c,
                                      error_dist = c("normal", "t"),
                                      df = NULL,
                                      burnin = 200,
                                      seed = NULL,
                                      standardize_t = TRUE,
                                      return_phi = TRUE) {
  # Basic checks
  if (!is.null(seed)) set.seed(seed)
  
  error_dist <- match.arg(error_dist)
  
  if (length(n) != 1 || !is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be one positive integer.")
  }
  
  if (length(dims) < 1 || any(!is.numeric(dims)) || any(dims <= 0) ||
      any(dims != as.integer(dims))) {
    stop("dims must be a vector of positive integers.")
  }
  
  if (length(c) != 1 || !is.numeric(c) || c < 0 || c >= 1) {
    stop("c must be a number in [0, 1).")
  }
  
  if (error_dist == "t") {
    if (is.null(df) || !is.numeric(df) || length(df) != 1 || df <= 0) {
      stop("For t errors, df must be one positive number.")
    }
    if (standardize_t && df <= 2) {
      stop("If standardize_t = TRUE, df must be greater than 2.")
    }
  }
  
  # Number of spatial locations
  p <- prod(dims)
  total_n <- n + burnin
  
  # Draw one AR coefficient for each location
  phi <- runif(p, min = -c, max = c)
  
  # Generate innovations matrix: total_n x p
  if (error_dist == "normal") {
    eps <- matrix(rnorm(total_n * p), nrow = total_n, ncol = p)
  } else {
    eps <- matrix(rt(total_n * p, df = df), nrow = total_n, ncol = p)
    if (standardize_t) {
      # Scale t innovations to have variance 1
      eps <- eps / sqrt(df / (df - 2))
    }
  }
  
  # Allocate matrix for all AR(1) paths
  # Each column is one location-specific AR(1) process
  Y_mat <- matrix(0, nrow = total_n, ncol = p)
  
  # Initialize
  Y_mat[1, ] <- eps[1, ]
  
  # Time recursion, vectorized over all locations
  for (tt in 2:total_n) {
    Y_mat[tt, ] <- phi * Y_mat[tt - 1, ] + eps[tt, ]
  }
  
  # Remove burn-in
  Y_mat <- Y_mat[(burnin + 1):total_n, , drop = FALSE]
  
  # Reshape to tensor: c(n, d1, ..., dm)
  Y_tensor <- array(Y_mat, dim = c(n, dims))
  
  if (return_phi) {
    phi_tensor <- array(phi, dim = dims)
    return(list(
      Y = Y_tensor,      # Tensor time series
      phi = phi_tensor   # AR coefficients at each location
    ))
  } else {
    return(Y_tensor)
  }
}


Mat.k = function (A, k, eps = 10^-6)
{
  ev = eigen(A)$values
  mark = which(ev > eps)
  ev = ev[mark]
  evc = as.matrix(eigen(A)$vectors)[, mark]
  Matk = evc %*% diag(ev^k) %*% t(evc)
  return(Matk)
}



sigma_k = function(Y,Y_mean,k,n){
  
  Y = t( scale(t(Y),scale = F) )
  
  return(Y[,(k+1):(n)]%*%t(Y[,(1):(n-k)])/(n-k))
  
}


rho2.loss = function(A_hat,A){
  
  max(apply(1-(t(A_hat)%*%A)^2,2,min))
  
}

rho2.f.loss = function(f_hat,f){
  
  max( apply(1- cor(f_hat,f)^2 ,2,min))
  
}


rho2.loss.list = function(A_hat,A){
  
  m = length(A_hat)
  
  rho2  = c()
  for (j in 1:m) {
    rho2[j] = max(apply(1-(t(A_hat[[j]])%*%A[[j]])^2,2,min))
  }
  return(rho2)
}


vecpsi.loss = function(A_hat,A){
  apply(1-(t(A_hat)%*%A)^2,2,min)
}

vecpsi.loss.list = function(A_hat,A){
  m = length(A_hat)
  rho2  = vector()
  for (j in 1:m) {
    rho2  = rbind(rho2, vecpsi.loss(A_hat[[j]],A[[j]]))
  }
  apply(rho2, 2, max)
}

psi2.loss.list = function(A_hat,A){
  
  m = length(A_hat)
  
  rho2  = c()
  for (j in 1:m) {
    rho2[j] = min(apply(1-(t(A_hat[[j]])%*%A[[j]])^2,2,min))
  }
  return(rho2)
}


D.loss = function(Q_hat,Q){
  
  r_breve = NCOL(as.matrix(Q_hat))
  
  
  rho2 = sqrt(1 - svd(t(Q_hat)%*%Q)$d[1])
  
  return(rho2)
}

D.loss.list = function(Q_hat,Q){
  
  m = length(Q_hat)
  r_breve = NCOL(as.matrix(Q_hat[[1]]))
  rho2  = c()
  for (j in 1:m) {
    rho2[j] = sqrt(1 - svd(t(Q_hat[[j]])%*%Q[[j]])$d[1])
  }
  return(rho2)
}




DGP.TCP = function(n,m,D,r,w,ar.coef,factor.loading = c("sparse-random","sparse-exact","sparse-exact-corr1","sparse-exact-corr2"), factor.corr = 0,alpha = 0.5,delta = 0.75,par_E = 1,heavytail = FALSE,error.ar = FALSE){
  
  if(alpha > 0){
    
    if(factor.loading == "sparse-random"){
      A = list()
      omega = rep(1,r)
      for (j in 1:m){
        tau = 0
        while (0 %in%  tau || 1 %in%  tau ) {
          Aj = matrix(runif(r*D[j],-1,1),D[j],r)
          Aj[which(abs(Aj) < alpha*1)] <- 0
          tau = apply(Aj,2, function(x){length(which(x != 0))})
        }
        A[[j]] = apply(Aj,2,l2s)
        
      }
    }
    
    # if(factor.loading == "sparse-random-corr1"){
    #   A = list()
    #   omega = rep(1,r)
    #   for (j in 1:m){
    #     tau = 0
    #     while (0 %in%  tau || 1 %in%  tau ) {
    #       Aj = matrix(runif(r*D[j],-1,1),D[j],r)
    #
    #       Aj[which(abs(Aj) < alpha*1)] <- 0
    #
    #       tau = apply(Aj,2, function(x){length(which(x != 0))})
    #     }
    #
    #     for (i in 2:r) {
    #       Aj[,i] = delta*Aj[,i-1] + Aj[,i]
    #     }
    #
    #     A[[j]] = apply(Aj,2,l2s)
    #
    #   }
    # }
    
    if(factor.loading == "sparse-random-corr1"){
      A = list()
      omega = rep(1,r)
      for (j in 1:m){
        
        Aj = matrix(runif(r*D[j],-1,1),D[j],r)
        
        for (i in 2:r) {
          Aj[,i] = delta*Aj[,i-1] + Aj[,i]
          
        }
        
        for (i in 1:r) {
          Aj[sample(D[j],alpha*D[j]),i] <- 0
          
        }
        
        A[[j]] = apply(Aj,2,l2s)
        
      }
    }
    
    if(factor.loading == "sparse-random-corr2"){
      A = list()
      omega = rep(1,r)
      for (j in 1:m){
        Aj = matrix(runif(r*D[j],-1,1),D[j],r)
        
        Aj[sample(D[j],alpha*D[j]), ] <- 0
        
        for (i in 2:r) {
          Aj[,i] = delta*Aj[,i-1] + Aj[,i]
        }
        
        A[[j]] = apply(Aj,2,l2s)
        
      }
    }
    
    if(factor.loading == "sparse-exact"){
      A = list()
      for (j in 1:m){
        Aj = matrix(1,D[j],r)
        for (i in 1:r) {
          Aj[sample(1:D[j],alpha*D[j]), i] <- 0
        }
        A[[j]] = apply(Aj,2,l2s)
      }
    }
    
    if(factor.loading == "sparse-exact-corr1"){
      A = list()
      
      for (j in 1:m){
        Aj = matrix(1,D[j],r)
        for (i in 1:r) {
          Aj[sample(1:D[j],alpha*D[j]), i] <- 0
          Aj = apply(Aj,2,l2s)
        }
        
        for (i in 2:r) {
          Aj[,i] = delta*Aj[,i-1] + Aj[,i]
        }
        
        A[[j]] = apply(Aj,2,l2s)
      }
    }
    
    if(factor.loading == "sparse-exact-corr2"){
      A = list()
      
      for (j in 1:m){
        Aj = matrix(0,D[j],r)
        
        Aj[sample(1:D[j],(1-alpha)*D[j]), 1] <- 1
        
        for (i in 2:r) {
          Aj[which(Aj[,i-1] != 0),i] = delta*Aj[which(Aj[,i-1] != 0),i-1] + (1-delta)*runif(length(which(Aj[,i-1] != 0)),-1,1)
        }
        
        A[[j]] = apply(Aj,2,l2s)
      }
    }
    
  }else{
    
    A = list()
    for (j in 1:m){
      Aj = matrix(runif(r*D[j],-3,3),D[j],r)
      
      for (i in 2:r) {
        Aj[,i] = delta*Aj[,i-1] + Aj[,i]
      }
      
      A[[j]] = apply(Aj,2,l2s)
    }
    
  }
  
  
  
  
  
  f_m = matrix(NA,n,r)
  S = 0
  for (ii in 1:r) {
    xd = arima.sim(model = list(ar = ar.coef[[ii]]),n = n)
    
    f_m[,ii] = xd
  }
  
  COV = diag(1-factor.corr,r,r) + matrix(factor.corr,r,r)
  
  COV.half = eigen(COV)$vectors %*%  diag(sqrt( eigen(COV)$values)) %*% t(eigen(COV)$vectors)
  
  f_m = f_m%*%COV.half
  
  for (ii in 1:r) {
    S_r = f_m[,ii]
    for (jj in 1:m) {
      S_r = (S_r %o% A[[jj]][,ii])
    }
    S = S + w[ii]*S_r
  }
  
  if(error.ar == FALSE){
    if(heavytail == FALSE){
      E = array(rnorm(prod(D)*n,0,par_E),c(n,D))
    }else{
      E = array(rt(prod(D)*n,heavytail),c(n,D))
    }
  }else{
    if(heavytail == FALSE){
      resE = generate_tensor_ar1_indep(
        n = n,
        dims = D,
        c = error.ar,
        error_dist = "normal"
      )
    }else{
      resE = generate_tensor_ar1_indep(
        n = n,
        dims = D,
        c = error.ar,
        error_dist = "t",
        df = 5
      )
    }
    E = resE$Y
  }
  
  
  
  generate_tensor_ar1_indep(
    n = 500,
    dims = c(4, 5),
    c = 0.8,
    error_dist = "normal",
    seed = 123
  )
  
  
  
  Y = S + E
  
  return(list(Y = Y,C = S,E = E,A = A,f = f_m, w = w))
  
}


Mat.tensor = function(Y,j){ # fold on j-th mode
  t(apply(Y,j,c))
}

Threshold.Tensor = function(SigmaY,n,sigma0,delta){
  d = prod(dim(SigmaY))
  SigmaY[which(abs(SigmaY) < delta*sigma0*sqrt(log(d)/n))] <- 0
  return(SigmaY)
}

tensor.est.xi = function(Y,d_max = 10,thresh_per = 0.99, random = F, seed = NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  n = dim(Y)[1];D = dim(Y)[-1]
  
  Y.mat = Mat.tensor(Y,1)
  
  if(n > prod(D)){
    eig_Y.mat = eigen(HDTSA:::MatMult(t(Y.mat),Y.mat))
    cfr =  cumsum(eig_Y.mat$values)/sum(eig_Y.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    
    w_inl = as.matrix(eig_Y.mat$vectors[,1:d_fin])
    
    #w_inl = rep(1,NROW(w_inl))%*%t(sign(w_inl[1,]))*w_inl #ensure 1st element of each column is positive
    
    if(random == T){
      w_hat = (w_inl)%*%pracma:::randortho(d_fin)
    }else{
      w_hat = (w_inl)
    }
    
    w_hat = make_first_row_positive(w_hat)
    
    xi.f = scale(Y.mat%*%w_hat)
    xi   = rowMeans(xi.f)
    
  }else{
    eig_Y.mat = eigen(HDTSA:::MatMult(Y.mat,t(Y.mat)))
    cfr = cumsum(eig_Y.mat$values)/sum(eig_Y.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    
    w_inl = as.matrix(eig_Y.mat$vectors[,1:d_fin])
    
    #w_inl = rep(1,NROW(w_inl))%*%t(sign(w_inl[1,]))*w_inl #ensure 1st element of each column is positive
    
    if(random == T){
      xi.f = (w_inl)%*%pracma:::randortho(d_fin)
    }else{
      xi.f = as.matrix(w_inl)
    }
    
    weight = sqrt(eig_Y.mat$values[1:d_fin])
    
    xi.f = xi.f   #%*%diag(weight)
    
    xi.f = make_first_row_positive(xi.f)
    
    xi = rowMeans(xi.f)
    
  }
  return(scale(xi))
}



Autocov_xi_Y = function(Y,eta,lag.k = k){
  
  if(length(dim(Y)) == 3){
    n = dim(Y)[1]
    k = lag.k
    
    Y_mean = 0
    eta_mean = 0
    for (ii in 1:n) {
      Y_mean = Y_mean + Y[ii,,]
      eta_mean = eta_mean + eta[ii]
    }
    Y_mean   = Y_mean/n
    eta_mean = eta_mean/n
    
    Sigma_Y_eta_k = 0
    for (ii in (k+1):n) {
      Sigma_Y_eta_k = Sigma_Y_eta_k + (Y[ii,,] - Y_mean)*(eta[ii-k] - eta_mean)
    }
    
  }
  
  if(length(dim(Y)) == 2){
    n = dim(Y)[1]
    k = lag.k
    
    Y_mean = 0
    eta_mean = 0
    for (ii in 1:n) {
      Y_mean = Y_mean + Y[ii,]
      eta_mean = eta_mean + eta[ii]
    }
    Y_mean   = Y_mean/n
    eta_mean = eta_mean/n
    
    Sigma_Y_eta_k = 0
    for (ii in (k+1):n) {
      Sigma_Y_eta_k = Sigma_Y_eta_k + (Y[ii,] - Y_mean)*(eta[ii-k] - eta_mean)
    }
    
  }
  
  return(Sigma_Y_eta_k/(n-k))
}



tensor.Autocov_xi_Y = function(Y,xi,lag.k = k){ # reuturn a autocovariance  tensor (d1 x d2 x ... x d_m)
  
  n = dim(Y)[1]
  k = lag.k
  
  Y_mean = apply(Y, c(2:length(dim(Y))), mean)
  xi_mean = mean(xi)
  
  Sigma_Y_xi_k = 0
  for (ii in (k+1):n) {
    Sigma_Y_xi_k = Sigma_Y_xi_k + (R.utils::extract(Y,"1" = ii,drop = T) - Y_mean)*(xi[ii-k] - xi_mean)
  }
  return(Sigma_Y_xi_k/(n-k))
}



HDTTS.CP.est =  function(Y, xi = NULL, Rank = NULL,K = 10, Threshold = FALSE, delta = NULL, Ratio.type = c("classical","log"), augmented = FALSE,grid_delta1 = 50,delta_max = 0.1){
  est.ABr =   function(Sigma.tensor.Y.k_list_orginal,Rank,delta,sigma0,...){
    
    Sigma.tensor.Y.k_list = list()
    
    for (kk in 1:K) {
      Sigma.tensor.Y.k = Sigma.tensor.Y.k_list_orginal[[kk]]
      
      if(Threshold == TRUE & kk == 1){
        Sigma.tensor.Y.k = Threshold.Tensor(Sigma.tensor.Y.k,n = n,sigma0 = sigma0,delta = delta)
        index = Sigma.tensor.Y.k
      }
      if(Threshold == TRUE & kk > 1){
        Sigma.tensor.Y.k[which(index == 0)] <- 0
      }
      
      Sigma.tensor.Y.k_list[[kk]] = Sigma.tensor.Y.k
    }
    
    r_tol  = vector()
    Mjk_list = MMjk_list = list()
    Sigma.Y.k_list_tol = list()
    for (j in 1:m){
      Mjk = 0
      MMjk = 0
      Sigma.Y.k_list = list()
      for (kk in 1:K){
        Sigma.Y.k = Mat.tensor(Sigma.tensor.Y.k_list[[kk]],j)
        Mjk = Mjk + Sigma.Y.k%*%t(Sigma.Y.k)
        MMjk = MMjk + t(Sigma.Y.k)%*%Sigma.Y.k
        Sigma.Y.k_list[[kk]] = Sigma.Y.k
      }
      
      Mjk_list[[j]] = Mjk
      MMjk_list[[j]] = MMjk
      Sigma.Y.k_list_tol[[j]] = Sigma.Y.k_list
      
      if(Ratio.type == "log"){
        eigenvalue_j =  log(eigen(Mjk)$values + 1)
      }
      if(Ratio.type == "classical"){
        eigenvalue_j =  eigen(Mjk)$values
      }
      
      ratio_j = (eigenvalue_j[-1] + sigma0/n)/(eigenvalue_j[-length(eigenvalue_j)] + sigma0/n)
      
      r_j = which.min(ratio_j[1:(0.5*length(ratio_j))])
      
      r_tol = c(r_tol,r_j)
    }
    
    length_enc <- rle(r_tol)
    
    # r = length_enc$values[which.max(length_enc$lengths)]
    # r.hat = length_enc$values[which.max(length_enc$lengths)]
    
    r =  max(r_tol)
    r.hat = max(r_tol)
    
    if(!is.null(Rank)){
      r = Rank
    }
    
    A_hat = A12_hat = A1K_hat = list()
    j = 1
    
    P_list = list()
    Q_list = list()
    K_12_tilde_list = list()
    
    for (j in 1:m) {
      Pj = as.matrix(eigen(Mjk_list[[j]])$vectors[,1:r])
      Qj = as.matrix(eigen(MMjk_list[[j]])$vectors[,1:r])
      
      P_list[[j]] = Pj
      Q_list[[j]] = Qj
      
      Sigma.Y.k_list = Sigma.Y.k_list_tol[[j]]
      if(augmented == T){
        
        ev.K12 = names = vector()
        ss = 1
        Aj_hat_list = list()
        for (kk in 1) {
          for (jj in (kk+1):K) {
            bb1 = (Sigma.Y.k_list[[kk]])%*%Qj
            bb2 = (Sigma.Y.k_list[[jj]])%*%Qj
            
            K12_j_tilde = bb1%*%MASS::ginv(t(bb2)%*%bb2)%*%t(bb2)
            
            evd.K12 = eigen(K12_j_tilde)
            
            if(sum(abs(Im(evd.K12$values)[1:r])) < 10^-10){
              ev.jk  =  Re(evd.K12$values)[1:r]
              ev.K12 =  rbind(ev.K12,ev.jk)
              
              Aj_hat_list[[ss]] = Re(evd.K12$vectors[,1:r])
              
              names = c(names,paste0(kk,jj))
              
              ss = ss+1
            }
            
          }
        }
        
        row.names(ev.K12) <- names
        
        ev.minl1 = abs(ev.K12)[apply(abs(ev.K12), 1, min) > 1,]
        
        if(length(ev.minl1) == 0){
          ev.K12 = names = vector()
          ss = 1
          Aj_hat_list = list()
          for (kk in 1) {
            for (jj in (kk+1):K) {
              
              bb1 = (Sigma.Y.k_list[[kk]])%*%Qj
              bb2 = (Sigma.Y.k_list[[jj]])%*%Qj
              
              K12_j_tilde = bb2%*%MASS::ginv(t(bb1)%*%bb1)%*%t(bb1)
              
              evd.K12 = eigen(K12_j_tilde)
              
              if(sum(abs(Im(evd.K12$values)[1:r])) < 10^-10){
                ev.jk  =  Re(evd.K12$values)[1:r]
                ev.K12 =  rbind(ev.K12,ev.jk)
                
                Aj_hat_list[[ss]] = Re(evd.K12$vectors[,1:r])
                
                names = c(names,paste0(kk,jj))
                
                ss = ss+1
              }
              
              
            }
          }
          
          row.names(ev.K12) <- names
          
          ev.minl1 =  abs(ev.K12)[apply(abs(ev.K12), 1, min) > 1,]
          
        }
        
        if(length(ev.minl1) == 0){
          Aj_hat = Aj_hat_list[[which(names %in% names(which.max(apply(abs(ev.K12),1,prod))))]]
          
          Aj12_hat = Aj_hat
          
          Aj1K_hat = Aj_hat
        }
        
        if(length(ev.minl1) == r){
          Aj_hat = Aj_hat_list[[which( abs(ev.K12)[,1] %in% ev.minl1[1])]]
          
          Aj12_hat = Aj_hat
          
          Aj1K_hat = Aj_hat
        }
        
        if(length(ev.minl1) > r){
          Aj_hat = Aj_hat_list[[which(names %in% names(which.max(apply(ev.minl1,1,prod))))]]
          
          Aj12_hat = Aj_hat
          
          Aj1K_hat = Aj_hat
          
        }
        
      }
      
      if(augmented == F){
        bb1 = Sigma.Y.k_list[[1]]%*%Qj
        bb2 = Sigma.Y.k_list[[2]]%*%Qj
        bbk = Sigma.Y.k_list[[K]]%*%Qj
        
        K21_j_tilde = bb2%*%MASS::ginv(t(bb1)%*%bb1)%*%t(bb1)
        
        K12_j_tilde = bb1%*%MASS::ginv(t(bb2)%*%bb2)%*%t(bb2)
        
        K1K_j_tilde = bb1%*%MASS::ginv(t(bbk)%*%bbk)%*%t(bbk)
        
        evd.K12 = eigen(K12_j_tilde)
        
        if(sum(abs(Im(evd.K12$values)[1:r])) < 10^-10){
          Aj_hat = Re(eigen(K21_j_tilde)$vectors[,1:r])
          
          Aj12_hat = Re(eigen(K12_j_tilde)$vectors[,1:r])
          
          Aj1K_hat = Re(eigen(K1K_j_tilde)$vectors[,1:r])
        }else{
          Aj_hat = Aj1K_hat = Aj12_hat =  apply(Complex2Real(eigen(K12_j_tilde)$vectors[,1:r]) ,2, l2s)
        }
        
      }
      
      if(augmented == "summation"){
        K12_j_tilde_tol = 0
        for (jj in 1:1){
          bb1 = Sigma.Y.k_list[[1]]%*%Qj
          bb2 = Sigma.Y.k_list[[2]]%*%Qj
          
          K12_j_tilde = bb1%*%MASS::ginv(t(bb2)%*%bb2)%*%t(bb2)
          K21_j_tilde = bb2%*%MASS::ginv(t(bb1)%*%bb1)%*%t(bb1)
          
          K12_j_tilde_tol =  K12_j_tilde_tol + K12_j_tilde%*%K12_j_tilde + K21_j_tilde%*%K21_j_tilde
        }
        
        Aj12_hat = Re(eigen(K12_j_tilde_tol)$vectors[,1:r])
        
        Aj_hat =  Aj12_hat
        
        Aj1K_hat =  Aj12_hat
      }
      
      K_12_tilde_list[[j]] = K12_j_tilde
      
      A_hat[[j]] = as.matrix(Aj_hat)
      
      A12_hat[[j]] = as.matrix(Aj12_hat)
      
      A1K_hat[[j]] = as.matrix(Aj1K_hat)
    }
    
    # rho2.loss.list(A_hat,data$A)
    # rho2.loss.list(A12_hat,data$A)
    # rho2.loss.list(A1K_hat,data$A)
    
    delta_sel = delta
    
    return(list(A.hat = A12_hat,A21.hat = A_hat, A1K.hat = A1K_hat, r.hat = r.hat, delta_sel = delta,Sigma.tensor.Y.k_list_threshold = Sigma.tensor.Y.k_list,P_list = P_list,Q_list = Q_list, K_12_tilde_list = K_12_tilde_list))
  }
  
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  if(is.null(xi)){
    xi = tensor.est.xi(Y,random = F)
  }
  
  #sigma0 = mean((Y - rep(1,n) %o% apply(Y, c(2:length(dim(Y))), mean))^2)
  sigma0 = sqrt(sum(Y^2)/(n*prod(D)))
  
  Sigma.tensor.Y.k_list_orginal = list()
  
  for (kk in 1:K) {
    Sigma.tensor.Y.k_list_orginal[[kk]]  = tensor.Autocov_xi_Y(Y,scale(xi),kk)
  }
  
  if(is.null(delta) & Threshold == TRUE){
    
    test_list = z_list = r_list  = vector()
    
    net_delta = seq(0,delta_max,length.out = grid_delta1)
    
    A_hat_list = list()
    
    for (delta in net_delta) {
      
      Sigma.tensor.Y.k_list = list()
      
      for (kk in 1:K) {
        Sigma.tensor.Y.k = Sigma.tensor.Y.k_list_orginal[[kk]]
        
        if(Threshold == TRUE & kk == 1){
          Sigma.tensor.Y.k = Threshold.Tensor(Sigma.tensor.Y.k,n = n,sigma0 = sigma0,delta = delta)
          index = Sigma.tensor.Y.k
        }
        if(Threshold == TRUE & kk > 1){
          Sigma.tensor.Y.k[which(index == 0)] <- 0
        }
        Sigma.tensor.Y.k_list[[kk]] = Sigma.tensor.Y.k
      }
      
      test_tol = z_tol = r_tol  = vector()
      
      A_hat = list()
      
      for (j in 1:m) {
        Mjk = 0
        MMjk = 0
        Sigma.Y.k_list = list()
        for (kk in 1:K){
          Sigma.Y.k = Mat.tensor(Sigma.tensor.Y.k_list[[kk]],j)
          Mjk = Mjk + Sigma.Y.k%*%t(Sigma.Y.k)
          MMjk = MMjk + t(Sigma.Y.k)%*%Sigma.Y.k
          Sigma.Y.k_list[[kk]] = Sigma.Y.k
        }
        
        if(Ratio.type == "log"){
          eigenvalue_j =  log(eigen(Mjk)$values + 1)
        }
        if(Ratio.type == "classical"){
          eigenvalue_j =  eigen(Mjk)$values
        }
        
        ratio_j = (eigenvalue_j[-1] + sigma0/n)/(eigenvalue_j[-length(eigenvalue_j)] + sigma0/n)
        
        ratio_j = ratio_j[1:(0.5*length(ratio_j))]
        
        r_j = which.min(ratio_j)
        
        test_j =   eigenvalue_j[r_j]
        
        z_j = min(ratio_j)
        
        test_tol = c(test_tol,test_j)
        z_tol = c(z_tol,z_j)
        
        r_tol = c(r_tol,r_j)
      }
      
      test_list = rbind(test_list, test_tol)
      z_list = rbind(z_list, z_tol)
      r_list = rbind(r_list, r_tol)
      
    }
    
    place1 = which.min(rowMeans(z_list))
    
    #place2 =  min(apply(z_list,2,function(x){which(diff(x) > -10^-8)[1]}))
    
    delta_sel_1 = net_delta[place1]
    
    
    res.1 = est.ABr(Sigma.tensor.Y.k_list_orginal,Rank = Rank,delta = delta_sel_1,sigma0)
    
    res.3 = est.ABr(Sigma.tensor.Y.k_list_orginal,Rank = Rank,delta = 0,sigma0)
    
    
  }else{ #no thresholding or given delta thresholding
    z_list = r_list =  test_list = NULL
    
    res.1 = est.ABr(Sigma.tensor.Y.k_list_orginal,Rank = Rank,delta = delta,sigma0,augmented=augmented)
    res.3 = est.ABr(Sigma.tensor.Y.k_list_orginal,Rank = Rank,delta = 0,sigma0,augmented=augmented)
  }
  
  return(list(results =  res.1, results.NT =  res.3, Sigma.tensor.Y.k_list_orginal = Sigma.tensor.Y.k_list_orginal, sigma0 = sigma0, z_list = z_list, r_list =  r_list, test_list = test_list))
  
}

#
# CP.iter =  function(A.hat, Sigma_orginal, n, delta2, sigma0, augmented = F, iter_max = 100,eps = 10^-3, print.eps = F){
#
#   K =  length(Sigma_orginal)
#   m =  length(A.hat)
#
#   A.tol = A.1 = A.0 = A.hat
#   r = NCOL(A.0[[1]])
#
#   z1_tol = z2_tol = z3_tol =vector()
#
#   for (ll in 1:iter_max) {
#
#     for (j in 1:m) {
#
#       delta_2j = delta2[j]
#
#
#       if(m >= 3){
#         B.tilde.j = rTensor::khatri_rao_list(A.0[c(m:1)[-(m-j+1)]])
#
#       }else{
#         B.tilde.j = A.0[c(1:m)[-j]][[1]]
#       }
#
#
#       Sigma.j.k_list = list()
#
#       for (kk in 1:K) {
#
#         Sigma.Yj.k = Mat.tensor(Sigma_orginal[[kk]],j)
#
#         if(kk == 1){
#           Sigma.j.k = Threshold.Tensor( Sigma.Yj.k %*% B.tilde.j , n , sigma0 , delta_2j)
#           index     = Sigma.j.k
#         }
#         if(kk > 1){
#           Sigma.j.k = Sigma.Yj.k %*% B.tilde.j
#           Sigma.j.k[which(index == 0)] <- 0
#         }
#         Sigma.j.k_list[[kk]] = Sigma.j.k
#       }
#
#       bb1 = 0
#       bb2 = 0
#       if(augmented == T){
#         for (kk in 2:K) {
#           bb1 = bb1 + Sigma.j.k_list[[kk-1]]
#           bb2 = bb2 + Sigma.j.k_list[[kk]]
#         }
#       }else{
#         bb1 = Sigma.j.k_list[[1]]
#         bb2 = Sigma.j.k_list[[2]]
#       }
#
#       if(ll == 1){
#         z1_tol[j] = sum((svd(bb1)$d))
#         z2_tol[j] = sum((svd(bb1)$d)^2) # trace
#         z3_tol[j] = (svd(bb1)$d)[r]
#       }
#
#       K21_j_tilde = bb2%*%MASS::ginv(t(bb1)%*%bb1)%*%t(bb1)
#
#       A.1[[j]] = Re(eigen(K21_j_tilde)$vectors[,1:r])
#
#       wps = abs(sum(rho2.loss.list(A.1,A.0)))
#
#       if(print.eps == T){
#         cat("\r",round(wps,4))
#       }
#       A.0 = A.1
#
#     }
#
#     if(abs(sum(rho2.loss.list(A.tol,A.0))) < eps){
#       break
#     }
#
#     A.tol = A.0
#
#   }
#
#
#   return(list(A.hat = A.tol, A.inl = A.hat, r = r, delta2_sel = delta2, Sigma_orginal = Sigma_orginal, iter_step = ll, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
#
# }

# CP.iter.BMP.xi.han =  function(A.hat, Sigma_orginal,Y, n, delta2, sigma0, augmented = F, iter_max = 100,eps = 10^-5, print.eps = T){
#
#   K =  length(Sigma_orginal)
#   m =  length(A.hat)
#
#   A.tol = A.tol.han = A.1 = A.1.han = A.0 = A.0.han = A.hat
#   r = NCOL(A.0[[1]])
#
#   z1_tol = z2_tol = z3_tol =vector()
#
#   wps.tol = vector()
#
#   for (ll in 1:iter_max) {
#
#     for(j in 1:m){
#       delta_2j = delta2[j]
#
#       Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j]))
#
#       if(m >= 3){
#         B.tilde.j = rTensor::khatri_rao_list(A.0[c(m:1)[-(m-j+1)]])
#       }else{
#         B.tilde.j = A.0[c(1:m)[-j]][[1]]
#       }
#
#       B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
#
#       if(m >= 3){
#         dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m+1):3])) # 将d3和d4合并为一个维度
#       }
#
#       Uj = rTensor::ttl(as.tensor(Yp),list_mat = list(t(B.MP.j)), ms = 3 )@data
#
#       for (k in 1:r) {
#
#         vjk = aperm(Y,c(j+1,c(2:(m+1))[-j],1))
#
#         for(jj in (1:m)[-j]){
#           A.tilde.j =  A.0.han[[jj]]
#           A.MPk = (A.tilde.j%*%MASS::ginv(t(A.tilde.j)%*%A.tilde.j))[,k]
#           vjk <- tensor(vjk,A.MPk,2,1)
#         }
#
#         vjk = t(vjk)
#         ujk = Uj[,,k]
#
#         sigma.k.vjk = sigma_k(t(vjk),as.matrix(colMeans(vjk)),1,n)
#
#         sigma.k.ujk = sigma_k(t(ujk),as.matrix(colMeans(ujk)),1,n)
#
#         sigma.k.vjk = Threshold.Tensor(sigma.k.vjk,n,sqrt(mean(vjk^2)),delta_2j)
#
#         sigma.k.ujk = Threshold.Tensor(sigma.k.ujk,n,sqrt(mean(ujk^2)),delta_2j)
#
#         M.vjk <- sigma.k.vjk/2+t(sigma.k.vjk)/2
#         M.ujk <- sigma.k.ujk/2+t(sigma.k.ujk)/2
#
#         eig.vjk <- eigen(M.vjk)
#         eig.ujk <- eigen(M.ujk)
#
#         A.1[[j]][,k]     <- eig.ujk$vectors[,1]
#         A.1.han[[j]][,k] <- eig.vjk$vectors[,1]
#
#       }
#
#       wps = c(sqrt(sum(rho2.loss.list(A.1.han,A.0.han))),sqrt(sum(rho2.loss.list(A.1,A.0))))
#
#       wps.tol = rbind(wps.tol,wps)
#
#       if(print.eps == T){
#         cat("\r",round(wps,6))
#       }
#       A.0 = A.1
#       A.0.han = A.1.han
#     }
#
#     if( wps[1] < eps){
#       break
#     }
#
#   }
#
#   A.tol = A.0
#   A.tol.han = A.0.han
#
#   return(list(A.hat = A.tol,A.hat.2 = A.tol.han, A.inl = A.hat, r = r, delta2_sel = delta2, iter_step = ll, wps.tol = wps.tol, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
#
# }



CP.iter.AMP.xi =  function(A.hat, K, Y, n, delta2, sigma0, augmented = F, iter_max = 20,eps = 10^-5, print.eps = T,A = NULL,iter_lag = 1,Componet = NULL,...){
  
  m =  length(A.hat)
  
  A.tol = A.tol.han = A.1 = A.1.han = A.0 = A.0.han = A.hat
  r = NCOL(A.0[[1]])
  
  A.tol.list = list()
  
  z1_tol = z2_tol = z3_tol =vector()
  
  wps.tol = vector()
  
  iter.error = vector() # it should be a m*iter_max vector: 1-m iter error of first mode;(m+1)-(2m) iter error of second mode;...
  
  iter.error.mat = matrix(0,iter_max+1,1)
  
  
  fnorm.resid <- vector()
  
  for (ll in 1:iter_max) {
    
    if(!is.null(A)){
      iter.error.mat[c(ll:(iter_max+1)),] = matrix( max(rho2.loss.list(A.0,A)), length(c(ll:(iter_max+1))), 1 , byrow = T )
    }
    
    A.tol = A.0
    
    for(j in 1:m){
      delta_2j = delta2[j]
      
      Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
      
      
      A.MP = list()
      for (jjj in 1:m) {
        A.MP[[jjj]] = A.0[[jjj]]%*%MASS::ginv(t(A.0[[jjj]])%*%A.0[[jjj]])
      }
      
      if(m >= 3){
        B.MP.j = rTensor::khatri_rao_list(A.MP[c(m:1)[-(m-j+1)]])
      }else{
        B.MP.j = A.MP[c(1:m)[-j]][[1]]
      }
      
      if(m >= 3){
        B.tilde.j = rTensor::khatri_rao_list(A.0[c(m:1)[-(m-j+1)]])
      }else{
        B.tilde.j = A.0[c(1:m)[-j]][[1]]
      }
      
      A.tilde.j  = A.0[[j]]
      
      #B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
      A.MP.j = A.tilde.j%*%MASS::ginv(t(A.tilde.j)%*%A.tilde.j)
      
      
      if(m >= 3){
        dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m+1):3])) # combine dj j \neq i into one mode, then Yp be a n x di x prod(dj, j \neq i)  3 dimensional tensor
      }
      
      Uj = rTensor::ttl(as.tensor(Yp),list_mat = list(t(A.MP.j),t(B.MP.j)), ms = c(2,3) )@data
      
      # if(j == 1 && ll == 1){
      #   Uj0 = Uj
      # }
      #
      # #print(abs(cor(Uj,Uj0)))
      #
      # if(ll > 20 && NCOL(Uj) == 1 && abs(cor(Uj,Uj0)) < 0.95){
      #   Uj = Uj0
      # }
      #
      # Uj0 = Uj
      
      if(NCOL(Uj) > 1){
        Uj.mat = t(apply(Uj, 1, diag))
      }else{
        Uj.mat = as.matrix(Uj)
      }
      
      
      if(j == 1){
        Yp_hat = 0
        for (k in 1:r) {
          Yp_hat = Yp_hat +  Uj.mat[,k] %o% A.tilde.j[,k] %o% B.tilde.j[,k]
        }
        fnorm.resid[ll] <- fnorm(Yp-Yp_hat)/fnorm(Yp)
      }
      
      if(ll == 1){
        Uj.mat.inl =  as.matrix(Uj.mat)
      }
      
      Uj.reg = as.matrix( Uj.mat )
      
      #sigma_breve = sqrt(mean(Uj.reg^2))
      
      for (k in 1:r) {
        A.fnorm = vector()
        A.breve.jk.list = vector()
        for (qq in 1:iter_lag) {
          
          if(r > 1){
            ujk = lm(Uj.reg[-c((n-qq+1):n),k] ~ Uj.reg[-c(1:qq),-k]-1)$residuals #Uj.reg[,k]
            sigma.k.ujk = tensor.Autocov_xi_Y(Y,c(scale(ujk),rep(0,qq)),qq)
          }else{
            ujk = Uj.reg
            sigma.k.ujk = tensor.Autocov_xi_Y(Y,scale(ujk),qq)
          }
          
          
          # ujk = lm(Uj.reg[,k] ~ Uj.reg[,-k]-1)$residuals #Uj.reg[,k]
          # sigma.k.ujk = tensor.Autocov_xi_Y(Y,scale(ujk),1)
          
          A.breve.jk = Mat.tensor(sigma.k.ujk,j)%*%B.MP.j[,k]
          
          A.fnorm[qq] = fnorm(A.breve.jk)
          
          A.breve.jk.list = as.matrix(cbind(A.breve.jk.list,A.breve.jk))
          
          
        }
        
        #print(apply(A.breve.jk.list, 2, l2s))
        
        #print(A.fnorm)
        
        A.breve.jk = as.matrix(A.breve.jk.list[,which.max(A.fnorm)])
        
        A.breve.jk.thres = Threshold.Tensor(A.breve.jk, n, sigma0, delta_2j)
        
        if(sum(A.breve.jk.thres) < 10^-8){
          A.breve.jk.thres = Threshold.Tensor(A.breve.jk, n, sigma0, 0.1*delta_2j)
        }
        
        if(sum(A.breve.jk.thres) < 10^-8){
          A.breve.jk.thres = Threshold.Tensor(A.breve.jk, n, sigma0, 0)
        }
        
        A.1[[j]][,k] <- apply(A.breve.jk.thres, 2, l2s)
        
      }
      
      A.0 = A.1
      
    }
    
    A.tol.new = A.1
    
    A.tol.list[[ll]] = A.tol.new
    
    wps = sum(sqrt(abs(rho2.loss.list(A.tol.new,A.tol))))
    
    wps.tol = rbind(wps.tol,wps)
    
    if(print.eps == T){
      cat("\n",round(wps,6))
    }
    
    if( wps < eps){
      break
    }
    
  }
  
  iter.error = c(iter.error.mat)
  
  A.tol = A.tol.new
  
  if(ll == iter_max){
    A.tol = A.tol.list[[which.min(fnorm.resid)]]
  }
  
  
  j = 1
  Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  if(m >= 3){
    B.tilde.j = rTensor::khatri_rao_list(A.tol[c(m:1)[-(m-j+1)]])
  }else{
    B.tilde.j = A.tol[c(1:m)[-j]][[1]]
  }
  A.tilde.j  = A.tol[[j]]
  
  B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
  A.MP.j = A.tilde.j%*%MASS::ginv(t(A.tilde.j)%*%A.tilde.j)
  if(m >= 3){
    dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m+1):3])) # combine dj j \neq i into one mode, then Yp be a n x di x prod(dj, j \neq i)  3 dimensional tensor
  }
  
  Uj = rTensor::ttl(as.tensor(Yp),list_mat = list(t(A.MP.j),t(B.MP.j)), ms = c(2,3) )@data
  
  if(r > 1){
    Uj.mat = t(apply(Uj, 1, diag))
  }else{
    Uj.mat = as.matrix(Uj)
  }
  
  
  
  return(list(A.hat = A.tol, A.inl = A.hat, A.tol.list = A.tol.list, r = r, delta2_sel = delta2, iter_step = ll,iter_error = iter.error, wps.tol = wps.tol, fnorm.resid = fnorm.resid,f_hat = Uj.mat, f_hat_inl = Uj.mat.inl, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
  
}


CP.iter.BMP.xi =  function(A.hat, K, Y, n, delta2, sigma0, augmented = F, iter_max = 20,eps = 10^-5, print.eps = T,A = NULL,iter_lag = 1, Componet = NULL,...){
  
  m =  length(A.hat)
  
  A.tol = A.tol.han = A.1 = A.1.han = A.0 = A.0.han = Sigma.yij.xii.1 = A.hat
  r = NCOL(A.0[[1]])
  
  A.tol.list = list()
  
  z1_tol = z2_tol = z3_tol =vector()
  
  wps.tol = vector()
  
  iter.error = vector() # it should be a m*iter_max vector: 1-m iter error of first mode;(m+1)-(2m) iter error of second mode;...
  
  iter.error.mat = matrix(0,iter_max+1,1)
  
  
  fnorm.resid <- vector()
  
  for (ll in 1:iter_max) {
    
    if(!is.null(A)){
      iter.error.mat[c(ll:(iter_max+1)),] = matrix( max(rho2.loss.list(A.0,A)), length(c(ll:(iter_max+1))), 1 , byrow = T )
    }
    
    A.tol = A.0
    
    for(j in 1:m){
      delta_2j = delta2[j]
      
      Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
      
      if(m >= 3){
        B.tilde.j = rTensor::khatri_rao_list(A.0[c(m:1)[-(m-j+1)]])
      }else{
        B.tilde.j = A.0[c(1:m)[-j]][[1]]
      }
      A.tilde.j  = A.0[[j]]
      
      B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
      A.MP.j = A.tilde.j%*%MASS::ginv(t(A.tilde.j)%*%A.tilde.j)
      
      
      if(m >= 3){
        dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m+1):3])) # combine dj j \neq i into one mode, then Yp be a n x di x prod(dj, j \neq i)  3 dimensional tensor
      }
      
      Uj = rTensor::ttl(as.tensor(Yp),list_mat = list(t(A.MP.j),t(B.MP.j)), ms = c(2,3) )@data
      
      # if(j == 1 && ll == 1){
      #   Uj0 = Uj
      # }
      #
      # #print(abs(cor(Uj,Uj0)))
      #
      # if(ll > 20 && NCOL(Uj) == 1 && abs(cor(Uj,Uj0)) < 0.95){
      #   Uj = Uj0
      # }
      #
      # Uj0 = Uj
      
      if(NCOL(Uj) > 1){
        Uj.mat = t(apply(Uj, 1, diag))
      }else{
        Uj.mat = as.matrix(Uj)
      }
      
      
      if(j == 1){
        Yp_hat = 0
        for (k in 1:r) {
          Yp_hat = Yp_hat +  Uj.mat[,k] %o% A.tilde.j[,k] %o% B.tilde.j[,k]
        }
        fnorm.resid[ll] <- fnorm(Yp-Yp_hat)/fnorm(Yp)
        
        if(ll == 1){
          Uj.mat.inl =  as.matrix(Uj.mat)
          A.tilde.j.inl = A.tilde.j
          B.tilde.j.inl = B.tilde.j
        }
        
      }
      
      
      
      Uj.reg = as.matrix( Uj.mat )
      
      #sigma_breve = sqrt(mean(Uj.reg^2))
      
      for (k in 1:r) {
        A.fnorm = vector()
        A.breve.jk.list = vector()
        for (qq in 1:iter_lag) {
          
          if(r > 1){
            ujk = lm(Uj.reg[-c((n-qq+1):n),k] ~ Uj.reg[-c(1:qq),-k]-1)$residuals #Uj.reg[,k]
            sigma.k.ujk = tensor.Autocov_xi_Y(Y,c(scale(ujk),rep(0,qq)),qq)
          }else{
            ujk = Uj.reg
            sigma.k.ujk = tensor.Autocov_xi_Y(Y,scale(ujk),qq)
          }
          
          
          # ujk = lm(Uj.reg[,k] ~ Uj.reg[,-k]-1)$residuals #Uj.reg[,k]
          # sigma.k.ujk = tensor.Autocov_xi_Y(Y,scale(ujk),1)
          
          A.breve.jk = Mat.tensor(sigma.k.ujk,j)%*%B.MP.j[,k]
          
          A.fnorm[qq] = fnorm(A.breve.jk)
          
          A.breve.jk.list = as.matrix(cbind(A.breve.jk.list,A.breve.jk))
          
          
        }
        
        #print(apply(A.breve.jk.list, 2, l2s))
        
        #print(A.fnorm)
        
        A.breve.jk = as.matrix(A.breve.jk.list[,which.max(A.fnorm)])
        
        A.breve.jk.thres = Threshold.Tensor(A.breve.jk, n, sigma0, delta_2j)
        
        if(sum(A.breve.jk.thres) < 0.01){
          A.breve.jk.thres = Threshold.Tensor(A.breve.jk, n, sigma0, 0)
        }
        
        Sigma.yij.xii.1[[j]][,k] <- A.breve.jk.thres
        
        A.1[[j]][,k] <- apply(A.breve.jk.thres, 2, l2s)
        
      }
      
      A.0 = A.1
      
    }
    
    A.tol.new = A.1
    
    A.tol.list[[ll]] = A.tol.new
    
    wps = sum(sqrt(abs(rho2.loss.list(A.tol.new,A.tol))))
    
    wps.tol = rbind(wps.tol,wps)
    
    if(print.eps == T){
      cat("\n",round(wps,6))
    }
    
    if( wps < eps){
      break
    }
    
  }
  
  iter.error = c(iter.error.mat)
  
  A.tol = A.tol.new
  
  if(ll == iter_max){
    A.tol = A.tol.list[[which.min(fnorm.resid)]]
  }
  
  
  j = 1
  Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  if(m >= 3){
    B.tilde.j = rTensor::khatri_rao_list(A.tol[c(m:1)[-(m-j+1)]])
  }else{
    B.tilde.j = A.tol[c(1:m)[-j]][[1]]
  }
  A.tilde.j  = A.tol[[j]]
  
  B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
  A.MP.j = A.tilde.j%*%MASS::ginv(t(A.tilde.j)%*%A.tilde.j)
  if(m >= 3){
    dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m+1):3])) # combine dj j \neq i into one mode, then Yp be a n x di x prod(dj, j \neq i)  3 dimensional tensor
  }
  
  Uj = rTensor::ttl(as.tensor(Yp),list_mat = list(t(A.MP.j),t(B.MP.j)), ms = c(2,3) )@data
  
  if(r > 1){
    Uj.mat = t(apply(Uj, 1, diag))
  }else{
    Uj.mat = as.matrix(Uj)
  }
  
  CP_loss_iter = CP_loss_inl = -999
  if(!is.null(Componet)){
    
    Cp <- aperm(Componet, c(1,j+1,c(2:(m+1))[-j]))
    Yp_hat = 0
    Yp_hat_inl = 0
    for (k in 1:r) {
      Yp_hat = Yp_hat +  Uj.mat[,k] %o% A.tilde.j[,k] %o% B.tilde.j[,k]
      Yp_hat_inl = Yp_hat_inl +  Uj.mat.inl[,k] %o% A.tilde.j.inl[,k] %o% B.tilde.j.inl[,k]
    }
    CP_loss_iter = fnorm(Yp_hat - Cp)/sqrt(n*prod(D))
    CP_loss_inl  = fnorm(Yp_hat_inl - Cp)/sqrt(n*prod(D))
  }
  
  
  
  
  return(list(A.hat = A.tol, A.inl = A.hat, A.tol.list = A.tol.list, Sigma.yij.xii.1 = Sigma.yij.xii.1, r = r, delta2_sel = delta2, iter_step = ll,iter_error = iter.error, wps.tol = wps.tol, fnorm.resid = fnorm.resid,f_hat = Uj.mat, f_hat_inl = Uj.mat.inl, CP_loss = c(CP_loss_inl,CP_loss_iter), Yp_hat = Yp_hat, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
  
}


CP.iter.BMP =  function(A.hat, Sigma_orginal, n, delta2, sigma0, augmented = F, iter_max = 100,eps = 10^-5, print.eps = F){
  
  K =  length(Sigma_orginal)
  m =  length(A.hat)
  
  A.tol = A.1 = A.0 = A.hat
  r = NCOL(A.0[[1]])
  
  z1_tol = z2_tol = z3_tol =vector()
  wps.tol = vector()
  for (ll in 1:iter_max) {
    
    for (j in 1:m) {
      
      delta_2j = delta2[j]
      
      if(m >= 3){
        B.tilde.j = rTensor::khatri_rao_list(A.0[c(m:1)[-(m-j+1)]])
      }else{
        B.tilde.j = A.0[c(1:m)[-j]][[1]]
      }
      
      Sigma.Yj.k = Mat.tensor(Sigma_orginal[[1]],j)
      
      B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
      
      Sigma.j.k = Threshold.Tensor( Sigma.Yj.k %*% B.MP.j , n , sigma0 , delta_2j)
      
      bb1 = Sigma.j.k
      
      if(0 %in% svd(bb1)$d){
        bb1[which(bb1 == 0)] <- 0.01
      }
      
      if(ll == 1){
        z1_tol[j] = sum((svd(bb1)$d))
        z2_tol[j] = sum((svd(bb1)$d)^2) # trace
        z3_tol[j] = (svd(bb1)$d)[r]
      }
      
      A.1[[j]] = apply(bb1,2,l2s)
      
      
      wps = sqrt(abs(sum(rho2.loss.list(A.1,A.0))) )
      wps.tol = rbind(wps.tol,wps)
      if(print.eps == T){
        cat("\r",round(wps,6))
      }
      A.0 = A.1
      
    }
    
    if(wps < eps){
      break
    }
    
  }
  
  A.tol = A.0
  
  return(list(A.hat = A.tol, A.inl = A.hat, r = r, delta2_sel = delta2, Sigma_orginal = Sigma_orginal, iter_step = ll,wps.tol=wps.tol, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
  
}

CP.iter.BMP.col =  function(A.hat, Sigma_orginal, n, delta2, sigma0, augmented = F, iter_max = 100,eps = 10^-5, print.eps = F){
  
  K =  length(Sigma_orginal)
  m =  length(A.hat)
  
  A.tol = A.1 = A.0 = A.hat
  r = NCOL(A.0[[1]])
  
  z1_tol = z2_tol = z3_tol =vector()
  
  for (ll in 1:iter_max) {
    
    for (j in 1:m){
      
      delta_2j = delta2[j]
      
      for (k in 1:r) {
        
        if(m >= 3){
          B.tilde.j = rTensor::khatri_rao_list(A.0[c(m:1)[-(m-j+1)]])
        }else{
          B.tilde.j = A.0[c(1:m)[-j]][[1]]
        }
        
        Sigma.Yj.k = Mat.tensor(Sigma_orginal[[1]],j)
        
        B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
        
        Sigma.j.k = Threshold.Tensor( Sigma.Yj.k %*% B.MP.j , n , sigma0 , delta_2j)
        
        bb1 = Sigma.j.k
        
        A.1[[j]][,k] = apply(bb1,2,l2s)[,k]
        
        wps = abs(sum(rho2.loss.list(A.1,A.0)))
        
        if(print.eps == T){
          cat("\n",round(wps,4))
        }
        A.0 = A.1
        
      }
      
    }
    
    if(abs(sum(rho2.loss.list(A.tol,A.0))) < eps){
      break
    }
    
    A.tol = A.0
    
  }
  
  
  
  return(list(A.hat = A.tol, A.inl = A.hat, r = r, delta2_sel = delta2, Sigma_orginal = Sigma_orginal, iter_step = ll, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
  
}


CP.iter.AMP =  function(A.hat, Sigma_orginal, n, delta2, sigma0, augmented = F, iter_max = 100,eps = 10^-5, print.eps = F){
  
  K =  length(Sigma_orginal)
  m =  length(A.hat)
  
  A.tol = A.1 = A.0 = A.hat
  r = NCOL(A.0[[1]])
  
  z1_tol = z2_tol = z3_tol =vector()
  
  for (ll in 1:iter_max) {
    
    A.MP = list()
    for (j in 1:m) {
      A.MP[[j]] = A.0[[j]]%*%MASS::ginv(t(A.0[[j]])%*%A.0[[j]])
    }
    
    for (j in 1:m) {
      
      delta_2j = delta2[j]
      
      if(m >= 3){
        B.MP.j = rTensor::khatri_rao_list(A.MP[c(m:1)[-(m-j+1)]])
      }else{
        B.MP.j = A.MP[c(1:m)[-j]][[1]]
      }
      
      Sigma.Yj.k = Mat.tensor(Sigma_orginal[[1]],j)
      
      Sigma.j.k = Threshold.Tensor( Sigma.Yj.k %*% B.MP.j , n , sigma0 , delta_2j)
      
      bb1 = Sigma.j.k
      
      if(0 %in% svd(bb1)$d){
        bb1[which(bb1 == 0)] <- 0.01
      }
      
      if(ll == 1){
        z1_tol[j] = sum((svd(bb1)$d))
        z2_tol[j] = sum((svd(bb1)$d)^2) # trace
        z3_tol[j] = (svd(bb1)$d)[r]
      }
      
      A.1[[j]] = apply(bb1,2,l2s)
      
      wps = abs(sum(rho2.loss.list(A.1,A.0)))
      
      if(print.eps == T){
        cat("\r",round(wps,4))
      }
      A.0 = A.1
      
    }
    
    if(abs(sum(rho2.loss.list(A.tol,A.0))) < eps){
      break
    }
    
    A.tol = A.0
    
  }
  
  
  return(list(A.hat = A.tol, A.inl = A.hat, r = r, delta2_sel = delta2, Sigma_orginal = Sigma_orginal, iter_step = ll, z_tol = list(z1_tol=z1_tol,z2_tol=z2_tol,z3_tol=z3_tol)))
  
}



HDTTS.CP.iter = function(Y,  xi = NULL, Rank = NULL,lag.k = 5, Threshold = FALSE, delta = NULL, delta2 = NULL, Ratio.type = c("classical","log"), augmented = FALSE, iter_max = 10, eps = 10^-4, grid.num = 50, delta_max = 3,all.out = F, ev.all = T,A = NULL,...){
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  mse_tuning_delta2 = NULL
  
  res.inl =  HDTTS.CP.est(Y = Y,xi = xi, Rank = Rank, K = lag.k, Threshold = Threshold, delta = delta, Ratio.type = Ratio.type, augmented = augmented)
  
  A.hat = res.inl$results$A.hat
  A.hat.NT = res.inl$results.NT$A.hat
  
  Sigma_orginal = res.inl$Sigma.tensor.Y.k_list_orginal
  
  sigma0 = res.inl$sigma0
  
  if(is.null(delta2)){
    
    if(all.out == F){ #only delta2 output
      
      delta2_net =  matrix(rep(seq(0,delta_max,length.out = grid.num),m),grid.num,m)
      
      vv.mat = vv.mat.alter = vv.mat.NT = matrix(NA,grid.num,m)
      
      delta2_vv1 = delta2_vv2 = delta2 = delta2.alter = delta2.NT = rep(0,m)
      
      for(j in order(D,decreasing = F) ){
        
        vv1 = vv2 = vv = vv.alter = vv.NT = vector()
        for (kk in seq(0,delta_max,length.out = grid.num) ) {
          
          delta2[j] =  kk
          
          res.iter        =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max,eps = eps, print.eps = F)
          
          z_tol           =   res.iter$z_tol
          
          vv              =   rbind(vv     , z_tol$z2_tol) #### use z2_tol : trace
          
          # if(ev.all == T){
          #   vv1              =   rbind(vv1 , z_tol$z1_tol)
          #   vv2              =   rbind(vv2 , z_tol$z2_tol)
          # }
          #
          
          # cat("\r",paste(j,"-",round(kk,2)))
        }
        
        delta2[j]               =  delta2_net[which.max(vv[,j]),j]
        
        # if(ev.all == T){
        #   delta2_vv1[j]           =  delta2_net[which.max(vv1[,j]),j]
        #   delta2_vv2[j]           =  delta2_net[which.max(vv2[,j]),j]
        # }
        
        
      }
      
      res.iter.fin              =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      res.iter.NT.fin           =   CP.iter(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      
      # if(ev.all == T){
      #   res.iter.fin_vv1        =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv1, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #   res.iter.NT.fin_vv1     =   CP.iter(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv1, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #
      #   res.iter.fin_vv2        =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #   res.iter.NT.fin_vv2     =   CP.iter(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #
      #   }
      
      
    }else{
      
      delta2_net =  matrix(rep(seq(0,delta_max,length.out = grid.num),m),grid.num,m)
      
      delta2  = delta2.NT = rep(0,m)
      
      for(j in order(D,decreasing = F) ){
        
        vv  = vv.NT = vector()
        for (kk in seq(0,delta_max,length.out = grid.num) ) {
          
          delta2[j] =  kk
          delta2.NT[j] =  kk
          
          res.iter        =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max,eps = eps, print.eps = F)
          
          res.iter.NT     =   CP.iter(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2.NT, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
          
          vv              =  rbind(vv      , res.iter$z_tol$z2_tol) ###use z2
          vv.NT           =  rbind(vv.NT   , res.iter.NT$z_tol$z2_tol)
          
          # cat("\r",paste(j,"-",round(kk,2)))
        }
        
        delta2[j]       =  delta2_net[which.max(vv[,j]),j]
        delta2.NT[j]    =  delta2_net[which.max(vv.NT[,j]),j]
      }
      
      res.iter.fin        =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      res.iter.NT.fin     =   CP.iter(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2.NT, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      
    }
    
  }else{
    
    
    res.iter.fin        =   CP.iter(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
    
    res.iter.NT.fin     =   CP.iter(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
    
    
  }
  
  if(!is.null(A)){
    rho.iter.fin            =  mean( rho2.loss.list(res.iter.fin$A.hat,A) )
    rho.iter.NT.fin         =  mean( rho2.loss.list(res.iter.NT.fin$A.hat,A) )
    rho.inl                 =  mean( rho2.loss.list(A.hat,A) )
    rho.inl.NT              =  mean( rho2.loss.list(A.hat.NT,A) )
    
    mse_tuning_delta2 = c(rho.iter.fin,
                          rho.iter.NT.fin,
                          rho.inl,
                          rho.inl.NT)
  }else{
    mse_tuning_delta2 = NULL
  }
  
  
  return(list(  res.iter   =  res.iter.fin,  res.iter.NT     =  res.iter.NT.fin, res.inl = res.inl, mse_tuning_delta2 = mse_tuning_delta2 ))
  
}

Boostrap.xi.sel = function(Y, r_breve = NULL,lag.k = 10, Threshold = FALSE, delta = NULL, delta2 = NULL,  Randomized.time = 20, Ratio.type = c("classical","log"), augmented = FALSE, all.out = F, ev.all = T,A = NULL,...){
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  Rank = r_breve
  
  mse_tuning_delta2 = NULL
  
  A.hat.NT_list = list()
  xi_list = list()
  
  for (ss in 1:Randomized.time) {
    xi = tensor.est.xi(Y,random = T)
    res.inl =  HDTTS.CP.est(Y = Y,xi =  xi, Rank = Rank, K = lag.k, Threshold = F, delta = 0, Ratio.type = Ratio.type, augmented = F)
    xi_list[[ss]] = xi
    A.hat.NT_list[[ss]] = res.inl$results.NT$A.hat
    
    #print(paste0("Randomized.iter-",ss))
  }
  
  G  =  matrix(1,Randomized.time,Randomized.time)
  H  =  vector()
  
  for (jj in 1:(Randomized.time-1)) {
    for (kk in (jj+1):Randomized.time){
      G[kk,jj] = G[jj,kk] = max(psi2.loss.list(A.hat.NT_list[[jj]],A.hat.NT_list[[kk]]))
      H = rbind(H,c(jj,kk,max(psi2.loss.list(A.hat.NT_list[[jj]],A.hat.NT_list[[kk]]))))
    }
  }
  
  
  
  diag(G) <- 0
  
  eps0 =  quantile(H[,3],0.1)
  
  res_graph = build_eps_graph(D = G, eps = eps0, weight = "gaussian")
  
  place = res_graph$idx_strength_max
  
  # res_graph$idx_deg_max
  #
  # res_graph$degree
  # res_graph$strength
  #
  # plot(res_graph$graph)
  #
  
  xi.sel =  xi_list[[place]]
  
  
  return(list(xi.sel = xi.sel,H = H, G = G,  res_graph = res_graph))
  
}

Boostrap.xi.sel.new = function(Y, r_breve = NULL,lag.k = 10, Threshold = FALSE, delta = NULL, delta2 = NULL,  Randomized.time = 20, Ratio.type = c("classical","log"), augmented = FALSE, all.out = F, ev.all = T,A = NULL,...){
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  Rank = r_breve
  
  mse_tuning_delta2 = NULL
  
  gg = 0
  iter = 1
  while (gg  == 0 && iter < 20) {
    
    
    A.hat.NT_list = list()
    xi_list = list()
    
    for (ss in 1:Randomized.time) {
      xi = tensor.est.xi(Y, random = T)
      res.inl =  HDTTS.CP.est(Y = Y,xi =  xi, Rank = Rank, K = lag.k, Threshold = F, delta = 0, Ratio.type = "log", augmented = F)
      xi_list[[ss]] = xi
      A.hat.NT_list[[ss]] = res.inl$results.NT$A.hat
      
      #print(paste0("Randomized.iter-",ss))
    }
    
    G  =  matrix(1,Randomized.time,Randomized.time)
    H  =  vector()
    
    diag(G) <- 0
    
    eps0 =  0.1
    # min(min(apply(H[,-c(1:2)],2, function(x){quantile(x,0.3)} )),0.05)
    res_graph = W_all = list()
    strength_all = 0
    for (vv in 1:r_breve) {
      G_i = matrix(0,Randomized.time,Randomized.time)
      
      for (jj in 1:(Randomized.time-1)){
        for (kk in (jj+1):Randomized.time){
          H = rbind(H,c(jj,kk, vecpsi.loss.list(A.hat.NT_list[[jj]],A.hat.NT_list[[kk]])) )
          G_i[kk,jj] = G_i[jj,kk] = vecpsi.loss.list(A.hat.NT_list[[jj]],A.hat.NT_list[[kk]])[vv]
        }
      }
      res_graph_i = build_eps_graph(D = G_i, eps = eps0, weight = "gaussian")
      
      W_i = res_graph_i$W
      
      W_all[[vv]] =  W_i
      
      res_graph[[vv]] = res_graph_i
      
    }
    
    gg = sum(find_max_pair_matrix(W_all)$count)
    iter = iter + 1
    #print(gg)
  }
  
  if(iter == 20){
    xi.sel = tensor.est.xi(Y,d_max = 1)
  }else{
    place = find_max_pair_matrix(W_all)$i[1]
    xi.sel =  xi_list[[place]]
  }
  
  return(list(xi.sel = xi.sel,H = H, G = G,  res_graph = res_graph))
  
}


Boostrap.xi.sel.v3 = function(Y, r_breve = NULL, eps = 0.1, lag.k = 10, Threshold = FALSE, delta = NULL, delta2 = NULL,  Randomized.time = 20, Ratio.type = c("classical","log"), augmented = FALSE, all.out = F, ev.all = T,A = NULL,...){
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  Rank = r_breve
  
  mse_tuning_delta2 = NULL
  
  gg = 0
  iter = 1
  while (gg  == 0 && iter < 20) {
    
    
    A.hat.NT_list = list()
    xi_list = list()
    
    for (ss in 1:Randomized.time) {
      xi = tensor.est.xi(Y, random = T)
      res.inl =  HDTTS.CP.est(Y = Y,xi =  xi, Rank = Rank, K = lag.k, Threshold = F, delta = 0, Ratio.type = "log", augmented = F)
      xi_list[[ss]] = xi
      A.hat.NT_list[[ss]] = res.inl$results.NT$A.hat
      
      #print(paste0("Randomized.iter-",ss))
    }
    
    G  =  matrix(1,Randomized.time,Randomized.time)
    H  =  vector()
    
    diag(G) <- 100
    
    eps0 =  eps
    # min(min(apply(H[,-c(1:2)],2, function(x){quantile(x,0.3)} )),0.05)
    res_graph = W_all = list()
    strength_all = 0
    DD_tol = vector()
    for (vv in 1:r_breve) {
      G_i = matrix(0,Randomized.time,Randomized.time)
      
      for (jj in 1:(Randomized.time)){
        for (kk in 1:Randomized.time){
          if(jj != kk){
            H = rbind(H,c(jj,kk, vecpsi.loss.list(A.hat.NT_list[[jj]],A.hat.NT_list[[kk]])) )
            D_ijk  = abs(vecpsi.loss.list(A.hat.NT_list[[jj]],A.hat.NT_list[[kk]])[vv])
            G_i[jj,kk] = ifelse(D_ijk < eps0, 1, 0)
          }
          
        }
      }
      
      
      DD_tol = rbind(DD_tol, colSums(G_i))
      
    }
    
    gg = max(colSums(DD_tol))
    iter = iter + 1
    #print(gg)
  }
  
  if(iter == 20){
    xi.sel = tensor.est.xi(Y,d_max = 1)
  }else{
    place = which.max(colSums(DD_tol))
    xi.sel =  xi_list[[place]]
  }
  
  return(list(xi.sel = xi.sel,H = H, G = G))
  
}


find_max_pair_matrix <- function(adj_list){
  C <- Reduce(`+`, adj_list)
  C[lower.tri(C, diag=TRUE)] <- -Inf
  mx <- max(C)
  idx <- which(C == mx, arr.ind = TRUE)
  data.frame(i = idx[,1], j = idx[,2], count = mx)
}

# Build an epsilon-graph from a pairwise distance matrix and pick central nodes
# D   : n x n symmetric distance matrix with zeros on the diagonal
# eps : epsilon threshold; connect i-j if distance(i,j) <= eps
# weight :
#   "linear"   -> similarity w_ij = 1 - D_ij/eps  (in [0,1], closer = larger)
#   "inverse"  -> similarity w_ij = 1/(D_ij + 1e-9)
#   "gaussian" -> similarity w_ij = exp(-D_ij^2 / (2*sigma^2))
# sigma : bandwidth for the Gaussian similarity; defaults to eps/3 if NULL
#
# Returns a list with:
#   graph              : igraph object (undirected, weighted by similarity)
#   degree             : unweighted degree of each node (number of neighbors)
#   idx_deg_max        : indices of nodes with maximal degree
#   strength           : weighted degree (sum of similarities) of each node
#   idx_strength_max   : indices of nodes with maximal strength
#   W                  : similarity (weight) matrix actually used

build_eps_graph <- function(D, eps, weight = c("linear","inverse","gaussian"), sigma = NULL) {
  require(igraph)
  stopifnot(is.matrix(D), nrow(D) == ncol(D))
  n <- nrow(D)
  weight <- match.arg(weight)
  
  # Prevent self-loops: set diagonal to Inf so it will never pass the eps filter
  D2 <- D
  diag(D2) <- Inf
  
  # Boolean mask of edges that satisfy the epsilon threshold
  A_mask <- (D2 <= eps)
  
  # Construct a nonnegative, symmetric similarity matrix W
  W <- matrix(0, n, n)
  if (weight == "linear") {
    # Scale to [0,1]; edges farther than eps get weight 0
    W[A_mask] <- pmax(0, 1 - D2[A_mask] / eps)
  } else if (weight == "inverse") {
    # Larger weight for smaller distance
    W[A_mask] <- 1 / (D2[A_mask] + 1e-9)
  } else { # "gaussian"
    if (is.null(sigma)) sigma <- eps / 3
    W[A_mask] <- exp(-(D2[A_mask]^2) / (2 * sigma^2))
  }
  # Force symmetry and remove self-weights
  W <- (W + t(W)) / 2
  diag(W) <- 0
  
  # Build an undirected, weighted graph from the similarity matrix
  g <- graph_from_adjacency_matrix(W, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  
  # Unweighted degree: number of neighbors (ignores weights)
  deg <- degree(g, mode = "all")
  idx_deg_max <- which(deg == max(deg))
  
  # Strength: weighted degree = sum of incident edge weights (similarities)
  strg <- strength(g, vids = V(g), mode = "all", weights = E(g)$weight)
  idx_strength_max <- which(strg == max(strg))
  
  list(graph = g,
       degree = deg,           idx_deg_max = idx_deg_max,
       strength = strg,        idx_strength_max = idx_strength_max,
       W = W, eps = eps)
}

# Boostrap.xi.sel.Q = function(Y,  xi = NULL, Rank = Rank,lag.k = 5, Threshold = FALSE, delta = NULL, delta2 = NULL, Randomized.xi = T, Randomized.time = 20, Ratio.type = c("classical","log"), augmented = FALSE, all.out = F, ev.all = T,A = NULL,print.eps=F,iter_max = 10,eps.ghost = NULL,cv.rule = 1,...){
#
#   n = dim(Y)[1]
#   D = dim(Y)[-1]
#   m = length(D)
#
#   r = Rank
#
#   mse_tuning_delta2 = NULL
#
#   A.hat.NT_list = list()
#   xi_list = list()
#   Q_list = list()
#
#   zz = 0
#
#
#
#      for (ss in 1:Randomized.time) {
#         xi = tensor.est.xi(Y,random = T)
#         res.inl =  HDTTS.CP.est(Y = Y,xi =  xi, Rank = Rank, K = lag.k, Threshold = F, delta = 0, Randomized.xi=Randomized.xi, Ratio.type = Ratio.type, augmented = augmented)
#
#         Q_list[[ss]] = res.inl$results.NT$Q_list
#         xi_list[[ss]] = xi
#         A.hat.NT_list[[ss]] = res.inl$results.NT$A.hat
#
#         #print(paste0("Randomized.iter-",ss))
#       }
#
#       G  =  matrix(1,Randomized.time,Randomized.time)
#       H  =  vector()
#
#       for (jj in 1:(Randomized.time-1)) {
#         for (kk in (jj+1):Randomized.time) {
#           G[kk,jj] = G[jj,kk] = max(D.loss.list(Q_list[[jj]],Q_list[[kk]]))
#           H = rbind(H,c(jj,kk,max(D.loss.list(Q_list[[jj]],Q_list[[kk]]))))
#         }
#       }
#       diag(G) <- 0
#
#       eps0 =  quantile(H[,3],0.05)
#
#       res_graph = build_eps_graph(D = G, eps = eps0, weight = "gaussian")
#
#       res_graph$idx_strength_max
#       res_graph$idx_deg_max
#
#       res_graph$degree
#       res_graph$strength
#
#       plot(res_graph$graph)
#
#
#       rho2.loss.list(A.hat.NT_list[[19]],A)
#
#       plot(res_graph$graph)
#
#       pos   =  res_graph$idx_strength_max
#
#
#
#   xi.sel =  xi_list[[pos[1]]]
#
#   if(is.null(xi.sel)){
#     xi.sel = xi_list[[ghost]]
#   }
#
#   return(list(xi.sel = xi.sel,H = H, G = G))
#
# }
#


HDTTS.CP.iter.BMP = function(Y,  xi = NULL, Rank = NULL,lag.k = 5, Threshold = FALSE, delta = NULL, delta2 = NULL, Ratio.type = c("classical","log"), Random.Project = F, augmented = FALSE, BMP = T,iter_max = 20, eps = 10^-4, grid.num = 50, delta_max = 3,all.out = F, ev.all = T,A = NULL,A.inl = NULL,print.eps=F,iter_lag=1,Componet = NULL,...){
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  mse_tuning_delta2 = NULL
  
  if(is.null(xi)){
    
    xi.chang  = tensor.est.xi(Y)
    
    if(Random.Project == T){
      res.only.used.rank = HDTTS.CP.est(Y = Y,xi = xi.chang, K = lag.k, Ratio.type = "log")
      
      r_hat_first =  res.only.used.rank$results$r.hat
      
      r_breve = 2*r_hat_first
      if(!is.null(Rank)){
        r_breve = 2*Rank
      }
      
      xi.res = Boostrap.xi.sel.v3(Y, r_breve = r_breve, eps = 0.1, lag.k = lag.k, Threshold = F, delta = NULL, Randomized.time = 50,
                                  BMP = T, Ratio.type = "log", augmented = F,grid_delta1 = 50, all.out = T, A = A,print.eps = T)
      xi = xi.res$xi.sel
    }else{
      
      xi = xi.chang
    }
    
    
  }
  
  if(is.null(A.inl)){
    res.inl =  HDTTS.CP.est(Y = Y,xi = xi, Rank = Rank, K = lag.k, Threshold = Threshold, delta = delta, Ratio.type = Ratio.type, augmented = augmented)
    
    A.hat    = res.inl$results$A.hat
    A.hat.NT = res.inl$results.NT$A.hat
    sigma0 =  res.inl$sigma0
  }else{
    res.inl = NULL
    A.hat    = A.inl
    A.hat.NT = A.inl
    sigma0 = sqrt(sum(Y^2)/(n*prod(D)))
    
  }
  
  vv = vv.NT = NULL
  
  if(Threshold == F){
    delta2 = rep(0,m)
  }
  
  if(is.null(delta2)){
    
    if(all.out == F){ #only delta2 output
      
      delta2_net =  matrix(rep(seq(0,delta_max,length.out = grid.num),m),grid.num,m)
      
      vv.mat = vv.mat.alter = vv.mat.NT = matrix(NA,grid.num,m)
      
      delta2_vv1 = delta2_vv2 = delta2 = delta2.alter = delta2.NT = rep(0,m)
      
      for(j in order(D,decreasing = F) ){
        
        vv1 = vv2 = vv = vv.alter = vv.NT = vector()
        for (kk in seq(0,delta_max,length.out = grid.num) ) {
          
          delta2[j] =  kk
          
          if(BMP == T){
            res.iter        =   CP.iter.BMP.xi(Y = Y, A.hat = A.hat, K = K , n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max,eps = eps, print.eps = F)
          }else{
            res.iter        =   CP.iter.BMP(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max,eps = eps, print.eps = F)
          }
          
          
          
          
          z_tol           =   res.iter$z_tol
          
          vv              =   rbind(vv  , z_tol$z2_tol) #### use z2_tol : trace
          
          # if(ev.all == T){
          #   vv1              =   rbind(vv1 , z_tol$z1_tol)
          #   vv2              =   rbind(vv2 , z_tol$z2_tol)
          # }
          #
          
          # cat("\r",paste(j,"-",round(kk,2)))
        }
        
        delta2[j]               =  delta2_net[which.max(vv[,j]),j]
        
        # if(ev.all == T){
        #   delta2_vv1[j]           =  delta2_net[which.max(vv1[,j]),j]
        #   delta2_vv2[j]           =  delta2_net[which.max(vv2[,j]),j]
        # }
        
        
      }
      
      
      if(BMP == T){
        res.iter.fin              =   CP.iter.BMP.xi(Y = Y,A.hat = A.hat, K=K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
        res.iter.NT.fin           =    res.iter.fin #CP.iter.BMP.xi(Y = Y,A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps)
      }else{
        res.iter.fin              =   CP.iter.AMP.xi(Y = Y,A.hat = A.hat, K=K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
        res.iter.NT.fin           =    res.iter.fin #CP.iter.BMP.xi(Y = Y,A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps)
      }
      
      
      # if(ev.all == T){
      #   res.iter.fin_vv1        =   CP.iter.BMP(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv1, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #   res.iter.NT.fin_vv1     =   CP.iter.BMP(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv1, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #
      #   res.iter.fin_vv2        =   CP.iter.BMP(A.hat = A.hat, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #   res.iter.NT.fin_vv2     =   CP.iter.BMP(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2_vv2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = F)
      #
      #   }
      
      
    }else{
      
      delta2_net =  matrix(rep(seq(0,delta_max,length.out = grid.num),m),grid.num,m)
      
      delta2  = delta2.NT = rep(0,m)
      
      for(j in order(D,decreasing = F) ){
        
        vv  = vv.NT = vector()
        for (kk in seq(0,delta_max,length.out = grid.num) ) {
          
          delta2[j] =  kk
          delta2.NT[j] =  kk
          
          if(BMP == T){
            res.iter.fin              =   CP.iter.BMP.xi(Y = Y, A.hat = A.hat, K=K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
            res.iter.NT.fin           =   CP.iter.BMP.xi(Y = Y, A.hat = A.hat.NT, K=K, n = n, delta2 = delta2.NT, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
          }else{
            res.iter.fin              =   CP.iter.AMP.xi(Y = Y, A.hat = A.hat, K=K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
            res.iter.NT.fin           =   CP.iter.AMP.xi(Y = Y, A.hat = A.hat.NT, K=K, n = n, delta2 = delta2.NT, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
          }
          
          vv              =  rbind(vv      , res.iter$z_tol$z2_tol) ###use z2
          vv.NT           =  rbind(vv.NT   , res.iter.NT$z_tol$z2_tol)
          
          # cat("\r",paste(j,"-",round(kk,2)))
        }
        
        delta2[j]       =  delta2_net[which.max(vv[,j]),j]
        delta2.NT[j]    =  delta2_net[which.max(vv.NT[,j]),j]
      }
      
      if(BMP == T){
        res.iter.fin              =   CP.iter.BMP.xi(Y = Y, A.hat = A.hat, K=K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
        res.iter.NT.fin           =   CP.iter.BMP.xi(Y = Y, A.hat = A.hat.NT, K=K, n = n, delta2 = delta2.NT, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
      }else{
        res.iter.fin              =   CP.iter.AMP.xi(Y = Y, A.hat = A.hat, K=K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
        res.iter.NT.fin           =   CP.iter.AMP.xi(Y = Y, A.hat = A.hat.NT, K=K, n = n, delta2 = delta2.NT, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,iter_lag = iter_lag)
      }
      
      
    }
    
  }else{
    
    if(BMP == T){
      res.iter.fin              =   CP.iter.BMP.xi(Y = Y, A.hat = A.hat, K = K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,A = A,iter_lag = iter_lag,Componet = Componet)
      res.iter.NT.fin           =    res.iter.fin #CP.iter.BMP.xi(Y = Y, A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,A = A)
    }else{
      res.iter.fin              =   CP.iter.AMP.xi(Y = Y, A.hat = A.hat, K = K, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps,A = A,iter_lag = iter_lag,Componet = Componet)
      res.iter.NT.fin           =    res.iter.fin #CP.iter.BMP(A.hat = A.hat.NT, Sigma_orginal = Sigma_orginal, n = n, delta2 = delta2, sigma0 = sigma0, augmented = augmented, iter_max = iter_max, eps = eps, print.eps = print.eps)
    }
    
  }
  
  if(!is.null(A)){
    if(BMP == T){
      rho.iter.fin            =  max( rho2.loss.list(res.iter.fin$A.hat,A) )
      rho.iter.NT.fin         =  max( rho2.loss.list(res.iter.NT.fin$A.hat,A) )
      
      rho.inl                 =  max( rho2.loss.list(A.hat,A) )
      rho.inl.NT              =  max( rho2.loss.list(A.hat.NT,A) )
    }else{
      rho.iter.fin            =  max( rho2.loss.list(res.iter.fin$A.hat,A) )
      rho.iter.NT.fin         =  max( rho2.loss.list(res.iter.NT.fin$A.hat,A) )
      
      rho.inl                 =  max( rho2.loss.list(A.hat,A) )
      rho.inl.NT              =  max( rho2.loss.list(A.hat.NT,A) )
    }
    
    
    mse_tuning_delta2 = c(rho.iter.fin,
                          rho.iter.NT.fin,
                          
                          rho.inl,
                          rho.inl.NT)
  }else{
    mse_tuning_delta2 = NULL
  }
  
  
  return(list(  res.iter   =  res.iter.fin,  res.iter.NT     =  res.iter.NT.fin, res.inl = res.inl, mse_tuning_delta2 = mse_tuning_delta2, vv = vv, vv.NT = vv.NT ))
  
}


R.j     = function(a,K12j){
  K12j%*%a - as.numeric(t(a)%*%K12j%*%a)*a
}

Theta.j.tilde = function(a,K12j,K12j_tilde,index){
  if(index == "0000"){
    return(K12j_tilde - as.numeric(t(a)%*%K12j_tilde%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j_tilde + t(K12j_tilde)))
  }
  if(index == "1000"){
    return(K12j  - as.numeric(t(a)%*%K12j_tilde%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j_tilde + t(K12j_tilde)))
  }
  if(index == "0100"){
    return(K12j_tilde - as.numeric(t(a)%*%K12j%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j_tilde + t(K12j_tilde)))
  }
  if(index == "0010"){
    return(K12j_tilde - as.numeric(t(a)%*%K12j_tilde%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j + t(K12j_tilde)))
  }
  if(index == "0001"){
    return(K12j_tilde - as.numeric(t(a)%*%K12j_tilde%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j_tilde + t(K12j)))
  }
  if(index == "1111"){
    return(K12j - as.numeric(t(a)%*%K12j%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j + t(K12j)))
  }
}

Theta.j = function(a,K12j){
  K12j - as.numeric(t(a)%*%K12j%*%a)*diag(1,dim(K12j)) - a%*%t(a)%*%(K12j + t(K12j))
}

K12j.breve = function(Sigma.Y.1,Sigma.Y.2,Sigma.Y.2.threshold,Qj){
  
  B12j = t(Qj)%*%t(Sigma.Y.2.threshold - Sigma.Y.2)%*%Sigma.Y.2.threshold%*%Qj
  
  K12j_breve = Sigma.Y.1%*%Qj%*%(MASS::ginv(t(Qj)%*%t(Sigma.Y.2.threshold)%*%Sigma.Y.2.threshold%*%Qj - B12j - t(B12j)))%*%t(Qj)%*%t(Sigma.Y.2)
  
  return(K12j_breve)
}




aij.debias.total = function(aij,A,j,xi,f,w,Sigma.Y.1,Sigma.Y.2,Sigma.Y.1.threshold,Sigma.Y.2.threshold,Qj){
  
  r = NCOL(Qj)
  
  G.1.xi = diag(w*Autocov_xi_Y(f,xi,1))
  G.2.xi = diag(w*Autocov_xi_Y(f,xi,2))
  
  Aj = A[[j]]
  
  K12j = Aj%*%G.1.xi%*%MASS::ginv(G.2.xi)%*%MASS::ginv(t(Aj)%*%Aj)%*%t(Aj)
  
  K12j_tilde = Sigma.Y.1.threshold%*%Qj%*%MASS::ginv(t(Qj)%*%Sigma.Y.2.threshold%*%Sigma.Y.2.threshold%*%Qj)%*%t(Qj)%*%t(Sigma.Y.2.threshold)
  
  K12j_breve = K12j.breve(Sigma.Y.1, Sigma.Y.2, Sigma.Y.2.threshold, Qj)
  
  #Theta0 R0
  Theta_0000     =  Theta.j.tilde(aij,K12j,K12j_tilde,"0000")
  Theta_1000     =  Theta.j.tilde(aij,K12j,K12j_tilde,"1000")
  Theta_0100     =  Theta.j.tilde(aij,K12j,K12j_tilde,"0100")
  Theta_0010     =  Theta.j.tilde(aij,K12j,K12j_tilde,"0010")
  Theta_0001     =  Theta.j.tilde(aij,K12j,K12j_tilde,"0001")
  Theta_1111     =  Theta.j.tilde(aij,K12j,K12j_tilde,"1111")
  
  R_breve     =  R.j(aij,K12j_breve)
  
  aij.0000 = aij - MASS::ginv(Theta_0000)%*%R_breve
  aij.1000 = aij - MASS::ginv(Theta_1000)%*%R_breve
  aij.0100 = aij - MASS::ginv(Theta_0100)%*%R_breve
  aij.0010 = aij - MASS::ginv(Theta_0010)%*%R_breve
  aij.0001 = aij - MASS::ginv(Theta_0001)%*%R_breve
  aij.1111 = aij - MASS::ginv(Theta_1111)%*%R_breve
  
  
  return(list(aij.0000 = aij.0000,
              aij.1000 = aij.1000,
              aij.0100 = aij.0100,
              aij.0010 = aij.0010,
              aij.0001 = aij.0001,
              aij.1111 = aij.1111,
              Theta_0000 = Theta_0000,
              Theta_1000 = Theta_1000,
              Theta_0100 = Theta_0100,
              Theta_0010 = Theta_0010,
              Theta_0001 = Theta_0001,
              Theta_1111 = Theta_1111,
              K12j = K12j,K12j_tilde=K12j_tilde, K12j_breve = K12j_breve ))
}


aij.debias.K2check = function(aij,Sigma.Y.1,Sigma.Y.2,Sigma.Y.1.threshold,Sigma.Y.2.threshold,Qj){
  
  r = NCOL(Qj)
  
  K12j_breve = K12j.breve(Sigma.Y.1, Sigma.Y.2, Sigma.Y.2.threshold, Qj)
  
  K12j_tilde = Sigma.Y.1.threshold%*%Qj%*%MASS::ginv(t(Qj)%*%t(Sigma.Y.2.threshold)%*%Sigma.Y.2.threshold%*%Qj)%*%t(Qj)%*%t(Sigma.Y.2.threshold)
  
  #Theta0 R0
  
  Theta_0100     =  Theta.j.tilde(aij,K12j_breve,K12j_tilde,"0100")
  
  R_breve        =  R.j(aij,K12j_breve)
  
  aij.0100 = aij - MASS::ginv(Theta_0100)%*%R_breve
  
  
  return(list(aij.0100 = aij.0100,
              Theta_0100 = Theta_0100,
              K12j_tilde=K12j_tilde,
              K12j_breve = K12j_breve ))
}


 
aij.debias.iter = function(aij,i,j,Y,A){
  m = length(A)
  r = NCOL(A[[1]])
  
  Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  
  
  if(m >= 3){
    B.tilde.j = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    B.tilde.j = A[c(1:m)[-j]][[1]]
  }
  
  A.tilde.j  = A[[j]]
  
  B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
  A.MP.j = A.tilde.j%*%MASS::ginv(t(A.tilde.j)%*%A.tilde.j)
  
  
  if(m >= 3){
    dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m+1):3])) # combine dj j \neq i into one mode, then Yp be a n x di x prod(dj, j \neq i)  3 dimensional tensor
  }
  
  Uj = rTensor::ttl(as.tensor(Yp),list_mat = list(t(A.MP.j),t(B.MP.j)), ms = c(2,3) )@data
  
  Uj.reg = t(apply(Uj, 1, diag))
  
  Uj.reg = as.matrix(scale(Uj.reg))
  
  ujk = lm(Uj.reg[-n,i] ~ Uj.reg[-1,-i]-1)$residuals #Uj.reg[,k]
  
  sigma.k.ujk = tensor.Autocov_xi_Y(Y,scale(c(ujk,0)),1)
  
  A.breve.jk = Mat.tensor(sigma.k.ujk,j)%*%B.MP.j[,i]
  
  aij.de  = A.breve.jk/as.numeric(t(aij)%*%A.breve.jk)
  
  return(list(aij.de = aij.de,A.breve.jk=A.breve.jk))
}


aij.debias.iter.v2 = function(aij,i,j,Sigma.yij.xii.1, A){
  m = length(A)
  r = NCOL(A[[1]])
  
  vartheta_ij = (as.numeric(t(aij)%*% Sigma.yij.xii.1[[j]][,i])*aij - Sigma.yij.xii.1[[j]][,i])/as.numeric(t(aij)%*% Sigma.yij.xii.1[[j]][,i])
  
  
  aij.de = aij - vartheta_ij
  
  return(list(aij.de = aij.de, vartheta_ij = vartheta_ij))
}


cov.aij.debias = function(h,A,j,aij,xi,C,f,COV.VECE,w,thres = 0.01){
  
  m = length(A)
  n = length(xi)
  
  Sigma.Y.xi.1.tensor = tensor.Autocov_xi_Y(C,xi,1)
  Sigma.Y.xi.2.tensor = tensor.Autocov_xi_Y(C,xi,2)
  
  Sigma.Y.xi.1  = Mat.tensor(Sigma.Y.xi.1.tensor,j)
  Sigma.Y.xi.2  = Mat.tensor(Sigma.Y.xi.2.tensor,j)
  
  G.1.xi = diag(w*Autocov_xi_Y(f,xi,1))
  G.2.xi = diag(w*Autocov_xi_Y(f,xi,2))
  
  Aj = A[[j]]
  
  i = which(1 - abs(t(aij)%*%Aj) < 10^-8)
  
  dj = NROW(Aj)
  
  if(m >= 3){
    Bj = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    Bj = A[c(1:m)[-j]][[1]]
  }
  
  Bj.MP = Bj%*%MASS::ginv(t(Bj)%*%Bj)
  Aj.MP = Aj%*%MASS::ginv(t(Aj)%*%Aj)
  
  aij.MP = Aj.MP[,i]
  bij.MP = Bj.MP[,i]
  
  K12j = Aj%*%G.1.xi%*%MASS::ginv(G.2.xi)%*%MASS::ginv(t(Aj)%*%Aj)%*%t(Aj)
  Thetaj = Theta.j(aij,K12j)
  
  G = t(h)%*%svd_inverse(Thetaj,thres)%*%(diag(dj) - aij%*%t(aij))
  
  L1 = t(bij.MP) %x% K12j
  
  L2 = t(bij.MP) %x% diag(dj)
  
  
  VAR1 =  sum(xi[1:(n-1)]^2)/(n-1)
  VAR2 =  sum(xi[1:(n-2)]^2)/(n-2)
  COV12 = sum(xi[1:(n-2)]*xi[2:(n-1)])/(n-2)
  
  cov.org = G%*%(VAR1 * L1%*%COV.VECE%*%t(L1) + VAR2 * L2%*%COV.VECE%*%t(L2) - COV12*L1%*%COV.VECE%*%t(L2) - COV12*L2%*%COV.VECE%*%t(L1) )%*%t(G)
  
  COV.TOL = (G.2.xi[i,i])^(-2)  * cov.org
  
  print(cov.org)
  print(G.2.xi[i,i])
  
  # beta1 = bij.MP %x% t( t(h)%*% svd_inverse(Theta.j(aij,K12j),thres)%*%(diag(dj) - aij%*%t(aij))%*%K12j)
  #
  # beta2 = bij.MP %x% t( t(h)%*% svd_inverse(Theta.j(aij,K12j),thres)%*%(diag(dj) - aij%*%t(aij)))
  #
  # Wj = rTensor::khatri_rao(Bj,Aj)
  # #
  # Ept <- aperm(E, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  # Ejt = Mat.tensor(Ept,1)
  #
  # q1j = Ejt%*%beta1
  # q2j = Ejt%*%beta2
  #
  # # Ypt <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  # # Yjt = Mat.tensor(Ypt,1)
  # # q1jy = Yjt%*%beta1
  # # q2jy = Yjt%*%beta2
  #
  # g_2ixi_hat = as.numeric(t(aij.MP)%*%SY2Tj%*%bij.MP)
  #
  # cov.vec = mean((xi[1:(n-2)] - mean(xi))^2)*((q1j[2:(n-1)])^2 + (q2j[2:(n-1)])^2)  - 2 * mean((xi[2:(n-1)] - mean(xi))*(xi[1:(n-2)] - mean(xi)))*((q1j*q2j)[3:(n)])
  #
  # print(mean(cov.vec))
  
  return(as.numeric(COV.TOL))
}


cov.aij.debias.est = function(h,place.inl,j,Y,xi,res_CP,aij.tilde.other = NULL,iter = F, A_fit_C = NULL,thres = 0.1){
  
  if(iter == T){
    AA = res_CP$res.iter$A.hat
  }else{
    AA = res_CP$res.inl$results$A.hat
  }
  
  con.inl = res_CP$res.inl
  
  Q_list  = con.inl$results$Q_list
  
  SY1.tensor  = con.inl$Sigma.tensor.Y.k_list_orginal[[1]]
  SY2.tensor  = con.inl$Sigma.tensor.Y.k_list_orginal[[2]]
  
  SY1T.tensor = con.inl$results$Sigma.tensor.Y.k_list_threshold[[1]]
  SY2T.tensor = con.inl$results$Sigma.tensor.Y.k_list_threshold[[2]]
  
  SY1j  = Mat.tensor(SY1.tensor,j)
  SY2j  = Mat.tensor(SY2.tensor,j)
  
  SY1Tj  = Mat.tensor(SY1T.tensor,j)
  SY2Tj  = Mat.tensor(SY2T.tensor,j)
  
  Qj =  Q_list[[j]]
  
  m = length(AA)
  n = length(xi)
  
  Aj.tilde = AA[[j]]
  
  dj = NROW(Aj.tilde)
  
  aij.tilde = AA[[j]][,place.inl]
  
  
  res_debias_breve = aij.debias.K2check(aij.tilde, Sigma.Y.1 = SY1j , Sigma.Y.2 = SY2j, Sigma.Y.1.threshold = SY1Tj, Sigma.Y.2.threshold = SY2Tj, Qj = Qj)
  
  aij_debias_breve = c(res_debias_breve$aij.0100)
  
  Thetaj.tilde = res_debias_breve$Theta_0100
  
  K12j.tilde = res_debias_breve$K12j_tilde
  K12j_breve = res_debias_breve$K12j_breve
  
  if(m >= 3){
    Bj.tilde = rTensor::khatri_rao_list(AA[c(m:1)[-(m-j+1)]])
  }else{
    Bj.tilde = AA[c(1:m)[-j]][[1]]
  }
  
  Bj.MP.tilde = Bj.tilde%*%MASS::ginv(t(Bj.tilde)%*%Bj.tilde)
  Aj.MP.tilde = Aj.tilde%*%MASS::ginv(t(Aj.tilde)%*%Aj.tilde)
  
  bij.MP.tilde = Bj.MP.tilde[,place.inl]
  aij.MP.tilde = Aj.MP.tilde[,place.inl]
  
  
  beta1_hat = bij.MP.tilde %x% t(t(h)%*% svd_inverse(Thetaj.tilde,thres)%*%(diag(dj) - aij.tilde%*%t(aij.tilde))%*%K12j.tilde)
  
  beta2_hat = bij.MP.tilde %x% t(t(h)%*% svd_inverse(Thetaj.tilde,thres)%*%(diag(dj) - aij.tilde%*%t(aij.tilde)))
  
  if(!is.null(A_fit_C)){
    Aj_fit_C = A_fit_C[[j]]
    if(m >= 3){
      Bj_fit_C = rTensor::khatri_rao_list(A_fit_C[c(m:1)[-(m-j+1)]])
    }else{
      Bj_fit_C = A_fit_C[c(1:m)[-j]][[1]]
    }
    Wj.tilde = rTensor::khatri_rao(Bj_fit_C,Aj_fit_C)
    Wj.MP.tilde = Wj.tilde%*%MASS::ginv(t(Wj.tilde)%*%Wj.tilde)
  }else{
    Wj.tilde = rTensor::khatri_rao(Bj.tilde,Aj.tilde)
    Wj.MP.tilde = Wj.tilde%*%MASS::ginv(t(Wj.tilde)%*%Wj.tilde)
  }
  
  Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  
  Yj = Mat.tensor(Yp,1)
  
  #Yj = Yj - rep(1,n)%*%t(colMeans(Yj))
  
  Ej.tilde = Yj - Yj%*%Wj.MP.tilde%*%t(Wj.tilde)
  
  q1j_hat = Ej.tilde%*%beta1_hat
  q2j_hat = Ej.tilde%*%beta2_hat
  
  g_2ixi_hat = as.numeric(t(aij.MP.tilde)%*%SY2Tj%*%bij.MP.tilde)
  
  cov.vec = mean((xi[1:(n-2)] - mean(xi))^2)*((q1j_hat[2:(n-1)])^2 + (q2j_hat[2:(n-1)])^2)  - 2 * mean((xi[2:(n-1)] - mean(xi))*(xi[1:(n-2)] - mean(xi)))*((q1j_hat*q2j_hat)[3:(n)])
  
  COV.TOL = mean(cov.vec)/abs(g_2ixi_hat)^2
  
  print(mean(cov.vec))
  print(g_2ixi_hat)
  
  tau_hat = mean(cov.vec)
  
  return(list(var.est = COV.TOL, tau_hat = tau_hat, q1j_hat = q1j_hat, q2j_hat = q2j_hat, g_2ixi_hat = g_2ixi_hat)  )
}

cov.aij.debias.iter = function(h,A,j,aij,f,COV.VECE,w){
  
  m = length(A)
  n = NROW(f)
  dj = length(aij)
  
  Aj = A[[j]]
  
  i = which(1 - abs(t(aij)%*%Aj) < 10^-8)
  
  if(m >= 3){
    B.tilde.j = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    B.tilde.j = A[c(1:m)[-j]][[1]]
  }
  
  B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
  
  Fjk = lm(f[-n,i] ~ f[-1,-i]-1)$residuals
  
  sigma.fi.xi = Autocov_xi_Y(f,c(Fjk,0),1)[i]/sqrt(var(Fjk))
  
  G = t(h)%*%(diag(dj) - aij%*%t(aij))%*%(t(B.MP.j[,i]) %x% diag(D[j]))
  
  COV.TOL = G%*%COV.VECE%*%t(G)/sigma.fi.xi^2/(w[i])^2
  
  return(as.numeric(COV.TOL))
}

cov.aij.debias.iter.est = function(h,A,j,aij,f,Y){
  
  m  = length(A)
  n  = NROW(f)
  dj = length(aij)
  
  Aj = A[[j]]
  
  i = which(1 - abs(t(aij)%*%Aj) < 10^-8)
  
  if(m >= 3){
    Bj = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    Bj = A[c(1:m)[-j]][[1]]
  }
  
  B.MP.j = Bj%*%MASS::ginv(t(Bj)%*%Bj)
  
  ff = scale(f)
  
  xi.vmax.j = lm(ff[-n,i] ~ ff[-1,-i]-1)$residuals
  
  sigma.fi.xi.wi = Autocov_xi_Y(f,c(xi.vmax.j,0),1)[i]
  
  beta = (B.MP.j[,i]) %x% t(t(h)%*%(diag(dj) - aij%*%t(aij)))
  
  Wj = rTensor::khatri_rao(Bj,Aj)
  
  Wj.MP = Wj%*%MASS::ginv(t(Wj)%*%Wj)
  
  Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  
  Yj = Mat.tensor(Y,1)
  
  Ej = Yj - Yj%*%Wj.MP%*%t(Wj)
  
  qj  = Ej%*%beta
  qj2 = Yj%*%beta
  
  cov.tol = mean(xi.vmax.j^2*(qj[2:n])^2)
  
  COV.TOL = cov.tol/sigma.fi.xi.wi^2
  
  return(as.numeric(COV.TOL))
}

Boostrap.sigma.fi.xi = function(i,j,rep,seed_num,n, m, D, r, w, ar.coef, factor.loading, factor.corr, alpha, delta){
  sig_tol = vector()
  set.seed(seed_num)
  for (ss in 1:rep) {
    data  = DGP.TCP(n=n, m=m, D=D, r=r, w=w, ar.coef = ar.coef, factor.loading = factor.loading, factor.corr = factor.corr, alpha = ALPHA, delta = delta.fl, par_E = 1)
    
    f  = data$f
    
    Fjk = lm(f[-n,i] ~ f[-1,-i]-1)$residuals
    
    sigma.fi.xi = Autocov_xi_Y(f,c(Fjk,0),1)[i]/sqrt(var(Fjk))
    
    sig_tol = c(sig_tol,sigma.fi.xi)
  }
  
  return(mean(sig_tol))
  
  
}

cov.aij.debias.iter.f = function(h,A,j,aij,f,COV.VECE,sigma.fi.xi,w){
  
  m = length(A)
  n = NROW(f)
  dj = length(aij)
  
  Aj = A[[j]]
  
  i = which(1 - abs(t(aij)%*%Aj) < 10^-8)
  
  if(m >= 3){
    B.tilde.j = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    B.tilde.j = A[c(1:m)[-j]][[1]]
  }
  
  B.MP.j = B.tilde.j%*%MASS::ginv(t(B.tilde.j)%*%B.tilde.j)
  
  
  G = t(h)%*%(diag(dj) - aij%*%t(aij))%*%(t(B.MP.j[,i]) %x% diag(D[j]))
  
  COV.TOL = G%*%COV.VECE%*%t(G)/sigma.fi.xi^2/(w[i])^2
  
  return(as.numeric(COV.TOL))
}

cov.aij.debias.f = function(h,A,j,aij,n,COV.VECE,J,beta,w){
  
  
  m = length(A)
  r = NCOL(A[[1]])
  
  VAR.org = as.numeric(t(rep(1/r,r)) %*% Mat.k(J,0.5) %*% diag(1/(1-beta^2)) %*% Mat.k(J,0.5) %*% rep(1/r,r))
  
  G.1.xi = VAR.org^(-1/2)*diag(c( Mat.k(J,0.5) %*% diag(w*beta/(1-beta^2)) %*% Mat.k(J,0.5) %*% rep(1/r,r)) )
  G.2.xi = VAR.org^(-1/2)*diag(c( Mat.k(J,0.5) %*% diag(w*beta/(1-beta^2)) %*% diag(beta) %*% Mat.k(J,0.5) %*% rep(1/r,r) ))
  
  Aj = A[[j]]
  
  i = which(1 - abs(t(aij)%*%Aj) < 10^-8)
  
  dj = NROW(Aj)
  
  if(m >= 3){
    Bj = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    Bj = A[c(1:m)[-j]][[1]]
  }
  
  Bj.MP = Bj%*%MASS::ginv(t(Bj)%*%Bj)
  
  bij.MP = Bj.MP[,i]
  
  K12j = Aj%*%G.1.xi%*%MASS::ginv(G.2.xi)%*%MASS::ginv(t(Aj)%*%Aj)%*%t(Aj)
  
  G = t(h)%*%MASS::ginv(Theta.j(aij,K12j))%*%(diag(dj) - aij%*%t(aij))
  
  L1 = t(bij.MP) %x% K12j
  
  L2 = t(bij.MP) %x% diag(dj)
  
  
  VAR1 = 1 #as.numeric(t(rep(1/r,r)) %*% Mat.k(J,0.5) %*% diag(1/(1-beta^2)) %*% Mat.k(J,0.5) %*% rep(1/r,r))
  
  VAR2 = 1 #as.numeric(t(rep(1/r,r)) %*% Mat.k(J,0.5) %*% diag(1/(1-beta^2)) %*% Mat.k(J,0.5) %*% rep(1/r,r))
  
  VAR.org = as.numeric(t(rep(1/r,r)) %*% Mat.k(J,0.5) %*% diag(1/(1-beta^2)) %*% Mat.k(J,0.5) %*% rep(1/r,r))
  
  
  COV12 =  VAR.org^(-1) * as.numeric(t(rep(1/r,r)) %*% Mat.k(J,0.5) %*% diag(beta/(1-beta^2)) %*% Mat.k(J,0.5) %*% rep(1/r,r))
  
  
  COV.TOL = (G.2.xi[i,i])^(-2) * G%*%(VAR1 * L1%*%COV.VECE%*%t(L1) + VAR2 * L2%*%COV.VECE%*%t(L2) -  COV12*L2%*%COV.VECE%*%t(L1) - COV12*L1%*%COV.VECE%*%t(L2) )%*%t(G) #
  
  return(as.numeric(COV.TOL))
}

forecast.factor.tensor = function(f,lag.max = 6,sp = 2,diff = FALSE){
  if(diff == F){
    if(NCOL(f) > 1){
      # modeling ft with VAR
      colnames(f) <- paste0("y",1:NCOL(f))
      var_f = try(vars::VAR(f, type = "const", lag.max = lag.max, ic = "HQ"),silent = T)
      
      if(class(var_f) == "try-error"){ #sometimes when ft are almostly perfect correlated, VAR report error, we need do AR for each ft
        f_tenstep = vector()
        colnames(f) <- paste0("y",1:NCOL(f))
        
        for (jj in 1:NCOL(f)) {
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
          
          for (jj in 1:NCOL(f)) {
            ar_f  = ar(f[,jj])
            pre_f = predict(ar_f,n.ahead = 10)
            f_tenstep_i = as.matrix(pre_f$pred)
            f_tenstep = cbind(f_tenstep,f_tenstep_i)
          }
          
          colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
          row.names(f_tenstep) <- paste(1:10,"step")
        }
      }
    }else{
      ar_f  = ar(f)
      pre_f = predict(ar_f,n.ahead = 10)
      f_tenstep = as.matrix(pre_f$pred)
      colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
      row.names(f_tenstep) <- paste(1:10,"step")
      f_tenstep = as.matrix(f_tenstep)
    }
  }else{
    
    f_tenstep = vector()
    colnames(f) <- paste0("y",1:NCOL(f))
    
    for (jj in 1:NCOL(f)) {
      fit <- forecast::auto.arima(
        as.matrix(f)[,jj],
        d = 1,
        seasonal = TRUE,
        D = 1,
        max.P = 1, max.Q = 1,
        max.p = lag.max, max.q = lag.max,
        stepwise = F,
        approximation = F,
        ic = "aic"
      )
      fc <- forecast::forecast(fit, h = 10)
      f_tenstep_i = as.matrix(fc$mean)
      f_tenstep = cbind(f_tenstep,f_tenstep_i)
    }
    
    colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
    row.names(f_tenstep) <- paste(1:10,"step")
    
  }
  
  return(f_forecast = f_tenstep[1:sp,])
}

forecast.factor.tensor.seasonal = function(f,lag.max = 6,sp = 2,seasonal = F,s = NULL,diff = F){
  
  if(NCOL(f) > 1){
    # modeling ft with VAR
    colnames(f) <- paste0("y",1:NCOL(f))
    var_f = try(vars::VAR(f, type = "const", lag.max = lag.max, ic = "HQ"),silent = T)
    
    if(class(var_f) == "try-error"){ #sometimes when ft are almostly perfect correlated, VAR report error, we need do AR for each ft
      f_tenstep = vector()
      colnames(f) <- paste0("y",1:NCOL(f))
      
      for (jj in 1:NCOL(f)) {
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
        
        for (jj in 1:NCOL(f)) {
          ar_f  = ar(f[,jj])
          pre_f = predict(ar_f,n.ahead = 10)
          f_tenstep_i = as.matrix(pre_f$pred)
          f_tenstep = cbind(f_tenstep,f_tenstep_i)
        }
        
        colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
        row.names(f_tenstep) <- paste(1:10,"step")
      }
    }
  }else{
    if(seasonal == T){
      s <- s
      if(diff == F){
        dd = 0
      }else{
        dd = 1
      }
      fit <- forecast::auto.arima(
        f,
        d = dd,
        seasonal = TRUE,
        D = 1,
        max.P = 1, max.Q = 1,
        stepwise = F,
        approximation = F,
        ic = "aic"
      )
      fc <- forecast::forecast(fit, h = 10)
      f_tenstep = as.matrix(fc$mean)
      colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
      row.names(f_tenstep) <- paste(1:10,"step")
      f_tenstep = as.matrix(f_tenstep)
    }else{
      ar_f  = ar(f)
      pre_f = predict(ar_f,n.ahead = 10)
      f_tenstep = as.matrix(pre_f$pred)
      colnames(f_tenstep)  <- paste0("f_pre",1:NCOL(f))
      row.names(f_tenstep) <- paste(1:10,"step")
      f_tenstep = as.matrix(f_tenstep)
    }
    
  }
  
  
  return(f_forecast = f_tenstep[1:sp,])
}



Pro.forecast =  function(Y, forecast.step = 2, xi = NULL, Rank = NULL, lag.k = 10, Threshold = F, delta = NULL, delta2 = NULL, A.inl =  NULL, iter_lag = 1, seasonal = F, diff = F, Random.Project = F,
                         BMP = T, Ratio.type = "log", augmented = F, iter_max = 100, eps = 10^-5, grid_delta1 = 50, grid.num = 50, delta_max = 0.1, all.out = T, A = NULL,print.eps = F){
  
  con.BMP = HDTTS.CP.iter.BMP(Y,  xi = xi, Rank = Rank,lag.k = lag.k , Threshold = Threshold, delta = delta, delta2 = delta2, Random.Project = Random.Project,
                              Ratio.type = Ratio.type, augmented = augmented, BMP = T,iter_max = iter_max, eps = eps, grid.num = grid.num, delta_max = delta_max,all.out = T, ev.all = T,A = A,A.inl = A.inl,print.eps=print.eps,iter_lag=iter_lag)
  
  
  A.inl = con.BMP$res.iter$A.hat
  A.iter = con.BMP$res.iter$A.inl
  
  f.inl = con.BMP$res.iter$f_hat_inl
  f.iter = con.BMP$res.iter$f_hat
  
  
  f.fore.inl  =  as.matrix(forecast.factor.tensor.seasonal(f.inl,lag.max = 6,sp = forecast.step,seasonal = seasonal,s = 365,diff =diff))
  f.fore.iter =  as.matrix(forecast.factor.tensor.seasonal(f.iter,lag.max = 6,sp = forecast.step,seasonal = seasonal,s = 365,diff =diff))
  
  Y.forecast.inl = Y.forecast.iter  = list()
  
  for (gg in 1:forecast.step) {
    Y.forecast.iter.i = Y.forecast.inl.i = 0
    for (i in 1:NCOL(f.inl)) {
      
      Coef.inl = A.inl[[1]][,i]
      Coef.iter = A.iter[[1]][,i]
      
      for (j in 2:length(A.inl)) {
        Coef.inl = Coef.inl %o% A.inl[[j]][,i]
        Coef.iter = Coef.iter %o% A.iter[[j]][,i]
      }
      Y.forecast.inl.i = Y.forecast.inl.i + f.fore.inl[gg,i]*Coef.inl
      Y.forecast.iter.i = Y.forecast.iter.i + f.fore.iter[gg,i]*Coef.iter
    }
    
    Y.forecast.inl[[gg]]  = Y.forecast.inl.i
    Y.forecast.iter[[gg]] = Y.forecast.iter.i
    
  }
  
  return(list(Y.forecast.inl=Y.forecast.inl,Y.forecast.iter=Y.forecast.iter,A.inl=A.inl,A.iter=A.iter,f.inl=f.inl,f.iter=f.iter,con.BMP=con.BMP))
  
}#

PCATS.forecast = function(Y,permatation = "max",forecast.step = c(1,2,6),...){
  p = dim(Y)[2];q = dim(Y)[3];
  Y_m =  Vec.tensor(Y)
  
  pcats = HDTSA::PCA_TS(Y_m,permutation = permatation,beta = 0.1)
  
  Z     = pcats$Z
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
  
  Y_vect = Mat.tensor(Y,1)
  
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


UniARMA.forecast.tensor = function(Y,forecast.step = c(1,2)){
  
  D = dim(Y)[-1]
  
  Y_vect = Mat.tensor(Y,1)
  
  Y_tenstep = apply(Y_vect, 2, function(x){
    res = ar(x)
    arma_fore = as.numeric(predict(res,n.ahead = 10)$pred)
    return(arma_fore)
  })
  
  
  Y.forecast = list()
  for (gg in forecast.step) {
    Y.forecast[[gg]] <- array(Y_tenstep[gg,], dim = D)
  }
  
  return(list(Y.forecast = Y.forecast))
}



UniARMA.forecast.tensor.diff = function(Y,forecast.step = c(1,2)){
  
  D = dim(Y)[-1]
  
  Y_vect = Mat.tensor(Y,1)
  
  Y_tenstep = apply(Y_vect, 2, function(x){
    fit <- forecast::auto.arima(
      x,
      d = 1,
      seasonal = TRUE,
      D = 1,
      max.P = 1, max.Q = 1,
      max.p = 6, max.q = 6,
      stepwise = F,
      approximation = F,
      ic = "aic"
    )
    fc <- forecast::forecast(fit, h = 10)
    arma_fore = as.numeric(fc$mean)
    return(arma_fore)
  })
  
  
  Y.forecast = list()
  for (gg in forecast.step) {
    Y.forecast[[gg]] <- array(Y_tenstep[gg,], dim = D)
  }
  
  return(list(Y.forecast = Y.forecast))
}

PCATS.forecast.tensor = function(Y,permatation = "max",forecast.step = c(1,2),...){
  D = dim(Y)[-1]
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
    Y.forecast[[jj]] <- array(Y_fore[,jj],dim = D)
  }
  
  return(list(Y.forecast = Y.forecast, Z = Z, B = B))
  
}
