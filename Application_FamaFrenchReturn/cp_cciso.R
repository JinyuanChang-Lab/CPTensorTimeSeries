### cp.init: Algorithm 1 without RP, c-pca 
### cp.init.rppca: Algorithm 1 with RP, rc-pca 
### cp.iso: Algorithm 2, CC-ISO

library(tensor)
library(rTensor)

#--------------------CC-ISO Algorithm-------------------------#
# initialization (CPCA) for any Dim tensor time series with CP structure 
cp.init <- function(x,r,r_est = F){
  # x: tensor of any dimension
  # CP initialization
  # x: d1 * d2 * d3 * ... * dk * n
  # use 2-th lag h cross moment
  
  dd <- dim(x)
  d <- length(dd) # d >= 2, d=k+1
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1
  
  x.matrix <- matrix(x,ncol=n)
  M.temp <- x.matrix %*% t(x.matrix) / n 
  
  ans.eig <- svd(x.matrix)
  #ans.eig <- eigen(M.temp)
  ans.lambda <- (ans.eig$d)^2
  #ans.lambda <- ans.eig$values
  if (r_est == T){
    rmax <- n
    # estimate the number of factors using eigenvalue-ratio method 
    # temp1 <- tail(ans.lambda, 1) 
    temp1 <- ans.lambda[1:(n-1)]
    temp2 <- ans.lambda[2:n] 
    r <- min(which.max(temp1 / temp2))
  } 
  #ans.Q2th <- ans.eig$vectors[,1:r^2,drop=FALSE]
  ans.Q2th <- ans.eig$u[,1:r^2,drop=FALSE]
  
  #ans.lambda[1:r^3]
  ans.Q0 <- matrix(0,sum(dd[-d]),r) 
  # ans.Q0 collects a_ik in a matrix, i.e. a_i1 is the top d1 rows and etc. 
  for(i in 1:r){
    ans.Q2 <- array(ans.Q2th[,i],dd[-d])
    for(l in 1:k){
      m.temp <- tensor(ans.Q2,ans.Q2,c(1:k)[-l],c(1:k)[-l])
      ans.eig.mtemp <- eigen(m.temp)
      ans.Q0[row(ans.Q0)>(sum(dd[1:(l-1)])*(l>1)) & row(ans.Q0)<=sum(dd[1:l]) & col(ans.Q0)==i] <- 
        ans.eig.mtemp$vectors[,1,drop=FALSE]
    }
  }
  ans.Qinit <- NULL
  # ans.Qinit takes the loading matrix for each dimension out so length = k, for each k, dim = dk X r 
  for(i in 1:k){
    qtemp <- matrix(ans.Q0[row(ans.Q0)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0)<=sum(dd[1:i])],ncol=r)
    ans.Qinit <- c(ans.Qinit,list(qtemp))
  }
  
  list("Q"=ans.Qinit,"Qmatrix"=ans.Q0,"lambda"=ans.lambda, "r"=r)
}


non_distinct_eigvals_index_func <- function(eigvals, r, c0){ 

  lambda_r <- eigvals[r] 
  temp1 <- eigvals[-r]
  temp2 <- eigvals[-1]

  set1 <- sort(which( (temp1 - temp2) < c0 * lambda_r))

  if (length(set1) == 0){
    set_distinct <- 1:r 
    I_list <- NULL  
  } else {
    I_list <- split(set1, cumsum(c(1, diff(set1)) != 1))
    
    for (i in seq_along(I_list)){
      I_list[[i]] <- union(I_list[[i]], I_list[[i]] + 1)
    } 
    set <- union(set1, set1 + 1) 
    set_distinct <- setdiff(1:r, set)
  } 

  list("distinct_index" = set_distinct, "non_distinct_index_list" = I_list)
}

# Same as cp.init, but with detection of close eigenvalues and randomized projection if close eigenvaluees are detected. 
cp.init.rppca <- function(x,r,c0, L = 100, v = 0.5, 
                           distinct_set = NULL, non_distinct_set = NULL){
  
  dd <- dim(x)
  d <- length(dd) # d >= 2, d=k+1
  d1 <- dd[1]
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1
  
  x.matrix <- matrix(x,ncol=n)
  M.temp <- x.matrix %*% t(x.matrix) / n 
  
  # ans.eig <- eigen(M.temp)
  ans.eig <- svd(M.temp)
  ans.lambda <- (ans.eig$d)^2
  
  # test the distinct eigenvalue conditions 
  if (is.null(distinct_set) && is.null(non_distinct_set)){
    index <- non_distinct_eigvals_index_func(ans.lambda[1:r], r, c0)
    distinct_set <- index$distinct_index 
    non_distinct_set <- index$non_distinct_index_list 
  }

  N <- length(non_distinct_set)
  
  ans.Q2th <- ans.eig$u[,1:r^2,drop=FALSE]
  
  #ans.lambda[1:r^3]
  ans.Q0 <- matrix(0,sum(dd[-d]),r) 
  if (length(distinct_set) > 0){ 
    for(i in distinct_set){
      ans.Q2 <- array(ans.Q2th[,i],dd[-d])
      for(l in 1:k){
        m.temp <- tensor(ans.Q2,ans.Q2,c(1:k)[-l],c(1:k)[-l])
        ans.eig.mtemp <- eigen(m.temp)
        ans.Q0[row(ans.Q0)>(sum(dd[1:(l-1)])*(l>1)) & row(ans.Q0)<=sum(dd[1:l]) & col(ans.Q0)==i] <- 
          ans.eig.mtemp$vectors[,1,drop=FALSE]
      }
    }
  } 
  
  #------------------------------------------------------------------------#
  # Algorithm 2: Randomized Projection. 
  if (!is.null(non_distinct_set)){
    print('close eigenvalues detected, using randomized projection...')
    for (j in 1:N){
      Ij <- non_distinct_set[[j]] 
      eigvals_j <- ans.lambda[Ij]
      eigvecs_j <- ans.Q2th[, Ij] 

      Sigma_j <- array( eigvecs_j %*% diag(eigvals_j) %*% t(eigvecs_j), c(dd[-d], dd[-d])) # d1 X ... X dk X d1 X ... X dk

      SL <- NULL 

      for (l in 1:L){ 
        theta <- matrix(rnorm(d1), d1, 1)
        temp <- aperm(tensor(Sigma_j, theta, 1,1), length(dim(Sigma_j)):1) 
        temp <- matrix(tensor(temp, theta, d,1), dd.prod/d1, dd.prod/d1) 

        eig.temp <- eigen(temp) 
        etal <- eig.temp$values[1] 
        ul <- array(eig.temp$vectors[,1], c(1,dd[-c(1,d)])) # 1 X d2 X ... X dk 

        a_list <- NULL 
        for (ik in 2:k){
          temp_k <- tensor(ul, ul, c(1:k)[-ik], c(1:k)[-ik]) 
          eig.temp_k <- eigen(temp_k) 
          a_list <- c(a_list, list(matrix(Re(eig.temp_k$vectors[,1]), nrow = 1)))
        }

        temp2 <- matrix(ttl(as.tensor(Sigma_j), c(a_list,a_list), c(2:k, (k+2):(2*k)))@data, d1, d1) 
        a_1 <- matrix(Re(eigen(temp2 %*% t(temp2))$vectors[,1]), nrow = 1) 
        a_list <- c(list(a_1), a_list) 
        SL <- c(SL, list(a_list)) 
      }

      for (i in Ij){
        norm_SL <- NULL 
        for (l in 1:L){ 
          norm_l <- rTensor::fnorm(ttl(as.tensor(Sigma_j), c(SL[[l]], SL[[l]]), c(1:k, (k+1):(2*k))))
          norm_SL <- append(norm_SL, norm_l) 
        } 
        ai_list <- SL[[which.max(norm_SL)]] 
        ### remove the elment in SL that are too close to ai_list 
        norm_SL_2 <- NULL 
        for (l in 1:L){ 
          norm_l_2 <- 0 
          for (ik in 2:k){
            for (jk in 2:k){
              norm_l_2 <- max(norm_l_2, abs(ai_list[[ik]] %*% t(SL[[l]][[jk]])))
            }
          }
          norm_SL_2 <- append(norm_SL_2, norm_l_2)
        }
        remove_idx <- which(norm_SL_2 > v) 
        if (length(remove_idx) != 0){
          SL <- SL[-which(norm_SL_2 > v)]
        }
        L <- length(SL)

        # add ai_list to ans.Q0 
        for (l in 1:k){
          ans.Q0[row(ans.Q0)>(sum(dd[1:(l-1)])*(l>1)) & row(ans.Q0)<=sum(dd[1:l]) & col(ans.Q0)==i] <- t(ai_list[[l]])
        }
      }
    }
  } 
  #------------------------------------------------------------------------#
  ans.Qinit <- NULL
  # ans.Qinit takes the loading matrix for each dimension out so length = k, for each k, dim = dk X r 
  for(i in 1:k){
    qtemp <- matrix(ans.Q0[row(ans.Q0)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0)<=sum(dd[1:i])],ncol=r)
    ans.Qinit <- c(ans.Qinit,list(qtemp))
  }
  
  list("Q"=ans.Qinit,"Qmatrix"=ans.Q0,"lambda"=ans.lambda, "r"=r, "distinct_index" = distinct_set, "non_distinct_index_list" = non_distinct_set)
}


eigratio_test <- function(x,kmax){
  dd <- dim(x)
  d <- length(dd) 
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1

  
  x.matrix <- matrix(x,ncol=n)
  ans.lambda <- svd(x.matrix)$d[1:(kmax+1)]^2

  temp1 <- ans.lambda[1:kmax]
  temp2 <- ans.lambda[-1] 
  r <- which.max(temp1 / temp2)
  
  return(r)
} 


# iteration algorithm (ISO) for any Dim tensor time series with CP structure
cp.iso <- function(x,r,cst=0.1,tol=1e-5,niter=100, A.real = NULL,ans.Qiter = NULL,
                   detect_close_eigenvals = FALSE, distinct_set = NULL, non_distinct_set = NULL,
                   c0 = 0.1, L = 1000, v = 0.5){
  # x: tensor of any dimension
  # iterative projection for CP tensor factor models
  # x: d1 * d2 * d3 * ... * dk * n
  # use 2-th lag h cross moment
  
  dd <- dim(x)
  d <- length(dd) # d >= 2, d=k+1
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1

  if(is.null(ans.Qiter)){
    
    if (detect_close_eigenvals){
      ans.init <- cp.init.rppca(x,r,c0,L,v, distinct_set = distinct_set, non_distinct_set = non_distinct_set)
    } else{ ans.init <- cp.init(x,r) }
    
    iiter <- 1
    dis <- 1
    ans.Qiter <- ans.init$Q
    ans.Q0iter <- ans.init$Qmatrix
  }else{
    iiter <- 1
    dis <- 1
    ans.init = list(Q = ans.Qiter)
    ans.Q0iter = vector()
    for (ss in 1:k) {
      ans.Q0iter = rbind(ans.Q0iter,ans.Qiter[[ss]])
    }
  }
  
  
  
  # calculate B^(0) with dim = d X r 
  B <- matrix(0,sum(dd[-d]),r)
  for(i in 1:k){
    A <- matrix(ans.Q0iter[row(ans.Q0iter) > (sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0iter)<=sum(dd[1:i])],ncol=r)
    Sig <- t(A)%*%A
    sig.eig <- eigen(Sig)
    sig.eigvalue <- sig.eig$values*(sig.eig$values>cst)+cst*(sig.eig$values<=cst)
    inv.Sig <- sig.eig$vectors%*%diag(1/sig.eigvalue,nrow=r,ncol=r)%*%t(sig.eig$vectors)
    B[row(B) > (sum(dd[1:(i-1)]) * (i>1)) & row(B)<=sum(dd[1:i])] <- A%*%inv.Sig
    
  }
  
  fnorm.resid <- rep(0,niter)  # collect RMSE
  x.tnsr <- as.tensor(x)
  x.hat.inl <- as.tensor(array(0,dd))
  for(i in 1:r){
    ans.P <- NULL
    for(j in 1:k){
      ans.B <- matrix(B[row(B)>(sum(dd[1:(j-1)])*(j>1)) & row(B)<=sum(dd[1:j])],ncol=r)
      ans.P <- c(ans.P,list(t(ans.B[,i]%*%t(ans.Qiter[[j]][,i]))))
    }
    x.hat.inl <- x.hat.inl + ttl(x.tnsr,ans.P,1:k)   # first X^(iso)
  }
  
  fnorm.resid[1] <- fnorm(x.tnsr@data-x.hat.inl@data)/fnorm(x.tnsr@data)   # error rate in X^(iso)
  x.hat.inl  = x.hat.inl@data
    
  ft.inl <- matrix(0,r,n)
  for(ii in 1:r){
    fx <- x
    for(i in 1:k){
      Bk <- matrix(B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])],ncol=r)
      fx <- tensor(fx,Bk[,ii],1,1)
    }
    ft.inl[ii,] <- fx
  } 
  
  
  iter.error.mat = matrix(0,niter + 1,1) 
  
  while ((dis > tol) && (iiter < niter)){

  if(!is.null(A.real)){
        iter.error.mat[c(iiter:(niter+1)),] = matrix(max(rho2.loss.list(ans.Qiter,A.real)), length(c(iiter:(niter+1))), 1,byrow = T)
  }
    
    # Store loading vectors a_ik
    ans.Q <- matrix(0,sum(dd[-d]),r)
    dis.temp <- 0
    
    # Big loop for Step 3 to 9
    for (i in 1:k){ # i: dimension of data
      Zt <- aperm(x,c(i,(1:d)[-i]))  # put the i^th dimension to the front for convenience
      
      for(ii in 1:r){ # ii: dimension of factors
        z <- Zt
        
        # Step 7 for Z_{_,ii,i}
        for(ll in (1:k)[-i]){
          Bk <- matrix(B[row(B)>(sum(dd[1:(ll-1)])*(ll>1)) & row(B)<=sum(dd[1:ll])],ncol=r)
          z <- tensor(z,Bk[,ii],2,1)
        }
        
        # Step 7 for Sigma 
        M.temp <- tensor(z,z,2,2) / n
  
        # Step 7 for eigenvector and a_{ii,i}
        ans.eig <- eigen(M.temp)
        ans.Q2th <- ans.eig$vectors[,1,drop=FALSE] 
        
        # Step 10: calculate the distance looped on ii and i
        ans.V <- matrix(ans.Q0iter[row(ans.Q0iter)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0iter)<=sum(dd[1:i])],ncol=r)
        dis.temp <- max(dis.temp, sqrt(abs(1-sum(ans.V[,ii]*ans.Q2th)^2)))
        ans.Q[row(ans.Q)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q)<=sum(dd[1:i]) & col(ans.Q)==ii] <- ans.Q2th
      }  
      
      # Back to Step 5: calculate B_i
      A <- matrix(ans.Q[row(ans.Q)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q)<=sum(dd[1:i])],ncol=r)
      
      Sig <- t(A)%*%A
      sig.eig <- eigen(Sig)
      sig.eigvalue <- sig.eig$values*(sig.eig$values>cst)+cst*(sig.eig$values<=cst)
      inv.Sig <- sig.eig$vectors%*%diag(1/sig.eigvalue,nrow=r,ncol=r)%*%t(sig.eig$vectors)
      B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])] <- A%*%inv.Sig
    }
    
    # Write the updated a to ans.Qiter
    for(i in 1:k){
      qtemp <- matrix(ans.Q[row(ans.Q)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q)<=sum(dd[1:i])],ncol=r)
      ans.Qiter[[i]] <- qtemp
    }
    
    # Calculate X^(iso)
    x.hat <- as.tensor(array(0,dd))
    for(i in 1:r){
      ans.P <- NULL
      
      for(j in 1:k){
        ans.B <- matrix(B[row(B)>(sum(dd[1:(j-1)])*(j>1)) & row(B)<=sum(dd[1:j])],ncol=r)
        ans.P <- c(ans.P,list(t(ans.B[,i]%*%t(ans.Qiter[[j]][,i]))))
      }
      x.hat <- x.hat + ttl(x.tnsr,ans.P,1:k) 
    }
    
    fnorm.resid[iiter+1] <- rTensor::fnorm(x.tnsr-x.hat)/rTensor::fnorm(x.tnsr)
    ans.Q0iter <- ans.Q
    
    # The estimation error is the min between error of a and error of X^iso
    dis <- min(dis,dis.temp,abs(fnorm.resid[iiter+1] - fnorm.resid[iiter]))
    iiter <- iiter + 1
    if(iiter==2){
      ans.Qfirst <- ans.Qiter
    }
  }
  
  #iter.error
  iter.error = c(iter.error.mat)
  
  #x.tnsr <- as.tensor(x)
  x0 <- matrix(x,prod(dd[-d]))
  x0 <- t(scale(t(x0),scale=FALSE) )
  x0 <- array(x0,dd)
  norm.percent <- fnorm.resid[iiter]  #fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
  rsquare = 1 -  rTensor::fnorm(x.tnsr-x.hat)^2/rTensor::fnorm(as.tensor(x0))^2
  x.hat <- x.hat@data
  
  f_x0 = rTensor::fnorm(as.tensor(x0))
  
  # Output 
  ft <- matrix(0,r,n)
  for(ii in 1:r){
    fx <- x
    for(i in 1:k){
      Bk <- matrix(B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])],ncol=r)
      fx <- tensor(fx,Bk[,ii],1,1)
    }
    ft[ii,] <- fx
  } 
  w <- sqrt(apply(ft^2, 1, sum)) 
  #ft <- sqrt(n) * diag(w^(-1)) %*% ft
  list("Q"=ans.Qiter,"Qfirst"=ans.Qfirst,"Qinit"=ans.init$Q,"ft"=ft,"ft.inl" = ft.inl,"dis"=dis,"niter"=iiter,
       "norm.percent"=norm.percent,"rsquare"=rsquare,"fnorm.resid"=fnorm.resid[1:iiter],"x.hat"=x.hat, 
       "w" = w, "x.hat.inl" = x.hat.inl, "iter_error" = iter.error, "f_x0" = f_x0)
} 

#-----------------------AC-ISO----------------------------#
cp.init.tensor <- function(x,r,h0=1){
  # x: tensor of any dimension
  # CP initialization
  # x: d1 * d2 * d3 * ... * dk * n
  # use 2-th lag h cross moment
  
  dd <- dim(x)
  d <- length(dd) # d >= 2, d=k+1
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1
  x.matrix <- matrix(x,ncol=n)
  M.temp <- matrix(0,dd.prod,dd.prod)
  #for(h in 1:h0){
  h <- h0
  x.left <- array(x.matrix[,1:(n-h)],c(dd[-d],n-h))
  x.right <- array(x.matrix[,(h+1):n],c(dd[-d],n-h))
  Omega <- tensor(x.left,x.right,d,d)/(n-h)
  Omega.matrix <- matrix(Omega,dd.prod,dd.prod)
  M.temp <- M.temp + Omega.matrix/2+t(Omega.matrix)/2   #Omega.matrix%*%t(Omega.matrix)
  #}
  ans.eig <- eigen(M.temp)
  ans.Q2th <- ans.eig$vectors[,1:r^2,drop=FALSE]
  ans.lambda <- ans.eig$values
  #ans.lambda[1:r^3]
  ans.Q0 <- matrix(0,sum(dd[-d]),r)
  for(i in 1:r){
    ans.Q2 <- array(ans.Q2th[,i],dd[-d])
    for(l in 1:k){
      m.temp <- tensor(ans.Q2,ans.Q2,c(1:k)[-l],c(1:k)[-l])
      ans.eig.mtemp <- eigen(m.temp)
      ans.Q0[row(ans.Q0)>(sum(dd[1:(l-1)])*(l>1)) & row(ans.Q0)<=sum(dd[1:l]) & col(ans.Q0)==i] <- 
        ans.eig.mtemp$vectors[,1,drop=FALSE]
    }
  }
  ans.Qinit <- NULL
  for(i in 1:k){
    qtemp <- matrix(ans.Q0[row(ans.Q0)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0)<=sum(dd[1:i])],ncol=r)
    ans.Qinit <- c(ans.Qinit,list(qtemp))
  }
  
  list("Q"=ans.Qinit,"Qmatrix"=ans.Q0,"lambda"=ans.lambda)
}  

# AC-ISO
cp.iso.han <- function(x,r,h0=1,A.real = NULL,ans.Qiter = NULL,cst=0.1,tol=1e-5,niter=100){
  # x: tensor of any dimension
  # iterative projection for CP tensor factor models
  # x: d1 * d2 * d3 * ... * dk * n
  # use 2-th lag h cross moment
  
  dd <- dim(x)
  d <- length(dd) # d >= 2, d=k+1
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1
  if(is.null(ans.Qiter)){
    ans.init <- cp.init.tensor(x,r,h0)
    iiter <- 1
    dis <- 1
    ans.Qiter <- ans.init$Q
    ans.Q0iter <- ans.init$Qmatrix
  }else{
    iiter <- 1
    dis <- 1
    ans.init = list(Q = ans.Qiter)
    ans.Q0iter = vector()
    for (ss in 1:k) {
      ans.Q0iter = rbind(ans.Q0iter,ans.Qiter[[ss]])
    }

  }
  
  B <- matrix(0,sum(dd[-d]),r)
  for(i in 1:k){
    A <- matrix(ans.Q0iter[row(ans.Q0iter)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0iter)<=sum(dd[1:i])],ncol=r)
    Sig <- t(A)%*%A
    sig.eig <- eigen(Sig)
    sig.eigvalue <- sig.eig$values*(sig.eig$values>cst)+cst*(sig.eig$values<=cst)
    inv.Sig <- sig.eig$vectors%*%diag(1/sig.eigvalue,nrow=r,ncol=r)%*%t(sig.eig$vectors)
    B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])] <- A%*%inv.Sig
  }
  fnorm.resid <- rep(0,niter)
  x.tnsr <- as.tensor(x)
  x.hat.inl <- as.tensor(array(0,dd))
  for(i in 1:r){
    ans.P <- NULL
    for(j in 1:k){
      ans.B <- matrix(B[row(B)>(sum(dd[1:(j-1)])*(j>1)) & row(B)<=sum(dd[1:j])],ncol=r)
      ans.P <- c(ans.P,list(t(ans.B[,i]%*%t(ans.Qiter[[j]][,i]) ) ))
    }
    x.hat.inl <- x.hat.inl + ttl(x.tnsr,ans.P,1:k) 
  }
  fnorm.resid[1] <- rTensor::fnorm(x.tnsr-x.hat.inl)/rTensor::fnorm(x.tnsr)

  x.hat.inl <- x.hat.inl@data
  ft.inl <- matrix(0,r,n)
  for(ii in 1:r){
    fx <- x
    for(i in 1:k){
      Bk <- matrix(B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])],ncol=r)
      fx <- tensor(fx,Bk[,ii],1,1)
    }
    ft.inl[ii,] <- fx
  }
  
  
  iter.error.mat = matrix(0,niter+1,1) 
  
 while((dis > tol) && (iiter < niter)){
    
   if(!is.null(A.real)){
     iter.error.mat[c(iiter:(niter+1)),] = matrix(max(rho2.loss.list(ans.Qiter,A.real)), length(c(iiter:(niter+1))), 1,byrow = T)
   }
   
    ans.Q <- matrix(0,sum(dd[-d]),r)
    dis.temp <- 0
    for(i in 1:k){
      Zt <- aperm(x,c(i,(1:d)[-i]))
      for(ii in 1:r){
        z <- Zt
        for(ll in (1:k)[-i]){
          Bk <- matrix(B[row(B)>(sum(dd[1:(ll-1)])*(ll>1)) & row(B)<=sum(dd[1:ll])],ncol=r)
          z <- tensor(z,Bk[,ii],2,1)
        }
        
        M.temp <- matrix(0,dd[i],dd[i])
        #for(h in 1:h0){
        h <- h0
        x.left <- z[,1:(n-h)]
        x.right <- z[,(h+1):n]
        Omega <- tensor(x.left,x.right,2,2)/(n-h)
        M.temp <- M.temp + Omega/2+t(Omega)/2   
        #}
        ans.eig <- eigen(M.temp)
        ans.Q2th <- ans.eig$vectors[,1,drop=FALSE]
        ans.V <- matrix(ans.Q0iter[row(ans.Q0iter)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q0iter)<=sum(dd[1:i])],ncol=r)
        dis.temp <- max(dis.temp, sqrt(abs(1-sum(ans.V[,ii]*ans.Q2th)^2)))
        ans.Q[row(ans.Q)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q)<=sum(dd[1:i]) & col(ans.Q)==ii] <- ans.Q2th
      }  
      A <- matrix(ans.Q[row(ans.Q)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q)<=sum(dd[1:i])],ncol=r)
      Sig <- t(A)%*%A
      sig.eig <- eigen(Sig)
      sig.eigvalue <- sig.eig$values*(sig.eig$values>cst)+cst*(sig.eig$values<=cst)
      inv.Sig <- sig.eig$vectors%*%diag(1/sig.eigvalue,nrow=r,ncol=r)%*%t(sig.eig$vectors)
      B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])] <- A%*%inv.Sig  
    }
    for(i in 1:k){
      qtemp <- matrix(ans.Q[row(ans.Q)>(sum(dd[1:(i-1)])*(i>1)) & row(ans.Q)<=sum(dd[1:i])],ncol=r)
      ans.Qiter[[i]] <- qtemp
    }
    x.hat <- as.tensor(array(0,dd))
    for(i in 1:r){
      ans.P <- NULL
      for(j in 1:k){
        ans.B <- matrix(B[row(B)>(sum(dd[1:(j-1)])*(j>1)) & row(B)<=sum(dd[1:j])],ncol=r)
        ans.P <- c(ans.P,list(t(ans.B[,i]%*%t(ans.Qiter[[j]][,i]) ) ))
      }
      x.hat <- x.hat + ttl(x.tnsr,ans.P,1:k) 
    }
    fnorm.resid[iiter+1] <- rTensor::fnorm(x.tnsr-x.hat)/rTensor::fnorm(x.tnsr)
    ans.Q0iter <- ans.Q
    dis <- min(dis,dis.temp,abs(fnorm.resid[iiter+1] - fnorm.resid[iiter]))
    iiter <- iiter + 1
    if(iiter==2){
      ans.Qfirst <- ans.Qiter
    }
 }
  
  #iter.error
  iter.error = c(iter.error.mat)
  #x.tnsr <- as.tensor(x)
  x0 <- matrix(x,prod(dd[-d]))
  x0 <- t(scale(t(x0),scale=FALSE) )
  x0 <- array(x0,dd)
  norm.percent <- fnorm.resid[iiter]  #fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
  rsquare = 1 - rTensor::fnorm(x.tnsr-x.hat)^2/rTensor::fnorm(as.tensor(x0))^2
  f_x0 = rTensor::fnorm(as.tensor(x0))
  x.hat <- x.hat@data
  ft <- matrix(0,r,n)
  for(ii in 1:r){
    fx <- x
    for(i in 1:k){
      Bk <- matrix(B[row(B)>(sum(dd[1:(i-1)])*(i>1)) & row(B)<=sum(dd[1:i])],ncol=r)
      fx <- tensor(fx,Bk[,ii],1,1)
    }
    ft[ii,] <- fx
  }
  w <- sqrt(apply(ft^2, 1, sum)) 
  #ft <- diag(w^(-1)) %*% ft
  list("Q"=ans.Qiter,"Qfirst"=ans.Qfirst,"Qinit"=ans.init$Q,"ft"=ft,"dis"=dis,"niter"=iiter, iter_error = iter.error,"ft.inl" = ft.inl,
       "norm.percent"=norm.percent,"rsquare"=rsquare,"fnorm.resid"=fnorm.resid[1:iiter],"x.hat"=x.hat, "x.hat.inl"=x.hat.inl, "w" = w)
}


#------------------theorem 4.3(i) CLT-------------------# 
Bi_list <- function(Alist,i,k,trans = FALSE){
    bi_list <- NULL 
    for (ik in 1:2){
        if (ik == k){
            bi_list[[ik]] <- 0
        }
        else{ 
          if (trans){
            bi_list[[ik]]<- matrix((Alist[[ik]] %*% solve((t(Alist[[ik]]) %*% Alist[[ik]])))[,i], nrow = 1 )
          }
          else{
            bi_list[[ik]] <- matrix((Alist[[ik]] %*% solve((t(Alist[[ik]]) %*% Alist[[ik]])))[,i], ncol = 1)
          }

        }
    }
    bi_list 
}

cltvar <- function(u, k, dk, K, lambda_i, Theta_ii, aik, bi_list, sigma_e){ 
  ### bi_list includes b_k. 
  Paik <- diag(x=1, nrow = dk) - aik %*% t(aik) 
  templist <- bi_list 
  templist[[k]] <- Paik %*% u 
  tempvec <- 1

  for (ik in 1:K){
    tempvec <- kronecker(tempvec, templist[[K - ik + 1]])
  }

  # for (ik in 1:K){
  #   tempvec <- kronecker(tempvec, templist[[ik]])
  # }

  tempvec <- matrix(tempvec, ncol = 1)

  Theta_ii * ( t(tempvec) %*%  sigma_e %*% tempvec ) / lambda_i^2
} 

cltstat <- function(n,u,aik_hat,aik){
  
  sqrt(n) * t(u) %*% (aik_hat - (2 * (sum(aik_hat * aik) > 0) - 1) * aik)
}

cltstat_leadingterm <- function(n, u, k, dk, K, lambda_i, aik, bi_list, fi, E){
  Paik <- diag(x=1, nrow = dk) - aik %*% t(aik) 
  bi_list[[k]] <- t(Paik %*% u)

  temp <- matrix(ttl(as.tensor(E), bi_list, 1:K)@data, ncol = 1) 
  
  t(fi) %*% temp / (lambda_i * n) 
}


#------------------sign permutation for max R2-------------------# 

# sign permutation vector, diagonalizing yields sign permutation matrix
perm <- function(n){
  temp <- matrix(0,nrow=n,ncol=1) 
  for (i in 0:n){
    temp <- cbind(temp,combn(n,i,function(x)replace(numeric(n)-1,x,1)))
  }
  temp[,-1]
} 

# get Xhat from loadings and factors 
get_X_ndim <- function(Alist,ft){
  r <- dim(ft)[1] 
  n <- dim(ft)[2] 
  K <- length(Alist) 
  dim_diagf <- append(rep(r,K),n)
  diagf <- array(0,dim_diagf)
  for (ir in 1:r){
    for (ii in 1:n){
      diag_idx <- append(rep(ir,K),ii)
      diagf[t(diag_idx)] <-ft[ir,ii]
    }
  }
  ttl(as.tensor(diagf),Alist,1:K)@data 
}

# R2 
get_r2 <- function(xhat,x,n){ 
  dd <- dim(x) 
  d <- length(dd) 
  x0 <- matrix(x,prod(dd[-d]))
  x0 <- t(scale(t(x0),scale=FALSE) )
  x0 <- array(x0,dd)
  1 - sum( (xhat - x)^2) / sum(x0^2)
}

# get max R2 by sign permuting
max_r2 <- function(Alist, ft,x){ 
  n <- dim(ft)[2] 
  r <- dim(ft)[1] 
  permidx <- perm(r) 
  
  r2vec <- vector(length = dim(permidx)[2])
  for (i in 1:dim(permidx)[2]){
    idx <- permidx[,i] 
    ftperm <-  diag(idx) %*% ft 
    xperm <- get_X_ndim(Alist, ft)
    r2vec[i] <- get_r2(xperm,x,n)
  }
  
  maxr2 <- max(r2vec) 
  maxr2_perm <- permidx[,which.max(r2vec)]

  list("max_r2" = maxr2, "max_r2_perm" = maxr2_perm)
}

#------------------rank estimation-------------------#

topup <- function(x,dd,K,n){
    Omega <- tensor(x,x,K+1,K+1) / n 

    M <- NULL
    for (k in 1:K){
        temp_k <- tensor(Omega,Omega,c(1:(2*K))[-k], c(1:(2*K))[-k])
        M <- c(M, list(temp_k)) 
    }
    M 
}

tipup <- function(x,dd,K,n){ 
    M <- NULL 
    for (k in 1:K){
        temp_k <- tensor(x,x,c(1:(K+1))[-k], c(1:(K+1))[-k]) / n
        M <- c(M, list(temp_k)) 
    }
    M
}

rank_er_est_modewide <- function(x, tucker = FALSE){
    dd <- dim(x) 
    K <- length(dd) - 1
    n <- dd[K + 1] 
    M_tipup <- tipup(x,dd,K,n) 
    M_topup <- topup(x,dd,K,n) 

    r_est <- matrix(0,nrow = K, ncol = 2) 
    for (k in 1:K){
        M_ip_K <- M_tipup[[k]] 
        M_up_K <- M_topup[[k]] 

        M_ip_K_eigen <- eigen(M_ip_K)$values 
        M_up_K_eigen <- eigen(M_up_K)$values 

        r_est[k,1] <- which.min(M_ip_K_eigen[2:floor(dd[k]/2)] / M_ip_K_eigen[1:(floor(dd[k]/2) - 1)])

        r_est[k,2] <- which.min(M_up_K_eigen[2:floor(dd[k]/2)] / M_up_K_eigen[1:(floor(dd[k]/2) - 1)])
    } 
    if (tucker){
        return(r_est) 
    } 
    else {
        res <- matrix(0,nrow = 4, ncol = 2)
        res[1,] <- apply(r_est,2,mean) 
        res[2,] <- apply(r_est,2,median) 
        res[3,] <- apply(r_est,2,max) 
        res[4,] <- apply(r_est,2,min) 

        res
    }
}

rank_er_est_agg <- function(x,kmax, kmin = 1){
  dd <- dim(x)
  d <- length(dd)
  n <- dd[d]
  dd.prod <- prod(dd) / n
  k <- d-1
  
  x.matrix <- matrix(x,ncol=n)
  ans.lambda <- svd(x.matrix)$d[kmin:(kmax+1)]
  
  temp1 <- ans.lambda[-length(ans.lambda)]
  temp2 <- ans.lambda[-1] 
  r <- which.max(temp1 / temp2) + (kmin - 1)
  
  r
} 




Han.forecast = function(Y, r, forecast.step = 2, ans.Qiter = NULL, A.real = NULL, niter = 100,diff = F){
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  Yh =  aperm(Y,c(2:(m+1),1)) 
  
  res.han  = cp.iso.han(x = Yh,r = r, ans.Qiter = NULL, A.real = NULL, niter = 100)
  
  
  A.iter = res.han$Q
  A.inl = res.han$Qinit 
  
  f.iter = t(res.han$ft)
  f.inl = t(res.han$ft.inl)
  
  f.fore.inl  =   as.matrix(forecast.factor.tensor(f.inl,lag.max = 6,sp = forecast.step,diff = diff))
  f.fore.iter =   as.matrix(forecast.factor.tensor(f.iter,lag.max = 6,sp = forecast.step,diff = diff))
  
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
  
  return(list(Y.forecast.inl=Y.forecast.inl,Y.forecast.iter=Y.forecast.iter,A.inl=A.inl,A.iter=A.iter,f.inl=f.inl,f.iter=f.iter,con.AMP=res.han))
}

Chen.forecast = function(Y, r, forecast.step = 2, A.real = NULL, niter = 100,diff = F){
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  Yh =  aperm(Y,c(2:(m+1),1)) 
  
  res  = cp.iso(x = Yh,r = r, A.real = NULL, niter = 100)
  
  
  A.iter = res$Q
  A.inl  = res$Qinit 
  
  f.iter = t(res$ft)
  f.inl = t(res$ft.inl)
  
  f.fore.inl  =  as.matrix(forecast.factor.tensor(f.inl,lag.max = 6,sp = forecast.step, diff = diff))
  f.fore.iter =  as.matrix(forecast.factor.tensor(f.iter,lag.max = 6,sp = forecast.step, diff = diff))
  
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
  
  return(list(Y.forecast.inl=Y.forecast.inl,Y.forecast.iter=Y.forecast.iter,A.inl=A.inl,A.iter=A.iter,f.inl=f.inl,f.iter=f.iter,con.AMP=res))
}



