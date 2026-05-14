library(pracma)
library(Rcpp)
library(RcppEigen)
library(foreach)
library(doParallel)
library(jointDiag)
library(HDTSA)



detectCores() 
cl <- makeCluster(120)
registerDoParallel(cl)


n_tol = c(400,800)           #sample size   
delta.fl_tol = c(0.25,0.75)  #correlation of factor loadings 
factor.loading_tol = c(1)    #1:random sparse; 2:structure sparse
factor.corr_tol = c(0,0.75)  #correlation of factors
ALPHA_tol = c(0,0.3,0.6)     #sparsity

D = c(20,20)
m = length(D)
r = 3
beta = c(0.8,0.75,0.7)
ar.coef = as.list(beta)
w = rep(15,r) 

K = 10

delta.fl       = 0.25
factor.loading = "sparse-random-corr1"
factor.corr    = 0.25
ALPHA          = 0
n              = 800


rho2.loss.list <- HDTSA:::rho2.loss.list
fnorm <- HDTSA:::fnorm
rho2.f.loss <- HDTSA:::rho2.f.loss
aij.debias.iter <- HDTSA:::aij.debias.iter
cov.aij.debias.iter <- HDTSA:::cov.aij.debias.iter
cov.aij.debias.iter.est <-  HDTSA:::cov.aij.debias.iter.est
DGP.TCP <- HDTSA:::DGP.TCP
tensor.est.xi  <- HDTSA:::tensor.est.xi 
HDTTS.CP.est <- HDTSA:::HDTTS.CP.est 
RP.xi.sel <- HDTSA:::RP.xi.sel


Rolling_fun2 <- function(x, n, ALPHA, delta.fl, factor.loading, factor.corr, seed_num) {
  
  source("cp_cciso.R")
  
  set.seed(seed_num + x)
  
  D <- c(20, 20)
  m <- length(D)
  r <- 3
  beta <- c(0.8, 0.75, 0.7)
  ar.coef <- as.list(beta)
  w <- rep(15, r)
  K <- 10
  
  data <- DGP.TCP(
    n = n, m = m, D = D, r = r, w = w, ar.coef = ar.coef,
    factor.loading = factor.loading, factor.corr = factor.corr,
    alpha = ALPHA, delta = delta.fl, par_E = 1, heavytail = 5,
    error.ar = FALSE
  )
  
  A <- data$A
  Y <- data$Y
  f <- data$f
  C <- data$C
  
  J <- diag(1 - factor.corr, r, r) + matrix(factor.corr, r, r)
  
  make_bmp_control <- function(ratio.type = "log",
                               threshold = TRUE,
                               delta2 = rep(1, m),
                               iter.max = 20,
                               eps = 1e-5,
                               print.eps = FALSE,
                               iter.lag = 1,
                               all.put = TRUE,
                               A = NULL,
                               component = NULL) {
    list(
      lag.k.dpi = K,
      threshold = threshold,
      delta = NULL,
      delta2 = delta2,
      ratio.type = ratio.type,
      random.projection = FALSE,
      iter.max = iter.max,
      eps = eps,
      grid.num = 50,
      delta.max = 0.1,
      print.eps = print.eps,
      iter.lag = iter.lag,
      all.put = all.put,
      A = A,
      component = component
    )
  }
  
  run_bmp <- function(xi, Rank = NULL, Ratio.type = "log", A.init = NULL,
                      Threshold = TRUE, delta2 = rep(1, m),
                      iter_max = 20, eps = 1e-5, print.eps = FALSE,
                      iter_lag = 1, all.put = TRUE,
                      A = NULL, Component = NULL) {
    
    rank.arg <- Rank
    if (is.list(rank.arg)) {
      rank.arg <- rank.arg$d
    }
    
    CP_TTS(
      Y = Y,
      xi = xi,
      r = rank.arg,
      A.init = A.init,
      control.DPI = make_bmp_control(
        ratio.type = Ratio.type,
        threshold = Threshold,
        delta2 = delta2,
        iter.max = iter_max,
        eps = eps,
        print.eps = print.eps,
        iter.lag = iter_lag,
        all.put = all.put,
        A = A,
        component = Component
      )
    )
  }
  
  get_rank <- function(fit) {
    if (!is.null(fit$r)) {
      return(fit$r)
    }
    if (!is.null(fit$A.hat)) {
      return(ncol(fit$A.hat[[1]]))
    }
    stop("Cannot determine the estimated rank from the CP_TTS output.")
  }
  
  iter.num <- 20
  gg <- 0
  
  while (iter.num > 10 && gg < 50) {
    xi.chang <- tensor.est.xi(Y)
    
    res.only.used.rank <- HDTTS.CP.est(
      Y = Y, xi = xi.chang, K = K, Ratio.type = "log"
    )
    
    r_hat_first <- res.only.used.rank$r.hat
    
    xi.res <- RP.xi.sel(
      Y,
      r_breve = 2 * r_hat_first,
      eps = 0.1,
      lag.k = K,
      Randomized.time = 50
    )
    
    con.BMP <- run_bmp(
      xi = xi.res$xi.sel,
      Rank = NULL,
      Ratio.type = "log",
      A.init = NULL,
      Threshold = TRUE,
      delta2 = rep(1, m),
      iter_max = 20,
      eps = 1e-5,
      print.eps = FALSE,
      iter_lag = 1,
      all.put = TRUE,
      A = A,
      Component = C
    )
    
    iter.num <- con.BMP$iter.step
    gg <- gg + 1
  }
  
  con_pro <- c(
    max(rho2.loss.list(con.BMP$A.hat, A)),
    max(rho2.loss.list(con.BMP$A.init, A))
  )
  
  CP_loss_pro <- con.BMP$CP.loss
  f.pro.inl <- con.BMP$f.hat.inl
  f.pro.iter <- con.BMP$f.hat
  
  con.BMP.classical <- run_bmp(
    xi = xi.res$xi.sel,
    Rank = NULL,
    Ratio.type = "classical",
    A.init = NULL,
    Threshold = TRUE,
    delta2 = rep(1, m),
    iter_max = 20,
    eps = 1e-5,
    print.eps = FALSE,
    iter_lag = 1,
    all.put = TRUE,
    A = A,
    Component = NULL
  )
  
  A.init.pro <- con.BMP$A.init
  
  rank.est <- get_rank(con.BMP)
  rank.est.NT <- get_rank(con.BMP.classical)
  
  ## Chang (2023)
  res.chang <- HDTSA::CP_MTS(
    Y = Y, xi = NULL, Rank = NULL, lag.k = 10, method = "CP.Refined"
  )
  
  A_chang <- res.chang$A
  B_chang <- res.chang$B
  f.chang <- res.chang$f
  
  C_chang <- 0
  for (k in 1:res.chang$Rank$d) {
    C_chang <- C_chang + f.chang[, k] %o% A_chang[, k] %o% B_chang[, k]
  }
  
  CP_loss_chang <- fnorm(C_chang - C) / sqrt(n * prod(D))
  rank.chang <- res.chang$Rank$d
  A.chang <- list(res.chang$A, res.chang$B)
  
  ## Han (2024)
  Y1 <- aperm(Y, c(2:(m + 1), 1))
  res.han <- cp.iso.han(x = Y1, r = rank.est, A.real = data$A, niter = 20)
  
  con_han <- c(
    max(rho2.loss.list(res.han$Q, data$A)),
    max(rho2.loss.list(res.han$Qfirst, data$A)),
    max(rho2.loss.list(res.han$Qinit, data$A))
  )
  
  f.han.inl <- t(res.han$ft.inl)
  f.han.iter <- t(res.han$ft)
  C1 <- aperm(C, c(2:(m + 1), 1))
  CP_loss_han_iter <- fnorm(C1 - res.han$x.hat) / sqrt(n * prod(D))
  CP_loss_han_inl <- fnorm(C1 - res.han$x.hat.inl) / sqrt(n * prod(D))
  
  ## Chen (2024)
  rank.chen <- eigratio_test(x = Y1, kmax = 5)
  res.chen <- cp.iso(
    x = Y1, r = rank.chen, A.real = A, niter = 20,
    detect_close_eigenvals = TRUE
  )
  
  con_chen <- c(
    max(rho2.loss.list(res.chen$Q, data$A)),
    max(rho2.loss.list(res.chen$Qfirst, data$A)),
    max(rho2.loss.list(res.chen$Qinit, data$A))
  )
  
  f.chen.inl <- t(res.chen$ft.inl)
  f.chen.iter <- t(res.chen$ft)
  CP_loss_chen_iter <- fnorm(C1 - res.chen$x.hat) / sqrt(n * prod(D))
  CP_loss_chen_inl <- fnorm(C1 - res.chen$x.hat.inl) / sqrt(n * prod(D))
  
  f_loss_inl <- c(
    rho2.f.loss(f.pro.inl, f),
    rho2.f.loss(f.han.inl, f),
    rho2.f.loss(f.chen.inl, f),
    rho2.f.loss(f.chang, f)
  )
  
  f_loss_iter <- c(
    rho2.f.loss(f.pro.iter, f),
    rho2.f.loss(f.han.iter, f),
    rho2.f.loss(f.chen.iter, f)
  )
  
  f_loss <- c(f_loss_iter, f_loss_inl)
  
  CP_loss_iter <- c(CP_loss_pro[2], CP_loss_han_iter, CP_loss_chen_iter)
  CP_loss_inl <- c(CP_loss_pro[1], CP_loss_han_inl, CP_loss_chen_inl, CP_loss_chang)
  CP_loss <- c(CP_loss_iter, CP_loss_inl)
  
  #### Iterative comparison (all methods know rank r)
  res.han.iter <- cp.iso.han(x = Y1, r = r, A.real = data$A, niter = 20)
  A.init.han <- res.han.iter$Qinit
  iter.error.han <- res.han.iter$iter_error
  iter.num.han <- res.han.iter$niter
  
  res.chen.iter <- cp.iso(
    x = Y1, r = r, A.real = A, niter = 20,
    detect_close_eigenvals = TRUE
  )
  A.init.chen <- res.chen.iter$Qinit
  iter.error.chen <- res.chen.iter$iter_error
  iter.num.chen <- res.chen.iter$niter
  
  res.pro.pro.iter <- run_bmp(
    xi = xi.res$xi.sel,
    Rank = r,
    Ratio.type = "log",
    A.init = NULL,
    Threshold = TRUE,
    delta2 = rep(1, m),
    iter_max = 20,
    eps = 1e-5,
    print.eps = FALSE,
    iter_lag = 1,
    all.put = TRUE,
    A = A,
    Component = NULL
  )
  
  res.pro.han.iter <- run_bmp(
    xi = xi.res$xi.sel,
    Rank = r,
    Ratio.type = "log",
    A.init = A.init.han,
    Threshold = TRUE,
    delta2 = rep(1, m),
    iter_max = 20,
    eps = 1e-5,
    print.eps = FALSE,
    iter_lag = 1,
    all.put = TRUE,
    A = A,
    Component = NULL
  )
  
  res.pro.chen.iter <- run_bmp(
    xi = xi.res$xi.sel,
    Rank = r,
    Ratio.type = "log",
    A.init = A.init.chen,
    Threshold = TRUE,
    delta2 = rep(1, m),
    iter_max = 20,
    eps = 1e-5,
    print.eps = FALSE,
    iter_lag = 1,
    all.put = TRUE,
    A = A,
    Component = NULL
  )
  
  iter.num.pro <- res.pro.pro.iter$iter.step
  iter.error.pro <- res.pro.pro.iter$iter.error
  iter.error.han_init.pro_iter <- res.pro.han.iter$iter.error
  iter.error.chen_init.pro_iter <- res.pro.chen.iter$iter.error
  
  res.han.pro <- cp.iso.han(
    x = Y1, r = r, A.real = A,
    ans.Qiter = res.pro.pro.iter$A.init,
    niter = 20
  )
  iter.error.han_iter.pro_init <- res.han.pro$iter_error
  
  res.chen.pro <- cp.iso(
    x = Y1, r = r, A.real = A,
    ans.Qiter = res.pro.pro.iter$A.init,
    niter = 20
  )
  iter.error.chen_iter.pro_init <- res.chen.pro$iter_error
  
  iter_tol <- c(
    iter.error.pro,
    iter.error.han,
    iter.error.chen,
    iter.error.han_init.pro_iter,
    iter.error.chen_init.pro_iter,
    iter.error.han_iter.pro_init,
    iter.error.chen_iter.pro_init
  )
  
  #### Normality (known r)
  A.hat <- res.pro.pro.iter$A.init
  A.iter <- res.pro.pro.iter$A.hat
  f.iter <- res.pro.pro.iter$f.hat
  Sigma.yij.xii.1 <- res.pro.pro.iter$Sigma.yij.xii.1
  
  if (max(rho2.loss.list(A.hat, A)) < 0.3) {
    j <- 1
    i <- 1
    
    place.iter <- which.max(abs(cor(A[[j]][, i], A.iter[[j]])))
    
    aij.iter <- A.iter[[j]][, place.iter]
    aij <- A[[j]][, i]
    
    res_debias_iter <- aij.debias.iter(
      A = A.iter,
      i = place.iter,
      j = j,
      Sigma.yij.xii.1 = Sigma.yij.xii.1
    )
    
    aij_debias_iter <- c(res_debias_iter$aij.de)
    
    h1 <- c(1, rep(0, length(aij_debias_iter) - 1))
    h2 <- rep(1 / sqrt(length(aij_debias_iter)), length(aij_debias_iter))
    
    v.iter.1 <- cov.aij.debias.iter(
      h = h1, A = A, i = i, j = j, f = f,
      COV.VECE = diag(prod(D)), w = w
    )
    
    v.iter.2 <- cov.aij.debias.iter(
      h = h2, A = A, i = i, j = j, f = f,
      COV.VECE = diag(prod(D)), w = w
    )
    
    v.iter.1.est <- cov.aij.debias.iter.est(
      h1, i = place.iter, j = j, A = A.iter, f = f.iter, Y = Y
    )
    
    v.iter.2.est <- cov.aij.debias.iter.est(
      h2, i = place.iter, j = j, A = A.iter, f = f.iter, Y = Y
    )
    
    mse1 <- sqrt(n) * as.numeric(
      t(h1) %*% (aij_debias_iter - sign(cor(aij, aij_debias_iter)) * aij)
    ) / sqrt(v.iter.1)
    
    mse2 <- sqrt(n) * as.numeric(
      t(h2) %*% (aij_debias_iter - sign(cor(aij, aij_debias_iter)) * aij)
    ) / sqrt(v.iter.2)
    
    mse3 <- sqrt(n) * as.numeric(
      t(h1) %*% (aij_debias_iter - sign(cor(aij, aij_debias_iter)) * aij)
    ) / sqrt(v.iter.1.est)
    
    mse4 <- sqrt(n) * as.numeric(
      t(h2) %*% (aij_debias_iter - sign(cor(aij, aij_debias_iter)) * aij)
    ) / sqrt(v.iter.2.est)
    
    mse.dis <- c(mse1, mse2, mse3, mse4, mse5 = 0, mse6 = 0)
  } else {
    mse.dis <- rep(-999, 6)
  }
  
  res <- c(
    con_pro,
    con_han,
    con_chen,
    max(rho2.loss.list(A.chang, A)),
    iter.num,
    iter.num.han,
    iter.num.chen,
    rank.est,
    rank.est.NT,
    rank.chang,
    rank.chen,
    mse.dis,
    iter_tol,
    f_loss,
    CP_loss
  )
  
  return(res)
}



 
############################## settings ##############################

rep <- 2000
seed_num <- 123456

con_list <- list()
meta_list <- list()

hgl <- 1

############################## 1. run simulation and only store con_list ##############################

for (factor.corr in factor.corr_tol) {
  
  for (factor.loading.num in factor.loading_tol) {
    
    factor.loading <- paste0("sparse-random-corr", factor.loading.num)
    
    for (delta.fl in delta.fl_tol) {
      
      for (ALPHA in ALPHA_tol) {
        
        for (n in n_tol) {
          
          source("cp_cciso.R")
          
          con <- foreach(
            x = 1:rep,
            .combine = "rbind",
            .packages = c("HDTSA","R.utils","tensor","rTensor","MASS",
                          "pracma","Rcpp","RcppArmadillo","RcppEigen",
                          "vars","base"),
            .errorhandling = "pass",
            .verbose = FALSE
          ) %dopar% Rolling_fun2(
            x = x,
            n = n,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading = factor.loading,
            factor.corr = factor.corr,
            seed_num = seed_num
          )
          
          con_list[[hgl]] <- con
          
          meta_list[[hgl]] <- data.frame(
            n = n,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading = factor.loading.num,
            factor.corr = factor.corr
          )
          
          name_i <- paste0(
            "(factor.corr)", factor.corr,
            "(random-sparse-corr)", factor.loading.num,
            "(delta.fl)", delta.fl,
            "(alpha)", ALPHA,
            "(n)", n,
            "(count)", NROW(con)
          )
          
          print(name_i)
          cat("\n")
          
          seed_num <- seed_num + rep
          hgl <- hgl + 1
        }
      }
    }
  }
}

meta_df <- do.call(rbind, meta_list)
rownames(meta_df) <- NULL

direction = paste0("simulation_estimation_heavytailed_",format(Sys.time(), "%Y%m%d-%H-%M-%S"),".RData")

save.image(direction)


stopImplicitCluster()
stopCluster(cl)                    


############################## 2. summarize rho(A_hat, A) from con_list ##############################

rho_method_name_old <- c(
  "iter", "inl",
  "han.iter", "han.1", "han.inl",
  "chen.iter", "chen.1", "chen.inl",
  "chang"
)

rho_method_keep_old <- c(
  "iter", "han.iter", "chen.iter",
  "inl", "han.inl", "chen.inl", "chang"
)

rho_method_keep_new <- c(
  "Pro.iter", "HOPE", "CC-ISO",
  "Pro.init", "cPCA", "RP-PCA", "RCP"
)

extract_rho_table_from_conlist <- function(con_list, meta_df) {
  
  result_tol <- vector("list", length(con_list))
  
  for (ii in seq_along(con_list)) {
    
    con <- con_list[[ii]]
    
    if (is.null(con) || length(con) == 0) next
    
    con <- as.matrix(con)
    storage.mode(con) <- "numeric"
    
    # Updated Rolling_fun2:
    # first 9 columns are
    # iter, inl, han.iter, han.1, han.inl, chen.iter, chen.1, chen.inl, chang
    con_error <- con[, 1:length(rho_method_name_old), drop = FALSE]
    colnames(con_error) <- rho_method_name_old
    
    # remove han.1 and chen.1
    con_error <- con_error[, rho_method_keep_old, drop = FALSE]
    colnames(con_error) <- rho_method_keep_new
    
    mean.k <- colMeans(con_error, na.rm = TRUE)
    sd.k   <- apply(con_error, 2, sd, na.rm = TRUE)
    
    names(mean.k) <- paste0(colnames(con_error), ".mean")
    names(sd.k)   <- paste0(colnames(con_error), ".sd")
    
    result_tol[[ii]] <- c(
      meta_df[ii, ],
      mean.k,
      sd.k
    )
  }
  
  result_tol <- do.call(rbind, result_tol)
  result_tol <- as.data.frame(result_tol, stringsAsFactors = FALSE)
  
  numeric_cols <- setdiff(colnames(result_tol),
                          c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr"))
  result_tol[numeric_cols] <- lapply(result_tol[numeric_cols], as.numeric)
  result_tol[c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr")] <-
    lapply(result_tol[c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr")], as.numeric)
  
  rownames(result_tol) <- NULL
  result_tol
}

result_tol_rho <- extract_rho_table_from_conlist(
  con_list = con_list,
  meta_df = meta_df
)
##############################  table of rho(A_hat, A) 

method_type <- c(
  "Pro.iter", "HOPE", "CC-ISO",
  "Pro.init", "cPCA", "RP-PCA", "RCP"
)

method_type_mean <- paste0(method_type, ".mean")
method_type_sd   <- paste0(method_type, ".sd")

delta.fl_tol <- c(0.25, 0.75)
factor.loading <- 1

con_tol <- list()
kk <- 1

for (factor.corr in factor.corr_tol) {
  
  for (delta.fl in delta.fl_tol) {
    
    for (ALPHA in ALPHA_tol) {
      
      for (n in n_tol) {
        
        place <- (result_tol_rho$n == n) &
          (result_tol_rho$ALPHA == ALPHA) &
          (result_tol_rho$delta.fl == delta.fl) &
          (result_tol_rho$factor.loading == factor.loading) &
          (result_tol_rho$factor.corr == factor.corr)
        
        result_i <- result_tol_rho[place, , drop = FALSE]
        
        if (nrow(result_i) == 0) next
        
        mse_A <- paste0(
          sprintf("%0.2f", as.numeric(result_i[1, method_type_mean]) * 100),
          "(",
          sprintf("%0.2f", as.numeric(result_i[1, method_type_sd]) * 100),
          ")"
        )
        
        names(mse_A) <- method_type
        
        con_i <- data.frame(
          result_i[1, c("factor.corr", "factor.loading", "delta.fl", "ALPHA", "n")],
          t(mse_A),
          check.names = FALSE
        )
        
        con_tol[[kk]] <- con_i
        kk <- kk + 1
      }
    }
  }
}

con_tol <- do.call(rbind, con_tol)
rownames(con_tol) <- NULL


write.csv(con_tol, "table_rho_Ahat_A_heavytailed.csv", row.names = FALSE)

