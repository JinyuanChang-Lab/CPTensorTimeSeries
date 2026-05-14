library(pracma)
library(ggplot2)
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

#helper functions
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


Rolling_fun2 <- function(x, n, ALPHA, delta.fl, factor.loading, factor.corr,
                         seed_num) {
  
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
    alpha = ALPHA, delta = delta.fl, par_E = 1
  )
  
  A <- data$A
  Y <- data$Y
  f <- data$f
  C <- data$C
  
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
  
  run_bmp <- function(Y_input, xi_input, Rank_input = NULL,
                      Ratio.type = "log", A.init = NULL,
                      Threshold = TRUE, delta2 = rep(1, m),
                      iter_max = 20, eps = 1e-5,
                      print.eps = FALSE, iter_lag = 1,
                      all.put = TRUE, A = NULL, Component = NULL) {
    
    CP_TTS(
      Y = Y_input,
      xi = xi_input,
      r = Rank_input,
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
  
  ############################
  ## First-stage proposed method
  ############################
  
  iter.num <- 20
  gg <- 0
  
  while (iter.num > 10 && gg < 50) {
    xi.chang <- tensor.est.xi(Y)
    
    res.only.used.rank <- HDTTS.CP.est(
      Y = Y,
      xi = xi.chang,
      K = K,
      Ratio.type = "log"
    )
    
    r_hat_first <- res.only.used.rank$r.hat
    
    xi.res <- RP.xi.sel(
      Y,
      r_breve = 2 * r_hat_first,
      eps = 0.1,
      lag.k = K,
      Randomized.time = 50
    )
    
    con.BMP.1 <- run_bmp(
      Y_input = Y,
      xi_input = xi.res$xi.sel,
      Rank_input = NULL,
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
    
    iter.num <- con.BMP.1$iter.step
    gg <- gg + 1
  }
  
  ############################
  ## First-stage results
  ############################
  
  f.pro.iter.1 <- con.BMP.1$f.hat
  f.pro.inl.1 <- con.BMP.1$f.hat.inl
  
  A.iter.pro.1 <- con.BMP.1$A.hat
  A.init.pro.1 <- con.BMP.1$A.init
  
  rank.est.1 <- get_rank(con.BMP.1)
  rank.est.1step <- rank.est.1
  
  error_iter_1step <- max(
    rho2.loss.list(A.iter.pro.1, A)
  )
  
  error_inl_1step <- max(
    rho2.loss.list(A.init.pro.1, A)
  )
  
  ############################
  ## Second-stage residual method
  ############################
  
  residual.step1 <- cp_residuals_general(
    Y,
    f.pro.iter.1,
    A.iter.pro.1
  )
  
  Y.resid.1 <- residual.step1$residual
  
  xi.chang.resid.1 <- tensor.est.xi(Y.resid.1)
  
  res.rank.1 <- HDTTS.CP.est(
    Y = Y.resid.1,
    xi = xi.chang.resid.1,
    K = K,
    Ratio.type = "log",
    Threshold = TRUE
  )
  
  con.BMP.2 <- run_bmp(
    Y_input = Y.resid.1,
    xi_input = xi.chang.resid.1,
    Rank_input = min(res.rank.1$r.hat, min(D) / 4),
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
  
  ############################
  ## Second-stage results
  ############################
  
  f.pro.iter.2 <- con.BMP.2$f.hat
  f.pro.inl.2 <- con.BMP.2$f.hat.inl
  
  A.iter.pro.2 <- con.BMP.2$A.hat
  A.init.pro.2 <- con.BMP.2$A.init
  
  rank.est.2 <- get_rank(con.BMP.2)
  
  ############################
  ## Combine first two stages
  ############################
  
  A.iter.pro.2step <- Map(
    cbind,
    A.iter.pro.1,
    A.iter.pro.2
  )
  
  A.init.pro.2step <- Map(
    cbind,
    A.init.pro.1,
    A.init.pro.2
  )
  
  error_iter_2step <- max(
    rho2.loss.list(A.iter.pro.2step, A)
  )
  
  error_inl_2step <- max(
    rho2.loss.list(A.init.pro.2step, A)
  )
  
  rank.est.2step <- rank.est.1 + rank.est.2
  
  ############################
  ## Final summary
  ############################
  
  res <- c(
    error_iter_1step,
    error_iter_2step,
    error_inl_1step,
    error_inl_2step,
    rank.est.1step,
    rank.est.2step
  )
  
  return(res)
}

#res = Rolling_fun2(1, n, ALPHA, delta.fl, factor.loading, factor.corr, seed_num = 1)

 


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
          
          #con <- rbind(res,res,res,res,res)  
          
          con_list[[hgl]] <- con
          
          
          
          
          meta_list[[hgl]] <- data.frame(
            n = n,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading = factor.loading.num,
            factor.corr = factor.corr
          )
          
          
          save.image("temp_estimation_WrongR_factornumber.RData")
          
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


direction = paste0("simulation_wrongR_factornumber_",format(Sys.time(), "%Y%m%d-%H-%M-%S"),".RData")

save.image(direction)

meta_df <- do.call(rbind, meta_list)
rownames(meta_df) <- NULL

############################## 2. table rank ##############################
########################################
## rank estimation: two-stage methods
########################################

## columns in con_list from Rolling_fun2:
## 5 = 1-stage estimated rank
## 6 = 2-stage estimated rank

r <- 3

rank_method_name <- c(
  "1-stage",
  "2-stage"
)

rank_col_id <- c(5, 6)

rank_tol_stage <- list()
kk <- 1

for (ggg in seq_along(con_list)) {
  
  con_i <- con_list[[ggg]]
  meta_i <- meta_list[[ggg]]
  
  if (is.null(con_i) || length(con_i) == 0) {
    next
  }
  
  con_i <- as.matrix(con_i)
  suppressWarnings(storage.mode(con_i) <- "numeric")
  
  for (jj in seq_along(rank_col_id)) {
    
    rank_i <- con_i[, rank_col_id[jj]]
    
    valid_n <- sum(!is.na(rank_i))
    
    if (valid_n == 0) {
      P_less_r <- NA
      P_geq_r <- NA
    } else {
      P_less_r <- 100 * sum(rank_i < r, na.rm = TRUE) / valid_n
      P_geq_r <- 100 * sum(rank_i >= r, na.rm = TRUE) / valid_n
    }
    
    rank_tol_stage[[kk]] <- data.frame(
      n = meta_i$n,
      ALPHA = meta_i$ALPHA,
      delta.fl = meta_i$delta.fl,
      factor.loading = meta_i$factor.loading,
      factor.corr = meta_i$factor.corr,
      method = rank_method_name[jj],
      P_less_r = round(P_less_r, 2),
      P_geq_r = round(P_geq_r, 2)
    )
    
    kk <- kk + 1
  }
}

rank_tol_stage <- do.call(rbind, rank_tol_stage)
rownames(rank_tol_stage) <- NULL


########################################
## reshape into wide table
########################################

rank_tol2 <- rank_tol_stage

rank_tol2$rho <- rank_tol2$factor.corr
rank_tol2$phi <- rank_tol2$delta.fl
rank_tol2$s   <- rank_tol2$ALPHA


########################################
## helper function
########################################

reshape_one_method <- function(df, method_name, prefix_name) {
  
  df_sub <- df[
    df$method == method_name,
    c("rho", "phi", "s", "n", "factor.loading", "P_less_r", "P_geq_r")
  ]
  
  colnames(df_sub) <- c(
    "rho",
    "phi",
    "s",
    "n",
    "factor.loading",
    paste0(prefix_name, "_P(<r)"),
    paste0(prefix_name, "_P(>=r)")
  )
  
  df_sub
}


########################################
## create separate tables
########################################

rank_1stage <- reshape_one_method(
  rank_tol2,
  "1-stage",
  "1-stage"
)

rank_2stage <- reshape_one_method(
  rank_tol2,
  "2-stage",
  "2-stage"
)


########################################
## merge side by side
########################################

rank_wide_stage <- merge(
  rank_1stage,
  rank_2stage,
  by = c("rho", "phi", "s", "n", "factor.loading"),
  all = TRUE
)


########################################
## sort rows
########################################

rank_wide_stage <- rank_wide_stage[
  order(
    rank_wide_stage$rho,
    rank_wide_stage$factor.loading,
    rank_wide_stage$phi,
    rank_wide_stage$s,
    rank_wide_stage$n
  ),
]

rownames(rank_wide_stage) <- NULL


########################################
## export
########################################

write.csv(
  rank_wide_stage,
  "table_rank_estimation_stage.csv",
  row.names = FALSE
)