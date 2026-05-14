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


n_tol = c(400,800,1600,3200)           #sample size   
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
factor.corr    = 0
ALPHA          = 0.6
n              = 1600

# helper functions
rho2.loss.list <- HDTSA:::rho2.loss.list
fnorm <- HDTSA:::fnorm
rho2.f.loss <- HDTSA:::rho2.f.loss
aij.debias.iter <- HDTSA:::aij.debias.iter
cov.aij.debias.iter <- HDTSA:::cov.aij.debias.iter
cov.aij.debias.iter.est <- HDTSA:::cov.aij.debias.iter.est
DGP.TCP <- HDTSA:::DGP.TCP
tensor.est.xi <- HDTSA:::tensor.est.xi
HDTTS.CP.est <- HDTSA:::HDTTS.CP.est
RP.xi.sel <- HDTSA:::RP.xi.sel

Rolling_fun2 <- function(x,
                         n,
                         ALPHA,
                         delta.fl,
                         factor.loading,
                         factor.corr,
                         seed_num,
                         i = 1,
                         j = 1,
                         K = 10,
                         max_try = 50) {
  
  set.seed(seed_num + x)
  
  D <- c(20, 20)
  m <- length(D)
  r <- 3
  beta <- c(0.8, 0.75, 0.7)
  ar.coef <- as.list(beta)
  w <- rep(15, r)
  
  data <- DGP.TCP(
    n = n,
    m = m,
    D = D,
    r = r,
    w = w,
    ar.coef = ar.coef,
    factor.loading = factor.loading,
    factor.corr = factor.corr,
    alpha = ALPHA,
    delta = delta.fl,
    par_E = 1
  )
  
  A <- data$A
  Y <- data$Y
  f <- data$f
  
  make_control <- function() {
    list(
      lag.k.dpi = K,
      threshold = TRUE,
      delta = NULL,
      delta2 = rep(1, m),
      ratio.type = "log",
      random.projection = FALSE,
      iter.max = 20,
      eps = 1e-5,
      grid.num = 50,
      delta.max = 0.1,
      print.eps = FALSE,
      iter.lag = 1,
      all.put = TRUE,
      A = A,
      component = NULL
    )
  }
  
  out_na <- list(
    true_var_iter = c(h1 = NA_real_, h2 = NA_real_),
    est_var_iter  = c(h1 = NA_real_, h2 = NA_real_)
  )
  
  for (trying_time in seq_len(max_try)) {
    
    fit <- tryCatch({
      
      xi_chang <- tensor.est.xi(Y)
      
      res_rank0 <- HDTTS.CP.est(
        Y = Y,
        xi = xi_chang,
        K = K,
        Ratio.type = "log"
      )
      
      r_hat_first <- res_rank0$r.hat
      
      if (
        length(r_hat_first) != 1 ||
        is.na(r_hat_first) ||
        !is.finite(r_hat_first) ||
        r_hat_first < 1
      ) {
        stop("Invalid r_hat_first.")
      }
      
      xi_res <- RP.xi.sel(
        Y = Y,
        r_breve = max(2L, 2L * as.integer(round(r_hat_first))),
        eps = 0.1,
        lag.k = K,
        Randomized.time = 50,
        A = A
      )
      
      CP_TTS(
        Y = Y,
        xi = xi_res$xi.sel,
        r = r,
        A.init = NULL,
        control.DPI = make_control()
      )
      
    }, error = function(e) {
      NULL
    })
    
    if (is.null(fit)) {
      next
    }
    
    A.init <- fit$A.init
    A.iter <- fit$A.hat
    f_hat <- as.matrix(fit$f.hat)
    
    if (is.null(A.init) || is.null(A.iter) || is.null(f_hat)) {
      next
    }
    
    if (max(rho2.loss.list(A.init, A)) >= 0.3) {
      next
    }
    
    place.iter <- tryCatch(
      which.max(abs(cor(A[[j]][, i], A.iter[[j]]))),
      error = function(e) NA_integer_
    )
    
    if (length(place.iter) != 1 || is.na(place.iter) || place.iter < 1) {
      next
    }
    
    aij.iter <- A.iter[[j]][, place.iter]
    d_j <- length(aij.iter)
    
    h1 <- c(1, rep(0, d_j - 1))
    h2 <- rep(1 / sqrt(d_j), d_j)
    
    true_var_h1 <- tryCatch(
      cov.aij.debias.iter(
        h = h1,
        A = A,
        i = i,
        j = j,
        f = f,
        COV.VECE = diag(prod(D)),
        w = w
      ),
      error = function(e) NA_real_
    )
    
    true_var_h2 <- tryCatch(
      cov.aij.debias.iter(
        h = h2,
        A = A,
        i = i,
        j = j,
        f = f,
        COV.VECE = diag(prod(D)),
        w = w
      ),
      error = function(e) NA_real_
    )
    
    est_var_h1 <- tryCatch(
      cov.aij.debias.iter.est(
        h = h1,
        i = place.iter,
        j = j,
        A = A.iter,
        f = f_hat,
        Y = Y
      ),
      error = function(e) NA_real_
    )
    
    est_var_h2 <- tryCatch(
      cov.aij.debias.iter.est(
        h = h2,
        i = place.iter,
        j = j,
        A = A.iter,
        f = f_hat,
        Y = Y
      ),
      error = function(e) NA_real_
    )
    
    out <- list(
      true_var_iter = c(
        h1 = true_var_h1,
        h2 = true_var_h2
      ),
      est_var_iter = c(
        h1 = est_var_h1,
        h2 = est_var_h2
      )
    )
    
    if (all(is.finite(unlist(out)))) {
      return(out)
    }
  }
  
  return(out_na)
}


Rolling_fun2_vec <- function(x,
                             n,
                             ALPHA,
                             delta.fl,
                             factor.loading,
                             factor.corr,
                             seed_num,
                             i = 1,
                             j = 1,
                             K = 10,
                             max_try = 50) {
  
  res <- Rolling_fun2(
    x = x,
    n = n,
    ALPHA = ALPHA,
    delta.fl = delta.fl,
    factor.loading = factor.loading,
    factor.corr = factor.corr,
    seed_num = seed_num,
    i = i,
    j = j,
    K = K,
    max_try = max_try
  )
  
  c(
    true_var_h1 = res$true_var_iter["h1"],
    true_var_h2 = res$true_var_iter["h2"],
    est_var_h1  = res$est_var_iter["h1"],
    est_var_h2  = res$est_var_iter["h2"]
  )
}

 

rep <- 2000
seed_num <- 123456

con_list <- list()
meta_list <- list()
hgl <- 1

for (factor.corr in factor.corr_tol) {
  
  for (factor.loading.num in factor.loading_tol) {
    
    factor.loading <- paste0("sparse-random-corr", factor.loading.num)
    
    for (delta.fl in delta.fl_tol) {
      
      for (ALPHA in ALPHA_tol) {
        
        for (n in n_tol) {
          
          con <- foreach(
            x = 1:rep,
            .combine = "rbind",
            .packages = c(
              "HDTSA", "R.utils", "tensor", "rTensor", "MASS",
              "pracma", "Rcpp", "RcppArmadillo", "RcppEigen",
              "vars", "base", "utils"
            ),
            .export = c(
              "CP_TTS",
              "Rolling_fun2",
              "Rolling_fun2_vec",
              "rho2.loss.list",
              "DGP.TCP",
              "tensor.est.xi",
              "HDTTS.CP.est",
              "RP.xi.sel",
              "cov.aij.debias.iter",
              "cov.aij.debias.iter.est"
            ),
            .errorhandling = "pass",
            .verbose = FALSE
          ) %dopar% Rolling_fun2_vec(
            x = x,
            n = n,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading = factor.loading,
            factor.corr = factor.corr,
            seed_num = seed_num,
            i = 1,
            j = 1,
            K = 10,
            max_try = 50
          )
          
          con <- as.data.frame(con)
          colnames(con) <- c(
            "true_var_h1",
            "true_var_h2",
            "est_var_h1",
            "est_var_h2"
          )
          
          con_list[[hgl]] <- con
          
          meta_list[[hgl]] <- data.frame(
            n = n,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading.num = factor.loading.num,
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
          print(colMeans(con, na.rm = TRUE))
          cat("\n")
          
          save.image("temp_variance_estimation.RData")
          
          seed_num <- seed_num + rep
          hgl <- hgl + 1
        }
      }
    }
  }
}

direction <- paste0(
  "simulation_variance_",
  format(Sys.time(), "%Y%m%d-%H-%M-%S"),
  ".RData"
)

save.image(direction)

if (exists("cl")) {
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
}


label.mat <- do.call(rbind, meta_list)
rownames(label.mat) <- NULL


rmse_iter_tol_name_fin <- NULL

for (factor.corr in factor.corr_tol) {
  
  for (factor.loading.num in factor.loading_tol) {
    
    for (delta.fl in delta.fl_tol) {
      
      for (ALPHA in ALPHA_tol) {
        
        rmse_h1_iter_tol <- c()
        rmse_h2_iter_tol <- c()
        
        for (n in n_tol) {
          
          place <- (label.mat$n == n) &
            (label.mat$ALPHA == ALPHA) &
            (label.mat$delta.fl == delta.fl) &
            (label.mat$factor.loading.num == factor.loading.num) &
            (label.mat$factor.corr == factor.corr)
          
          idx <- which(place)
          
          if (length(idx) != 1) {
            rmse_h1_iter_tol <- c(rmse_h1_iter_tol, NA_real_)
            rmse_h2_iter_tol <- c(rmse_h2_iter_tol, NA_real_)
            next
          }
          
          result_i <- as.data.frame(con_list[[idx]])
          
          colnames(result_i) <- c(
            "true_var_h1",
            "true_var_h2",
            "est_var_h1",
            "est_var_h2"
          )
          
          rmse_h1_iter <- sqrt(
            mean(
              (result_i$true_var_h1 - result_i$est_var_h1)^2,
              na.rm = TRUE
            )
          )
          
          rmse_h2_iter <- sqrt(
            mean(
              (result_i$true_var_h2 - result_i$est_var_h2)^2,
              na.rm = TRUE
            )
          )
          
          rmse_h1_iter_tol <- c(rmse_h1_iter_tol, rmse_h1_iter)
          rmse_h2_iter_tol <- c(rmse_h2_iter_tol, rmse_h2_iter)
        }
        
        rmse_iter_tol <- rbind(
          rmse_h1_iter_tol,
          rmse_h2_iter_tol
        )
        
        colnames(rmse_iter_tol) <- paste0("n = ", n_tol)
        
        rmse_iter_tol_name <- cbind(
          factor.corr = factor.corr,
          factor.loading.num = factor.loading.num,
          delta.fl = delta.fl,
          ALPHA = ALPHA,
          h = c("h1", "h2"),
          rmse_iter_tol
        )
        
        rmse_iter_tol_name_fin <- rbind(
          rmse_iter_tol_name_fin,
          rmse_iter_tol_name
        )
      }
    }
  }
}

rmse_iter_tol_name_fin <- as.data.frame(
  rmse_iter_tol_name_fin,
  stringsAsFactors = FALSE
)

rownames(rmse_iter_tol_name_fin) <- NULL

write.csv(
  rmse_iter_tol_name_fin,
  "rmse_iter_variance.csv",
  row.names = FALSE
)
