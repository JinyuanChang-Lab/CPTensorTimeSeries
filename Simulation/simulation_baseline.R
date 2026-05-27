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


#helper functions
rho2.loss.list <- HDTSA:::rho2.loss.list
fnorm <- HDTSA:::fnorm
rho2.f.loss <- HDTSA:::rho2.f.loss
aij.debias.iter <- HDTSA:::aij.debias.iter
cov.aij.debias.iter <- HDTSA:::cov.aij.debias.iter
cov.aij.debias.iter.est <-  HDTSA:::cov.aij.debias.iter.est
cov.aij.debias.iter.longrun.est <- HDTSA:::cov.aij.debias.iter.longrun.est
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
    alpha = ALPHA, delta = delta.fl, par_E = 1
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

  rank.est <- con.BMP$r
  rank.est.NT <- con.BMP.classical$r

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

    v.iter.1 <- cov.aij.debias.iter.longrun.est(
      h1, i = place.iter, j = j, A = A.iter, f = f.iter, Y = Y
    )

    v.iter.2 <- cov.aij.debias.iter.longrun.est(
      h2, i = place.iter, j = j, A = A.iter, f = f.iter, Y = Y
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

#res = Rolling_fun2(x = 1,n = n, ALPHA = ALPHA,delta.fl = delta.fl,factor.loading = factor.loading, factor.corr = factor.corr,seed_num = 1)


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

          #con <- rbind(res,res,res,res,res)

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

direction = paste0("simulation_estimation_corr1_log_",format(Sys.time(), "%Y%m%d-%H-%M-%S"),".RData")

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


write.csv(con_tol, "table_rho_Ahat_A.csv", row.names = FALSE)




######################################## 3. rank estimation ########################################

delta.fl_tol <- c(0.25, 0.75)
factor.loading <- 1

rank_method_name <- c("log-ER", "ER", "Unfold-ER")
rank_col_id <- c(13,14,16)

rank_tol <- list()
kk <- 1
ggg <- 1

for (factor.corr in factor.corr_tol) {

  for (delta.fl in delta.fl_tol) {

    for (ALPHA in ALPHA_tol) {

      for (n in n_tol) {

        con_i <- con_list[[ggg]]

        if (is.null(con_i) || length(con_i) == 0) {
          ggg <- ggg + 1
          next
        }

        con_i <- as.matrix(con_i)
        storage.mode(con_i) <- "numeric"

        for (jj in seq_along(rank_col_id)) {

          rank_i <- con_i[, rank_col_id[jj]]

          res_i <- 100 * c(
            "Under"   = sum(rank_i < r, na.rm = TRUE),
            "Correct" = sum(rank_i == r, na.rm = TRUE),
            "Over"    = sum(rank_i > r, na.rm = TRUE)
          ) / sum(!is.na(rank_i))

          rank_k <- data.frame(
            factor.corr = factor.corr,
            factor.loading = factor.loading,
            delta.fl = delta.fl,
            ALPHA = ALPHA,
            n = n,
            method = rank_method_name[jj],
            Under = round(res_i["Under"], 2),
            Correct = round(res_i["Correct"], 2),
            Over = round(res_i["Over"], 2)
          )

          rank_tol[[kk]] <- rank_k
          kk <- kk + 1
        }

        ggg <- ggg + 1
      }
    }
  }
}

rank_tol <- do.call(rbind, rank_tol)
rownames(rank_tol) <- NULL


# rank_tol is assumed to be the long-format table created previously
# with columns:
# factor.corr, factor.loading, delta.fl, ALPHA, n, method, Under, Correct, Over

rank_tol2 <- rank_tol

# keep only the three methods
rank_tol2 <- rank_tol2[rank_tol2$method %in% c("log-ER", "ER", "Unfold-ER"), ]

# rename columns for display
rank_tol2$rho <- rank_tol2$factor.corr
rank_tol2$phi <- rank_tol2$delta.fl
rank_tol2$s   <- rank_tol2$ALPHA

# helper function: reshape one method into wide columns
reshape_one_method <- function(df, method_name, prefix_name) {
  df_sub <- df[df$method == method_name,
               c("rho", "phi", "s", "n", "Under", "Correct", "Over")]

  colnames(df_sub) <- c(
    "rho", "phi", "s", "n",
    paste0(prefix_name, "_P(<r)"),
    paste0(prefix_name, "_P(=r)"),
    paste0(prefix_name, "_P(>r)")
  )

  df_sub
}

rank_logER   <- reshape_one_method(rank_tol2, "log-ER", "Log-ER")
rank_ER      <- reshape_one_method(rank_tol2, "ER", "ER")
rank_unfold  <- reshape_one_method(rank_tol2, "Unfold-ER", "Unfold-ER")

# merge side by side
rank_wide <- merge(rank_logER, rank_ER, by = c("rho", "phi", "s", "n"), all = TRUE)
rank_wide <- merge(rank_wide, rank_unfold, by = c("rho", "phi", "s", "n"), all = TRUE)

# sort rows
rank_wide <- rank_wide[order(rank_wide$rho, rank_wide$phi, rank_wide$s, rank_wide$n), ]

rownames(rank_wide) <- NULL


write.csv(rank_wide, "table_rank_estimation.csv", row.names = FALSE)

######################################## 4. iterative steps ########################################

# The updated Rolling_fun2 returns:
# [1:22]     scalar summaries
# [23:169]   iter_tol = 21 * 7
# [170:176]  f_loss
# [177:183]  CP_loss

iter_start <- 23
iter_block_len <- 21
iter_group_num_total <- 7

# Only keep the first five groups for plotting
# 1: Pro.iter
# 2: HOPE
# 3: CC-ISO
# 4: Pro.iter + cPCA
# 5: Pro.iter + RP-PCA
iter_group_num_keep <- 5

iter_meta_tol <- vector()
result_iter_tol <- vector()

ggg <- 1

for (factor.corr in factor.corr_tol) {
  for (factor.loading.num in factor.loading_tol) {
    for (delta.fl in delta.fl_tol) {
      for (ALPHA in ALPHA_tol) {
        for (n in n_tol) {

          con_i <- con_list[[ggg]]

          if (is.null(con_i) || length(con_i) == 0) {
            ggg <- ggg + 1
            next
          }

          con_i <- as.matrix(con_i)
          storage.mode(con_i) <- "numeric"

          iter_end <- iter_start + iter_block_len * iter_group_num_total - 1
          con_iter_error <- con_i[, iter_start:iter_end, drop = FALSE]

          mean.iter.error.k <- colMeans(con_iter_error, na.rm = TRUE)

          # keep only the first 5 groups = 105 columns
          mean.iter.error.k <- mean.iter.error.k[1:(iter_block_len * iter_group_num_keep)]

          result_iter_tol <- rbind(result_iter_tol, mean.iter.error.k)

          iter_meta_k <- c(
            n = n,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading = factor.loading.num,
            factor.corr = factor.corr
          )

          iter_meta_tol <- rbind(iter_meta_tol, iter_meta_k)

          ggg <- ggg + 1
        }
      }
    }
  }
}

result_iter_tol <- as.data.frame(result_iter_tol)
iter_meta_tol <- as.data.frame(iter_meta_tol)

colnames(iter_meta_tol) <- c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr")
colnames(result_iter_tol) <- 1:ncol(result_iter_tol)

######################################## split five iterative groups
result_iter_max_pro      <- result_iter_tol[,   1:21]
result_iter_max_han      <- result_iter_tol[,  22:42]
result_iter_max_chen     <- result_iter_tol[,  43:63]
result_iter_max_pro_han  <- result_iter_tol[,  64:84]
result_iter_max_pro_chen <- result_iter_tol[,  85:105]

colnames(result_iter_max_pro)      <- 0:20
colnames(result_iter_max_han)      <- 0:20
colnames(result_iter_max_chen)     <- 0:20
colnames(result_iter_max_pro_han)  <- 0:20
colnames(result_iter_max_pro_chen) <- 0:20

######################################## build long-format plotting data
con_iter_pro <- con_iter_han <- con_iter_chen <- con_iter_pro_han <- con_iter_pro_chen <- vector()

for (jj in 1:NROW(iter_meta_tol)) {

  con_iter_pro_i <- data.frame(iter_meta_tol[jj, ], step = 0:20, mse = as.numeric(result_iter_max_pro[jj, ]))
  con_iter_han_i <- data.frame(iter_meta_tol[jj, ], step = 0:20, mse = as.numeric(result_iter_max_han[jj, ]))
  con_iter_chen_i <- data.frame(iter_meta_tol[jj, ], step = 0:20, mse = as.numeric(result_iter_max_chen[jj, ]))

  con_iter_pro_han_i <- data.frame(iter_meta_tol[jj, ], step = 0:20, mse = as.numeric(result_iter_max_pro_han[jj, ]))
  con_iter_pro_chen_i <- data.frame(iter_meta_tol[jj, ], step = 0:20, mse = as.numeric(result_iter_max_pro_chen[jj, ]))

  con_iter_pro <- rbind(con_iter_pro, con_iter_pro_i)
  con_iter_han <- rbind(con_iter_han, con_iter_han_i)
  con_iter_chen <- rbind(con_iter_chen, con_iter_chen_i)
  con_iter_pro_han <- rbind(con_iter_pro_han, con_iter_pro_han_i)
  con_iter_pro_chen <- rbind(con_iter_pro_chen, con_iter_pro_chen_i)
}

con_iter <- rbind(
  cbind(method = "Pro.iter", con_iter_pro),
  cbind(method = "Pro.iter+cPCA", con_iter_pro_han),
  cbind(method = "Pro.iter+RP-PCA", con_iter_pro_chen),
  cbind(method = "HOPE", con_iter_han),
  cbind(method = "CC-ISO", con_iter_chen)
)

con_iter$method <- factor(
  con_iter$method,
  levels = c("Pro.iter", "Pro.iter+cPCA", "Pro.iter+RP-PCA", "HOPE", "CC-ISO")
)

######################################## plot iterative-step figures
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

color_type <- c( "black", "darkgreen", "darkorange", "red", "blue")
line_type  <- c("solid", "solid", "solid", "solid", "solid")
shape_type <- c( 15, 2, 1, 17, 16)

plot_list <- list()
kk <- 1

delta.fl_tol <- c(0.25, 0.75)
n <- 400
factor.loading <- 1

for (factor.corr in factor.corr_tol) {
  for (delta.fl in delta.fl_tol) {
    for (ALPHA in ALPHA_tol) {

      place <- (con_iter$n == n) &
        (con_iter$ALPHA == ALPHA) &
        (con_iter$delta.fl == delta.fl) &
        (con_iter$factor.loading == factor.loading) &
        (con_iter$factor.corr == factor.corr)

      result_i <- con_iter[place, ]

      n1 <- substitute(
        expression(rho == value1 * "," ~ phi == value2 * "," ~ s == value3),
        list(value1 = factor.corr, value2 = delta.fl, value3 = ALPHA)
      )

      a1A <- ggplot(result_i, aes(x = step, y = mse, col = method, linetype = method, shape = method)) +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_shape_manual(values = shape_type) +
        ggtitle(eval(n1)) +
        xlab("iterative step") +
        ylab("estimation error") +
        scale_x_continuous(breaks = c(0, 5, 10, 15, 20)) +
        theme(legend.position = "none", legend.direction = "horizontal") +
        scale_color_manual(values = color_type) +
        scale_linetype_manual(values = line_type) +
        theme(
          axis.line = element_line(linewidth = 0.5, linetype = "solid"),
          legend.text = element_text(family = "serif"),
          panel.background = element_rect(fill = "gray100", colour = alpha("black", 0.5)),
          legend.key = element_rect(fill = alpha("white", 1)),
          panel.grid.major = element_line(linewidth = 0.4, linetype = "dashed", colour = "gray92"),
          legend.background = element_rect(colour = "black"),
          legend.title = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 15, face = "bold")
        )

      plot_list[[kk]] <- a1A
      kk <- kk + 1
    }
  }
}

######################################## save figure
write.csv(con_iter, "iterative_steps_long.csv", row.names = FALSE)

graphics.off()
while (!is.null(dev.list())) dev.off()

pdf("iterative_steps_selected_methods.pdf", width = 12, height = 10)
pig <- arrangeGrob(grobs = plot_list, ncol = 3)
grid.draw(pig)
dev.off()

######################################## 5. normality ########################################

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

## con_list comes from the updated Rolling_fun2
## columns 17:22 are mse.dis
## here we only keep:
## 17 = Pro.iter long-run variance h1
## 18 = Pro.iter long-run variance h2
## 19 = Pro.iter estimated variance h1
## 20 = Pro.iter estimated variance h2

normality_names <- c(
  "Pro.iter_longrun_h1",
  "Pro.iter_longrun_h2",
  "Pro.iter_estimated_h1",
  "Pro.iter_estimated_h2"
)

############################## 4.1 extract normality results

con_normality_list <- vector("list", length(con_list))

for (ii in seq_along(con_list)) {
  con_i <- as.matrix(con_list[[ii]])
  storage.mode(con_i) <- "numeric"

  con_normality_i <- con_i[, 17:20, drop = FALSE]
  colnames(con_normality_i) <- normality_names

  con_normality_list[[ii]] <- con_normality_i
}

############################## 4.2 build label table
label.mat <- NULL
ww <- 1

for (factor.corr in factor.corr_tol) {
  for (factor.loading.num in factor.loading_tol) {
    for (delta.fl in delta.fl_tol) {
      for (ALPHA in ALPHA_tol) {
        for (n in n_tol) {
          label.mat <- rbind(
            label.mat,
            data.frame(
              n = n,
              ALPHA = ALPHA,
              delta.fl = delta.fl,
              factor.loading.num = factor.loading.num,
              factor.corr = factor.corr,
              label = ww
            )
          )
          ww <- ww + 1
        }
      }
    }
  }
}

rownames(label.mat) <- NULL

##############################   plotting function

make_hist_plot_normality <- function(sample_vec, factor.corr, delta.fl, ALPHA) {

  sample_vec <- sample_vec[is.finite(sample_vec)]
  sample_vec <- sample_vec[!is.na(sample_vec)]

  ## keep the same trimming rule as your previous code
  sample_vec <- sample_vec[abs(sample_vec) < qnorm(0.9995)]

  df <- data.frame(sample = sample_vec)

  title_expr <- substitute(
    expression(rho == value1 * "," ~ phi == value2 * "," ~ s == value3),
    list(value1 = factor.corr, value2 = delta.fl, value3 = ALPHA)
  )

  p <- ggplot(df, aes(x = sample)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = 0.2,
      color = "black",
      fill = "lightblue",
      alpha = 0.7
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = 1),
      color = "red",
      linewidth = 1
    ) +
    ggtitle(eval(title_expr)) +
    labs(x = "Value", y = "Density") +
    theme(
      axis.line = element_line(linewidth = 0.5, linetype = "solid"),
      legend.text = element_text(family = "serif"),
      panel.background = element_rect(fill = "gray100", colour = alpha("black", 0.5)),
      legend.key = element_rect(fill = alpha("white", 1)),
      panel.grid.major = element_line(linewidth = 0.4, linetype = "dashed", colour = "gray92"),
      legend.background = element_rect(colour = "black"),
      legend.title = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 15, face = "bold")
    )

  return(p)
}

##############################  draw the four big figures

delta.fl_tol <- c(0.25, 0.75)
factor.loading.num <- 1
n_selected <- 400

for (aq in 1:4) {

  plot_list_hist <- list()
  kk <- 1

  for (factor.corr in factor.corr_tol) {
    for (delta.fl in delta.fl_tol) {
      for (ALPHA in ALPHA_tol) {

        place <- (label.mat$n == n_selected) &
          (label.mat$ALPHA == ALPHA) &
          (label.mat$delta.fl == delta.fl) &
          (label.mat$factor.loading.num == factor.loading.num) &
          (label.mat$factor.corr == factor.corr)

        idx <- which(place)

        if (length(idx) != 1) next

        result_i <-  con_normality_list[[idx]][, aq]
        result_i <- result_i[!is.na(result_i)]

        ph <- make_hist_plot_normality(
          sample_vec = result_i,
          factor.corr = factor.corr,
          delta.fl = delta.fl,
          ALPHA = ALPHA
        )

        plot_list_hist[[kk]] <- ph
        kk <- kk + 1
      }
    }
  }

  graphics.off()
  while (!is.null(dev.list())) dev.off()

  pdf(
    paste0(normality_names[aq], "_hist.pdf"),
    width = 12,
    height = 10
  )

  pig <- arrangeGrob(grobs = plot_list_hist, ncol = 3)
  grid.draw(pig)
  dev.off()

}


############################## factor/common component summary from con_list

method_type <- c(
  "Pro.iter", "HOPE", "CC-ISO",
  "Pro.init", "cPCA", "RP-PCA", "RCP"
)

############################## helper function

extract_summary_table_from_conlist <- function(con_list, meta_df, col_range, method_type) {

  result_tol <- vector("list", length(con_list))

  for (ii in seq_along(con_list)) {

    con <- con_list[[ii]]

    if (is.null(con) || length(con) == 0) next

    con <- as.matrix(con)
    storage.mode(con) <- "numeric"

    con_error <- con[, col_range, drop = FALSE]
    colnames(con_error) <- method_type

    mean.k <- colMeans(con_error, na.rm = TRUE)
    sd.k   <- apply(con_error, 2, sd, na.rm = TRUE)

    names(mean.k) <- paste0(method_type, ".mean")
    names(sd.k)   <- paste0(method_type, ".sd")

    result_tol[[ii]] <- c(
      meta_df[ii, ],
      mean.k,
      sd.k
    )
  }

  result_tol <- do.call(rbind, result_tol)
  result_tol <- as.data.frame(result_tol, stringsAsFactors = FALSE)

  numeric_cols <- setdiff(
    colnames(result_tol),
    c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr")
  )

  result_tol[numeric_cols] <- lapply(result_tol[numeric_cols], as.numeric)
  result_tol[c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr")] <-
    lapply(result_tol[c("n", "ALPHA", "delta.fl", "factor.loading", "factor.corr")], as.numeric)

  rownames(result_tol) <- NULL
  result_tol
}

############################## helper function for final formatted table

build_display_table <- function(result_tol_obj,
                                method_type,
                                factor.corr_tol,
                                delta.fl_tol = c(0.25, 0.75),
                                ALPHA_tol,
                                n_tol,
                                factor.loading = 1,
                                scale_factor = 1) {

  method_type_mean <- paste0(method_type, ".mean")
  method_type_sd   <- paste0(method_type, ".sd")

  con_tol <- list()
  kk <- 1

  for (factor.corr in factor.corr_tol) {
    for (delta.fl in delta.fl_tol) {
      for (ALPHA in ALPHA_tol) {
        for (n in n_tol) {

          place <- (result_tol_obj$n == n) &
            (result_tol_obj$ALPHA == ALPHA) &
            (result_tol_obj$delta.fl == delta.fl) &
            (result_tol_obj$factor.loading == factor.loading) &
            (result_tol_obj$factor.corr == factor.corr)

          result_i <- result_tol_obj[place, , drop = FALSE]

          if (nrow(result_i) == 0) next

          mse_A <- paste0(
            sprintf("%0.2f", as.numeric(result_i[1, method_type_mean]) * scale_factor),
            "(",
            sprintf("%0.2f", as.numeric(result_i[1, method_type_sd]) * scale_factor),
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
  con_tol
}


############################## 6. common component estimation error ##############################

## Updated Rolling_fun2:
## CP_loss is in columns 177:183

result_tol_cp_new <- extract_summary_table_from_conlist(
  con_list = con_list,
  meta_df = meta_df,
  col_range = 177:183,
  method_type = method_type
)

con_tol_cp <- build_display_table(
  result_tol_obj = result_tol_cp_new,
  method_type = method_type,
  factor.corr_tol = factor.corr_tol,
  delta.fl_tol = c(0.25, 0.75),
  ALPHA_tol = ALPHA_tol,
  n_tol = n_tol,
  factor.loading = 1,
  scale_factor = 1
)


write.csv(con_tol_cp, "table_common_component_error.csv", row.names = FALSE)
