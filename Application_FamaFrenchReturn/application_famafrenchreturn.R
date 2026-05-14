


############################################################
## Rolling forecast comparison for tensor time series
############################################################

rm(list = ls())
gc()

## =========================
## 0. Packages
## =========================
library(pracma)
library(foreach)
library(doParallel)
library(vars)
library(HDTSA)

## =========================
## 1. User inputs
## =========================
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  
  frame_files <- vapply(sys.frames(), function(x) {
    if (!is.null(x$ofile)) x$ofile else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[length(frame_files)])))
  }
  
  getwd()
}

script_dir <- get_script_dir()
local_file <- function(...) file.path(script_dir, ...)

source_files <- c(
  local_file("tensor_cp_functions.R"),
  local_file("cp_cciso.R"),
  local_file("CP_functions.R")
)
cpp_source_file <- local_file("tools4cp.cpp")

load_custom_functions <- function() {
  for (f in source_files) source(f)
  if (!exists("sigmak", mode = "function")) {
    Rcpp::sourceCpp(cpp_source_file)
  }
  invisible(TRUE)
}

# ---- Forecast settings ----
t0 <- 456
sp <- 2
r  <- 1
lag_k <- 10
K_main <- 20
Ktilde_main <- 10

# ---- Parallel settings ----
ncore <- 120

# ---- Whether to report SD or SE in parentheses ----
dispersion_type <- "sd"   # choose from "sd" or "se"

# ---- Load data ----
data <- read.csv(local_file("100_Portfolios.CSV"), row.names = 1)
data <- as.matrix(data[1:696, ])

data[data == -99.99] <- 0

data_fit <- as.matrix(data[, -1]) -  as.matrix(data[, 1]) %*% t(rep(1, 100))

Y <- array(NA, dim = c(NROW(data_fit), 10, 10))
for (tt in 1:NROW(data_fit)) {
  for (ii in 1:10) {
    Y[tt, ii, ] <- data_fit[tt, (1 + 10 * (ii - 1)):(10 * ii)]
  }
}

stopifnot(exists("Y"))
stopifnot(length(dim(Y)) == 3)

n <- dim(Y)[1]
p <- dim(Y)[2]
q <- dim(Y)[3]


######plot Figure F10#######

pdf(local_file("Y_time_series_100.pdf"), width = 12, height = 6)

par(mfrow = c(10, 10))

par(mar = c(0.2, 0.2, 0.2, 0.2) + 0.1)

for (i in 1:10) {
  for (j in 1:10) {
    plot(
      Y[, i, j],
      xlab = NULL,
      ylab = NULL,
      main = NULL,
      yaxt = "n",
      xaxt = "n",
      ann = TRUE,
      type = "l"
    )
  }
} 

dev.off()




if (n <= t0) stop("Need n > t0 for rolling forecast.")

# ---- Utility functions ----
ffnorm <- function(X) sqrt(sum(X^2))
MAE_mat <- function(E) mean(abs(E))

safe_run <- function(expr, fallback = NULL) {
  out <- try(eval.parent(substitute(expr)), silent = TRUE)
  if (inherits(out, "try-error")) fallback else out
}

get_forecast_mat <- function(obj, step, slot_name) {
  if (is.null(obj)) return(NULL)
  x <- obj[[slot_name]]
  if (is.null(x)) return(NULL)
  if (length(x) < step) return(NULL)
  x[[step]]
}

calc_err <- function(pred, truth, scale = 100) {
  if (is.null(pred)) {
    return(c(mse = NA_real_, mae = NA_real_))
  }
  c(
    mse = ffnorm(pred - truth)^2 / scale,
    mae = MAE_mat(pred - truth)
  )
}

# ---- Methods to report ----
method_order <- c(
  "pro_iter",
  "pro_inl",
  "chen_iter",
  "chen_inl",
  "han_inl",
  "han_iter",
  "cp221",
  "topup22",
  "tipup22",
  "fac11",
  "chang1",
  "mar",
  "pca",
  "arma"
)

method_labels <- c(
  pro_iter  = "Pro.iter",
  pro_inl   = "Pro.init",
  chen_iter = "CC-ISO",
  chen_inl  = "RP-PCA",
  han_inl   = "cPCA",
  han_iter  = "HOPE",
  cp221     = "CP-JAD",
  topup22   = "TOPUP",
  tipup22   = "TIPUP",
  fac11     = "FAC",
  chang1    = "RCP",
  mar       = "MAR",
  pca       = "TS-PCA",
  arma      = "UniAR"
)

# ---- One block forecast evaluator ----
evaluate_block <- function(Y_train, Y_target, step_to_use, r, p, q,
                           lag_k, K_main, Ktilde_main, tt) {
  
  ny <- dim(Y_train)[1]
  m <- 2
  
  sigma0 <- sqrt(mean((Y_train - rep(1, ny) %o% apply(Y_train, c(2, 3), mean))^2))
  xi <- est.xi(Y_train, d_max = 1)$xi
  delta <- sigma0 * sqrt(log(p * q) / ny)
  Y1_train <- aperm(Y_train, c(2:(m + 1), 1))
  
  set.seed(2222 + tt)
  
  # Pro.iter / Pro.init
  # xi.res = Boostrap.xi.sel.v3(Y, r_breve = 2, eps = 0.1, lag.k = lag_k, Threshold = F, delta = NULL, Randomized.time = 50,
  #                             BMP = T, Ratio.type = "log", augmented = F,grid_delta1 = 50, all.out = T, A = NULL,print.eps = F)
  # xi.rp = xi.res$xi.sel
  
  fore_pro <- safe_run(
    Pro.forecast(
      Y = Y_train,
      forecast.step = sp,
      xi = NULL,
      Rank = r,
      lag.k = lag_k,
      Threshold = TRUE,
      delta = NULL,
      delta2 = rep(1, 2),
      A.inl = NULL,
      iter_lag = 1,
      seasonal = FALSE,
      diff = FALSE,
      BMP = TRUE,
      Ratio.type = "log",
      augmented = FALSE,
      iter_max = 100,
      eps = 1e-5,
      grid_delta1 = 50,
      grid.num = 50,
      delta_max = 0.1,
      all.out = TRUE,
      A = NULL,
      print.eps = FALSE,
      Random.Project = T
    ),
    fallback = NULL
  )
  
  # HOPE / cPCA with r = 1
  fore_han <- safe_run(
    Han.forecast(Y = Y_train, r = r, forecast.step = sp),
    fallback = NULL
  )
  
  # CC-ISO / RP-PCA with r = 1
  fore_chen <- safe_run(
    Chen.forecast(Y = Y_train, r = r, forecast.step = sp),
    fallback = NULL
  )
  
  # CP-JAD
  fore_cp221 <- safe_run(
    HDMTS.CP.Unified.forecast(
      Y_train,
      xi,
      Rank = list(d=2,d1=2,d2=1),
      K = K_main,
      Ktilde = Ktilde_main,
      solve.UV = NULL,
      thresh1 = FALSE,
      thresh2 = FALSE,
      thresh3 = FALSE,
      delta1 = delta,
      delta2 = delta,
      delta3 = delta,
      c1 = sigma0 / ny,
      c2 = sigma0 / ny,
      c3 = sigma0 / ny,
      forecast.step = 1:sp,
      All_out = FALSE
    ),
    fallback = NULL
  )
  
  # RCP
  fore_chang1 <- safe_run(
    Chang2022.forecast(
      Y_train,
      xi,
      Kmax = 10,
      Rank = 1,
      forecast.step = 1:sp
    ),
    fallback = NULL
  )
  
  # TOPUP
  fore_TOPUP22 <- safe_run(
    FAC.forecast(
      Y_train,
      r = c(2,2),
      h0 = 1,
      method = "TOPUP",
      forecast.step = 1:sp
    ),
    fallback = NULL
  )
  
  # TIPUP
  fore_TIPUP22 <- safe_run(
    FAC.forecast(
      Y_train,
      r = c(2,2),
      h0 = 1,
      method = "TIPUP",
      forecast.step = 1:sp
    ),
    fallback = NULL
  )
  
  # FAC
  fore_fac11 <- safe_run(
    Wang.forecast(
      Y_train,
      r = 1,
      h0 = 1,
      forecast.step = 1:sp
    ),
    fallback = NULL
  )
  
  # TS-PCA
  fore_pca <- safe_run(
    PCATS.forecast(Y_train, forecast.step = 1:sp),
    fallback = NULL
  )
  
  # UniAR
  fore_arma <- safe_run(
    UniARMA.forecast(Y_train, forecast.step = 1:sp),
    fallback = NULL
  )
  
  # Fallbacks
  if (is.null(fore_TOPUP22)) fore_TOPUP22 <- fore_arma
  if (is.null(fore_TIPUP22)) fore_TIPUP22 <- fore_arma
  if (is.null(fore_fac11))   fore_fac11   <- fore_arma
  if (is.null(fore_cp221))   fore_cp221   <- fore_arma
  
  # MAR(10,1)
  fore_mar <- safe_run(
    MAR.forecast(
      Y_train,
      method = "RRLSE",
      k1 = 10,
      k2 = 1,
      tol = 1e-2,
      forecast.step = 1:sp
    ),
    fallback = fore_arma
  )
  
  out <- list()
  
  out$pro_iter <- calc_err(
    get_forecast_mat(fore_pro, step_to_use, "Y.forecast.iter"),
    Y_target
  )
  
  out$pro_inl <- calc_err(
    get_forecast_mat(fore_pro, step_to_use, "Y.forecast.inl"),
    Y_target
  )
  
  out$chen_iter <- calc_err(
    get_forecast_mat(fore_chen, step_to_use, "Y.forecast.iter"),
    Y_target
  )
  
  out$chen_inl <- calc_err(
    get_forecast_mat(fore_chen, step_to_use, "Y.forecast.inl"),
    Y_target
  )
  
  out$han_inl <- calc_err(
    get_forecast_mat(fore_han, step_to_use, "Y.forecast.inl"),
    Y_target
  )
  
  out$han_iter <- calc_err(
    get_forecast_mat(fore_han, step_to_use, "Y.forecast.iter"),
    Y_target
  )
  
  out$cp221 <- calc_err(
    get_forecast_mat(fore_cp221, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$topup22 <- calc_err(
    get_forecast_mat(fore_TOPUP22, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$tipup22 <- calc_err(
    get_forecast_mat(fore_TIPUP22, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$fac11 <- calc_err(
    get_forecast_mat(fore_fac11, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$chang1 <- calc_err(
    get_forecast_mat(fore_chang1, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$mar <- calc_err(
    get_forecast_mat(fore_mar, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$pca <- calc_err(
    get_forecast_mat(fore_pca, step_to_use, "Y.forecast"),
    Y_target
  )
  
  out$arma <- calc_err(
    get_forecast_mat(fore_arma, step_to_use, "Y.forecast"),
    Y_target
  )
  
  mse_vec <- sapply(method_order, function(m) out[[m]]["mse"])
  mae_vec <- sapply(method_order, function(m) out[[m]]["mae"])
  
  c(mse_vec, mae_vec)
}

# ---- Rolling evaluator for one tt ----
Rolling_app_new <- function(tt, Y, t0, sp, r, p, q, lag_k, K_main, Ktilde_main) {
  
  load_custom_functions()
  
  # One-step
  Y_t1 <- Y[tt:(tt + t0 - 1), , ]
  Y_target <- Y[tt + t0, , ]
  
  res1 <- evaluate_block(
    Y_train = Y_t1,
    Y_target = Y_target,
    step_to_use = 1,
    r = r,
    p = p,
    q = q,
    lag_k = lag_k,
    K_main = K_main,
    Ktilde_main = Ktilde_main,
    tt = tt
  )
  
  names(res1) <- c(
    paste0("mse_", method_order, "_onestep"),
    paste0("mae_", method_order, "_onestep")
  )
  
  gc()
  
  # Two-step
  Y_t2 <- Y[tt:(tt + t0 - 2), , ]
  
  res2 <- evaluate_block(
    Y_train = Y_t2,
    Y_target = Y_target,
    step_to_use = 2,
    r = r,
    p = p,
    q = q,
    lag_k = lag_k,
    K_main = K_main,
    Ktilde_main = Ktilde_main,
    tt = tt
  )
  
  names(res2) <- c(
    paste0("mse_", method_order, "_twostep"),
    paste0("mae_", method_order, "_twostep")
  )
  
  c(res1, res2)
}

# ---- Source custom functions on master and workers ----
load_custom_functions()

cl <- parallel::makeCluster(ncore)
doParallel::registerDoParallel(cl)

parallel::clusterEvalQ(cl, {
  library(pracma)
  library(vars)
})

parallel::clusterExport(
  cl,
  varlist = c(
    "Y", "t0", "sp", "r", "p", "q", "lag_k", "K_main", "Ktilde_main",
    "source_files", "cpp_source_file", "load_custom_functions",
    "ffnorm", "MAE_mat", "safe_run", "get_forecast_mat", "calc_err",
    "method_order", "method_labels", "evaluate_block", "Rolling_app_new"
  ),
  envir = environment()
)

parallel::clusterEvalQ(cl, {
  load_custom_functions()
})

# ---- Run rolling forecasts in parallel ----
n_roll <- n - t0
cat("Number of rolling windows =", n_roll, "\n")

con <- foreach::foreach(
  tt = 1:n_roll,
  .combine = "rbind",
  .packages = c(
    "jointDiag", "HDTSA", "tensor", "rTensor", "tensorTS", "base",
    "pracma", "Rcpp", "RcppArmadillo", "RcppEigen", "vars",
    "forecast", "MASS"
  ),
  .errorhandling = "pass",
  .verbose = FALSE
) %dopar% {
  Rolling_app_new(
    tt = tt,
    Y = Y,
    t0 = t0,
    sp = sp,
    r = r,
    p = p,
    q = q,
    lag_k = lag_k,
    K_main = K_main,
    Ktilde_main = Ktilde_main
  )
}

save.image(local_file("temp_application_port.RData"))

parallel::stopCluster(cl)

if (is.list(con)) {
  stop("Some iterations returned non-standard outputs. Check workers.")
}

con <- as.data.frame(con)

# ---- Summaries ----
mse_cols <- grep("^mse_", names(con), value = TRUE)
mae_cols <- grep("^mae_", names(con), value = TRUE)

con[mse_cols] <- sqrt(con[mse_cols])

get_dispersion <- function(x, type = c("sd", "se")) {
  type <- match.arg(type)
  if (type == "sd") {
    return(sd(x, na.rm = TRUE))
  }
  if (type == "se") {
    return(sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  }
}

summarize_block <- function(prefix, horizon, dat, methods, disp_type = "sd") {
  cols <- paste0(prefix, "_", methods, "_", horizon)
  means <- sapply(cols, function(cc) mean(dat[[cc]], na.rm = TRUE))
  disps <- sapply(cols, function(cc) get_dispersion(dat[[cc]], type = disp_type))
  names(means) <- methods
  names(disps) <- methods
  list(mean = means, disp = disps)
}

rmse_one <- summarize_block("mse", "onestep", con, method_order, dispersion_type)
rmae_one <- summarize_block("mae", "onestep", con, method_order, dispersion_type)
rmse_two <- summarize_block("mse", "twostep", con, method_order, dispersion_type)
rmae_two <- summarize_block("mae", "twostep", con, method_order, dispersion_type)

summary_mean <- rbind(
  "rRMSE(one-step)" = rmse_one$mean,
  "rMAE(one-step)"  = rmae_one$mean,
  "rRMSE(two-step)" = rmse_two$mean,
  "rMAE(two-step)"  = rmae_two$mean
)

colnames(summary_mean) <- method_labels[colnames(summary_mean)]

summary_disp <- rbind(
  rmse_one$disp,
  rmae_one$disp,
  rmse_two$disp,
  rmae_two$disp
)

rownames(summary_disp) <- c(
  paste0(dispersion_type, "(rRMSE one-step)"),
  paste0(dispersion_type, "(rMAE one-step)"),
  paste0(dispersion_type, "(rRMSE two-step)"),
  paste0(dispersion_type, "(rMAE two-step)")
)

colnames(summary_disp) <- method_labels[colnames(summary_disp)]

print(round(summary_mean, 4))
print(round(summary_disp, 4))

# ---- Display table ----
make_display_row <- function(mean_vec, disp_vec, digits = 4) {
  paste0(
    sprintf(paste0("%.", digits, "f"), mean_vec),
    "\n(",
    sprintf(paste0("%.", digits, "f"), disp_vec),
    ")"
  )
}

display_table <- rbind(
  "rRMSE (one-step)" = make_display_row(rmse_one$mean, rmse_one$disp),
  "rMAE (one-step)"  = make_display_row(rmae_one$mean, rmae_one$disp),
  "rRMSE (two-step)" = make_display_row(rmse_two$mean, rmse_two$disp),
  "rMAE (two-step)"  = make_display_row(rmae_two$mean, rmae_two$disp)
)

colnames(display_table) <- method_labels 

print(display_table, quote = FALSE)

# ---- Final selected table ----

tab <- display_table[, unname(method_labels), drop = FALSE]

extract_mean_sd <- function(x) {
  x <- trimws(as.character(x))
  sp <- strsplit(x, "\n")[[1]]
  mean_val <- as.numeric(trimws(sp[1]))
  sd_val <- as.numeric(gsub("[()]", "", trimws(sp[2])))
  c(mean = mean_val, sd = sd_val)
}

mean_mat <- matrix(
  NA_real_,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = list(rownames(tab), colnames(tab))
)

sd_mat <- matrix(
  NA_real_,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = list(rownames(tab), colnames(tab))
)

for (i in 1:nrow(tab)) {
  for (j in 1:ncol(tab)) {
    tmp <- extract_mean_sd(tab[i, j])
    mean_mat[i, j] <- tmp["mean"]
    sd_mat[i, j] <- tmp["sd"]
  }
}

fmt_cell <- function(mean_val, sd_val, digits = 4) {
  paste0(
    sprintf(paste0("%.", digits, "f"), mean_val),
    "\n(",
    sprintf(paste0("%.", digits, "f"), sd_val),
    ")"
  )
}

tab_final <- matrix(
  NA_character_,
  nrow = nrow(mean_mat),
  ncol = ncol(mean_mat),
  dimnames = list(rownames(mean_mat), colnames(mean_mat))
)

for (i in 1:nrow(mean_mat)) {
  for (j in 1:ncol(mean_mat)) {
    tab_final[i, j] <- fmt_cell(mean_mat[i, j], sd_mat[i, j], digits = 4)
  }
}

tab_final_df <- data.frame(
  Metric = rownames(tab_final),
  tab_final,
  check.names = FALSE,
  row.names = NULL
)

mean_df <- data.frame(
  Metric = rownames(mean_mat),
  mean_mat,
  check.names = FALSE,
  row.names = NULL
)

sd_df <- data.frame(
  Metric = rownames(sd_mat),
  sd_mat,
  check.names = FALSE,
  row.names = NULL
)

print(tab_final_df, row.names = FALSE)

# ---- Save results ----
write.csv(tab_final_df, file = local_file("rolling_final_display_table.csv"), row.names = FALSE)
