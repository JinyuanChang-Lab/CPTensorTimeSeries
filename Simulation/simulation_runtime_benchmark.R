library(pracma)
library(Rcpp)
library(RcppEigen)
library(foreach)
library(doParallel)
library(jointDiag)
library(HDTSA) 


detectCores() 
cl <- makeCluster(100)
registerDoParallel(cl)


n_tol = c(400)
delta.fl_tol = c(0.25,0.75)
factor.loading_tol = c(1) #1:random sparse; 2:structure sparse
factor.corr_tol = c(0,0.75)
ALPHA_tol = c(0,0.3,0.6)

D_tol = list(c(20,20),c(40,40),c(60,60),c(80,80))
D = c(20,20)
m = length(D)
r = 3
beta = c(0.8,0.75,0.7)
ar.coef = as.list(beta)
w = rep(15,r) 

K = 10

delta.fl       = 0.25
factor.loading = "sparse-random-corr1"
factor.corr    = 0.75
ALPHA          = 0.6
n              = 400


 
source("cp_cciso.R")


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

#################### helper: build simulation data ####################

simulate_tensor_data <- function(n, D, ALPHA, delta.fl, factor.loading, factor.corr, seed_num, x) {
  m <- length(D)
  r <- 3
  beta <- c(0.8, 0.75, 0.7)
  ar.coef <- as.list(beta)
  w <- rep(15, r)
  
  set.seed(seed_num + x)
  
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
  
  list(
    Y = data$Y,
    A = data$A,
    Y1 = aperm(data$Y, c(2:(m + 1), 1))
  )
}

#################### helper: control list for the proposed method ####################

make_bmp_control_runtime <- function(K, m, A = NULL, random.projection = TRUE) {
  list(
    lag.k.dpi = K,
    threshold = TRUE,
    delta = NULL,
    delta2 = rep(1, m),
    ratio.type = "log",
    random.projection = random.projection,
    iter.max = 20,
    eps = 1e-4,
    grid.num = 50,
    delta.max = 0.1,
    print.eps = FALSE,
    iter.lag = 1,
    all.put = TRUE,
    A = A,
    component = NULL
  )
}
#################### methods to be benchmarked ####################

run_pro_method <- function(Y, A, K = 10) {
  m <- length(dim(Y)) - 1
  
  CP_TTS(
    Y = Y,
    xi = NULL,
    r  = NULL,
    A.init = NULL,
    control.DPI = make_bmp_control_runtime(
      K = K,
      m = m,
      A = A,
      random.projection = TRUE
    )
  )
}

run_pro_method_nosel <- function(Y, A, K = 10) {
  m <- length(dim(Y)) - 1
  
  CP_TTS(
    Y = Y,
    xi = NULL,
    r = NULL,
    A.init = NULL,
    control.DPI = make_bmp_control_runtime(
      K = K,
      m = m,
      A = A,
      random.projection = FALSE
    )
  )
}

run_han_method <- function(Y1) {
  rank_est <- eigratio_test(x = Y1, kmax = 5)
  cp.iso.han(x = Y1, r = rank_est, A.real = NULL, niter = 100)
}

run_chen_method <- function(Y1) {
  rank_est <- eigratio_test(x = Y1, kmax = 5)
  cp.iso(x = Y1, r = rank_est, A.real = NULL, niter = 100, detect_close_eigenvals = TRUE)
}

#################### one replication ####################

benchmark_one_rep <- function(x, n, D, ALPHA, delta.fl, factor.loading,
                              factor.corr, seed_num, K = 10) {
  
  dat <- simulate_tensor_data(
    n = n,
    D = D,
    ALPHA = ALPHA,
    delta.fl = delta.fl,
    factor.loading = factor.loading,
    factor.corr = factor.corr,
    seed_num = seed_num,
    x = x
  )
  
  out <- peakRAM::peakRAM(
    run_pro_method(dat$Y, dat$A, K = K),
    run_pro_method_nosel(dat$Y, dat$A, K = K),
    run_han_method(dat$Y1),
    run_chen_method(dat$Y1)
  )
  
  if (nrow(out) == 4) {
    out[, 1] <- c("Pro.sel", "Pro.nosel", "HOPE", "CC-ISO")
  }
  
  out
}

#################### main loop ####################

rep <- 100
seed_num <- 123456
n <- 400

con_list <- list()
hgl <- 1

for (factor.corr in factor.corr_tol) {
  for (factor.loading.num in factor.loading_tol) {
    factor.loading <- paste0("sparse-random-corr", factor.loading.num)
    
    for (delta.fl in delta.fl_tol) {
      for (ALPHA in ALPHA_tol) {
        for (D in D_tol) {
          
          con <- foreach(
            x = 1:rep,
            .packages = c(
              "HDTSA", "R.utils", "tensor", "rTensor", "MASS", "pracma",
              "Rcpp", "RcppArmadillo", "RcppEigen", "vars", "base",
              "peakRAM", "utils"
            ),
            .export = c(
              "CP_TTS",
              "simulate_tensor_data",
              "make_bmp_control_runtime",
              "run_pro_method",
              "run_pro_method_nosel",
              "run_han_method",
              "run_chen_method",
              "benchmark_one_rep"
            ),
            .errorhandling = "pass",
            .verbose = FALSE
          ) %dopar% benchmark_one_rep(
            x = x,
            n = n,
            D = D,
            ALPHA = ALPHA,
            delta.fl = delta.fl,
            factor.loading = factor.loading,
            factor.corr = factor.corr,
            seed_num = seed_num,
            K = 10
          )
          
          con_list[[hgl]] <- con
          
          save.image("temp_runtime_updated.RData")
          
          valid_con <- Filter(function(z) is.data.frame(z), con)
          
          name_i <- paste0(
            "(D)", D[1],
            "(factor.corr)", factor.corr,
            "(random-sparse-corr)", factor.loading.num,
            "(delta.fl)", delta.fl,
            "(alpha)", ALPHA,
            "(n)", n,
            "(count)", length(valid_con)
          )
          
          print(name_i)
          
          if (length(valid_con) > 0) {
            res_tol <- Reduce(`+`, lapply(valid_con, function(z) z[, -1, drop = FALSE]))
            row.names(res_tol) <- valid_con[[1]][, 1]
            print(round(res_tol / length(valid_con), 4))
          }
          
          cat("\n")
          
          seed_num <- seed_num + rep
          hgl <- hgl + 1
        }
      }
    }
  }
}

direction <- paste0(
  "simulation_runtime_",
  format(Sys.time(), "%Y%m%d-%H-%M-%S"),
  ".RData"
)

save.image(direction)


########################################## analysis ############################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

############################## helper functions ##############################

normalize_peakram_one <- function(obj) {
  df <- as.data.frame(obj, stringsAsFactors = FALSE)
  
  if (ncol(df) < 4) {
    stop("Each peakRAM result must have at least 4 columns.")
  }
  
  df <- df[, 1:4, drop = FALSE]
  colnames(df) <- c("Method", "Runtime", "RAM_used", "RAM_peak")
  
  df$Method <- as.character(df$Method)
  df$Runtime <- as.numeric(as.character(df$Runtime))
  df$RAM_used <- as.numeric(as.character(df$RAM_used))
  df$RAM_peak <- as.numeric(as.character(df$RAM_peak))
  
  bad_method <- is.na(df$Method) | df$Method == "" | grepl("^[0-9.]+$", df$Method)
  if (all(bad_method) && !is.null(rownames(df))) {
    rn <- rownames(df)
    if (length(rn) == nrow(df)) {
      df$Method <- rn
    }
  }
  
  df$Method <- trimws(df$Method)
  df$Method[df$Method %in% c("Rolling_pro", "Pro.sel", "Pro.iter1")] <- "Pro.sel"
  df$Method[df$Method %in% c("Rolling_pro_nosel", "Pro.nosel", "Pro.iter")] <- "Pro.nosel"
  df$Method[df$Method %in% c("Rolling_han", "HOPE")] <- "HOPE"
  df$Method[df$Method %in% c("Rolling_chen", "CC-ISO")] <- "CC-ISO"
  
  if (nrow(df) == 4 && !all(df$Method %in% c("Pro.sel", "Pro.nosel", "HOPE", "CC-ISO"))) {
    df$Method <- c("Pro.sel", "Pro.nosel", "HOPE", "CC-ISO")
  }
  
  df
}

build_runtime_long_df <- function(con_list,
                                  factor.corr_tol,
                                  factor.loading_tol,
                                  delta.fl_tol,
                                  ALPHA_tol,
                                  D_tol,
                                  n_selected = 400) {
  out_list <- list()
  hgl <- 1
  kk <- 1
  
  for (factor.corr in factor.corr_tol) {
    for (factor.loading.num in factor.loading_tol) {
      for (delta.fl in delta.fl_tol) {
        for (ALPHA in ALPHA_tol) {
          for (D in D_tol) {
            
            con_i <- con_list[[hgl]]
            
            if (is.null(con_i) || length(con_i) == 0) {
              hgl <- hgl + 1
              next
            }
            
            for (ss in seq_along(con_i)) {
              one_rep <- con_i[[ss]]
              if (is.null(one_rep)) next
              
              one_rep <- tryCatch(
                normalize_peakram_one(one_rep),
                error = function(e) NULL
              )
              if (is.null(one_rep)) next
              
              one_rep$D <- D[1]
              one_rep$n <- n_selected
              one_rep$factor.corr <- factor.corr
              one_rep$factor.loading <- factor.loading.num
              one_rep$delta.fl <- delta.fl
              one_rep$ALPHA <- ALPHA
              one_rep$rep_id <- ss
              
              out_list[[kk]] <- one_rep
              kk <- kk + 1
            }
            
            hgl <- hgl + 1
          }
        }
      }
    }
  }
  
  out_df <- do.call(rbind, out_list)
  rownames(out_df) <- NULL
  out_df
}

make_runtime_panel <- function(df_summary, factor.corr, delta.fl, ALPHA, ylab_text) {
  color_type <- c(
    "Pro.iter" = "black",
    "HOPE" = "red",
    "CC-ISO" = "blue"
  )
  
  line_type <- c(
    "Pro.iter" = "solid",
    "HOPE" = "solid",
    "CC-ISO" = "solid"
  )
  
  shape_type <- c(
    "Pro.iter" = 15,
    "HOPE" = 17,
    "CC-ISO" = 16
  )
  
  n1 <- substitute(
    expression(rho == value1 * "," ~ phi == value2 * "," ~ s == value3),
    list(value1 = factor.corr, value2 = delta.fl, value3 = ALPHA)
  )
  
  ggplot(df_summary,
         aes(x = D, y = mean_value, color = Method, fill = Method, group = Method, shape = Method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = mean_value - se, ymax = mean_value + se), alpha = 0.2, linewidth = 0) +
    scale_shape_manual(values = shape_type) +
    scale_fill_manual(values = color_type) +
    scale_color_manual(values = color_type) +
    ggtitle(eval(n1)) +
    xlab("dimension") +
    ylab(ylab_text) +
    theme(legend.position = "none", legend.direction = "horizontal") +
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
}

############################## build long-format data ##############################

n_selected <- 400
factor.loading_selected <- 1

runtime_df <- build_runtime_long_df(
  con_list = con_list,
  factor.corr_tol = factor.corr_tol,
  factor.loading_tol = factor.loading_tol,
  delta.fl_tol = delta.fl_tol,
  ALPHA_tol = ALPHA_tol,
  D_tol = D_tol,
  n_selected = n_selected
)

############################## keep exactly 3 methods as in the old figure ##############################

runtime_df <- runtime_df %>%
  filter(Method %in% c("Pro.sel", "HOPE", "CC-ISO"))

runtime_df$Method <- dplyr::recode(
  runtime_df$Method,
  "Pro.sel" = "Pro.iter"
)

runtime_df$Method <- factor(runtime_df$Method, levels = c("Pro.iter", "HOPE", "CC-ISO"))

############################## summary ##############################

runtime_summary <- runtime_df %>%
  group_by(Method, D, n, factor.corr, factor.loading, delta.fl, ALPHA) %>%
  summarize(
    mean_value = mean(Runtime, na.rm = TRUE),
    se = sd(Runtime, na.rm = TRUE),
    .groups = "drop"
  )

ram_summary <- runtime_df %>%
  group_by(Method, D, n, factor.corr, factor.loading, delta.fl, ALPHA) %>%
  summarize(
    mean_value = mean(RAM_peak, na.rm = TRUE),
    se = sd(RAM_peak, na.rm = TRUE),
    .groups = "drop"
  )

############################## make the two figures ##############################

plot_list_runtime <- list()
plot_list_ram <- list()
kk <- 1

for (factor.corr in factor.corr_tol) {
  for (delta.fl in delta.fl_tol) {
    for (ALPHA in ALPHA_tol) {
      
      df_runtime_i <- runtime_summary %>%
        filter(
          n == n_selected,
          ALPHA == !!ALPHA,
          delta.fl == !!delta.fl,
          factor.loading == factor.loading_selected,
          factor.corr == !!factor.corr
        )
      
      df_ram_i <- ram_summary %>%
        filter(
          n == n_selected,
          ALPHA == !!ALPHA,
          delta.fl == !!delta.fl,
          factor.loading == factor.loading_selected,
          factor.corr == !!factor.corr
        )
      
      plot_list_runtime[[kk]] <- make_runtime_panel(
        df_summary = df_runtime_i,
        factor.corr = factor.corr,
        delta.fl = delta.fl,
        ALPHA = ALPHA,
        ylab_text = "runtime(seconds)"
      )
      
      plot_list_ram[[kk]] <- make_runtime_panel(
        df_summary = df_ram_i,
        factor.corr = factor.corr,
        delta.fl = delta.fl,
        ALPHA = ALPHA,
        ylab_text = "peak RAM(MiB)"
      )
      
      kk <- kk + 1
    }
  }
}
 

############################## save ##############################

graphics.off()
while (!is.null(dev.list())) dev.off()

pdf("runtime_benchmark.pdf", width = 14, height = 10)
grid.draw(arrangeGrob(grobs = plot_list_runtime, ncol = 3))
dev.off()

pdf("ram_benchmark.pdf", width = 14, height = 10)
grid.draw(arrangeGrob(grobs = plot_list_ram, ncol = 3))
dev.off()
 
 