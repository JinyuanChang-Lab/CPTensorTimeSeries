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


n_tol = c(400)           #sample size   
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
DGP.TCP <- HDTSA:::DGP.TCP
tensor.est.xi <- HDTSA:::tensor.est.xi
HDTTS.CP.est <- HDTSA:::HDTTS.CP.est
RP.xi.sel <- HDTSA:::RP.xi.sel

Rolling_fun2 <- function(x, n, ALPHA, delta.fl, factor.loading, factor.corr,
                         seed_num, Rank_set = NULL) {
  
  set.seed(seed_num + x)
  
  D <- c(20, 20)
  m <- length(D)
  r <- 3
  beta <- c(0.8, 0.75, 0.7)
  ar.coef <- as.list(beta)
  w <- rep(15, r)
  K <- 10
  
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
  
  make_control <- function(ratio.type = "log",
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
  
  xi.chang <- tensor.est.xi(Y)
  
  res.only.used.rank <- HDTTS.CP.est(
    Y = Y,
    xi = xi.chang,
    K = K,
    Ratio.type = "log"
  )
  
  r_hat_first <- min(res.only.used.rank$r.hat, min(D) / 4)
  
  xi.res <- RP.xi.sel(
    Y = Y,
    r_breve = 2 * r_hat_first,
    eps = 0.1,
    lag.k = K,
    Randomized.time = 50,
    A = A
  )
  
  fit <- CP_TTS(
    Y = Y,
    xi = xi.res$xi.sel,
    r = Rank_set,
    A.init = NULL,
    control.DPI = make_control(
      ratio.type = "log",
      threshold = TRUE,
      delta2 = rep(1, m),
      iter.max = 20,
      eps = 1e-5,
      print.eps = FALSE,
      iter.lag = 1,
      all.put = TRUE,
      A = A,
      component = data$C
    )
  )
  
  error_iter <- max(rho2.loss.list(fit$A.hat, A))
  error_init <- max(rho2.loss.list(fit$A.init, A))
  
  c(
    Pro.iter = error_iter,
    Pro.init = error_init
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
          for (Rank_set in 1:10) {
            
            con <- foreach(
              x = 1:rep,
              .combine = "rbind",
              .packages = c(
                "HDTSA", "R.utils", "tensor", "rTensor", "MASS",
                "pracma", "Rcpp", "RcppArmadillo", "RcppEigen",
                "vars", "base", "tensorTS", "utils"
              ),
              .export = c(
                "Rolling_fun2", "CP_TTS", "rho2.loss.list",
                "DGP.TCP", "tensor.est.xi", "HDTTS.CP.est", "RP.xi.sel"
              ),
              .errorhandling = "remove",
              .verbose = FALSE
            ) %dopar% Rolling_fun2(
              x = x,
              n = n,
              ALPHA = ALPHA,
              delta.fl = delta.fl,
              factor.loading = factor.loading,
              factor.corr = factor.corr,
              seed_num = seed_num,
              Rank_set = Rank_set
            )
            
            con_list[[hgl]] <- con
            
            meta_list[[hgl]] <- data.frame(
              factor.corr = factor.corr,
              factor.loading.num = factor.loading.num,
              delta.fl = delta.fl,
              ALPHA = ALPHA,
              n = n,
              R = Rank_set
            )
            
            save.image("temp_estimation_WrongR.RData")
            
            name_i <- paste0(
              "(factor.corr)", factor.corr,
              "(random-sparse-corr)", factor.loading.num,
              "(delta.fl)", delta.fl,
              "(alpha)", ALPHA,
              "(n)", n,
              "(rtilde)", Rank_set,
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
}

direction <- paste0(
  "simulation_wrongR_",
  format(Sys.time(), "%Y%m%d-%H-%M-%S"),
  ".RData"
)

save.image(direction)

build_con_df <- function(con_list, meta_list) {
  
  out_list <- vector("list", length(con_list))
  
  for (h in seq_along(con_list)) {
    
    con_i <- con_list[[h]]
    meta_i <- meta_list[[h]]
    
    if (is.null(con_i) || length(con_i) == 0) {
      next
    }
    
    con_i <- as.data.frame(con_i)
    
    if (ncol(con_i) < 2) {
      next
    }
    
    con_i <- con_i[, 1:2, drop = FALSE]
    colnames(con_i) <- c("Pro.iter", "Pro.init")
    
    con_i$Pro.iter <- as.numeric(as.character(con_i$Pro.iter))
    con_i$Pro.init <- as.numeric(as.character(con_i$Pro.init))
    
    con_i$factor.corr <- meta_i$factor.corr
    con_i$factor.loading.num <- meta_i$factor.loading.num
    con_i$delta.fl <- meta_i$delta.fl
    con_i$ALPHA <- meta_i$ALPHA
    con_i$n <- meta_i$n
    con_i$R <- meta_i$R
    
    out_list[[h]] <- con_i
  }
  
  out_df <- do.call(rbind, out_list)
  rownames(out_df) <- NULL
  out_df
}

con_df_all <- build_con_df(con_list, meta_list)

selected_methods <- c("Pro.iter", "Pro.init")

make_long_df <- function(df, selected_methods) {
  
  long_list <- lapply(selected_methods, function(mth) {
    data.frame(
      factor.corr = df$factor.corr,
      factor.loading.num = df$factor.loading.num,
      delta.fl = df$delta.fl,
      ALPHA = df$ALPHA,
      n = df$n,
      R = df$R,
      method = mth,
      mse = df[[mth]]
    )
  })
  
  long_df <- do.call(rbind, long_list)
  rownames(long_df) <- NULL
  long_df$method <- factor(long_df$method, levels = selected_methods)
  long_df
}

con_long <- make_long_df(con_df_all, selected_methods)

con_mean <- aggregate(
  mse ~ factor.corr + factor.loading.num + delta.fl + ALPHA + n + R + method,
  data = con_long,
  FUN = function(x) mean(x, na.rm = TRUE)
)

con_sd <- aggregate(
  mse ~ factor.corr + factor.loading.num + delta.fl + ALPHA + n + R + method,
  data = con_long,
  FUN = function(x) sd(x, na.rm = TRUE)
)

colnames(con_mean)[colnames(con_mean) == "mse"] <- "mse_mean"
colnames(con_sd)[colnames(con_sd) == "mse"] <- "mse_sd"

con_plot <- merge(
  con_mean,
  con_sd,
  by = c("factor.corr", "factor.loading.num", "delta.fl", "ALPHA", "n", "R", "method"),
  all = TRUE
)

con_plot$method <- factor(con_plot$method, levels = selected_methods)
con_plot$mse_display <- sprintf("%.2f (%.2f)", con_plot$mse_mean, con_plot$mse_sd)

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(scales)

color_type <- c("black", "red")
line_type <- c("solid", "solid")
shape_type <- c(15, 18)

plot_list <- list()
kk <- 1

delta.fl_show <- c(0.25, 0.75)
n_show <- 400
factor.loading.show <- 1

for (factor.corr in factor.corr_tol) {
  for (delta.fl in delta.fl_show) {
    for (ALPHA in ALPHA_tol) {
      
      place <- (con_plot$n == n_show) &
        (con_plot$ALPHA == ALPHA) &
        (con_plot$delta.fl == delta.fl) &
        (con_plot$factor.loading.num == factor.loading.show) &
        (con_plot$factor.corr == factor.corr)
      
      result_i <- con_plot[place, ]
      result_i$method <- factor(result_i$method, levels = selected_methods)
      
      title_i <- substitute(
        expression(rho == value1 * "," ~ phi == value2 * "," ~ s == value3),
        list(value1 = factor.corr, value2 = delta.fl, value3 = ALPHA)
      )
      
      p_i <- ggplot(
        result_i,
        aes(
          x = R,
          y = mse_mean,
          col = method,
          linetype = method,
          shape = method,
          group = method
        )
      ) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        scale_shape_manual(
          values = shape_type,
          breaks = selected_methods,
          limits = selected_methods
        ) +
        scale_color_manual(
          values = color_type,
          breaks = selected_methods,
          limits = selected_methods
        ) +
        scale_linetype_manual(
          values = line_type,
          breaks = selected_methods,
          limits = selected_methods
        ) +
        ggtitle(eval(title_i)) +
        xlab(expression(tilde(r))) +
        ylab("estimation error") +
        scale_x_continuous(breaks = c(1, 5, 10)) +
        theme(
          legend.position = "none",
          legend.direction = "horizontal",
          axis.line = element_line(size = 0.5, linetype = "solid"),
          legend.text = element_text(family = "serif"),
          panel.background = element_rect(fill = "gray100", colour = alpha("black", 0.5)),
          legend.key = element_rect(fill = alpha("white", 1)),
          panel.grid.major = element_line(size = 0.4, linetype = "dashed", colour = "gray92"),
          legend.background = element_rect(colour = "black"),
          legend.title = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 15, face = "bold")
        )
      
      plot_list[[kk]] <- p_i
      kk <- kk + 1
    }
  }
}

graphics.off()
while (!is.null(dev.list())) dev.off()

pdf("rtilde_selected_methods.pdf", width = 10 * (950 / 880), height = 10)
pig <- arrangeGrob(grobs = plot_list, ncol = 3)
grid.draw(pig)
dev.off()