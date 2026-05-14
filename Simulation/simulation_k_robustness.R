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

delta.fl       = 0.75
factor.loading = "sparse-random-corr1"
factor.corr    = 0.75
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

 
Rolling_fun2_Krobust <- function(x, n, ALPHA, delta.fl, factor.loading, factor.corr,
                                 seed_num, K) {
  
  set.seed(seed_num + x)
  
  D <- c(20, 20)
  m <- length(D)
  r <- 3
  beta <- c(0.8, 0.75, 0.7)
  ar.coef <- as.list(beta)
  w <- rep(15, r)
  
  out_na <- c(
    Pro.iter = NA_real_,
    Pro.init = NA_real_,
    rank.est = NA_real_
  )
  
  res <- tryCatch({
    
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
    
    make_bmp_control <- function() {
      list(
        lag.k.dpi = K,
        threshold = TRUE,
        delta = NULL,
        delta2 = rep(1, m),
        ratio.type = "log",
        random.projection = FALSE,
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
    
    xi.chang <- tensor.est.xi(Y)
    
    res.only.used.rank <- HDTTS.CP.est(
      Y = Y,
      xi = xi.chang,
      K = K,
      Ratio.type = "log"
    )
    
    r_hat_first <- res.only.used.rank$r.hat
    
    if (
      length(r_hat_first) != 1 ||
      is.na(r_hat_first) ||
      !is.finite(r_hat_first) ||
      r_hat_first < 1
    ) {
      return(out_na)
    }
    
    r_breve_use <- max(2L, 2L * as.integer(round(r_hat_first)))
    
    xi.res <- RP.xi.sel(
      Y = Y,
      r_breve = r_breve_use,
      eps = 0.1,
      lag.k = K,
      Randomized.time = 50,
      A = A
    )
    
    con.BMP <- CP_TTS(
      Y = Y,
      xi = xi.res$xi.sel,
      r = NULL,
      A.init = NULL,
      control.DPI = make_bmp_control()
    )
    
    if (
      is.null(con.BMP$A.hat) ||
      is.null(con.BMP$A.init) ||
      is.null(con.BMP$r)
    ) {
      return(out_na)
    }
    
    A.iter.pro <- con.BMP$A.hat
    A.init.pro <- con.BMP$A.init
    rank.est <- con.BMP$r
    
    if (
      length(rank.est) != 1 ||
      is.na(rank.est) ||
      !is.finite(rank.est) ||
      rank.est < 1
    ) {
      return(out_na)
    }
    
    c(
      Pro.iter = max(rho2.loss.list(A.iter.pro, A)),
      Pro.init = max(rho2.loss.list(A.init.pro, A)),
      rank.est = as.numeric(rank.est)
    )
    
  }, error = function(e) {
    out_na
  })
  
  return(res)
}
 
############################## run simulation ##############################

rep <- 2000
seed_num <- 123456

K_tol <- 2:15
con_list <- list()

hgl <- 1

for (factor.corr in factor.corr_tol) {
  
  for (factor.loading.num in factor.loading_tol) {
    
    factor.loading <- paste0("sparse-random-corr", factor.loading.num)
    
    for (delta.fl in delta.fl_tol) {
      
      for (ALPHA in ALPHA_tol) {
        
        for (n in n_tol) {
          
          for (K in K_tol) {
            
            con <- foreach(
              x = 1:rep,
              .combine = "rbind",
              .packages = c(
                "HDTSA", "R.utils", "tensor", "rTensor", "MASS",
                "pracma", "Rcpp", "RcppArmadillo", "RcppEigen",
                "vars", "base", "utils"
              ),
              .errorhandling = "stop",
              .verbose = FALSE
            ) %dopar% Rolling_fun2_Krobust(
              x = x,
              n = n,
              ALPHA = ALPHA,
              delta.fl = delta.fl,
              factor.loading = factor.loading,
              factor.corr = factor.corr,
              seed_num = seed_num,
              K = K
            )
            
             
            
            con_list[[hgl]] <- con
            
            names(con_list)[hgl] <- paste0(
              "rho_", factor.corr,
              "_loading_", factor.loading.num,
              "_delta_", delta.fl,
              "_alpha_", ALPHA,
              "_n_", n,
              "_K_", K
            )
            
            print(names(con_list)[hgl])
            cat("\n")
            
            seed_num <- seed_num + rep
            hgl <- hgl + 1
          }
        }
      }
    }
  }
}

############################## save ##############################

direction = paste0("simulation_Krobust",format(Sys.time(), "%Y%m%d-%H-%M-%S"),".RData")

save.image(direction)

stopImplicitCluster()
stopCluster(cl)                    


 
########################################## analysis ##########################################

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

############################## 1. rebuild meta information ##############################

# If names(con_list) were created in the simulation loop, we parse them.
# Otherwise, we rebuild the metadata using the same loop order.

if (!is.null(names(con_list)) && all(nzchar(names(con_list)))) {
  
  parse_one_name <- function(x) {
    parts <- strsplit(x, "_")[[1]]
    data.frame(
      factor.corr = as.numeric(parts[2]),
      factor.loading = as.numeric(parts[4]),
      delta.fl = as.numeric(parts[6]),
      ALPHA = as.numeric(parts[8]),
      n = as.numeric(parts[10]),
      K = as.numeric(parts[12])
    )
  }
  
  meta_df <- do.call(rbind, lapply(names(con_list), parse_one_name))
  
} else {
  
  meta_list <- list()
  hgl <- 1
  
  for (factor.corr in factor.corr_tol) {
    for (factor.loading.num in factor.loading_tol) {
      for (delta.fl in delta.fl_tol) {
        for (ALPHA in ALPHA_tol) {
          for (n in n_tol) {
            for (K in 2:15) {
              
              meta_list[[hgl]] <- data.frame(
                factor.corr = factor.corr,
                factor.loading = factor.loading.num,
                delta.fl = delta.fl,
                ALPHA = ALPHA,
                n = n,
                K = K
              )
              
              hgl <- hgl + 1
            }
          }
        }
      }
    }
  }
  
  meta_df <- do.call(rbind, meta_list)
}

rownames(meta_df) <- NULL

############################## 2. summarize con_list ##############################

summary_list <- vector("list", length(con_list))

for (ii in seq_along(con_list)) {
  
  con_i <- con_list[[ii]]
  con_i <- as.matrix(con_i)
  storage.mode(con_i) <- "numeric"
  
  if (ncol(con_i) < 3) {
    stop("Each element of con_list must have at least 3 columns: Pro.iter, Pro.init, rank.est.")
  }
  
  colnames(con_i) <- c("Pro.iter", "Pro.init", "rank.est")
  
  summary_list[[ii]] <- data.frame(
    meta_df[ii, ],
    Pro.iter.mean = mean(con_i[, "Pro.iter"], na.rm = TRUE),
    Pro.iter.sd   = sd(con_i[, "Pro.iter"], na.rm = TRUE),
    Pro.init.mean = mean(con_i[, "Pro.init"], na.rm = TRUE),
    Pro.init.sd   = sd(con_i[, "Pro.init"], na.rm = TRUE),
    rank.acc      = 100 * mean(con_i[, "rank.est"] == 3, na.rm = TRUE)
  )
}

result_krobust <- do.call(rbind, summary_list)
rownames(result_krobust) <- NULL

############################## 3. long-format data for plotting ##############################

## estimation error
con_mean_iter <- data.frame(
  result_krobust[, c("factor.corr", "factor.loading", "delta.fl", "ALPHA", "n", "K")],
  mse = result_krobust$Pro.iter.mean,
  method = "Pro.iter"
)

con_mean_inl <- data.frame(
  result_krobust[, c("factor.corr", "factor.loading", "delta.fl", "ALPHA", "n", "K")],
  mse = result_krobust$Pro.init.mean,
  method = "Pro.init"
)

con_error_plot <- rbind(con_mean_iter, con_mean_inl)
con_error_plot$method <- factor(con_error_plot$method, levels = c("Pro.iter", "Pro.init"))

## rank accuracy
con_acc_plot <- data.frame(
  result_krobust[, c("factor.corr", "factor.loading", "delta.fl", "ALPHA", "n", "K")],
  mse = result_krobust$rank.acc,
  method = "Log-ER"
)

############################## 4. plotting style ##############################

color_type_error <- c("Pro.iter" = "blue", "Pro.init" = "red")
line_type_error  <- c("Pro.iter" = "solid", "Pro.init" = "solid")
shape_type_error <- c("Pro.iter" = 16, "Pro.init" = 17)

color_type_acc <- c("Log-ER" = "blue")
line_type_acc  <- c("Log-ER" = "solid")
shape_type_acc <- c("Log-ER" = 16)

delta.fl_tol <- c(0.25, 0.75)
n_selected <- 400
factor.loading_selected <- 1

############################## 5. plot rank accuracy vs K ##############################

plot_list_acc <- list()
kk <- 1

for (factor.corr in factor.corr_tol) {
  for (delta.fl in delta.fl_tol) {
    for (ALPHA in ALPHA_tol) {
      
      place <- (con_acc_plot$n == n_selected) &
        (con_acc_plot$ALPHA == ALPHA) &
        (con_acc_plot$delta.fl == delta.fl) &
        (con_acc_plot$factor.loading == factor.loading_selected) &
        (con_acc_plot$factor.corr == factor.corr)
      
      result_i <- con_acc_plot[place, ]
      
      n1 <- substitute(
        expression(rho == value1 * "," ~ phi == value2 * "," ~ s == value3),
        list(value1 = factor.corr, value2 = delta.fl, value3 = ALPHA)
      )
      
      p_acc <- ggplot(result_i, aes(x = K, y = mse, col = method, linetype = method, shape = method)) +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_shape_manual(values = shape_type_acc) +
        ggtitle(eval(n1)) +
        xlab("K") +
        ylab("accuracy") +
        scale_x_continuous(breaks = c(2, 5, 10, 15)) +
        coord_cartesian(ylim = c(80, 100)) +
        theme(
          legend.position = "none",
          legend.direction = "horizontal",
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
        ) +
        scale_color_manual(values = color_type_acc) +
        scale_linetype_manual(values = line_type_acc)
      
      plot_list_acc[[kk]] <- p_acc
      kk <- kk + 1
    }
  }
}

############################## 6. plot estimation error vs K ##############################

plot_list_error <- list()
kk <- 1

for (factor.corr in factor.corr_tol) {
  for (delta.fl in delta.fl_tol) {
    for (ALPHA in ALPHA_tol) {
      
      place <- (con_error_plot$n == n_selected) &
        (con_error_plot$ALPHA == ALPHA) &
        (con_error_plot$delta.fl == delta.fl) &
        (con_error_plot$factor.loading == factor.loading_selected) &
        (con_error_plot$factor.corr == factor.corr)
      
      result_i <- con_error_plot[place, ]
      
      n1 <- substitute(
        expression(rho == value1 * "," ~ phi == value2 * "," ~ s == value3),
        list(value1 = factor.corr, value2 = delta.fl, value3 = ALPHA)
      )
      
      p_err <- ggplot(result_i, aes(x = K, y = mse, col = method, linetype = method, shape = method)) +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_shape_manual(values = shape_type_error) +
        ggtitle(eval(n1)) +
        xlab("K") +
        ylab("estimation error") +
        scale_x_continuous(breaks = c(2, 5, 10, 15)) +
        coord_cartesian(ylim = c(0, 0.25)) +
        theme(
          legend.position = "none",
          legend.direction = "horizontal",
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
        ) +
        scale_color_manual(values = color_type_error) +
        scale_linetype_manual(values = line_type_error)
      
      plot_list_error[[kk]] <- p_err
      kk <- kk + 1
    }
  }
}

############################## 7. save plots ##############################
 
graphics.off()
while (!is.null(dev.list())) dev.off()

pdf("Krobust_accuracy.pdf", width = 14, height = 10)
pig_acc <- arrangeGrob(grobs = plot_list_acc, ncol = 3)
grid.draw(pig_acc)
dev.off()

pdf("Krobust_error.pdf", width = 14, height = 10)
pig_err <- arrangeGrob(grobs = plot_list_error, ncol = 3)
grid.draw(pig_err)
dev.off()
 













