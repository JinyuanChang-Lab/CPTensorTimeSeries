library(HDTSA)
library(tensorTS)
library(OLCPM)
library(ggplot2)
library(ggrepel)
library(grid)
library(tibble)
library(zoo)
library(xts)
library(scales)
library(dplyr)
library(sf)
library(tidyr)

data("BeijingAir")

Y1 <- BeijingAir

source("cp_cciso.R")


####################helper functions#################

winsorize_vec <- function(x, prob = 0.01) {
  stopifnot(is.numeric(x))
  stopifnot(length(prob) == 1, prob >= 0, prob < 0.5)
  qs <- quantile(x, probs = c(prob, 1 - prob), na.rm = TRUE, type = 7)
  x[x < qs[1]] <- qs[1]
  x[x > qs[2]] <- qs[2]
  x
}

winsorize_tensor_by_series <- function(X, prob = 0.01) {
  stopifnot(length(dim(X)) == 4)
  dims <- dim(X)
  X_win <- X
  for (i in seq_len(dims[2])) {
    for (j in seq_len(dims[3])) {
      for (k in seq_len(dims[4])) {
        X_win[, i, j, k] <- winsorize_vec(X[, i, j, k], prob = prob)
      }
    }
  }
  X_win
}

cp_residuals <- function(Y, f, A_list) {
  if (!is.array(Y) || length(dim(Y)) != 4) {
    stop("Y must be a 4-way array with dimensions n x d1 x d2 x d3.")
  }
  if (!is.matrix(f)) {
    stop("f must be an n x r matrix.")
  }
  if (!is.list(A_list) || length(A_list) != 3) {
    stop("A_list must be a list of length 3.")
  }
  
  dims <- dim(Y)
  n <- dims[1]
  d1 <- dims[2]
  d2 <- dims[3]
  d3 <- dims[4]
  r <- ncol(f)
  
  A1 <- A_list[[1]]
  A2 <- A_list[[2]]
  A3 <- A_list[[3]]
  
  if (nrow(f) != n) {
    stop("nrow(f) must equal dim(Y)[1].")
  }
  if (nrow(A1) != d1 || nrow(A2) != d2 || nrow(A3) != d3) {
    stop("Dimensions of A_list do not match Y.")
  }
  if (ncol(A1) != r || ncol(A2) != r || ncol(A3) != r) {
    stop("All loading matrices must have ncol = ncol(f).")
  }
  
  Y_hat <- array(0, dim = dims)
  for (t in seq_len(n)) {
    for (i in seq_len(r)) {
      comp12 <- outer(A1[, i], A2[, i])
      comp123 <- outer(as.vector(comp12), A3[, i])
      comp123 <- array(comp123, dim = c(d1, d2, d3))
      Y_hat[t, , , ] <- Y_hat[t, , , ] + f[t, i] * comp123
    }
  }
  
  list(residual = Y - Y_hat, fitted = Y_hat)
}

run_empirical_moment_test <- function(Y_array, alpha = 0.05, R = 1000) {
  dims <- dim(Y_array)
  out <- vector()
  for (i in seq_len(dims[2])) {
    for (j in seq_len(dims[3])) {
      for (k in seq_len(dims[4])) {
        out <- c(out, OLCPM::moment.determine(Y_array[, i, j, k], alpha = alpha, R = R))
      }
    }
  }
  out
}

moment_table_to_df <- function(moment_vec) {
  tb <- round(100 * table(moment_vec) / length(moment_vec), 2)
  data.frame(
    highest_supported_moment = names(tb),
    proportion_percent = as.numeric(tb),
    row.names = NULL
  )
}

####################alignment helper functions

get_row_index <- function(A_mode2, target_name = NULL, default_index = NULL) {
  rn <- rownames(A_mode2)
  if (!is.null(target_name) && !is.null(rn)) {
    hit <- which(rn == target_name)
    if (length(hit) == 1) {
      return(hit)
    }
  }
  return(default_index)
}

sign_to_target <- function(x, target = c("positive", "negative")) {
  target <- match.arg(target)
  if (is.na(x) || x == 0) {
    return(1)
  }
  if (target == "positive") {
    return(ifelse(x > 0, 1, -1))
  } else {
    return(ifelse(x < 0, 1, -1))
  }
}

detect_factor_roles <- function(A_list,
                                pm25_index = 1,
                                o3_index = 6) {
  # Use the pollution-variable mode (j = 2) to determine
  # which factor is ozone-related and which is general pollution.
  A2 <- A_list[[2]]
  r_hat <- ncol(A2)
  
  if (r_hat < 2) {
    stop("At least two factors are required for role detection.")
  }
  
  first5_idx <- setdiff(seq_len(nrow(A2)), o3_index)
  
  # Ozone-related factor:
  # large absolute loading on O3 relative to the first five pollutants
  ozone_score <- abs(A2[o3_index, ]) / (colMeans(abs(A2[first5_idx, , drop = FALSE])) + 1e-8)
  ozone_idx <- which.max(ozone_score)
  
  # General pollution factor:
  # relatively balanced loadings on the first five pollutants and small O3 loading
  remain_idx <- setdiff(seq_len(r_hat), ozone_idx)
  
  if (length(remain_idx) == 1) {
    general_idx <- remain_idx
  } else {
    general_score <- colMeans(abs(A2[first5_idx, remain_idx, drop = FALSE])) /
      (abs(A2[o3_index, remain_idx]) + 1e-8)
    general_idx <- remain_idx[which.max(general_score)]
  }
  
  # Put ozone-related factor first and general pollution factor second
  perm <- c(ozone_idx, general_idx, setdiff(seq_len(r_hat), c(ozone_idx, general_idx)))
  
  list(
    ozone_idx = ozone_idx,
    general_idx = general_idx,
    perm = perm
  )
}
apply_alignment_transform <- function(A_list,
                                      f_mat = NULL,
                                      Sigma_list = NULL,
                                      perm,
                                      sign_mode1,
                                      sign_mode2,
                                      sign_mode3) {
  r_hat <- ncol(A_list[[1]])
  
  # Reorder columns first
  A_list <- lapply(A_list, function(Aj) Aj[, perm, drop = FALSE])
  
  if (!is.null(Sigma_list)) {
    Sigma_list <- lapply(Sigma_list, function(Sj) Sj[, perm, drop = FALSE])
  }
  
  if (!is.null(f_mat)) {
    f_mat <- as.matrix(f_mat)[, perm, drop = FALSE]
  }
  
  # Apply sign conventions mode by mode
  S1 <- diag(sign_mode1, r_hat, r_hat)
  S2 <- diag(sign_mode2, r_hat, r_hat)
  S3 <- diag(sign_mode3, r_hat, r_hat)
  
  A_list[[1]] <- A_list[[1]] %*% S1
  A_list[[2]] <- A_list[[2]] %*% S2
  A_list[[3]] <- A_list[[3]] %*% S3
  
  if (!is.null(Sigma_list)) {
    Sigma_list[[1]] <- Sigma_list[[1]] %*% S1
    Sigma_list[[2]] <- Sigma_list[[2]] %*% S2
    Sigma_list[[3]] <- Sigma_list[[3]] %*% S3
  }
  
  # Adjust factors so that the CP reconstruction remains unchanged
  if (!is.null(f_mat)) {
    factor_sign <- sign_mode1 * sign_mode2 * sign_mode3
    f_mat <- f_mat %*% diag(factor_sign, r_hat, r_hat)
  }
  
  list(
    A_list = A_list,
    f_mat = f_mat,
    Sigma_list = Sigma_list
  )
}

align_one_estimate_main <- function(A_list,
                                    f_mat = NULL,
                                    Sigma_list = NULL,
                                    pm25_index = 1,
                                    o3_index = 6,
                                    station_index = 1,
                                    diurnal_index = 1) {
  role_info <- detect_factor_roles(
    A_list = A_list,
    pm25_index = pm25_index,
    o3_index = o3_index
  )
  
  perm <- role_info$perm
  A_tmp <- lapply(A_list, function(Aj) Aj[, perm, drop = FALSE])
  
  r_hat <- ncol(A_tmp[[1]])
  
  # Sign conventions from the manuscript:
  # (j = 2) ozone-related factor positive on O3,
  #         general pollution factor positive on PM2.5
  sign_mode2 <- rep(1, r_hat)
  sign_mode2[1] <- sign_to_target(A_tmp[[2]][o3_index, 1], "positive")
  sign_mode2[2] <- sign_to_target(A_tmp[[2]][pm25_index, 2], "positive")
  
  # (j = 1) both factors positive at the first station
  sign_mode1 <- rep(1, r_hat)
  sign_mode1[1] <- sign_to_target(A_tmp[[1]][station_index, 1], "positive")
  sign_mode1[2] <- sign_to_target(A_tmp[[1]][station_index, 2], "positive")
  
  # (j = 3) first elements of both factors negative
  sign_mode3 <- rep(1, r_hat)
  sign_mode3[1] <- sign_to_target(A_tmp[[3]][diurnal_index, 1], "negative")
  sign_mode3[2] <- sign_to_target(A_tmp[[3]][diurnal_index, 2], "negative")
  
  apply_alignment_transform(
    A_list = A_list,
    f_mat = f_mat,
    Sigma_list = Sigma_list,
    perm = perm,
    sign_mode1 = sign_mode1,
    sign_mode2 = sign_mode2,
    sign_mode3 = sign_mode3
  )
}


align_application_estimates_main <- function(A_iter, A.init, A_han, A_chen,
                                             f_iter, f_inl, f_han, f_chen,
                                             Sigma_list,
                                             pm25_index = 1,
                                             o3_index = 6,
                                             station_index = 1,
                                             diurnal_index = 1) {
  
  res_iter <- align_one_estimate_main(
    A_list = A_iter,
    f_mat = f_iter,
    Sigma_list = Sigma_list,
    pm25_index = pm25_index,
    o3_index = o3_index,
    station_index = station_index,
    diurnal_index = diurnal_index
  )
  
  res_inl <- align_one_estimate_main(
    A_list = A.init,
    f_mat = f_inl,
    Sigma_list = NULL,
    pm25_index = pm25_index,
    o3_index = o3_index,
    station_index = station_index,
    diurnal_index = diurnal_index
  )
  
  res_han <- align_one_estimate_main(
    A_list = A_han,
    f_mat = f_han,
    Sigma_list = NULL,
    pm25_index = pm25_index,
    o3_index = o3_index,
    station_index = station_index,
    diurnal_index = diurnal_index
  )
  
  res_chen <- align_one_estimate_main(
    A_list = A_chen,
    f_mat = f_chen,
    Sigma_list = NULL,
    pm25_index = pm25_index,
    o3_index = o3_index,
    station_index = station_index,
    diurnal_index = diurnal_index
  )
  
  list(
    A_iter = res_iter$A_list,
    A.init = res_inl$A_list,
    A_han = res_han$A_list,
    A_chen = res_chen$A_list,
    f_iter = res_iter$f_mat,
    f_inl = res_inl$f_mat,
    f_han = res_han$f_mat,
    f_chen = res_chen$f_mat,
    Sigma_list = res_iter$Sigma_list
  )
}

####################winsorized alignment
align_application_estimates_winsorized <- function(A_iter, f_iter, Sigma_list,
                                                   pm25_index = 1,
                                                   o3_index = 6,
                                                   station_index = 1,
                                                   diurnal_index = 1) {
  res_iter <- align_one_estimate_main(
    A_list = A_iter,
    f_mat = f_iter,
    Sigma_list = Sigma_list,
    pm25_index = pm25_index,
    o3_index = o3_index,
    station_index = station_index,
    diurnal_index = diurnal_index
  )
  
  list(
    A_iter = res_iter$A_list,
    f_iter = res_iter$f_mat,
    Sigma_list = res_iter$Sigma_list
  )
}

build_iter_inference_list <- function(A_iter, Sigma_list, f_iter, Y_array, D_vec, m_dim, r_hat, n_obs) {
  out_mode <- vector("list", m_dim)
  for (j in seq_len(m_dim)) {
    out_factor <- vector("list", r_hat)
    for (i in seq_len(r_hat)) {
      con <- vector()
      for (vv in seq_len(D_vec[j])) {
        aij_iter <- A_iter[[j]][, i]
        debias_res <- aij.debias.iter(
          A = A_iter,
          i = i,
          j = j,
          Sigma.yij.xii.1 = Sigma_list
        )
        aij_debias <- c(debias_res$aij.de)
        h <- rep(0, D_vec[j])
        h[vv] <- 1
        var_est <- cov.aij.debias.iter.est(h = h, A = A_iter, j = j, i = i, f = f_iter, Y = Y_array)
        std <- sqrt(var_est / n_obs)
        tstat <- as.numeric(t(h) %*% aij_debias / std)
        con_vv <- c(as.numeric(t(h) %*% aij_debias), std, tstat)
        con <- rbind(con, con_vv)
      }
      out_factor[[i]] <- con
    }
    out_mode[[j]] <- out_factor
  }
  out_mode
}

build_point_estimate_list <- function(A_list, r_hat) {
  m_dim <- length(A_list)
  out_mode <- vector("list", m_dim)
  for (j in seq_len(m_dim)) {
    out_factor <- vector("list", r_hat)
    for (i in seq_len(r_hat)) {
      out_factor[[i]] <- cbind(A_list[[j]][, i], NA_real_, NA_real_)
    }
    out_mode[[j]] <- out_factor
  }
  out_mode
}

build_mode_table <- function(con_list_iter, con_list_init, A_han, A_chen, mode_index, row_labels, method_han = "HOPE", method_chen = "CC-ISO") {
  n_factor <- length(con_list_iter[[mode_index]])
  out <- list()
  for (i in seq_len(n_factor)) {
    out[[i]] <- data.frame(
      factor = paste0("i=", i),
      item = row_labels,
      Pro.iter = con_list_iter[[mode_index]][[i]][, 1],
      Pro.iter.se = con_list_iter[[mode_index]][[i]][, 2],
      Pro.iter.t = con_list_iter[[mode_index]][[i]][, 3],
      Pro.init = con_list_init[[mode_index]][[i]][, 1],
      HOPE = A_han[[mode_index]][, i],
      CC.ISO = A_chen[[mode_index]][, i],
      row.names = NULL
    )
  }
  do.call(rbind, out)
}

build_mode_table_single_method <- function(con_list_iter, mode_index, row_labels) {
  n_factor <- length(con_list_iter[[mode_index]])
  out <- list()
  for (i in seq_len(n_factor)) {
    out[[i]] <- data.frame(
      factor = paste0("i=", i),
      item = row_labels,
      Pro.iter = con_list_iter[[mode_index]][[i]][, 1],
      Pro.iter.se = con_list_iter[[mode_index]][[i]][, 2],
      Pro.iter.t = con_list_iter[[mode_index]][[i]][, 3],
      row.names = NULL
    )
  }
  do.call(rbind, out)
}

save_mode3_plot <- function(i_factor,
                            con_list_iter,
                            con_list_init = NULL,
                            A_han = NULL,
                            A_chen = NULL,
                            file_stub,
                            only_proiter = FALSE) {
  main_col <- "#0072B2"
    zero_col <- "#7A7A7A"
      init_col <- "#D55E00"
        hope_col <- "#009E73"
          chen_col <- "#CC79A7"
            
          vals_iter <- con_list_iter[[3]][[i_factor]][, 1]
          ses_iter  <- con_list_iter[[3]][[i_factor]][, 2]
          
          dat_iter <- data.frame(
            hour = seq_along(vals_iter),
            value = vals_iter,
            se = ses_iter
          )
          dat_iter$lower <- dat_iter$value - qnorm(0.975) * dat_iter$se
          dat_iter$upper <- dat_iter$value + qnorm(0.975) * dat_iter$se
          
          if (only_proiter) {
            p <- ggplot(dat_iter, aes(x = hour, y = value)) +
              geom_ribbon(aes(ymin = lower, ymax = upper),
                          fill = "grey70", alpha = 0.5) +
              geom_hline(yintercept = 0, linetype = "dashed",
                         linewidth = 0.6, color = zero_col) +
              geom_line(linewidth = 0.8, color = main_col) +
              geom_point(size = 1.8, color = main_col) +
              scale_x_continuous(breaks = seq_along(vals_iter)) +
              labs(title = NULL, x = NULL, y = NULL) +
              theme_minimal(base_size = 12) +
              theme(
                axis.line.x = element_line(color = "#333333", linewidth = 0.6),
                axis.line.y = element_line(color = "#333333", linewidth = 0.6),
                axis.ticks = element_line(color = "#333333", linewidth = 0.5),
                axis.ticks.length = unit(3, "pt"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_line(color = "#ECECEC"),
                panel.grid.major.y = element_line(color = "#ECECEC"),
                axis.title = element_text(color = "#333333"),
                axis.text = element_text(color = "#333333")
              )
          } else {
            if (is.null(con_list_init) || is.null(A_han) || is.null(A_chen)) {
              stop("con_list_init, A_han, and A_chen must be provided when only_proiter = FALSE.")
            }
            
            vals_init <- con_list_init[[3]][[i_factor]][, 1]
            vals_hope <- A_han[[3]][, i_factor]
            vals_chen <- A_chen[[3]][, i_factor]
            
            dat_all <- rbind(
              data.frame(hour = seq_along(vals_iter), value = vals_iter, method = "Pro.iter"),
              data.frame(hour = seq_along(vals_init), value = vals_init, method = "Pro.init"),
              data.frame(hour = seq_along(vals_hope), value = vals_hope, method = "HOPE"),
              data.frame(hour = seq_along(vals_chen), value = vals_chen, method = "CC-ISO")
            )
            dat_all$method <- factor(dat_all$method, levels = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO"))
            
            p <- ggplot() +
              geom_ribbon(data = dat_iter, aes(x = hour, ymin = lower, ymax = upper),
                          fill = "grey70", alpha = 0.5) +
              geom_hline(yintercept = 0, linetype = "dashed",
                         linewidth = 0.6, color = zero_col) +
              geom_line(data = subset(dat_all, method != "Pro.iter"),
                        aes(x = hour, y = value, color = method, linetype = method),
                        linewidth = 0.8, alpha = 0.9) +
              geom_line(data = subset(dat_all, method == "Pro.iter"),
                        aes(x = hour, y = value, color = method, linetype = method),
                        linewidth = 1.2) +
              geom_point(data = subset(dat_all, method == "Pro.iter"),
                         aes(x = hour, y = value, color = method), size = 2.0) +
              scale_x_continuous(breaks = seq_along(vals_iter)) +
              scale_color_manual(
                values = c("Pro.iter" = main_col,
                           "Pro.init" = init_col,
                           "HOPE" = hope_col,
                           "CC-ISO" = chen_col),
                breaks = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO")
              ) +
              scale_linetype_manual(
                values = c("Pro.iter" = "solid",
                           "Pro.init" = "longdash",
                           "HOPE" = "dotdash",
                           "CC-ISO" = "twodash"),
                breaks = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO")
              ) +
              labs(title = NULL, x = NULL, y = NULL, color = NULL, linetype = NULL) +
              theme_minimal(base_size = 12) +
              theme(
                axis.line.x = element_line(color = "#333333", linewidth = 0.6),
                axis.line.y = element_line(color = "#333333", linewidth = 0.6),
                axis.ticks = element_line(color = "#333333", linewidth = 0.5),
                axis.ticks.length = unit(3, "pt"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_line(color = "#ECECEC"),
                panel.grid.major.y = element_line(color = "#ECECEC"),
                legend.position = "bottom",
                legend.background = element_rect(fill = "white", colour = "black"),
                legend.box.background = element_rect(fill = "white", colour = "black"),
                legend.key = element_rect(fill = "white", colour = NA)
              ) +
              guides(
                color = guide_legend(order = 1, nrow = 1, byrow = TRUE),
                linetype = guide_legend(order = 1, nrow = 1, byrow = TRUE)
              )
          }
          
          ggsave(paste0(file_stub, ".pdf"), p, width = 6.2, height = 4.8)
          
          invisible(p)
}

save_factor_ts_plot <- function(i_factor, f_iter, f_inl, f_han, f_chen, idx, file_stub) {
  main_col <- "#0072B2"
    init_col <- alpha("#D55E00", 0.80)
    hope_col <- "#009E73"
      chen_col <- "#CC79A7"
        col_grid <- "#ECECEC"
          col_axis <- "#333333"
            
          df <- tibble(
            date = idx,
            Pro.iter = as.numeric(f_iter[, i_factor]),
            Pro.init = as.numeric(f_inl[, i_factor]),
            HOPE = as.numeric(f_han[, i_factor]),
            `CC-ISO` = as.numeric(f_chen[, i_factor])
          ) |>
            pivot_longer(cols = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO"), names_to = "method", values_to = "value") |>
            group_by(method) |>
            arrange(date, .by_group = TRUE) |>
            mutate(ma30 = zoo::rollmean(value, k = 30, align = "right", fill = NA)) |>
            ungroup()
          
          df$method <- factor(df$method, levels = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO"))
          
          p <- ggplot() +
            geom_line(data = subset(df, method != "Pro.iter"), aes(x = date, y = ma30, color = method, linetype = method), linewidth = 0.95, na.rm = TRUE) +
            geom_line(data = subset(df, method == "Pro.iter"), aes(x = date, y = ma30, color = method, linetype = method), linewidth = 1.25, na.rm = TRUE) +
            geom_vline(xintercept = as.Date(c("2014-01-01", "2015-01-01", "2016-01-01", "2017-01-01")), linetype = "dotted", color = "#BEBEBE", linewidth = 0.5) +
            scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m", expand = expansion(mult = c(0.01, 0.02))) +
            scale_color_manual(values = c("Pro.iter" = main_col, "Pro.init" = init_col, "HOPE" = hope_col, "CC-ISO" = chen_col), breaks = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO")) +
            scale_linetype_manual(values = c("Pro.iter" = "solid", "Pro.init" = "longdash", "HOPE" = "dotdash", "CC-ISO" = "twodash"), breaks = c("Pro.iter", "Pro.init", "HOPE", "CC-ISO")) +
            labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = NULL) +
            theme_classic(base_size = 12) +
            theme(
              axis.text = element_text(color = col_axis),
              axis.title = element_text(color = col_axis),
              panel.grid.major.y = element_line(color = col_grid),
              panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.background = element_rect(fill = "white", colour = "black"),
              legend.box.background = element_rect(fill = "white", colour = "black"),
              legend.key = element_rect(fill = "white", colour = NA)
            ) +
            guides(color = guide_legend(order = 1, nrow = 1, byrow = TRUE), linetype = guide_legend(order = 1, nrow = 1, byrow = TRUE))
          
          ggsave(paste0(file_stub, ".pdf"), p, width = 6.8, height = 4.6)
          
          invisible(p)
}

save_beijing_map <- function(weights, file_stub) {
  adm2 <- st_read("gadm41_CHN_3.json", quiet = TRUE)
  bj_districts <- adm2[which(adm2$NAME_1 == "Beijing"), ]
  bj_outline <- bj_districts %>% summarise(geometry = st_union(geometry))
  bbox <- st_bbox(bj_districts)
  
  stations_df <- data.frame(
    station_en = c("Aotizhongxin", "Changping", "Dingling", "Dongsi", "Guanyuan",
                   "Gucheng", "Huairou", "Nongzhanguan", "Shunyi", "Tiantan", "Wanliu", "Wanshouxigong"),
    lon = c(116.407, 116.230, 116.170, 116.434, 116.361, 116.223, 116.643, 116.473, 116.720, 116.434, 116.315, 116.366),
    lat = c(40.003, 40.195, 40.287, 39.952, 39.942, 39.928, 40.394, 39.971, 40.144, 39.874, 39.994, 39.867)
  )
  
  weights_df <- data.frame(station_en = stations_df$station_en, weight = round(weights, 2))
  stations_map <- stations_df %>%
    left_join(weights_df, by = "station_en") %>%
    mutate(label_txt = sprintf("%.2f", weight), lab_col = "black")
  
  st_pts_w <- st_as_sf(stations_map, coords = c("lon", "lat"), crs = 4326)
  
  p <- ggplot() +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_sf(data = bj_districts, aes(geometry = geometry), fill = NA, color = "#7a92b2", linewidth = 0.5, linetype = "22") +
    geom_sf(data = bj_outline, fill = NA, color = "#1b3a57", linewidth = 1.1) +
    geom_sf(data = st_pts_w, aes(fill = weight), shape = 21, size = 5, stroke = 0) +
    geom_text_repel(data = stations_map, aes(x = lon, y = lat, label = station_en), size = 3.5, box.padding = 0.25, point.padding = 1, min.segment.length = 0.3, segment.color = "black", segment.size = 0.5, seed = 12345) +
    scale_fill_gradientn(colours = c("#2ecc71", "#f39c12", "#e74c3c"), limits = c(0, 0.5), oob = scales::squish, name = "Weight") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = 0) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "right") +
    labs(title = NULL, subtitle = NULL, x = NULL, y = NULL)
  
  ggsave(paste0(file_stub, ".pdf"), p, width = 6.0, height = 6.0)
  
  invisible(p)
}

save_factor_ts_plot_proiter <- function(i_factor,
                                        f_iter,
                                        file_stub,
                                        start_date = as.Date("2013-03-01"),
                                        end_date = as.Date("2017-02-28"),
                                        ma_window = 30) {
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(scales)
  
  x <- as.numeric(f_iter[, i_factor])
  idx <- seq.Date(start_date, end_date, by = "day")
  
  if (length(x) != length(idx)) {
    stop(sprintf("Length of factor series = %d does not match number of dates = %d.",
                 length(x), length(idx)))
  }
  
  df <- data.frame(
    date = idx,
    value = x
  )
  df$ma30 <- zoo::rollmean(df$value, k = ma_window, align = "right", fill = NA)
  
  col_raw  <- "#9A9A9A"
    col_ma   <- "#0072B2"
      col_grid <- "#ECECEC"
        col_axis <- "#333333"
          
        p <- ggplot(df, aes(date, value)) +
          geom_line(color = col_raw, linewidth = 0.45, alpha = 0.9, na.rm = TRUE) +
          geom_line(aes(y = ma30), color = col_ma, linewidth = 1.1, na.rm = TRUE) +
          geom_vline(
            xintercept = as.Date(c("2014-01-01", "2015-01-01", "2016-01-01", "2017-01-01")),
            linetype = "dotted",
            color = "#BEBEBE",
            linewidth = 0.5
          ) +
          scale_x_date(
            date_breaks = "6 months",
            date_labels = "%Y-%m",
            expand = expansion(mult = c(0.01, 0.02))
          ) +
          labs(title = NULL, x = NULL, y = NULL) +
          theme_classic(base_size = 12) +
          theme(
            axis.text = element_text(color = col_axis),
            axis.title = element_text(color = col_axis),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "#555555"),
            panel.grid.major.y = element_line(color = col_grid),
            panel.grid.minor = element_blank()
          )
        
        ggsave(paste0(file_stub, ".pdf"), p, width = 6.8, height = 4.6)
        
        invisible(p)
}


####################helper functions for figure data#################

build_iter_diurnal_df <- function(con_list_iter, factor_index) {
  data.frame(
    hour = seq_len(nrow(con_list_iter[[3]][[factor_index]])),
    value = con_list_iter[[3]][[factor_index]][, 1],
    se = con_list_iter[[3]][[factor_index]][, 2],
    tstat = con_list_iter[[3]][[factor_index]][, 3]
  )
}

build_iter_factor_ts_df <- function(f_iter, factor_index,
                                    start_date = as.Date("2013-03-01"),
                                    end_date = as.Date("2017-02-28")) {
  idx <- seq.Date(start_date, end_date, by = "day")
  if (length(idx) != nrow(f_iter)) {
    stop("The date index length does not match nrow(f_iter).")
  }
  data.frame(
    date = idx,
    value = as.numeric(f_iter[, factor_index])
  )
}


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

####################main estimation#################

set.seed(123456)

m <- length(dim(Y1)) - 1
n <- dim(Y1)[1]
D <- dim(Y1)[-1]
K <- 10

xi_chang <- tensor.est.xi(Y1)

res_rank0 <- HDTTS.CP.est(
  Y = Y1,
  xi = xi_chang,
  K = K,
  Ratio.type = "log"
)

r_hat_first <- res_rank0$r.hat
r_breve <- 2 * r_hat_first

xi_res <- RP.xi.sel(
  Y = Y1,
  r_breve = r_breve,
  eps = 0.1,
  lag.k = K,
  Randomized.time = 50,
  A = NULL
)

xi_sel <- xi_res$xi.sel

res_inl <- HDTTS.CP.est(
  Y = Y1,
  xi = xi_sel,
  Rank = NULL,
  K = K,
  Threshold = TRUE,
  delta = NULL,
  Ratio.type = "log",
  grid_delta1 = 50,
  delta_max = 0.1
)

r_hat <- res_inl$r.hat

con_BMP <- CP_TTS(
  Y = Y1,
  xi = xi_sel,
  r = r_hat,
  A.init = NULL,
  control.DPI = list(
    lag.k.dpi = 10,
    threshold = TRUE,
    delta = NULL,
    delta2 = rep(1, m),
    ratio.type = "log",
    random.projection = FALSE,
    iter.max = 100,
    eps = 1e-05,
    grid.num = 50,
    delta.max = 0.1,
    print.eps = TRUE,
    iter.lag = 1,
    all.put = TRUE,
    A = NULL,
    component = NULL
  )
)

YY1 <- aperm(Y1, c(2:(m + 1), 1))

res_han <- cp.iso.han(
  x = YY1,
  r = r_hat,
  A.real = NULL,
  niter = 100
)

res_chen <- cp.iso(
  x = YY1,
  r = r_hat,
  A.real = NULL,
  niter = 100,
  detect_close_eigenvals = TRUE
)

A_iter <- con_BMP$A.hat
A.init <- con_BMP$A.init
A_han <- res_han$Q
A_chen <- res_chen$Q

f_iter <- as.matrix(con_BMP$f.hat)
f_inl <- as.matrix(con_BMP$f.hat.inl)
f_han <- t(res_han$ft)
f_chen <- t(res_chen$ft)

Sigma_yij_xii_1 <- con_BMP$Sigma.yij.xii.1

aligned_main <- align_application_estimates_main(
  A_iter = A_iter,
  A.init = A.init,
  A_han = A_han,
  A_chen = A_chen,
  f_iter = f_iter,
  f_inl = f_inl,
  f_han = f_han,
  f_chen = f_chen,
  Sigma_list = Sigma_yij_xii_1
)

A_iter <- aligned_main$A_iter
A.init <- aligned_main$A.init
A_han <- aligned_main$A_han
A_chen <- aligned_main$A_chen

f_iter <- aligned_main$f_iter
f_inl <- aligned_main$f_inl
f_han <- aligned_main$f_han
f_chen <- aligned_main$f_chen

Sigma_yij_xii_1 <- aligned_main$Sigma_list

####################iterative inference and Pro.init point estimates#################

con_list_j_iter <- build_iter_inference_list(
  A_iter = A_iter,
  Sigma_list = Sigma_yij_xii_1,
  f_iter = f_iter,
  Y_array = Y1,
  D_vec = D,
  m_dim = m,
  r_hat = r_hat,
  n_obs = n
)

con_list_j_inl <- build_point_estimate_list(
  A.init,
  r_hat = r_hat
)



####################acf data#################

residual_obj <- cp_residuals(
  Y = Y1,
  f = f_iter,
  A_list = A_iter
)

E11 <- residual_obj$residual[, 1, , 1]
colnames(E11) <- c("PM2.5", "PM10", "SO2", "NO2", "CO", "O3")
colnames(f_iter) <- c("Factor 1", "Factor 2")

pdf("acf_factor.pdf", width = 10, height = 6)
acf(f_iter, lag.max = 200)
dev.off()

pdf("acf_residual.pdf", width = 10, height = 6)
acf(E11, lag.max = 200)
dev.off()

####################output tables and figures for the main analysis#################

pollutant_labels <- c("PM2.5", "PM10", "SO2", "NO2", "CO", "O3")

station_labels <- c(
  "Aotizhongxin", "Changping", "Dingling", "Dongsi", "Guanyuan",
  "Gucheng", "Huairou", "Nongzhanguan", "Shunyi", "Tiantan",
  "Wanliu", "Wanshouxigong"
)

diurnal_labels <- paste0("Hour_", 1:24)

table_a2_all <- build_mode_table(
  con_list_iter = con_list_j_iter,
  con_list_init = con_list_j_inl,
  A_han = A_han,
  A_chen = A_chen,
  mode_index = 2,
  row_labels = pollutant_labels
)

table_a1_all <- build_mode_table(
  con_list_iter = con_list_j_iter,
  con_list_init = con_list_j_inl,
  A_han = A_han,
  A_chen = A_chen,
  mode_index = 1,
  row_labels = station_labels
)

## Main Table: Pro.iter loading estimates for pollution-variable mode
table_a2_main <- table_a2_all[
  ,
  c("factor", "item", "Pro.iter", "Pro.iter.se", "Pro.iter.t"),
  drop = FALSE
]

write.csv(table_a2_main, "table_loading_a2_main.csv", row.names = FALSE)

## Supplement Tables T1 and T2
write.csv(table_a2_all, "table_loading_a2_all.csv", row.names = FALSE)
write.csv(table_a1_all, "table_loading_a1_all.csv", row.names = FALSE)

## Main Figure: station-mode maps
save_beijing_map(
  weights = A_iter[[1]][, 1],
  file_stub = "beijing_map_a1"
)

save_beijing_map(
  weights = A_iter[[1]][, 2],
  file_stub = "beijing_map_a2"
)

## Main Figure: diurnal-mode loading curves
save_mode3_plot(
  i_factor = 1,
  con_list_iter = con_list_j_iter,
  con_list_init = con_list_j_inl,
  A_han = A_han,
  A_chen = A_chen,
  file_stub = "mode3_a1_Pro.iter",
  only_proiter = TRUE
)

save_mode3_plot(
  i_factor = 2,
  con_list_iter = con_list_j_iter,
  con_list_init = con_list_j_inl,
  A_han = A_han,
  A_chen = A_chen,
  file_stub = "mode3_a2_Pro.iter",
  only_proiter = TRUE
)

## Main Figure: latent factor time series
save_factor_ts_plot_proiter(
  i_factor = 1,
  f_iter = f_iter,
  file_stub = "timeseries_factor1"
)

save_factor_ts_plot_proiter(
  i_factor = 2,
  f_iter = f_iter,
  file_stub = "timeseries_factor2"
)

## Supplement Figure F5: comparison of diurnal-mode loadings
save_mode3_plot(
  i_factor = 1,
  con_list_iter = con_list_j_iter,
  con_list_init = con_list_j_inl,
  A_han = A_han,
  A_chen = A_chen,
  file_stub = "mode3_a1_compare",
  only_proiter = FALSE
)

save_mode3_plot(
  i_factor = 2,
  con_list_iter = con_list_j_iter,
  con_list_init = con_list_j_inl,
  A_han = A_han,
  A_chen = A_chen,
  file_stub = "mode3_a2_compare",
  only_proiter = FALSE
)

## Supplement Figure F6: comparison of estimated latent factors
idx <- seq.Date(
  as.Date("2013-03-01"),
  as.Date("2017-02-28"),
  by = "day"
)

save_factor_ts_plot(
  i_factor = 1,
  f_iter = f_iter,
  f_inl = f_inl,
  f_han = f_han,
  f_chen = f_chen,
  idx = idx,
  file_stub = "factor1_timeseries_compare"
)

save_factor_ts_plot(
  i_factor = 2,
  f_iter = f_iter,
  f_inl = f_inl,
  f_han = f_han,
  f_chen = f_chen,
  idx = idx,
  file_stub = "factor2_timeseries_compare"
)




####################moment test: Supplement Table T3#################

set.seed(123456)

Y_win <- winsorize_tensor_by_series(
  X = Y1,
  prob = 0.005
)

moment_raw <- run_empirical_moment_test(
  Y_array = Y1,
  alpha = 0.05,
  R = 1000
)

moment_win <- run_empirical_moment_test(
  Y_array = Y_win,
  alpha = 0.05,
  R = 1000
)

moment_raw_chr <- as.character(moment_raw)
moment_win_chr <- as.character(moment_win)

suppressWarnings({
  moment_raw_num <- as.numeric(moment_raw_chr)
  moment_win_num <- as.numeric(moment_win_chr)
})

moment_raw_chr[!is.na(moment_raw_num) & moment_raw_num < 4] <- "below 4th moment"
moment_win_chr[!is.na(moment_win_num) & moment_win_num < 4] <- "below 4th moment"

moment_raw_df <- moment_table_to_df(moment_raw_chr)
moment_raw_df$data_type <- "raw"

moment_win_df <- moment_table_to_df(moment_win_chr)
moment_win_df$data_type <- "winsorized"

table_moment_test <- rbind(
  moment_raw_df,
  moment_win_df
)

table_moment_test <- table_moment_test[
  ,
  c("data_type", "highest_supported_moment", "proportion_percent"),
  drop = FALSE
]

write.csv(table_moment_test, "table_moment_test.csv", row.names = FALSE)

####################robust analysis with winsorized data####################

Y_robust <- winsorize_tensor_by_series(
  X = Y1,
  prob = 0.005
)

K_robust <- 10

set.seed(123456)

xi_chang_robust <- tensor.est.xi(Y_robust)

res_rank0_robust <- HDTTS.CP.est(
  Y = Y_robust,
  xi = xi_chang_robust,
  K = K_robust,
  Ratio.type = "log"
)

r_hat_first_robust <- res_rank0_robust$r.hat

xi_res_robust <- RP.xi.sel(
  Y = Y_robust,
  r_breve = 2 * r_hat_first_robust,
  eps = 0.1,
  lag.k = K_robust,
  Randomized.time = 50,
  A = NULL
)

res_inl_robust <- HDTTS.CP.est(
  Y = Y_robust,
  xi = xi_res_robust$xi.sel,
  Rank = NULL,
  K = K_robust,
  Threshold = TRUE,
  delta = NULL,
  Ratio.type = "log",
  grid_delta1 = 50,
  delta_max = 0.1
)

r_hat_robust <- res_inl_robust$r.hat

con_BMP_robust <- CP_TTS(
  Y = Y_robust,
  xi = xi_res_robust$xi.sel,
  r = r_hat_robust,
  A.init = NULL,
  control.DPI = list(
    lag.k.dpi = 10,
    threshold = TRUE,
    delta = NULL,
    delta2 = rep(1, m),
    ratio.type = "log",
    random.projection = FALSE,
    iter.max = 100,
    eps = 1e-05,
    grid.num = 50,
    delta.max = 0.1,
    print.eps = TRUE,
    iter.lag = 1,
    all.put = TRUE,
    A = NULL,
    component = NULL
  )
)

A_iter_robust <- con_BMP_robust$A.hat
f_iter_robust <- as.matrix(con_BMP_robust$f.hat)
Sigma_yij_xii_1_robust <- con_BMP_robust$Sigma.yij.xii.1

aligned_robust <- align_application_estimates_winsorized(
  A_iter = A_iter_robust,
  f_iter = f_iter_robust,
  Sigma_list = Sigma_yij_xii_1_robust
)

A_iter_robust <- aligned_robust$A_iter
f_iter_robust <- aligned_robust$f_iter
Sigma_yij_xii_1_robust <- aligned_robust$Sigma_list

con_list_j_iter_robust <- build_iter_inference_list(
  A_iter = A_iter_robust,
  Sigma_list = Sigma_yij_xii_1_robust,
  f_iter = f_iter_robust,
  Y_array = Y_robust,
  D_vec = D,
  m_dim = m,
  r_hat = r_hat_robust,
  n_obs = n
)

####################robust output: Supplement Table T4 and Figures F7--F9####################

table_a2_robust <- build_mode_table_single_method(
  con_list_iter = con_list_j_iter_robust,
  mode_index = 2,
  row_labels = pollutant_labels
)

write.csv(table_a2_robust, "table_loading_a2_robust.csv", row.names = FALSE)

## Supplement Figure F7: robust station-mode maps
save_beijing_map(
  weights = A_iter_robust[[1]][, 1],
  file_stub = "beijing_map_a1_robust"
)

save_beijing_map(
  weights = A_iter_robust[[1]][, 2],
  file_stub = "beijing_map_a2_robust"
)

## Supplement Figure F8: robust diurnal-mode loading curves
save_mode3_plot(
  i_factor = 1,
  con_list_iter = con_list_j_iter_robust,
  file_stub = "mode3_a1_robust",
  only_proiter = TRUE
)

save_mode3_plot(
  i_factor = 2,
  con_list_iter = con_list_j_iter_robust,
  file_stub = "mode3_a2_robust",
  only_proiter = TRUE
)

## Supplement Figure F9: robust latent factor time series
save_factor_ts_plot_proiter(
  i_factor = 1,
  f_iter = f_iter_robust,
  file_stub = "timeseries_factor1_robust"
)

save_factor_ts_plot_proiter(
  i_factor = 2,
  f_iter = f_iter_robust,
  file_stub = "timeseries_factor2_robust"
)

############time series plot##############


data_fit <- Y1[, , , 1]

n <- dim(data_fit)[1]
p <- dim(data_fit)[2]
q <- dim(data_fit)[3]

vec_data <- as.vector(data_fit)

date_tol <- seq.Date(
  from = as.Date("1964-01-01"),
  to = as.Date("2021-12-31"),
  by = "day"
)[1:n]

# Build a long-format data frame for all p x q marginal time series
data <- data.frame(
  Time = rep(date_tol, p * q),
  Value = vec_data,
  S = factor(
    rep(paste0("S", 1:q), each = n * p),
    levels = paste0("S", 1:q)
  ),
  BM = factor(
    rep(rep(paste0("BM", 1:p), each = n), q),
    levels = paste0("BM", 1:p)
  )
)

# Combine the two mode labels into one group variable
data_long <- data %>%
  unite("Group", S, BM, sep = "-", remove = FALSE)

# Draw the faceted time-series plot
pp <- ggplot(data_long, aes(x = Time, y = Value)) +
  geom_line(linewidth = 0.25) +
  facet_wrap(~ BM + S, scales = "free_y", ncol = q) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    # Remove facet labels and their background
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    
    # Remove axis tick labels
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    
    # Remove axis titles
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Adjust spacing between panels
    panel.spacing = unit(0.5, "lines"),
    
    # Add a border to each small panel
    panel.border = element_rect(color = "black", fill = NA),
    
    # Adjust plot margins
    plot.margin = margin(t = 20, r = 10, b = 0, l = 20)
  )

pdf("Y_time_series_airpollution.pdf", width = 8, height = 6)

print(pp)

dev.off()
