library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus); 
library(surveillance);library(tidyr); library(lubridate); library(magrittr); library(patchwork)
rm=ls()
days <- 0:150  
prob_density <- dnorm(days, mean = 75, sd = 25)  
case_counts <- prob_density / max(prob_density) * 50 
case_counts_scaled <- round(case_counts)

df_samples <- data.frame(
  time_onset = days,
  case_count = case_counts_scaled)

last_time_onset <- max(df_samples$time_onset)

temp_pre_days <- data.frame(
  time_onset = -400:-1,
  case_count = 0
)

# 새로운 10일치 데이터 생성
temp_lastdays <- data.frame(
  time_onset = (last_time_onset + 1):(last_time_onset + 10),
  case_count = 0  # 신규 날짜의 case_count는 0으로 설정
)

df_samples_full <- bind_rows(temp_pre_days, df_samples, temp_lastdays)

df_samples <- df_samples_full %>%
  mutate(time_onset = 0:(n() - 1))
#--------------------------------------------------------------------------------------------
sample_ip_normal <- function(n) {rnorm(n, mean = 5, sd = 1)}

# Short-term (감마 분포)
sample_ip_gamma <- function(n) rgamma(n, shape = 1.2, scale = 22.2)

# Bimodal (감마 + 정규 혼합 분포) 
sample_ip_bimodal <- function(n) {
  p <- 83 / (142 + 83)  # normal 쪽 비율
  n_normal <- rbinom(1, n, prob = p)
  n_gamma <- n - n_normal
  
  c(rgamma(n_gamma, shape = 1.2, scale = 22.2),
    rnorm(n_normal, mean = 337.4, sd = 38.5))
}
#--------------------------------------------------------------------------------------------
simulate_infection <- function(df_onset, ip_sampler, n = 100000, start_date = "2025-04-01") {
  infection_records <- vector("list", n)
  
  for (sim in 1:n) {
    infection_list <- list()
    
    for (i in 1:nrow(df_onset)) {
      n_cases <- df_onset$case_count[i]
      onset_day <- df_onset$time_onset[i]
      ip_samples <- ip_sampler(n_cases)
      infect_days <- floor(onset_day - ip_samples)
      infection_list[[i]] <- infect_days
    }
    
    all_infect_days <- unlist(infection_list)
    sim_counts <- table(factor(all_infect_days, levels = 0:max(df_onset$time_onset)))
    infection_records[[sim]] <- as.numeric(sim_counts)
  }
  
  infections_mat <- do.call(rbind, infection_records)
  mean_infections <- colMeans(infections_mat)
  
  # 랜덤 trajectory 선택
  set.seed(123)
  random_index <- sample(1:n, 1)
  selected_trajectory <- infection_records[[random_index]]
  
  # 날짜는 trajectory 기준으로 맞추기
  infection_days <- 0:(length(selected_trajectory) - 1)
  first_infection_day <- min(infection_days[selected_trajectory > 0])
  start_date <- as.Date(start_date)
  date_seq <- start_date + (infection_days - first_infection_day)
  
  summary_df <- data.frame(
    date = date_seq,
    mean = mean_infections,
    trajectory = selected_trajectory
  )
  
  return(summary_df)
}

#--------------------------------------------------------------------------------------------
result_normal <- simulate_infection(df_samples, sample_ip_normal)
sum(result_normal$mean)

# result_gamma <- simulate_infection(df_samples, sample_ip_gamma)
# sum(result_gamma$mean)
# 
# result_bimodal <- simulate_infection(df_samples, sample_ip_bimodal)
# sum(result_bimodal$mean)
# head(result_normal)
#--------------------------------------------------------------------------------------------
plot_backproj <- function(backproj_result, df_samples, title) {
  
  y_max <- max(max(df_samples$case_count), max(backproj_result$trajectory, na.rm = TRUE))
  
  p1 <- ggplot(df_samples, aes(x = time_onset, y = case_count)) +
    geom_col(fill = "gray50") +
    coord_cartesian(ylim = c(0, y_max)) +  
    labs(title = "Observed onset cases",
         x = "Symptom onset day", y = "Case count") +
    theme_minimal()
  
  p2 <- ggplot(backproj_result, aes(x = date)) +
    #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightcoral", alpha = 0.4) +
    geom_line(aes(y = trajectory), color = "firebrick", size = 1) +
    labs(title = title,
         x = "Estimated infection day", y = "Infection count (simulated)") +
    theme_minimal()
  
  return(p1 / p2 + plot_layout(heights = c(1, 1.2)))
}
#--------------------------------------------------------------------------------------------
options(repr.plot.width = 12, repr.plot.height = 10)
plot_backproj(result_normal, df_samples, "Back-projection Result(normal)")
# plot_backproj(result_gamma, df_samples, "Back-projection Result(gamma)")
# plot_backproj(result_bimodal, df_samples, "Back-projection Result(bimodal)")
#--------------------------------------------------------------------------------------------
result_normal<-result_normal%>%filter(date>=350)
result_normal$date <- 0:(nrow(result_normal) - 1)
# 
# result_gamma<-result_gamma%>%filter(time_infection>=350)
# result_gamma$time_infection <- 0:(nrow(result_gamma) - 1)
#--------------------------------------------------------------------------------------------
uni<-function(delay) {
  ifelse(delay >= 4 & delay <= 7, 1 / (7 - 4 + 1), 0)
}

# generation time distribution
gi_fit=list(shape=2.305, scale=5.452)
generation<-function(t){pweibull(t, shape=gi_fit$shape, scale=gi_fit$scale)-
    pweibull(t-1, shape=gi_fit$shape, scale=gi_fit$scale)}

# type 1) normal
incubation_normal<-function(n){pnorm(n, mean=5, sd=1)-pnorm(n-1, mean=5, sd=1)}
# type 2) short(gamma)
incubation_short<-function(r){pgamma(r, shape=1.2, scale=22.2)-pgamma(r-1, shape=1.2, scale=22.2)}

# type 3) Bimodal
inc_fit=list(p= 83/(142+83), shape=1.2, scale=22.2, mean=337.4, sd=38.5)
incubation_bimodal<-function(r){(1 - inc_fit$p) * pgamma(r, shape=inc_fit$shape, scale=inc_fit$scale) + 
    inc_fit$p * pnorm(r, mean=inc_fit$mean, sd=inc_fit$sd)-((1 - inc_fit$p) * pgamma(r-1, shape=inc_fit$shape, scale=inc_fit$scale) + 
                                                              inc_fit$p * pnorm(r-1, mean=inc_fit$mean, sd=inc_fit$sd))}

# 책 dist -> check용
# inc_fit=list(meanlog=1.519, sdlog=0.615)
# incubation<-function(t){plnorm(t, inc_fit$meanlog, inc_fit$sdlog)-
#     plnorm(t-1, inc_fit$meanlog, inc_fit$sdlog)}
#--------------------------------------------------------------------------------------------
prepare_matrix<- function(df_infection, incubation_func, Start.T){
  uni_vals <- uni(1:1000)
  incub_vals <- incubation_func(1:1000)
  
  conv <- convolve(uni_vals, rev(incub_vals), type = "open")
  
  today <- max(df_infection$date) - min(df_infection$date)
  
  max_time <- max(df_infection$date)
  precal <- matrix(0, nrow = max_time - Start.T + 3, ncol = max_time - 1)
  
  for (m in (Start.T - 2):max_time) {
    for (n in 1:(m - 1)) {
      idx <- today - m + n
      if (idx > 0 && idx <= length(conv)) {
        precal[m - Start.T + 3, n] <- sum(conv[1:idx])
      } else {
        precal[m - Start.T + 3, n] <- 0  # 인덱스 벗어나면 0
      }
    }
  }
  return(precal)
}
#--------------------------------------------------------------------------------------------
precal_normal<-prepare_matrix(result_normal, incubation_normal,10)
result_normal <- result_normal %>% dplyr::select(date, trajectory)

# precal_gamma<-prepare_matrix(result_gamma, incubation_short,10)
# result_gamma <- result_gamma %>% dplyr::select(time_infection, trajectory)
# 
# precal_bimodal <- prepare_matrix(result_bimodal, incubation_bimodal,1)
# result_bimodal <- result_bimodal %>% dplyr::select(time_infection, trajectory)
#--------------------------------------------------------------------------------------------
estimate_and_plot_Rt <- function(result_df, precal_mat, Start.T = 10, max_R = 20, by = 0.01, num) {
  max_tau <- max(result_df$date)
  gen_cache <- sapply(1:max_tau, generation)
  
  est.t <- list()
  est.CI <- list()
  
  for (TT in Start.T:max(result_df$date)) {
    dt.backproj.T <- result_df %>% filter(date <= TT)
    
    llk <- function(R) {
      t <- TT
      tau_seq <- 1:(t - 1)
      
      mean_part <- dt.backproj.T$trajectory[(t - 1):1]
      gen_part <- gen_cache[tau_seq]
      precal_part <- precal_mat[t - Start.T + 3, tau_seq]
      
      Css <- mean_part * gen_part / precal_part
      Cs <- sum(Css) * R
      Cs[Cs <= 0] <- 1e-5
      
      obs <- dt.backproj.T$trajectory[t + 1]
      return(-(-Cs + obs * log(Cs) - lgamma(obs + 1)))
    }
    
    opt_est <- optim(par = 0.7, fn = llk, method = "L-BFGS-B", lower = 0, control = list(maxit = 1000))
    est.t[[TT]] <- opt_est$par
    
    CI <- function(par_CI) {
      2 * (-llk(par_CI) + opt_est$value)
    }
    
    par_CI_seq <- seq(0, max_R, by = by)
    logLik_vals <- sapply(par_CI_seq, CI)
    data_CI <- data.frame(par_CI = par_CI_seq, logLik = logLik_vals)
    data_CI$logLik[data_CI$logLik < (max(data_CI$logLik, na.rm = TRUE) - 3.84)] <- NA
    data_CI <- na.omit(data_CI)
    
    ci_pro <- matrix(NA, ncol = 2, nrow = 1)
    if (nrow(data_CI) > 0) {
      ci_pro[1, 1] <- min(data_CI$par_CI)
      ci_pro[1, 2] <- max(data_CI$par_CI)
    }
    
    est.CI[[TT]] <- as.data.frame(ci_pro)
    colnames(est.CI[[TT]]) <- c("lower", "upper")
  }
  
  est_mat <- matrix(unlist(est.t), ncol = 1, byrow = TRUE)
  ci_mat <- matrix(unlist(est.CI), ncol = 2, byrow = TRUE)
  colnames(ci_mat) <- c("lower", "upper")
  est_df <- cbind(est_mat, ci_mat)
  colnames(est_df)[1] <- "Rt"
  est_df <- as.data.frame(est_df)
  est_df$date <- Start.T:max(result_df$date)
  
  result <- merge(result_df, est_df, by = "date", all.x = TRUE)
  result <- result %>% rename(case_delay = trajectory)
  result_withcut <- result %>% filter(date>=num &date <= max(date) - 10)
  
  write.csv(result_withcut,file = "estimated_Rt.csv",row.names = FALSE)
  
  # 시각화
  scaling_parameter <- max(result_withcut$case_delay) / max(result_withcut$upper, na.rm = TRUE)
  p <- result_withcut %>%
    ggplot() +
    geom_bar(aes(x = date, y = case_delay), stat = 'identity', fill = "gray35", width = 0.7) +
    geom_line(aes(x = date, y = Rt * scaling_parameter), color = "#1380A1", size = 1) +
    geom_ribbon(aes(x = date,
                    ymin = lower * scaling_parameter,
                    ymax = upper * scaling_parameter),
                fill = "#1380A1", alpha = 0.4) +
    labs(x = "\n 추정 감염일", y = "감염자 수\n") +
    theme_bw() +
    theme(text = element_text(size = 12, family = "sans", color = "black"),
          axis.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    scale_y_continuous(
      limits = c(0, 1.1 * max(result_withcut$case_delay)),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . / scaling_parameter,
                          breaks = c(0, 1, 2, 3),
                          name = "실질감염재생산수\n")
    ) +
    geom_hline(yintercept = 1 * scaling_parameter, linetype = "dashed", color = "#1380A1", size = 0.7)
  
  print(p)
  return(invisible(result_withcut))
}
#--------------------------------------------------------------------------------------------
estimate_and_plot_Rt(result_normal, precal_normal, Start.T = 10, num=400)
# estimate_and_plot_Rt(result_gamma, precal_gamma, Start.T = 10, num=300)
# estimate_and_plot_Rt(result_bimodal, precal_bimodal, Start.T = 1, num=380)
