library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus); library(readxl)
library(surveillance);library(tidyr); library(lubridate); library(magrittr); library(patchwork); library(lubridate) 
setwd("C:\\Users\\USER\\Desktop\\감염병 모델링")
rm=ls()
load(file = "mal_complt.rda")
#-------------------------------------- 한국 말라리아 케이스 --------------------------------------#
mal_complt <- mal_complt %>% 
  mutate(case_total = zero + ten + twenty + thirty + fourty + fifty + sixty + seventy)

mal <- mal_complt %>%
  filter(year(Date) >=2018 & year(Date)<=2019 ) %>%
  group_by(Date) %>% 
  summarise(case_total = sum(case_total, na.rm = TRUE))
mal%<>% mutate(Date=as.Date(Date))

temp_lastdays<-matrix(NA, ncol=1, nrow=10)
temp_lastdays%<>% as.data.frame()%>% mutate(Date =as.Date((max(mal$Date )+1):(max(mal$Date )+10),
                                                          origin="1970-01-01"), case_total=0) #%>%select(-V1)
temp_lastdays$V1<-NULL

mal<-rbind(mal, temp_lastdays)
mal%<>%mutate(time_onset=0:(nrow(mal)-1)) 
mal%<>%dplyr::select(Date, case_total, time_onset)

origin_date <- as.Date("2018-01-01")
mal <- mal %>%
  mutate(Date = seq.Date(from = origin_date, by = "day", length.out = nrow(.)))

temp_pre_days <- data.frame(
  Date = seq.Date(from = origin_date - 1000, by = "day", length.out = 1000),
  case_total = 0
)

mal <- bind_rows(temp_pre_days, mal)

mal <- mal %>%
  mutate(time_onset = as.integer(difftime(Date, origin_date, units = "days"))) %>%
  dplyr::select(Date, case_total, time_onset)
#head(mal)
#--------------------------------------------------------------------------------------------
# random number로 backprojection
# Bimodal (감마 + 정규 혼합 분포) 
sample_ip_bimodal <- function(n) {
  p <- 83 / (142 + 83)  # normal 쪽 비율
  n_normal <- rbinom(1, n, prob = p)
  n_gamma <- n - n_normal
  
  c(rgamma(n_gamma, shape = 1.2, scale = 22.2),
    rnorm(n_normal, mean = 337.4, sd = 38.5))
}

sample_ip_short <- function(n) {
  rgamma(n, shape = 1.2, scale = 22.2)
}

simulate_infection <- function(df_onset, ip_sampler, n = 100000) {
  # 기준일: time_onset == 0인 날짜 
  origin_date <- df_onset$Date[df_onset$time_onset == 0]
  infection_records <- vector("list", n)
  
  for (sim in 1:n) {
    infection_list <- list()
    
    for (i in 1:nrow(df_onset)) {
      n_cases <- df_onset$case_total[i]
      onset_day <- df_onset$time_onset[i]
      ip_samples <- ip_sampler(n_cases)
      infect_days <- floor(onset_day - ip_samples)
      infection_list[[i]] <- infect_days
    }
    
    all_infect_days <- unlist(infection_list)
    sim_counts <- table(factor(all_infect_days, levels = min(df_onset$time_onset):max(df_onset$time_onset)))
    infection_records[[sim]] <- as.numeric(sim_counts)
  }
  
  infections_mat <- do.call(rbind, infection_records)
  mean_infections <- colMeans(infections_mat)
  
  # 랜덤 시뮬레이션 하나 선택
  set.seed(123)
  random_index <- sample(1:n, 1)
  selected_trajectory <- infection_records[[random_index]]
  
  # 감염일(time_onset 기준)
  infection_days <- min(df_onset$time_onset):max(df_onset$time_onset)
  
  # 해당 감염일에 대한 날짜 생성
  date_seq <- origin_date + infection_days
  
  summary_df <- data.frame(
    date = date_seq,
    time_onset = infection_days,
    mean = mean_infections,
    trajectory = selected_trajectory
  )
  
  return(summary_df)
}
result_bimodal <- simulate_infection(mal, sample_ip_bimodal)
# sum(mal$case_total)
# sum(result_bimodal$mean)
# plot(result_bimodal$mean)
#---------------------------------------------------------------------------------
uni<-function(delay) {
  ifelse(delay >= 4 & delay <= 7, 1 / (7 - 4 + 1), 0)
}

# generation time distribution
generation <- function(r) {
  ifelse(r < 31, 0, incubation_short(r - 30))
}

# short term만 잘라서
inc=list(shape=1.2, scale=22.2)
incubation_short<-function(r){pgamma(r, shape=inc$shape, scale=inc$scale)-(pgamma(r-1, shape=inc$shape, scale=inc$scale))}

# type 3) Bimodal
inc_fit=list(p= 83/(142+83), shape=1.2, scale=22.2, mean=337.4, sd=38.5)
incubation_bimodal<-function(r){(1 - inc_fit$p) * pgamma(r, shape=inc_fit$shape, scale=inc_fit$scale) +
    inc_fit$p * pnorm(r, mean=inc_fit$mean, sd=inc_fit$sd)-((1 - inc_fit$p) * pgamma(r-1, shape=inc_fit$shape, scale=inc_fit$scale) +
                                                              inc_fit$p * pnorm(r-1, mean=inc_fit$mean, sd=inc_fit$sd))}

#--------------------------------------------------------------------------------------------
prepare_matrix<- function(df_infection, incubation_func, Start.T){
  uni_vals <- uni(1:1000)
  incub_vals <- incubation_func(1:1000)
  
  conv <- convolve(uni_vals, rev(incub_vals), type = "open")
  
  today <- max(df_infection$t) - min(df_infection$t)
  
  max_time <- max(df_infection$t)
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
precal_normal<-prepare_matrix(result_bimodal, incubation_bimodal,1)
head(result_bimodal)

result_bimodal%<>%mutate(t=0:(nrow(result_bimodal)-1)) 
#--------------------------------------------------------------------------------------------
max_tau <- max(result_bimodal$t)
gen_cache <- sapply(1:max_tau, generation)
est.t <- list(); est.CI <- list()
Start.T<-1
for (TT in Start.T:max(result_bimodal$t)) {
  dt.backproj.T <- result_bimodal %>% filter(t <= TT)
  
  llk <- function(R) {
    t <- TT
    tau_seq <- 1:(t - 1)
    
    mean_part <- dt.backproj.T$mean[(t):2]
    gen_part <- gen_cache[tau_seq]
    #precal_part <- precal_normal[t - Start.T + 3, tau_seq]
    
    Css <- mean_part * gen_part #/ precal_part
    Cs <- sum(Css) * R
    Cs[Cs <= 0] <- 1e-5
    
    obs <- dt.backproj.T$mean[t + 1]
    return(-(-Cs + obs * log(Cs) - lgamma(obs + 1)))
  }
  
  opt_est <- optim(par = 0.7, fn = llk, method = "L-BFGS-B", lower = 0, control = list(maxit = 1000))
  est.t[[TT]] <- opt_est$par
  
  CI <- function(par_CI) {
    2 * (-llk(par_CI) + opt_est$value)
  }
  
  par_CI_seq <- seq(0, 10, by = 0.01)
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
est_df$t <- Start.T:max(result_bimodal$t)

result <- merge(result_bimodal, est_df, by = "t", all.x = TRUE)
result <- result %>% rename(case_delay = mean)
result_withcut <- result %>% filter(t>=0)
result_withcut$date <- as.Date(result_withcut$date)
range = c(as.Date("2017-01-25"), (as.Date("2019-12-25")-(2)))
scaling_parameter <- max(result_withcut$case_delay) / max(result_withcut$upper, na.rm = TRUE)
result_withcut %>%
  ggplot() +
  geom_bar(aes(x = date, y = case_delay), stat = 'identity', fill = "gray35", width = 0.7, alpha=0.6) +
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
  scale_x_date(date_labels="%m/%d", date_breaks="60 day",
               limits=range, expand=c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1.1 * max(result_withcut$case_delay)),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . / scaling_parameter,
                        breaks = c(0, 1, 2, 3,4),
                        name = "실질감염재생산수\n")
  ) +
  geom_hline(yintercept = 1 * scaling_parameter, linetype = "dashed", color = "#1380A1", size = 0.7)
