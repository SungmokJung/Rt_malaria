library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus); 
library(surveillance);library(tidyr); library(lubridate); library(magrittr); library(patchwork)  

# 1. Toy data 
# --------------------------------------------------------------------------------------------
days <- 0:350  

prob_density <- dnorm(days, mean = 250, sd = 25)  # 평균 250, 표준편차 100
case_counts <- prob_density / max(prob_density) * 100  # 최댓값을 100으로 정규화

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
#plot(df_samples$case_count)
sum(df_samples$case_count)

# 2. IP 분포 함수들 
# --------------------------------------------------------------------------------------------
# r~ 쓰면 해당 분포에서 랜덤 샘플링 가능!
# normal 잠복기 (5일)
sample_ip_normal <- function(n) {
  rnorm(n, mean = 5, sd = 1)
}

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
#sample_ip_bimodal(10)

# 3. Back-projection 함수
# --------------------------------------------------------------------------------------------
# 확인 후 bimodal distribution으로 backprojection한 논문 있는지 문헌 고찰 필요 

simulate_infection <- function(df_onset, ip_sampler, n = 1000) {
  infection_records <- vector("list", n)
  
  for (sim in 1:n) {
    infection_list <- list()
    
    for (i in 1:nrow(df_samples)) {
      n_cases <- df_samples$case_count[i]
      onset_day <- df_samples$time_onset[i]
      ip_samples <- ip_sampler(n_cases) #해당 환자 수만큼 잠복기 샘플 생성
      infect_days <- floor(onset_day - ip_samples)
      infection_list[[i]] <- infect_days
      # if (n_cases >= 0) {
      #   ip_samples <- ip_sampler(n_cases) #해당 환자 수만큼 잠복기 샘플 생성
      #   infect_days <- floor(onset_day - ip_samples)
      #   infection_list[[i]] <- infect_days
      # }
    }
    
    all_infect_days <- unlist(infection_list)
    sim_counts <- table(factor(all_infect_days, levels = 0:max(df_samples$time_onset)))
    infection_records[[sim]] <- as.numeric(sim_counts)
  }
  
  infections_mat <- do.call(rbind, infection_records)
  #View(infections_mat)

  summary_df <- data.frame(
    time_infection = 0:max(df_samples$time_onset),
    mean = colMeans(infections_mat)
  )
  
  return(summary_df)
}

# --------------------------------------------------------------------------------------------
#sum(df_samples$case_count) #24754
# Normal IP
result_normal <- simulate_infection(df_samples, sample_ip_normal)
sum(result_normal$mean)

# short-term IP
result_gamma <- simulate_infection(df_samples, sample_ip_gamma)
sum(result_gamma$mean)

# Bimodal IP
result_bimodal <- simulate_infection(df_samples, sample_ip_bimodal)
sum(result_bimodal$mean)

# --------------------------------------------------------------------------------------------
# 원래 케이스 분포랑 같이 보여줄 것.
# 위에는 원래 케이스분포, 아래에는 Backprojection한 분포 

plot_backproj <- function(backproj_result, df_samples, title) {
  
  y_max <- max(max(df_samples$case_count), max(backproj_result$mean, na.rm = TRUE))
  
  p1 <- ggplot(df_samples, aes(x = time_onset, y = case_count)) +
    geom_col(fill = "gray50") +
    coord_cartesian(ylim = c(0, y_max)) +  
    labs(title = "Observed onset cases",
         x = "Symptom onset day", y = "Case count") +
    theme_minimal()
  
  p2 <- ggplot(backproj_result, aes(x = time_infection, y = mean)) +
    geom_col(fill = "firebrick") +
    coord_cartesian(ylim = c(0, y_max)) +  
    labs(title = title,
         x = "Estimated infection day", y = "Infection count (simulated)") +
    theme_minimal()
  
  return(p1 / p2 + plot_layout(heights = c(1, 1.2)))
}


# --------------------------------------------------------------------------------------------
plot_backproj(result_normal, df_samples, "Back-projection Result(normal)")
plot_backproj(result_gamma, df_samples, "Back-projection Result(gamma(short-term))")
plot_backproj(result_bimodal, df_samples, "Back-projection Result(bimodal)")

# --------------------------------------------------------------------------------------------
# Rt 추정 

uni<-function(delay) {
  ifelse(delay >= 4 & delay <= 7, 1 / (7 - 4 + 1), 0)
}

# generation time distribution
gi_fit=list(shape=2.305, scale=5.452)
generation<-function(t){pweibull(t, shape=gi_fit$shape, scale=gi_fit$scale)-
    pweibull(t-1, shape=gi_fit$shape, scale=gi_fit$scale)}

# incubation distribution
# type 1) normal
incubation_normal <- function(n) {dnorm(n, mean = 5, sd = 1)}

# type 2) short(gamma)
incubation_short<-function(r){dgamma(r, shape=1.2, scale=22.2)}

# type 3) Bimodal
inc_fit=list(p= 83/(142+83), shape=1.2, scale=22.2, mean=337.4, sd=38.5)
incubation_bimodal<-function(r){(1 - inc_fit$p) * dgamma(r, shape=inc_fit$shape, scale=inc_fit$scale) + 
    inc_fit$p * dnorm(r, mean=inc_fit$mean, sd=inc_fit$sd)}

#head(result_bimodal)

delay_precalculation <- function(t){sum(convolve(uni(1:(today-t)),
                                                 rev(incubation(1:(today-t))), type = c("open")))}

prepare_matrix <- function(df_infection, incubation_func, Start.T) {
  today <- 760 - min(df_infection$time_infection)
  
  conv <- function(t, tau) {
    sum(convolve(uni(1:1000), rev(incubation_func(1:1000)), type = "open")[1:(today - t + tau)])
  }
  
  max_time <- max(df_infection$time_infection)
  precal <- matrix(0, nrow = max_time - Start.T + 3, ncol = max_time - 1)
  
  for (m in (Start.T - 2):max_time) {
    for (n in 1:(m - 1)) {
      message(sprintf("m: %d, n: %d", m, n))
      precal[m - Start.T + 3, n] <- conv(m, n)
    }
  }
  return(precal)
}

precal_normal <- prepare_matrix(result_normal, incubation_normal,300)
precal_short <- prepare_matrix(result_gamma, incubation_short,300)
precal_bimodal <- prepare_matrix(result_bimodal, incubation_bimodal,200)
head(result_normal,6)

# --------------------------------------------------------------------------------------------
#result_normal,result_gamma,result_bimodal
# 🤔🤔🤔🤔
Start.T<-200 
est.t<-list(); est.CI<-list();
for (TT in Start.T:max(result_bimodal$time_infection)){
  result_bimodal %>% filter(time_infection <= TT) -> dt.backproj.T
  
  llk <- function(param){t=TT; R = param
  llk <- 0; Cs <- 0 ;　Css <- rep(0, t)
  for (tau in 1: t - 1){
    print(paste("TT:", TT, "tau:", tau))
    Css[tau] = (dt.backproj.T$mean[t-tau+1]) * generation(tau) /precal_bimodal[t - Start.T + 3, tau]
  }
  Cs = sum(Css) * R
  Cs[Cs <= 0] <- 1e-5
  return(-(-Cs+dt.backproj.T$mean[t+1] * log(Cs) - lgamma(dt.backproj.T$mean[t + 1] + 1)))
  }
  #param0 = c(0.7)
  param0=c(1e-5)
  opt_est <- optim(param0, fn=llk, method=c("L-BFGS-B"), lower=c(0), control=list(maxit = 10000))
  opt_est$par -> est.t[[TT]]
  
  ## 프로파일 우도를 활용한 95% 신뢰구간 추정 
  ci_pro <- matrix(NA, ncol = 2, nrow = 1)
  CI <- function(par_CI){return(2 * (-llk(par_CI)+opt_est$value))} #llk()
  par_CI <- seq(0, 20, by = 0.01)
  logLik <- sapply(par_CI, FUN = CI)
  as.data.frame(par_CI) -> par_CI; as.data.frame
  (logLik) -> logLik
  cbind(par_CI, logLik) -> data_CI
  data_CI$logLik[data_CI$logLik < (max(data_CI$logLik) - 3.84)] <- NA
  data_CI %<>% na.omit()
  min(data_CI$par_CI) -> ci_pro[1,1]; max(data_CI$par_CI) -> ci_pro[1,2]
  as.data.frame(ci_pro) -> ci_pro
  colnames(ci_pro) <- c("lower","upper")
  ci_pro -> est.CI[[TT]]
}
matrix(unlist(est.t),ncol=1,byrow=T) -> est.t
matrix(unlist(est.CI),ncol=2,byrow=T) -> est.CI
cbind(est.t, est.CI) -> est

est %>% as.data.frame() %>% mutate(time_infection = Start.T:max(result_bimodal$time_infection)) -> result
merge(result_bimodal, result, by='time_infection', all.x=TRUE) -> result
colnames(result) <- c("time_onset","case_delay","Rt","lower", "upper")
result %>% filter(time_onset>=200 &time_onset <= max(time_onset)-50) -> result_withcut
plot(result$Rt)
options(repr.plot.width = 7, repr.plot.height = 5)
theme_set(theme_bw())
scaling_parameter = max(result_withcut$case_delay) / max(result_withcut$upper[!is.na(result_withcut$upper)])

options(warn=-1) 
result_withcut %>% 
  ggplot() +
  geom_bar(aes(x=time_onset, y=case_delay), stat='identity', fill="gray35", width=0.7) + 
  geom_line(data=result_withcut[!is.na(result_withcut$Rt),],aes(x=time_onset ,y=Rt * scaling_parameter), color="#1380A1", size=1) +
  geom_ribbon(data=result_withcut, aes(ymax=result_withcut$upper * scaling_parameter, 
                                       ymin=result_withcut$lower * scaling_parameter, x=time_onset ), 
              fill="#1380A1", alpha = 0.4) + labs(x="\n 추정감염일", y=" 감염자 수\n") +
  theme (text = element_text ( size = 12 , family=
                                 "sans",color = "black"),
         axis.text = element_text (size = 12 , family="sans",color = "black"),
         panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
         legend.position="none") +
  scale_y_continuous(limit=c(0,60), expand = c(0, 0),
                     sec.axis=sec_axis(~./(scaling_parameter), breaks=c(0,2,4), name =" 실질감염재생산수\n")) + 
  geom_hline(yintercept=1 * scaling_parameter,
             linetype="dashed", color = "#1380A1", size =0.7)