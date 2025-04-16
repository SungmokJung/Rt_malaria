library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus); 
library(tidyr); library(lubridate); library(magrittr); 
#-------------------------------------------------------------------------------
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

sample_ip_normal <- function(n) {
  rnorm(n, mean = 5, sd = 1)
}

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
n<-10000
df_onset<-df_samples
ip_sampler<-sample_ip_bimodal
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

#--------------------------------------------------------------------------------------------------
# trajectory 임의로 선택
# 계속해서 모양 바뀌는 것까지 확인 완 
#set.seed(123)  # 고정을 원하면 주석 제거~
random_index <- sample(1:n, 1)
selected_trajectory <- infection_records[[random_index]]

#  첫 번째 대상자 날짜를 2025-04-01로 맞춤 -> 확인 완 
infection_days <- 0:(length(selected_trajectory) - 1)
first_infection_day <- min(infection_days[selected_trajectory > 0])  # 감염자 발생 첫 날
start_date <- as.Date("2025-04-01")
date_seq <- start_date + (infection_days - first_infection_day)

result_df <- data.frame(
  date = date_seq,
  infections = selected_trajectory
)

ggplot(result_df, aes(x = date, y = infections)) +
  geom_line(color = "blue") +
  labs(title = "Back-projected Infections (One Trajectory)",
       x = "Date",
       y = "Estimated Infections") +
  theme_minimal()
