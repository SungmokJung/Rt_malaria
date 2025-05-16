library(data.table); library(readr); library(dlnm) ; library(gnm) ; library(splines)
library(sf) ; library(terra);library(exactextractr); library(dplyr) ; library(tidyr); library(magrittr)
library(ggplot2) ; library(patchwork); library(tidyterra) ; library(zoo); library(lubridate); library(zoo); library(purrr)
rm=ls()
setwd("C:\\Users\\USER\\Desktop\\감염병 모델링\\CTS(malaria)")
load(file = "mal_complt.rda")
#-------------------------------------- 한국 말라리아 케이스 --------------------------------------#
mal_complt <- mal_complt %>% 
  mutate(case_total = zero + ten + twenty + thirty + fourty + fifty + sixty + seventy)

selected_regions <- c(23010, 23060, 23070, 23080, 23310, 
                      31030, 31080, 31101, 31200, 31230, 
                      31270, 31350, 32010, 32360, 32370, 
                      32380, 32390, 32400)

Divide_region <- function(df, region_code, origin_date = as.Date("2018-01-01")) {
  region_data <- df %>%
    filter(sgg_h == region_code) %>%
    filter(year(Date) >= 2018 & year(Date) <= 2019) %>%
    ungroup() %>%
    group_by(Date) %>%
    summarise(
      case_total = sum(case_total),
      tmean = first(tmean),
      rel_humid = first(rel_humid),
      total_prec = first(total_prec),
      .groups = "drop"
    ) %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(Date) %>%
    mutate(time_onset = 0:(n() - 1),
           Date = seq.Date(from = origin_date, by = "day", length.out = n()))
  
  # 400일 전 빈 데이터 추가
  temp_pre_days <- data.frame(
    Date = seq.Date(from = origin_date - 400, by = "day", length.out = 400),
    case_total = 0,
    tmean = NA,
    rel_humid = NA,
    total_prec = NA
  )
  
  region_data <- bind_rows(temp_pre_days, region_data) %>%
    mutate(time_onset = as.integer(difftime(Date, origin_date, units = "days")))
  
  return(region_data)
}

mal_region_list <- lapply(selected_regions, function(region) {
  Divide_region(mal_complt, region)
})
names(mal_region_list) <- paste0("mal_", selected_regions)
#----------------------------------- back projection ---------------------------------------
sample_ip_bimodal <- function(n) {
  p <- 83 / (142 + 83)
  n_normal <- rbinom(1, n, prob = p)
  n_gamma <- n - n_normal
  
  c(rgamma(n_gamma, shape = 1.2, scale = 22.2),
    rnorm(n_normal, mean = 337.4, sd = 38.5))
}

# 박사님 코드 보고 수정할 것....
simulate_infection <- function(df_onset, ip_sampler, n = 100000) {
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
  
  set.seed(123)
  random_index <- sample(1:n, 1)
  selected_trajectory <- infection_records[[random_index]]
  
  infection_days <- min(df_onset$time_onset):max(df_onset$time_onset)
  
  date_seq <- origin_date + infection_days
  
  summary_df <- data.frame(
    date = date_seq,
    time_onset = infection_days,
    mean = mean_infections,
    trajectory = selected_trajectory
  )
  
  summary_df%<>%mutate(mean_lag = lag(mean)) %>%  # 현재 값 제외 -> 맞겠지...?
    mutate(lag_7 = rollmean(mean_lag, k = 7, align = "right", fill = NA)) %>%
    select(-mean_lag)
  
  weather_vars <- df_onset %>%
    select(Date, tmean, rel_humid, total_prec) %>%
    rename(date = Date) 
  
  summary_df <- left_join(summary_df, weather_vars, by = "date")
  
  return(summary_df)
}

#-----------------------------------------------------------------------------------
simulation_results <- lapply(mal_region_list, function(df) {
  simulate_infection(df_onset = df, ip_sampler = sample_ip_bimodal, n = 100000)
})
names(simulation_results) <- names(mal_region_list)

# 하나로 합치기~
total_df <- bind_rows(
  map2(simulation_results, names(simulation_results), ~ mutate(.x, region = .y))
)
#-----------------------------------------------------------------------------------
