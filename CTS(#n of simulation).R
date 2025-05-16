library(data.table); library(readr); library(dlnm) ; library(gnm) ; library(splines)
library(sf) ; library(terra);library(exactextractr); library(dplyr) ; library(tidyr); library(magrittr)
library(ggplot2) ; library(patchwork); library(tidyterra) ; library(zoo); library(lubridate); library(zoo); library(purrr)
rm=ls()
setwd("C:\\Users\\USER\\Desktop\\감염병 모델링\\CTS(malaria)")
load(file = "mal_complt.rda")
#-------------------------------------- 한국 말라리아 케이스 --------------------------------------#
mal_complt <- mal_complt %>% 
  mutate(case_total = zero + ten + twenty + thirty + fourty + fifty + sixty + seventy)
mal_complt$total_prec_log <- log(mal_complt$total_prec + 0.001)


selected_regions <- c(23010, 23060, 23070, 23080, 23310, 
                      31030, 31080, 31101, 31200, 31230, 
                      31270, 31350, 32010, 32360, 32370, 
                      32380, 32390, 32400)

Divide_region <- function(df, region_code, origin_date = as.Date("2018-01-01")) {
  
  # join 용
  semi_region <- df %>%
    filter(sgg_h == region_code) %>%
    ungroup() %>%
    select(Date, tmean, rel_humid, total_prec, total_prec_log) %>%
    distinct()
  
  # @@@일 전 환경변수 있는 날짜 추출
  temp_pre_days <- data.frame(
    Date = seq.Date(from = as.Date("2018-01-01") - 300, by = "day", length.out = 300)
  ) %>%
    left_join(semi_region, by = "Date") %>%
    mutate(case_total = 0) %>%
    select(Date,case_total, tmean, rel_humid, total_prec, total_prec_log)
  
  
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
      total_prec_log=first(total_prec_log),
      .groups = "drop"
    ) %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(Date) %>%
    mutate(time_onset = 0:(n() - 1),
           Date = seq.Date(from = origin_date, by = "day", length.out = n()))
  
  
  region_data <- bind_rows(temp_pre_days, region_data) %>%
    mutate(time_onset = as.integer(difftime(Date, origin_date, units = "days")))
  
  return(region_data)
}

mal_region_list2 <- lapply(selected_regions, function(region) {
  Divide_region(mal_complt, region)
})
names(mal_region_list2) <- paste0(as.character(selected_regions))
#----------------------------------- back projection ---------------------------------------
# Bimodal 잠복기 분포
sample_ip_bimodal <- function(n) {
  p <- 83 / (142 + 83)
  n_normal <- rbinom(1, n, prob = p)
  n_gamma <- n - n_normal
  c(rgamma(n_gamma, shape = 1.2, scale = 22.2),
    rnorm(n_normal, mean = 337.4, sd = 38.5))
}

simulate_all_regions <- function(region_list, ip_sampler, n = 10000) {
  region_names <- names(region_list)
  
  results_tmean <- list()
  results_prec <- list()
  
  for (sim in 1:n) {
    message("Simulation ", sim, " / ", n)
    
    simulation_results <- map2(region_list, region_names, function(df, region_name) {
      origin_date <- df$Date[df$time_onset == 0]
      infection_list <- vector("list", length = nrow(df))
      
      for (i in seq_len(nrow(df))) {
        n_cases <- df$case_total[i]
        onset_day <- df$time_onset[i]
        ip_samples <- ip_sampler(n_cases)
        infect_days <- floor(onset_day - ip_samples)
        infection_list[[i]] <- infect_days
      }
      
      all_infect_days <- unlist(infection_list)
      sim_counts <- table(factor(all_infect_days, levels = min(df$time_onset):max(df$time_onset)))
      sim_counts <- as.numeric(sim_counts)
      
      infection_days <- min(df$time_onset):max(df$time_onset)
      date_seq <- df$Date[df$time_onset == 0] + infection_days
      
      summary_df <- data.frame(
        date = date_seq,
        time_onset = infection_days,
        trajectory = sim_counts
      ) %>%
        mutate(mean_lag = lag(trajectory)) %>%
        mutate(
          lag_7 = zoo::rollmean(mean_lag, k = 7, align = "right", fill = NA),
          lag_1week = lag(mean_lag, 7)
        ) %>%
        select(-mean_lag)
      
      weather_vars <- df %>%
        select(Date, tmean, rel_humid, total_prec, total_prec_log) %>%
        rename(date = Date)
      
      summary_df <- left_join(summary_df, weather_vars, by = "date") %>%
        mutate(region = region_name)
      
      return(summary_df)
    })
    
    total_df <- bind_rows(simulation_results)
    View(total_df)
    total_df %<>%
      mutate(year = year(date), month = month(date),
             stratum = factor(paste(region, year, month, sep = "-")))
    
    # DLNM 설정
    # temperature
    pct_temp <- quantile(total_df$tmean, prob = c(.01, .10, .33, .66, .90, .99), na.rm = TRUE)
    varknot_temp <- pct_temp[c(3, 4)]
    pct_temp_lag <- quantile(0:6, prob = c(.01, .10, .33, .66, .90, .99), na.rm = TRUE)
    varknot_temp_lag <- pct_temp_lag[c(3, 4)]
    
    argvar_temp <- list(fun = "ns", knots = varknot_temp)
    arglag_temp <- list(fun = "ns", knots = varknot_temp_lag)
    cb_tmean <- crossbasis(total_df$tmean, lag = c(0,6), argvar = argvar_temp, arglag = arglag_temp)
    
    # precipitation
    pct_preci <- quantile(total_df$total_prec_log, prob = c(.01, .10, .33, .66, .90, .99), na.rm = TRUE)
    varknot_preci <- pct_preci[c(5,6)]
    pct_preci_lag <- quantile(0:6, prob = c(.01, .10, .33, .50, .66, .90, .99), na.rm = TRUE)
    varknot_preci_lag <- pct_preci_lag[c(4)]
    
    argvar_preci <- list(fun = "ns", knots = varknot_preci)
    arglag_preci <- list(fun = "ns", knots = varknot_preci_lag)
    cb_prec <- crossbasis(total_df$total_prec_log, lag = c(0,6), argvar = argvar_preci, arglag = arglag_preci)
    
    model <- gnm(trajectory ~ cb_tmean + cb_prec + offset(lag_7),eliminate = stratum, data = total_df, family = quasipoisson)
    
    red <- crossreduce(cb_tmean, model, model.link = "log")
    plot(red)
    red_p <- crossreduce(cb_prec, model, model.link = "log")
    plot(red_p)
      
    results_tmean[[sim]] <- red$fit
    results_prec[[sim]] <- red_p$fit
    }
  
  # 데이터프레임으로 정리
  beta_df_tmean <- do.call(rbind, results_tmean) %>% as.data.frame()
  beta_df_prec <- do.call(rbind, results_prec) %>% as.data.frame()
  
  beta_df_tmean$sim <- 1:nrow(beta_df_tmean)
  beta_df_prec$sim <- 1:nrow(beta_df_prec)
  
  return(list(beta_tmean = beta_df_tmean, beta_prec = beta_df_prec))
}

result <- simulate_all_regions(mal_region_list2, sample_ip_bimodal, n = 1)
result_tmean<-result$beta_tmean
result_prec<-result$beta_prec

#------------------------------------- tmean 결과 
summary_tmean <- apply(result_tmean[, 1:(ncol(result_tmean) - 1)], 2, function(x) {
  c(median = median(x, na.rm = TRUE),
    lower = quantile(x, 0.025, na.rm = TRUE),
    upper = quantile(x, 0.975, na.rm = TRUE))
}) %>%
  t() %>%                 # 전치해서 행: 변수, 열: 통계량
  as.data.frame() %>%     
  tibble::rownames_to_column(var = "tmean")  
colnames(summary_tmean)

summary_tmean$tmean <- as.numeric(as.character(summary_tmean$tmean))

ggplot(summary_tmean, aes(x = tmean, y = median, group = 1)) +
  geom_line(color = "blue", size = 1.2) +          
  geom_point(color = "blue", size = 0.6) +           
  geom_ribbon(aes(ymin = `lower.2.5%`, ymax = `upper.97.5%`), 
              fill = "lightblue", alpha = 0.4) +    
  labs(title = "Simulation result of DLNMs",
       x = "Tmean", 
       y = "RR") +
  theme_minimal()

#------------------------------------- prec 결과 
summary_prec <- apply(result_prec[, 1:(ncol(result_prec) - 1)], 2, function(x) {
  c(median = median(x, na.rm = TRUE),
    lower = quantile(x, 0.025, na.rm = TRUE),
    upper = quantile(x, 0.975, na.rm = TRUE))
}) %>%
  t() %>%               
  as.data.frame() %>%     
  tibble::rownames_to_column(var = "prec")  
colnames(summary_prec)

summary_prec$prec <- as.numeric(as.character(summary_prec$prec))

ggplot(summary_prec, aes(x = prec, y = median, group = 1)) +
  geom_line(color = "blue", size = 1.2) +          
  geom_point(color = "blue", size = 0.6) +           
  geom_ribbon(aes(ymin = `lower.2.5%`, ymax = `upper.97.5%`), 
              fill = "lightblue", alpha = 0.4) +      
  labs(title = "Simulation result of DLNMs",
       x = "Precipitation", 
       y = "RR") +
  theme_minimal()