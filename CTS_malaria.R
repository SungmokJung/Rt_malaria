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

mal_region_list <- lapply(selected_regions, function(region) {
  Divide_region(mal_complt, region)
})
names(mal_region_list) <- paste0(as.character(selected_regions))
#View(mal_region_list[["23010"]])
#----------------------------------- back projection ---------------------------------------
sample_ip_bimodal <- function(n) {
  p <- 83 / (142 + 83)
  n_normal <- rbinom(1, n, prob = p)
  n_gamma <- n - n_normal
  
  c(rgamma(n_gamma, shape = 1.2, scale = 22.2),
    rnorm(n_normal, mean = 337.4, sd = 38.5))
}

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
  
  summary_df %<>%
    mutate(mean_lag = lag(trajectory)) %>%  
    mutate(
      lag_7 = rollmean(mean_lag, k = 7, align = "right", fill = NA), 
      lag_1week = lag(trajectory, 7) 
    ) %>%
    select(-mean_lag)
  
  weather_vars <- df_onset %>%
    select(Date, tmean, rel_humid, total_prec, total_prec_log) %>%
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
########################## Case-Time Series ########################## 
# 새봄쌤 코드-> tmean: lag 1~6, prec: lag 1~9
# lag_7을 offset으로 (현재는 7일동안의 값을 평균으로. (or 일주일 전 값만 사용할 지...)
# dlnm? cross-basis??
# 아니 이거 어떻게 하는거야!!!!!
# 논문 찾아보니 cross-basis한 결과들을 gnm(y~crossbasis한 x들, eliminate=stratum, data=data, family=binomial)
# 형식으로 사용하던데, offset은 어떻게 처리하고 crossbasis해야 하나?(한다면 tmean, total_prec만 ?)

#-------------------------------------------- 전체 기간 --------------------------------------------#
total_df %<>%mutate(year = year(date),month = month(date))
total_df$stratum <- with(total_df, factor(paste(region, year, month, sep="-")))

# temperature
summary(total_df$tmean)
pct_temp<-quantile(total_df$tmean, prob=c(.01, .10, .33, .66, .90, .99), na.rm=TRUE)
varknot_temp <- pct_temp[c(3, 4)]  # knots at 33rd and 66th

pct_temp_lag_1_6 <- quantile(0:6, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp_lag_1_6 <- pct_temp_lag_1_6[c(3,4)] # knots at 33rd and 66th

argvar_temp <- list(fun="ns", knots=varknot_temp)
arglag_temp_1_6 <- list(fun="ns", knots=varknot_temp_lag_1_6)

cb_tmean <- crossbasis(total_df$tmean, lag=c(0,6), argvar=argvar_temp, arglag=arglag_temp_1_6)

# precipitation -> 0에서부터 331.95로 범위가 너무 넓어 log취함 
summary(total_df$total_prec_log)
pct_preci <- quantile(total_df$total_prec_log, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(5,6)]

pct_preci_lag_1_9 <- quantile(0:6, prob=c(.01,.10,.33,.50,.66,.90,.99),na.rm=T)
varknot_preci_lag_1_9 <- pct_preci_lag_1_9[c(4)]

argvar_preci <- list(fun="ns", knots=varknot_preci)
arglag_preci_1_9 <- list(fun="ns", knots=varknot_preci_lag_1_9)

cb_prec <- crossbasis(total_df$total_prec_log, lag=c(0,6), argvar=argvar_preci, arglag=arglag_preci_1_9) # 얘로 인해 생기는 coef는 무시

model<-gnm(trajectory~cb_tmean+cb_prec+offset(lag_7),eliminate = stratum, data=total_df, family=quasipoisson)
#summary(model)

# Cross-prediction
cr_temp <- crosspred(cb_prec, model)
cpfull <- crosspred(cb_tmean, model)

red<-crossreduce(cb_tmean,model,model.link = "log")
plot(red,ylab="RR", xlab="tmean")


red_p<-crossreduce(cb_prec,model,model.link = "log",cen=4)
plot(red_p,ylab="RR", xlab="log(precipitation)")



#-------------------------------------------- 주요 기간(May-Oct) --------------------------------------------#
# 국내 말라리아는 5~10월 휴전선 접경지역에서 많이 발생 (https://www.kdca.go.kr/contents.es?mid=a20301050303)
monthly_df<-total_df%>%filter(month(date)>=5 & month(date)<=10)
monthly_df %<>%mutate(year = year(date),month = month(date))

monthly_df$stratum <- with(monthly_df, factor(paste(region, year, month, sep="-")))

# temperature
summary(monthly_df$tmean)
pct_temp<-quantile(monthly_df$tmean, prob=c(.01, .10, .33, .66, .90, .99), na.rm=TRUE)
varknot_temp <- pct_temp[c(3, 4)]  # knots at 33rd and 66th

pct_temp_lag_1_6 <- quantile(0:6, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_temp_lag_1_6 <- pct_temp_lag_1_6[c(3,4)] # knots at 33rd and 66th

argvar_temp <- list(fun="ns", knots=varknot_temp)
arglag_temp_1_6 <- list(fun="ns", knots=varknot_temp_lag_1_6)

cb_tmean <- crossbasis(monthly_df$tmean, lag=c(0,6), argvar=argvar_temp, arglag=arglag_temp_1_6)

# precipitation
summary(monthly_df$total_prec_log)
pct_preci <- quantile(monthly_df$total_prec_log, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci <- pct_preci[c(4,5)]

pct_preci_lag_1_9 <- quantile(0:6, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
varknot_preci_lag_1_9 <- pct_preci_lag_1_9[c(4,5)]

argvar_preci <- list(fun="ns", knots=varknot_preci)
arglag_preci_1_9 <- list(fun="ns", knots=varknot_preci_lag_1_9)

cb_prec <- crossbasis(monthly_df$total_prec_log, lag=c(0,6), argvar=argvar_preci, arglag=arglag_preci_1_9) 

## Modei fitting
model_5_10<-gnm(trajectory~cb_tmean+cb_prec+offset(lag_7),eliminate = stratum, data=monthly_df, family=quasipoisson)

## Cross-prediction
cpfull <- crosspred(cb_tmean, model_5_10)
cr_temp <- crosspred(cb_prec, model_5_10)

plot(cpfull)
plot(cr_temp)


red<-crossreduce(cb_tmean,model_5_10,model.link = "log",cen=20)
#beta_series1<-red$fit


plot(red,ylab="RR", xlab="tmean(May to Oct)")


red_p<-crossreduce(cb_prec,model_5_10,model.link = "log",cen=4)
plot(red_p,ylab="RR",xlab="log(precipitation)_(May to Oct)")

