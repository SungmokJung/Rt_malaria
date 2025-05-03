# Sys.which("make")
# #writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# pkgbuild::check_build_tools(debug = TRUE)
# 
# dotR <- file.path(Sys.getenv("HOME"), ".R")
# if (!file.exists(dotR)) dir.create(dotR)
# M <- file.path(dotR, "Makevars.win")
# if (!file.exists(M)) file.create(M)
# cat("CXX14FLAGS=-O3 -march=native -mtune=native -fPIC",
#     "CXX14=g++",
#     file = M, sep = "\n")
# 
# install.packages("RcppEigen", type = "source")
# install.packages("StanHeaders", type = "source")
# install.packages("loo", type = "source")
# install.packages("posterior", type = "source")
# install.packages("rstan", type = "source")
# install.packages("EpiNow2", repos = "https://epiforecasts.r-universe.dev", type = "source")
library(EpiNow2);library(distributional)
library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus); 
library(surveillance);library(tidyr); library(lubridate); library(magrittr); library(patchwork)
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
mal<-mal%>%rename(confirm=case_total, date=Date)

# temp_lastdays<-matrix(NA, ncol=1, nrow=10)
# temp_lastdays%<>% as.data.frame()%>% mutate(date =as.Date((max(mal$date )+1):(max(mal$date )+10),origin="1970-01-01"), confirm=0)
# temp_lastdays$V1<-NULL
# 
# origin_date <- as.Date("2018-01-01")
# mal <- mal %>%
#   mutate(date = seq.Date(from = origin_date, by = "day", length.out = nrow(.)))
# 
# temp_pre_days <- data.frame(
#   date = seq.Date(from = origin_date - 1000, by = "day", length.out = 1000),
#   confirm = 0
# )
# 
# mal <- bind_rows(temp_pre_days, mal,temp_lastdays)
# mal <- mal %>%
#   mutate(time_onset = as.integer(difftime(date, origin_date, units = "days"))) %>%
#   dplyr::select(date, confirm, time_onset)
# 
# 
# ### backprojection
# sample_ip_bimodal <- function(n) {
#   p <- 83 / (142 + 83)  # normal 쪽 비율
#   n_normal <- rbinom(1, n, prob = p)
#   n_gamma <- n - n_normal
#   
#   c(rgamma(n_gamma, shape = 1.2, scale = 22.2),
#     rnorm(n_normal, mean = 337.4, sd = 38.5))
# }
# 
# sample_ip_short <- function(n) {
#   rgamma(n, shape = 1.2, scale = 22.2)
# }
# 
# simulate_infection <- function(df_onset, ip_sampler, n = 100000) {
#   # 기준일: time_onset == 0인 날짜 
#   origin_date <- df_onset$date[df_onset$time_onset == 0]
#   infection_records <- vector("list", n)
#   
#   for (sim in 1:n) {
#     infection_list <- list()
#     
#     for (i in 1:nrow(df_onset)) {
#       n_cases <- df_onset$confirm[i]
#       onset_day <- df_onset$time_onset[i]
#       ip_samples <- ip_sampler(n_cases)
#       infect_days <- floor(onset_day - ip_samples)
#       infection_list[[i]] <- infect_days
#     }
#     
#     all_infect_days <- unlist(infection_list)
#     sim_counts <- table(factor(all_infect_days, levels = min(df_onset$time_onset):max(df_onset$time_onset)))
#     infection_records[[sim]] <- as.numeric(sim_counts)
#   }
#   
#   infections_mat <- do.call(rbind, infection_records)
#   mean_infections <- colMeans(infections_mat)
#   
#   # 랜덤 시뮬레이션 하나 선택
#   set.seed(123)
#   random_index <- sample(1:n, 1)
#   selected_trajectory <- infection_records[[random_index]]
#   
#   # 감염일(time_onset 기준)
#   infection_days <- min(df_onset$time_onset):max(df_onset$time_onset)
#   
#   # 해당 감염일에 대한 날짜 생성
#   date_seq <- origin_date + infection_days
#   
#   summary_df <- data.frame(
#     date = date_seq,
#     time_onset = infection_days,
#     mean = as.integer(round(mean_infections)),
#     trajectory = as.integer(round(selected_trajectory))
#   )
#   
#   return(summary_df)
# }
# 
# result_short<-simulate_infection(mal, sample_ip_short)
# short_mean<-result_short%>%dplyr::select(date, mean)
# short_mean<-short_mean%>%rename(confirm=mean, date=date)
# short_trajectory<-result_short%>%dplyr::select(date, trajectory)
# short_trajectory<-short_trajectory%>%rename(confirm=trajectory, date=date)
# 
# result_bimodal <- simulate_infection(mal, sample_ip_bimodal)
# bimodal_mean<-result_bimodal%>%dplyr::select(date, mean)
# bimodal_mean<-bimodal_mean%>%rename(confirm=mean, date=date)
# 
# bimodal_trajectory<-result_bimodal%>%dplyr::select(date, trajectory)
# bimodal_trajectory<-bimodal_trajectory%>%rename(confirm=trajectory, date=date)
#--------------------------------------------------------------------------------------------
delay_pmf <- rep(0, 11)  
delay_pmf[4:7] <- 1 / 4  
#plot(delay_pmf)
example_non_parametric_delay <- NonParametric(pmf = delay_pmf)

incubation_period <- Gamma(mean = 26.64, sd = 24.30, max = 150)

inc_fit=list(p= 83/(142+83), shape=1.2, scale=22.2, mean=337.4, sd=38.5)
incubation_bimodal<-function(r){(1 - inc_fit$p) * pgamma(r, shape=inc_fit$shape, scale=inc_fit$scale) +
    inc_fit$p * pnorm(r, mean=inc_fit$mean, sd=inc_fit$sd)-((1 - inc_fit$p) * pgamma(r-1, shape=inc_fit$shape, scale=inc_fit$scale) +
                                                              inc_fit$p * pnorm(r-1, mean=inc_fit$mean, sd=inc_fit$sd))}
days <- 0:500
pmf <- sapply(days, incubation_bimodal)
pmf[pmf < 0] <- 0
pmf <- pmf / sum(pmf)
plot(pmf)
incubation_bimodal_period <- NonParametric(pmf = pmf)

generation_time <- Gamma(
  mean = 26.64 + 30,  sd = 24.30,max = 180)

shift_days <- 30
gen_pmf <- c(rep(0, shift_days), pmf)[1:length(pmf)]
gen_pmf <- gen_pmf / sum(gen_pmf)
generation_bimodal <- NonParametric(pmf = gen_pmf)

delay_pmf <- convolve(pmf, rev(example_non_parametric_delay$pmf), type = "open")
delay_pmf[delay_pmf < 0] <- 0
delay_pmf <- delay_pmf / sum(delay_pmf)

delay_bimodal <- bound_dist(NonParametric(pmf = delay_pmf))
#plot(delay_pmf)

#--------------------------------------------------------------------------------------------
### short만 잘라서 
estimates <- epinow(
  data = mal,
  generation_time = generation_time_opts(generation_time),
  delays = delay_opts(example_non_parametric_delay + incubation_period),
  rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
  stan = stan_opts(cores = 4),
  verbose = interactive()
)
plot(estimates)

### bimodal 
estimates_bimodal <- epinow(
  data = mal,
  generation_time = generation_time_opts(generation_bimodal),
  delays = delay_opts(delay_bimodal),
  rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
  stan = stan_opts(cores = 4),
  verbose = interactive()
)
plot(estimates_bimodal)
#--------------------------------------------------------------------------------------------
# # 미리 backprojection한 결과 사용 -> X
# short_mean<-short_mean%>%filter(date>='2018-01-01')
# estimates_short_mean <- epinow(
#   data = short_mean,
#   generation_time = generation_time_opts(generation_time),
#   delays = delay_opts(example_non_parametric_delay + incubation_period),
#   rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
#   stan = stan_opts(cores = 4, control = list(adapt_delta = 0.95)),
#   verbose = interactive()
# )
# plot(estimates_short_mean)
# 
# short_trajectory<-short_trajectory%>%filter(date>='2018-01-01')
# estimates_short_trajectoty <- epinow(
#   data = short_trajectory,
#   generation_time = generation_time_opts(generation_time),
#   delays = delay_opts(example_non_parametric_delay + incubation_period),
#   rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
#   stan = stan_opts(cores = 4, control = list(adapt_delta = 0.95)),
#   verbose = interactive()
# )
# plot(estimates_short_trajectoty)
# 
# estimates_bimodal_mean <- epinow(
#   data = bimodal_mean,
#   generation_time = generation_time_opts(generation_bimodal),
#   delays = delay_opts(delay_bimodal),
#   rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
#   stan = stan_opts(cores = 4, control = list(adapt_delta = 0.95)),
#   verbose = interactive()
# )
# plot(estimates_bimodal_mean)
# 
# estimates_bimodal_trajectoty <- epinow(
#   data = bimodal_trajectory,
#   generation_time = generation_time_opts(generation_bimodal),
#   delays = delay_opts(delay_bimodal),
#   rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
#   stan = stan_opts(cores = 4, control = list(adapt_delta = 0.95)),
#   verbose = interactive()
# )
# plot(estimates_short_trajectoty)
# 
