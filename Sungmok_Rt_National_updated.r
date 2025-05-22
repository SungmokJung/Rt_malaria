libraries = c("dplyr", "magrittr", "tidyr", "ggplot2", "RColorBrewer", "zoo", "lubridate", "tidyverse")
for(x in libraries) {library(x,character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)}
theme_set(theme_bw())

#### The code is written by Sung-mok Jung

## load data
load(file = "../data/mal_complt.rda")
mal_complt %>% mutate(case_total = zero + ten + twenty + thirty + fourty + fifty + sixty + seventy) %>%
filter(year(Date) >= 2018 & year(Date) <= 2019 ) %>% group_by(Date) %>% 
summarise(case_total = sum(case_total, na.rm = TRUE)) -> mal_example

## sanity check (comparing the annual cumulative cases with the KDCA report)
mal_example %>% mutate(year_month = format(as.Date(Date), "%Y-%m")) %>%
group_by(year_month) %>% summarise(case_month = sum(case_total, na.rm = TRUE)) -> check_temp

## incubation period
p <- 83/(142+83)
incubation_pdf = (1-p)*dgamma(0:501, shape = 1.2, scale = 22.2) + p*dnorm(0:501, mean = 337.4, sd = 38.5)

options(repr.plot.width = 8,repr.plot.height = 6)
data.frame(prop = incubation_pdf, t = 0:501) %>%
ggplot() + 
geom_line(aes(x = t, y = prop), linewidth = 1) +
labs(x = "\n Incubation period (days)", y = "Probability \n") +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_continuous(expand = c(0, 0)) + 
scale_y_continuous(expand = c(0, 0)) +
coord_cartesian(ylim = c(0, 0.020)) 

## backprojection for a single trajectory
set.seed(123)

niter <- 1
min_date <- as.Date("2016-01-01")
backproj <- list()

for(i in 1: niter){
    mal_example %>% rowwise() %>%
    do({
        n_cases <- .$case_total
        if (n_cases > 0) {
            sampled_days <- sample(0:501, size = n_cases, replace = TRUE, prob = incubation_pdf)
            data.frame(infection_date = .$Date - days(sampled_days))
        } else {
            data.frame(infection_date = as.Date(character(0)))
        }
    }) %>% ungroup() %>% 
    count(infection_date, name = "case_total") %>%
    arrange(infection_date) -> temp
    data.frame(infection_date = as.Date(min_date:max(mal_example$Date))) -> cal_temp
    merge(cal_temp, temp, by = "infection_date", all.x = TRUE) %>% 
    mutate(t = 0:(nrow(cal_temp)-1)) -> backproj[[i]]
    backproj[[i]]$case_total[is.na(backproj[[i]]$case_total)] <- 0
}

## plotting a single trajectory
options(repr.plot.width = 10,repr.plot.height = 5, warn = -1)
range = c(as.Date("2016-09-01"), as.Date("2020-01-01"))

mal_example %>% 
ggplot(aes(x = Date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of reporting", y = "Incidence \n")

backproj[[1]] %>% 
ggplot(aes(x = infection_date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of infection", y = "Incidence \n")

## what if the generation time distribution is identical with the short-term incubation period distribution
generation <- function(t){pgamma(t, shape = 1.2, scale = 22.2) -
                          pgamma(t-1, shape = 1.2, scale = 22.2)}

## what if we borrow the generation time distribution of P. falciparum
t_grid <- 0:500

f1 <- dlnorm(t_grid, meanlog = 2.38, sdlog = 0.254)
f2 <- dgamma(t_grid, shape = 1.19, scale = 0.016)

f1p <- approx(t_grid, f1, xout = t_grid - 9, yleft = 0, yright = 0)$y
f2n <- approx(t_grid, f2, xout = t_grid - 10, yleft = 0, yright = 0)$y

conv12 <- fft(f1p) * fft(f2n)
fV <- Re(fft(conv12, inverse = TRUE)) / length(t_grid)

fV[fV < 0] <- 0
fV <- fV / sum(fV)

generation_fal <- function(t) {t <- floor(t)
                               ifelse(t < 0 | t > length(fV) - 1, 0, fV[t + 1])}

## comparison between two different generation time distribution assumptions
data.frame(prop = generation(0:101), t = 0:101) %>%
ggplot() + 
geom_line(aes(x = t, y = prop), linewidth = 1, color= "steelblue") +
geom_line(data = data.frame(prop = generation_fal(0:101), t = 0:101),
          aes(x = t, y = prop), linewidth = 1, color = "indianred2") +
labs(x = "\n Generatoion time (days)", y = "Probability \n") +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_continuous(expand = c(0, 0)) + 
scale_y_continuous(expand = c(0, 0)) +
coord_cartesian(ylim = c(0, 0.2)) 

## setting for Rt estimation
backproj[[1]] -> data
data %>% filter(case_total > 0) %>% arrange(infection_date) %>% slice(1) %>% pull(t) -> T_min
max(data$infection_date) -> T_max

## Rt estimation with a single trajectory
gen_vec <- sapply(1:(T_max + 1), generation_fal)

Start.T <- T_min-10

est_mod <-　list(); est_CI <- list()

for (TT in Start.T:T_max) {
    data_T <- data[data$t <= TT, ]
    t <- TT

    if ((t + 1) > nrow(data_T)) next
    Css <- rep(0, t)
    for (tau in 1:(t - 1)) {
        Css[tau] <- data_T$case_total[t - tau + 1] * gen_vec[tau]
    }

    llk <- function(R) {
        Cs <- sum(Css) * R
        if (Cs <= 0) Cs <- 1e-5
        obs <- data_T$case_total[t + 1]
        -(-Cs + obs * log(Cs) - lgamma(obs + 1))
    }

    opt_est <- optim(0.7, llk, method = "L-BFGS-B", lower = 0, control = list(maxit = 10000))
    est_mod[[TT]] <- opt_est$par
    
    par_CI_seq <- seq(0, 10, by = 0.01)
    CI <- function(par_CI) {2 * (-llk(par_CI) + opt_est$value)}
    logLik <- sapply(par_CI_seq, CI)

    data_CI <- data.frame(par_CI = par_CI_seq, logLik = logLik)
#     data_CI$logLik[data_CI$logLik < (max(data_CI$logLik, na.rm = TRUE) - 3.84)] <- NA
    data_CI$logLik[data_CI$logLik < (max(data_CI$logLik, na.rm = TRUE) - 2.71)] <- NA ## 90% confidence interval
    data_CI <- na.omit(data_CI)

    ci_pro <- data.frame(lower = min(data_CI$par_CI, na.rm = TRUE),
                         upper = max(data_CI$par_CI, na.rm = TRUE))

    est_CI[[TT]] <- ci_pro
}

Rt_vals <- unlist(est_mod[!sapply(est_mod, is.null)])
CI_vals <- do.call(rbind, est_CI[!sapply(est_CI, is.null)])
t_vals <- which(!sapply(est_mod, is.null))

result_mod <- data.frame(t = t_vals, Rt = Rt_vals, lower = CI_vals[, "lower"], upper = CI_vals[, "upper"])
merge(data, result_mod, by = 't', all.x = TRUE) -> result_mod

## smoothing the estimated Rt values
result_mod %>%
mutate(Rt_ma = zoo::rollmean(Rt, k = 14, fill = NA, align = "right"),
       lower_ma = rollmean(lower, k = 14, fill = NA, align = "right"),
       upper_ma = rollmean(upper, k = 14, fill = NA, align = "right")) %>% 
rename(date = infection_date) -> result_MA_mod

## plotting the estimated Rt
options(repr.plot.width = 10, repr.plot.height = 7)

scaling_parameter = max(result_MA_mod$case_total)/max(result_MA_mod$upper[!is.na(result_MA_mod$upper)])

adj = 2.9

result_MA_mod %>% filter(date >= as.Date("2018-01-01") & date <= as.Date("2018-12-31")) %>%
ggplot() + 
geom_bar(aes(x = date, y = case_total), stat = 'identity', fill = "#FAAB18") +  
geom_line(aes(x = date, y = Rt_ma*scaling_parameter*adj), color = "steelblue", size = 1) +
geom_ribbon(aes(x = date, ymax = upper_ma*scaling_parameter*adj, ymin = lower_ma*scaling_parameter*adj), 
            fill = "steelblue", alpha = 0.4) + 
scale_x_date(date_breaks = "1 months", labels = scales::date_format("%b"), expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0), breaks = c(2, 4, 6, 8, 10),
                   sec.axis = sec_axis(~./(scaling_parameter*adj), breaks = c(0, 1, 2, 3, 4, 5), 
                                       name = "Effective reproduction number \n")) + 
coord_cartesian(ylim = c(0, 10)) +
geom_hline(yintercept = scaling_parameter*adj, linetype = "dashed", color = "steelblue", size = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 17, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
labs(x = "\n Date of infection", y = "Incidence \n")

## what if we smooth the backprojected epicurve (using a single trajectory)
backproj[[1]] %>%
mutate(case_total_MA = rollmean(case_total, k = 14, fill = 0, align = "right")) %>%
dplyr::select(infection_date, case_total_MA, t) %>% rename(case_total = case_total_MA) -> data

data %>% filter(case_total > 0) %>% arrange(infection_date) %>% slice(1) %>% pull(t) -> T_min

options(repr.plot.width = 10,repr.plot.height = 5, warn = -1)

backproj[[1]] %>% 
ggplot(aes(x = infection_date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of infection", y = "Incidence \n")
            
data %>%
ggplot(aes(x = infection_date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of infection", y = "Incidence \n")

## Rt estimation with the smoothened epicurve
gen_vec <- sapply(1:(T_max + 1), generation_fal)

Start.T <- T_min-10


est_mod <-　list(); est_CI <- list()

for (TT in Start.T:T_max) {
    data_T <- data[data$t <= TT, ]
    t <- TT

    if ((t + 1) > nrow(data_T)) next
    Css <- rep(0, t)
    for (tau in 1:(t - 1)) {
        Css[tau] <- data_T$case_total[t - tau + 1] * gen_vec[tau]
    }

    llk <- function(R) {
        Cs <- sum(Css) * R
        if (Cs <= 0) Cs <- 1e-5
        obs <- data_T$case_total[t + 1]
        -(-Cs + obs * log(Cs) - lgamma(obs + 1))
    }

    opt_est <- optim(0.7, llk, method = "L-BFGS-B", lower = 0, control = list(maxit = 10000))
    est_mod[[TT]] <- opt_est$par
    
    par_CI_seq <- seq(0, 10, by = 0.005)
    CI <- function(par_CI) {2 * (-llk(par_CI) + opt_est$value)}
    logLik <- sapply(par_CI_seq, CI)

    data_CI <- data.frame(par_CI = par_CI_seq, logLik = logLik)
#     data_CI$logLik[data_CI$logLik < (max(data_CI$logLik, na.rm = TRUE) - 3.84)] <- NA
    data_CI$logLik[data_CI$logLik < (max(data_CI$logLik, na.rm = TRUE) - 2.71)] <- NA ## 90% confidence interval
    data_CI <- na.omit(data_CI)

    ci_pro <- data.frame(lower = min(data_CI$par_CI, na.rm = TRUE),
                         upper = max(data_CI$par_CI, na.rm = TRUE))

    est_CI[[TT]] <- ci_pro
}

Rt_vals <- unlist(est_mod[!sapply(est_mod, is.null)])
CI_vals <- do.call(rbind, est_CI[!sapply(est_CI, is.null)])
t_vals <- which(!sapply(est_mod, is.null))

result_mod <- data.frame(t = t_vals, Rt = Rt_vals, lower = CI_vals[, "lower"], upper = CI_vals[, "upper"])
merge(data, result_mod, by = 't', all.x = TRUE) %>% rename(date = infection_date) -> result_mod

## smoothing Rt again (do we really need this? not sure)
result_mod %>%
mutate(Rt_ma = zoo::rollmean(Rt, k = 14, fill = NA, align = "right"),
       lower_ma = rollmean(lower, k = 14, fill = NA, align = "right"),
       upper_ma = rollmean(upper, k = 14, fill = NA, align = "right")) -> result_MA_mod

## plotting the new version of Rt again
options(repr.plot.width = 10, repr.plot.height = 7)

scaling_parameter = max(result_mod$case_total)/max(result_mod$upper[!is.na(result_mod$upper)])

adj = 2.9

result_mod %>% filter(date >= as.Date("2018-01-01") & date <= as.Date("2018-12-31")) %>%
ggplot() + 
geom_bar(aes(x = date, y = case_total), stat = 'identity', fill = "#FAAB18") +  
geom_line(aes(x = date, y = Rt*scaling_parameter*adj), color = "steelblue", size = 1) +
geom_ribbon(aes(x = date, ymax = upper*scaling_parameter*adj, ymin = lower*scaling_parameter*adj), 
            fill = "steelblue", alpha = 0.4) + 
scale_x_date(date_breaks = "1 months", labels = scales::date_format("%b"), expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0), breaks = c(2, 4, 6, 8, 10),
                   sec.axis = sec_axis(~./(scaling_parameter*adj), breaks = c(0, 1, 2, 3, 4, 5), 
                                       name = "Effective reproduction number \n")) + 
coord_cartesian(ylim = c(0, 10)) +
geom_hline(yintercept = scaling_parameter*adj, linetype = "dashed", color = "steelblue", size = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 17, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
labs(x = "\n Date of infection", y = "Incidence \n")

## what if we run multiple iterations and calculate the mean for backprojections
niter <- 100
min_date <- as.Date("2016-01-01")

backproj_list <- vector("list", niter)

for(i in 1: niter){
    mal_example %>% rowwise() %>%
    do({
        n_cases <- .$case_total
        if (n_cases > 0) {
            sampled_days <- sample(0:501, size = n_cases, replace = TRUE, prob = incubation_pdf)
            data.frame(infection_date = .$Date - days(sampled_days))
        } else {
            data.frame(infection_date = as.Date(character(0)))
        }
    }) %>% ungroup() %>% 
    count(infection_date, name = "case_total") %>%
    arrange(infection_date) -> temp
    data.frame(infection_date = as.Date(min_date:max(mal_example$Date))) -> cal_temp
    merge(cal_temp, temp, by = "infection_date", all.x = TRUE) %>% 
    mutate(t = 0:(nrow(cal_temp)-1)) -> backproj_list[[i]]
    backproj_list[[i]]$case_total[is.na(backproj_list[[i]]$case_total)] <- 0
}


backproj_mean <- Reduce("+", lapply(backproj_list, function(x) x$case_total)) / niter

infection_dates <- backproj_list[[1]]$infection_date
backproj <- data.frame(infection_date = infection_dates, case_total = backproj_mean) %>% mutate(t = 0:(n() - 1))  
                                    
options(repr.plot.width = 10,repr.plot.height = 5, warn = -1)
range = c(as.Date("2016-09-01"), as.Date("2020-01-01"))

mal_example %>% 
ggplot(aes(x = Date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of reporting", y = "Incidence \n")

backproj %>% 
ggplot(aes(x = infection_date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of infection", y = "Incidence \n")

## estimating Rt with the multiple iterations
gen_vec <- sapply(1:(T_max + 1), generation_fal)

Start.T <- T_min-10

est_mod <-　list(); est_CI <- list()

for (TT in Start.T:T_max) {
    data_T <- data[data$t <= TT, ]
    t <- TT

    if ((t + 1) > nrow(data_T)) next
    Css <- rep(0, t)
    for (tau in 1:(t - 1)) {
        Css[tau] <- data_T$case_total[t - tau + 1] * gen_vec[tau]
    }

    llk <- function(R) {
        Cs <- sum(Css) * R
        if (Cs <= 0) Cs <- 1e-5
        obs <- data_T$case_total[t + 1]
        -(-Cs + obs * log(Cs) - lgamma(obs + 1))
    }

    opt_est <- optim(0.7, llk, method = "L-BFGS-B", lower = 0, control = list(maxit = 10000))
    est_mod[[TT]] <- opt_est$par
    
    par_CI_seq <- seq(0, 10, by = 0.01)
    CI <- function(par_CI) {2 * (-llk(par_CI) + opt_est$value)}
    logLik <- sapply(par_CI_seq, CI)

    data_CI <- data.frame(par_CI = par_CI_seq, logLik = logLik)
    data_CI$logLik[data_CI$logLik < (max(data_CI$logLik, na.rm = TRUE) - 3.84)] <- NA
    data_CI <- na.omit(data_CI)

    ci_pro <- data.frame(lower = min(data_CI$par_CI, na.rm = TRUE),
                         upper = max(data_CI$par_CI, na.rm = TRUE))

    est_CI[[TT]] <- ci_pro
}

Rt_vals <- unlist(est_mod[!sapply(est_mod, is.null)])
CI_vals <- do.call(rbind, est_CI[!sapply(est_CI, is.null)])
t_vals <- which(!sapply(est_mod, is.null))

result_mod <- data.frame(t = t_vals, Rt = Rt_vals, lower = CI_vals[, "lower"], upper = CI_vals[, "upper"])
merge(data, result_mod, by = 't', all.x = TRUE) %>% rename(date = infection_date) -> result_mod

## smoothing Rt
result_mod %>%
mutate(Rt_ma = zoo::rollmean(Rt, k = 14, fill = NA, align = "right"),
       lower_ma = rollmean(lower, k = 14, fill = NA, align = "right"),
       upper_ma = rollmean(upper, k = 14, fill = NA, align = "right")) -> result_MA_mod

## plotting the estimated Rt with multiple iterations
options(repr.plot.width = 10, repr.plot.height = 7)

scaling_parameter = max(result_mod$case_total)/max(result_mod$upper[!is.na(result_mod$upper)])

adj = 2.9

result_mod %>% filter(date >= as.Date("2018-01-01") & date <= as.Date("2018-12-31")) %>%
ggplot() + 
geom_bar(aes(x = date, y = case_total), stat = 'identity', fill = "#FAAB18") +  
geom_line(aes(x = date, y = Rt*scaling_parameter*adj), color = "steelblue", size = 1) +
geom_ribbon(aes(x = date, ymax = upper*scaling_parameter*adj, ymin = lower*scaling_parameter*adj), 
            fill = "steelblue", alpha = 0.4) + 
scale_x_date(date_breaks = "1 months", labels = scales::date_format("%b"), expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0), breaks = c(2, 4, 6, 8, 10),
                   sec.axis = sec_axis(~./(scaling_parameter*adj), breaks = c(0, 1, 2, 3, 4, 5), 
                                       name = "Effective reproduction number \n")) + 
coord_cartesian(ylim = c(0, 10)) +
geom_hline(yintercept = scaling_parameter*adj, linetype = "dashed", color = "steelblue", size = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 17, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
labs(x = "\n Date of infection", y = "Incidence \n")

## what if we backproject cases based on the distribution
cal_temp <- data.frame(infection_date = seq(min_date, max(mal_example$Date), by = "day")) %>%
mutate(case_total = 0)

mal_example %>% rowwise() %>%
do({
    n_cases <- .$case_total
    report_date <- .$Date
    
    if (!is.na(n_cases) && n_cases > 0) {
        delays <- 0:(length(incubation_pdf) - 1)
        inf_dates <- report_date - delays
        probs <- incubation_pdf[delays + 1]

        tibble(infection_date = inf_dates,
               case_contrib = n_cases * probs)} else {
        tibble(infection_date = as.Date(character()), case_contrib = numeric(0))
    }}) %>% ungroup() %>%
group_by(infection_date) %>%
summarise(case_total = sum(case_contrib, na.rm = TRUE), .groups = "drop") %>%
right_join(cal_temp, by = "infection_date") %>%
mutate(case_total = replace_na(case_total.x, 0), t = 0:(n() - 1)) %>%
select(infection_date, case_total, t) -> backproj

options(repr.plot.width = 10,repr.plot.height = 5, warn = -1)
range = c(as.Date("2016-09-01"), as.Date("2020-01-01"))

mal_example %>% 
ggplot(aes(x = Date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of reporting", y = "Incidence \n")

backproj %>% 
ggplot(aes(x = infection_date, y = case_total)) +
geom_bar(stat = "identity", fill = "gray35", width = 1) +
theme(text = element_text(size = 17, family = "sans", color = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
scale_x_date(date_breaks = "3 months", expand = c(0, 0),  limits = range, 
             labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
scale_y_continuous(limit = c(0, 14), expand = c(0, 0), breaks = c(2, 4, 6, 8, 10, 12, 14)) +
labs(x = "\n Date of infection", y = "Incidence \n")
