library(stats4); library(dplyr); library(ggplot2); library(MASS); library(fitdistrplus); 
library(surveillance);library(tidyr); library(lubridate); library(magrittr); 
rm=ls()
#-----------------------------------------------------------------------------------------------#
# 1. random으로 bell shape 만들어(report data)
set.seed(42)  # 재현성을 위해 시드 설정
days <- 0:100  

# 정규 분포 확률 밀도를 이용하여 날짜별 case 수 설정
prob_density <- dnorm(days, mean = 50, sd = 15)  # 평균 250, 표준편차 100
case_counts <- prob_density / max(prob_density) * 100  # 최댓값을 100으로 정규화

# 반올림하여 정수로 변환
case_counts_scaled <- round(case_counts)



df_samples <- data.frame(
  time_onset = days,
  case_count = case_counts_scaled)

last_time_onset <- max(df_samples$time_onset)

# 새로운 10일치 데이터 생성
temp_lastdays <- data.frame(
  time_onset = (last_time_onset + 1):(last_time_onset + 10),
  case_count = 0  # 신규 날짜의 case_count는 0으로 설정
)

# 기존 df_samples와 합치기
df_samples <- rbind(df_samples, temp_lastdays)


# head(df_samples)
# # 히스토그램 그리기 (날짜별 case_count를 고려하여)
# barplot(df_samples$case_count, names.arg = df_samples$time_onset,
#         col = "lightblue", border = "white",
#         main = "Bell Shape Case Distribution",
#         xlab = "Days", ylab = "Case Count")


# 2. backprojection
norm_fit <- list(mean = 5, sd = 1)
normal_dist <- function(t) {
  pnorm(t, norm_fit$mean, norm_fit$sd)-pnorm(t-1, norm_fit$mean, norm_fit$sd)
}

K=nrow(df_samples) 
normal_pmf = pnorm(1:K, norm_fit$mean, norm_fit$sd)-pnorm(1:K-1, norm_fit$mean, norm_fit$sd)
nor_pmf=c(0, normal_dist(1:K))
sum(nor_pmf)
#plot(nor_pmf)

sts=new("sts", epoch=df_samples$time_onset, observed=df_samples$case_count )
bpnp.control=list(k=2, eps=rep(1e-4,2), iter.max=rep(1000,2),
                  Tmark=nrow(sts), B=-1, alpha=0.01, verbose=FALSE, lambda0=NULL, eq3a.method=c("R","C"))
sts_bp=backprojNP(sts, incu.pmf=normal_pmf,
                  control=modifyList(bpnp.control, list(eq3a.method="C")))
df_samples$backproj=upperbound(sts_bp)

df_samples%>%mutate(
  total=backproj/sum(backproj)*sum(case_count))%>%
  filter(time_onset<=(101-1))->dt.backproj

ggplot(data = df_samples) +
  geom_bar(aes(x = time_onset, y = case_count, fill = "Case Count"), 
           stat = 'identity', alpha = 0.6, position = "identity") +  # 투명도 추가
  geom_bar(aes(x = time_onset, y = backproj, fill = "Back Projection"), 
           stat = 'identity', alpha = 0.6, position = "identity") +  # 투명도 추가
  labs(x = "Days", y = "Count", title = "Case Count and Back Projection") +
  scale_fill_manual(values = c("Case Count" = "lightblue", "Back Projection" = "lightcoral")) +  # 색상 설정
  theme_minimal() +
  theme(legend.title = element_blank())  # 범례 제목 없애기

# sum(df_samples$case_count) #3688
# sum(df_samples$backproj) #3688.116

#----------------------------------------------------------------------------------------------------------------------------
# 3. incubation * uniform -> 그려서 확인 (아마 원래 incubation과 비슷하게 나올 것)
norm_fit <- list(mean = 5, sd = 1)
normal_dist <- function(t) {
  pnorm(t, norm_fit$mean, norm_fit$sd)-pnorm(t-1, norm_fit$mean, norm_fit$sd)
}

uni<-function(delay) {
  ifelse(delay >= 4 & delay <= 7, 1 / (7 - 4 + 1), 0)
}

gi_fit=list(shape=2.305, scale=5.452)
generation<-function(t){pweibull(t, shape=gi_fit$shape, scale=gi_fit$scale)-
    pweibull(t-1, shape=gi_fit$shape, scale=gi_fit$scale)}

# generation<-function(r){
#   ifelse(r<30,0, normal_dist(r-30))
# }

# convolution 결과 확인 
# t_range <- 1:10
# normal_vec <- sapply(t_range, normal_dist)
# uni_vec <- sapply(t_range, uni)
# conv_result <- convolve(uni_vec, rev(normal_vec), type = "open")
# plot(conv_result, type = "h", lwd = 2, col = "darkblue",
#      main = "Convolution: Uniform(4~7) * Normal(mean=5, sd=1)",
#      xlab = "Time (days)", ylab = "Probability")
# plot(normal_dist(t_range), type='h')


today=101-min(dt.backproj$time_onset)
Start.T=1
conv<-function(t,tau){sum(convolve(uni(1:1000), rev(normal_dist(1:1000)), type=c("open"))[1:(today-t+tau)])}
precal<-matrix(0, nrow=max((dt.backproj$time_onset )-Start.T+2+1),ncol=(max(dt.backproj$time_onset )-1))
for (m in (Start.T-2):max(dt.backproj$time_onset )){
  for (n in 1:m-1){
    print(paste("m:", m, "n:", n)) 
    precal[m-Start.T+2+1, n]<-conv(m,n)}}

delay_precalculation <- function(t){sum(convolve(uni(1:(today-t)),
                                                 rev(normal_dist(1:(today-t))), type = c("open")))}
dt.backproj %<>% rowwise() %>%
  mutate(
    case_delay = total / delay_precalculation(time_onset),
    case_delay = ifelse(is.nan(case_delay) | is.infinite(case_delay), 0, case_delay)
  )

ggplot(data = dt.backproj) +
  geom_bar(aes(x = time_onset, y = case_count, fill = "Case Count"), 
           stat = 'identity', alpha = 0.6, position = "identity") +  # 투명도 추가
  geom_bar(aes(x = time_onset, y = case_delay, fill = "Case delay"), 
           stat = 'identity', alpha = 0.6, position = "identity") +  # 투명도 추가
  labs(x = "Days", y = "Count", title = "Case Count and Case delay") +
  scale_fill_manual(values = c("Case Count" = "lightblue", "Case delay" = "lightcoral")) +  # 색상 설정
  theme_minimal() +
  theme(legend.title = element_blank())  # 범례 제목 없애기
#---------------------------------------------------------------------------------------------------------
#head(df_samples)
est.t<-list(); est.CI<-list();
for (TT in Start.T:max(dt.backproj$time_onset)){
  dt.backproj %>% filter(time_onset <= TT) -> dt.backproj.T
  
  # R_t = I_t/sigma_tau=1~t-1(I_(t-tau)*g(tau)/F
  llk <- function(param){t=TT; R = param
  llk <- 0; Cs <- 0 ;　Css <- rep(0, t)
  for (tau in 1: t - 1){
    print(paste("TT:", TT, "tau:", tau))
    Css[tau] = (dt.backproj.T$case_delay[t-tau+1]) * generation(tau)/precal[t - Start.T + 3, tau]
  }
  Cs = sum(Css) * R
  Cs[Cs <= 0] <- 1e-5
  return(-(-Cs+dt.backproj.T$case_delay[t+1] * log(Cs) - lgamma(dt.backproj.T$case_delay[t + 1] + 1)))
  }
  param0 = c(0.7)
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
#head(dt.backproj)
est %>% as.data.frame() %>% mutate(time_onset = Start.T:max(dt.backproj$time_onset)) -> result
merge(dt.backproj, result, by='time_onset', all.x=TRUE) -> result
colnames(result) <- c("time_onset","case_count","backproj","total","case_delay","Rt","lower", "upper")
#result %>% filter(time_onset <= as.Date(as.Date("2023-01-01")-(17))) -> result_withcut

options(repr.plot.width = 7, repr.plot.height = 5)
theme_set(theme_bw())
scaling_parameter = max(result$case_delay) / max(result$upper[!is.na(result$upper)])
range = c(2, 100)
options(warn=-1) 
result %>% 
  ggplot() +
  #geom_line(data=result[!is.na(result$Rt),],aes(x=Date,y=case_total), color="red", size=1) +
  geom_bar(aes(x=time_onset, y=case_delay), stat='identity', fill="gray35", width=0.7) + 
  geom_line(data=result[!is.na(result$Rt),],aes(x=time_onset ,y=Rt * scaling_parameter), color="#1380A1", size=1) +
  geom_ribbon(data=result, aes(ymax=result$upper * scaling_parameter, 
                               ymin=result$lower * scaling_parameter, x=time_onset ), 
              fill="#1380A1", alpha = 0.4) + labs(x="\n 추정감염일", y=" 감염자 수\n") +
  theme (text = element_text ( size = 12 , family="sans",color = "black"),
         axis.text = element_text (size = 12 , family="sans",color = "black"),
         panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
         legend.position="none") +
  # scale_x_date(date_labels="%m/%d", date_breaks="10 day",
  #              limits=range, expand=c(0, 0)) +
  scale_y_continuous(limit=c(0,120), expand = c(0, 0),
                     sec.axis=sec_axis(~./(scaling_parameter), breaks=c(0,2,4), name =" 실질감염재생산수\n")) + 
  geom_hline(yintercept=1 * scaling_parameter,
             linetype="dashed", color = "#1380A1", size =0.7)
#----------------------------------------------------------------------------------------------------------------
