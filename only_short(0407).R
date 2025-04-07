
##########################################################################################################
days <- 0:100  

prob_density <- dnorm(days, mean = 50, sd = 15)  
case_counts <- prob_density / max(prob_density) * 50  


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

df_samples <- rbind(df_samples, temp_lastdays)



# inc_fit=list(p= 83/(142+83), shape=1.2, scale=22.2, mean=337.4, sd=38.5)
# incubation<-function(r){(1 - inc_fit$p) * dgamma(r, shape=inc_fit$shape, scale=inc_fit$scale) + 
#     inc_fit$p * dnorm(r, mean=inc_fit$mean, sd=inc_fit$sd)}

incubation<-function(r){dgamma(r, shape=1.2, scale=22.2)}
plot(incubation(1:100))
sum(incubation(1:100))

g_r <- function(delay) {
  ifelse(delay >= 4 & delay <= 7, 1 / (7 - 4 + 1), 0)
}

## Generation time distribution -> Endo and Nishiura 2015
# generation <- function(r) {
#   ifelse(r < 30, 0, incubation(r - 30))
# }
gi_fit=list(shape=2.305, scale=5.452)
generation<-function(t){pweibull(t, shape=gi_fit$shape, scale=gi_fit$scale)-
    pweibull(t-1, shape=gi_fit$shape, scale=gi_fit$scale)}


K=nrow(df_samples) 
#normal_pmf = pnorm(1:K, norm_fit$mean, norm_fit$sd)-pnorm(1:K-1, norm_fit$mean, norm_fit$sd)
nor_pmf=sapply(1:K, incubation)
sum(nor_pmf)
plot(nor_pmf)

sts=new("sts", epoch=df_samples$time_onset, observed=df_samples$case_count )
bpnp.control=list(k=2, eps=rep(1e-4,2), iter.max=rep(1000,2),
                  Tmark=nrow(sts), B=-1, alpha=0.01, verbose=FALSE, lambda0=NULL, eq3a.method=c("R","C"))
sts_bp=backprojNP(sts, incu.pmf=nor_pmf,
                  control=modifyList(bpnp.control, list(eq3a.method="C")))
df_samples$backproj=upperbound(sts_bp)

df_samples%>%mutate(
  total=backproj/sum(backproj)*sum(case_count))%>%
  filter(time_onset<=(101-1))->dt.backproj


plot_data <- dt.backproj %>%
  dplyr::select(time_onset, case_count, total) %>%
  rename(`Reported (onset)` = case_count,
         `Infection (backproj)` = total) %>%
  pivot_longer(cols = c(`Reported (onset)`, `Infection (backproj)`),
               names_to = "Type", values_to = "Count")

ggplot(plot_data, aes(x = time_onset, y = Count, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.6) +
  facet_grid(rows = vars(Type), scales = "fixed") +  
  scale_fill_manual(values = c("Reported (onset)" = "lightblue",
                               "Infection (backproj)" = "lightcoral")) +
  labs(x = "Days", y = "Count", title = "Reported vs BackProjection") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none", 
    panel.spacing = unit(1, "lines")
  )
##################################################################################################################
sum(df_samples$case_count)

today=101-min(dt.backproj$time_onset)
Start.T=1
conv<-function(t,tau){sum(convolve(uni(1:1000), rev(incubation(1:1000)), type=c("open"))[1:(today-t+tau)])}
precal<-matrix(0, nrow=max((dt.backproj$time_onset )-Start.T+2+1),ncol=(max(dt.backproj$time_onset )-1))
for (m in (Start.T-2):max(dt.backproj$time_onset )){
  for (n in 1:m-1){
    print(paste("m:", m, "n:", n)) 
    precal[m-Start.T+2+1, n]<-conv(m,n)}}

delay_precalculation <- function(t){sum(convolve(uni(1:(today-t)),
                                                 rev(incubation(1:(today-t))), type = c("open")))}
dt.backproj %<>% rowwise() %>%
  mutate(
    case_delay = total / delay_precalculation(time_onset),
    case_delay = ifelse(is.nan(case_delay) | is.infinite(case_delay), 0, case_delay)
  )


# dt.backproj$test <-dt.backproj$case_count*delay_precalculation(dt.backproj$time_onset)
# 
# 
# ggplot(data = dt.backproj) +
#   geom_bar(aes(x = time_onset, y = case_count, fill = "Case Count"), 
#            stat = 'identity', alpha = 0.6, position = "identity") +  # 투명도 추가
#   geom_bar(aes(x = time_onset, y = case_delay, fill = "Case delay"), 
#            stat = 'identity', alpha = 0.6, position = "identity") +  # 투명도 추가
#   labs(x = "Days", y = "Count", title = "Case Count and Case delay") +
#   scale_fill_manual(values = c("Case Count" = "lightblue", "Case delay" = "lightcoral")) +  # 색상 설정
#   theme_minimal() +
#   theme(legend.title = element_blank())  # 범례 제목 없애기
# sum(dt.backproj$case_delay)
# 
# ############################################################################################
# est.t<-list(); est.CI<-list();
# for (TT in Start.T:max(dt.backproj$time_onset)){
#   dt.backproj %>% filter(time_onset <= TT) -> dt.backproj.T
# 
#   # R_t = I_t/sigma_tau=1~t-1(I_(t-tau)*g(tau)/F
#   llk <- function(param){t=TT; R = param
#   llk <- 0; Cs <- 0 ;　Css <- rep(0, t)
#   for (tau in 1: t - 1){
#     print(paste("TT:", TT, "tau:", tau))
#     Css[tau] = (dt.backproj.T$case_delay[t-tau+1]) * generation(tau)/precal[t - Start.T + 3, tau]
#   }
#   Cs = sum(Css) * R
#   Cs[Cs <= 0] <- 1e-5
#   return(-(-Cs+dt.backproj.T$case_delay[t+1] * log(Cs) - lgamma(dt.backproj.T$case_delay[t + 1] + 1)))
#   }
#   param0 = c(0.7)
#   opt_est <- optim(param0, fn=llk, method=c("L-BFGS-B"), lower=c(0), control=list(maxit = 10000))
#   opt_est$par -> est.t[[TT]]
# 
#   ## 프로파일 우도를 활용한 95% 신뢰구간 추정
#   ci_pro <- matrix(NA, ncol = 2, nrow = 1)
#   CI <- function(par_CI){return(2 * (-llk(par_CI)+opt_est$value))} #llk()
#   par_CI <- seq(0, 20, by = 0.01)
#   logLik <- sapply(par_CI, FUN = CI)
#   as.data.frame(par_CI) -> par_CI; as.data.frame
#   (logLik) -> logLik
#   cbind(par_CI, logLik) -> data_CI
#   data_CI$logLik[data_CI$logLik < (max(data_CI$logLik) - 3.84)] <- NA
#   data_CI %<>% na.omit()
#   min(data_CI$par_CI) -> ci_pro[1,1]; max(data_CI$par_CI) -> ci_pro[1,2]
#   as.data.frame(ci_pro) -> ci_pro
#   colnames(ci_pro) <- c("lower","upper")
#   ci_pro -> est.CI[[TT]]
# }
# 
# matrix(unlist(est.t),ncol=1,byrow=T) -> est.t
# matrix(unlist(est.CI),ncol=2,byrow=T) -> est.CI
# cbind(est.t, est.CI) -> est
# #head(dt.backproj)
# est %>% as.data.frame() %>% mutate(time_onset = Start.T:max(dt.backproj$time_onset)) -> result
# merge(dt.backproj, result, by='time_onset', all.x=TRUE) -> result
# colnames(result) <- c("time_onset","case_count","backproj","total","case_delay","Rt","lower", "upper")
# #result %>% filter(time_onset <= as.Date(as.Date("2023-01-01")-(17))) -> result_withcut
# result %>% filter(time_onset>=5 &time_onset <= 60)-> result_withcut
# options(repr.plot.width = 7, repr.plot.height = 5)
# theme_set(theme_bw())
# scaling_parameter = max(result_withcut$case_delay) / max(result_withcut$upper[!is.na(result_withcut$upper)])
# range = c(2, 50)
# options(warn=-1)
# result_withcut %>%
#   ggplot() +
#   #geom_line(data=result[!is.na(result$Rt),],aes(x=Date,y=case_total), color="red", size=1) +
#   geom_bar(aes(x=time_onset, y=case_delay), stat='identity', fill="gray35", width=0.7) +
#   geom_line(data=result_withcut[!is.na(result_withcut$Rt),],aes(x=time_onset ,y=Rt * scaling_parameter), color="#1380A1", size=1) +
#   geom_ribbon(data=result_withcut, aes(ymax=result_withcut$upper * scaling_parameter,
#                                ymin=result_withcut$lower * scaling_parameter, x=time_onset ),
#               fill="#1380A1", alpha = 0.4) + labs(x="\n 추정감염일", y=" 감염자 수\n") +
#   theme (text = element_text ( size = 12 , family="sans",color = "black"),
#          axis.text = element_text (size = 12 , family="sans",color = "black"),
#          panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
#          legend.position="none") +
#   # scale_x_date(date_labels="%m/%d", date_breaks="10 day",
#   #              limits=range, expand=c(0, 0)) +
#   scale_y_continuous(limit=c(0,200), expand = c(0, 0),
#                      sec.axis=sec_axis(~./(scaling_parameter), breaks=c(0,2,4), name =" 실질감염재생산수\n")) +
#   geom_hline(yintercept=1 * scaling_parameter,
#              linetype="dashed", color = "#1380A1", size =0.7)