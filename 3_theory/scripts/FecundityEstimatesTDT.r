#+++++++++++++++++++++++++++++
#Estimate fecundity incorporating TDT data
#++++++++++++++++++++++++++++

# ============================================================
# Increasing-Variance Log-Logistic TDT Model
# Direct adaptation of Bullard et al. (2026, PLOS Biology) code
# https://github.com/gbullard87/Beyond_Mean_TDT_2025
# ============================================================
#bayesTLS https://daniel1noble.github.io/bayesTLS/#sec-aphids

library(tidyverse)
library(flexsurv)
library(survival)

#read in data and output
temps.all <-read.csv("3_theory/data/buckley_data/TempTimeSeries.csv")
PerfDat <- read.csv("3_theory/data/buckley_data/PerformanceData.csv")

#TDT Data from Li et al. 2023
aphid<-read.table("3_theory/data/buckley_data/surv_Lietal.txt", header = T)
aphid<-aphid%>%mutate(dur_h=dur/60)%>%within(.,spec<-factor(spec,levels = c("MD","SA","RP")),age<-factor(age,levels=c("2","6","12")))

#rename treatment in expt 3
t3<- c("22_0","22_5","22_9","22_13")
t3.lab<- c("00","05","09","13")
temps.all$treatment[temps.all$expt==3]<- t3.lab[match(temps.all$treatment[temps.all$expt==3],t3)]
PerfDat$treatment[PerfDat$expt==3]<- t3.lab[match(PerfDat$treatment[PerfDat$expt==3],t3)]

#====================
#Constant rate TPCs
#English grain aphid, Sitobion avenae
#Zhao et al 2013
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12196

#https://journals.biologists.com/jeb/article/218/14/2289/14375/Daily-temperature-extremes-play-an-important-role

#development rate
#use Ma et al 2015
dr= function(T, Tmax=32.947, a=0.137, b=4.514, c=7.267){ 
  d=exp(a*T)-exp(b-(Tmax-T)/c)
  d[d<0]<- 0
  return(d)
}

#survival
sur= function(T, a=345.51, b=0.35) {
  sur=1-T/(a-b*T^2)
  sur[T>30]<- 0
  return(sur)
}

long= function(T, a=32.73, b= -0.91) 
{
  lon= a + b*T
  lon[lon<0]<- 0
  return(lon)
}

#fedundity TPC (nymphs/adult)
fec= function(T, a= -69.1, b=12.49, c= -0.34){
  fec=a +b*T +c*T^2
  fec[fec<0]<- 0
  return(fec)
}

#find fecundity Topt
ts=seq(0,40,0.1)
topt<- ts[which.max(fec(ts) )]
fecmax= fec(topt)

# ------------------------------------------------------------
# 1. MODEL FUNCTION -- increasing_var_mod_functions.Rmd
# ------------------------------------------------------------
llogis.mod <- function(TDT.dat, temp.dat, shape = "linear") {
  # step 1: fit log-logistic model separately at each temperature
  shape.vec <- c()
  shape.vec.var <- c()
  scale.vec <- c()
  scale.vec.var <- c()
  for (j in unique(TDT.dat$Tref)) {
    mod.j <- flexsurvreg(Surv(time, event) ~ 1,
                          data = TDT.dat[(TDT.dat$Tref == j), ],
                          dist = "llogis")
    shape.vec.var <- c(shape.vec.var, mod.j$cov[[1, 1]])
    scale.vec.var <- c(scale.vec.var, exp(mod.j$cov[[2, 2]]))
    shape.vec <- c(shape.vec, mod.j$coefficients[[1]])
    scale.vec <- c(scale.vec, exp(mod.j$coefficients[[2]]))
  }

  # step 2: regress shape and log10(scale) against temperature
  shape.lm <- lm(shape.vec ~ unique(TDT.dat$Tref))
  shape.slope <- shape.lm$coefficients[[2]]
  shape.int <- shape.lm$coefficients[[1]]
  cons.shape.lm <- lm(shape.vec ~ 1)
  cons.shape.int <- cons.shape.lm$coefficients[[1]]
  scale.lm <- lm(log10(scale.vec) ~ unique(TDT.dat$Tref))
  scale.slope <- scale.lm$coefficients[[2]]
  scale.int <- scale.lm$coefficients[[1]]

  shape.pars <- c(shape.slope, shape.int, cons.shape.int)
  logscale.pars <- c(scale.slope, scale.int)

  if (shape == "linear") {
    cum.surv <- rep(NA, nrow(temp.dat))
    cum.surv[1] <- 1
    fail.rates <- rep(0, nrow(temp.dat))
    mean.list <- c(); shape.list <- c(); log.list <- c()

    for (i in 2:nrow(temp.dat)) {
      T_i <- temp.dat$Temperature[i]
      T_prev <- temp.dat$Temperature[i - 1]
      t_i <- temp.dat$Time_min[i]
      t_prev <- temp.dat$Time_min[i - 1]
      t_diff <- t_i - t_prev
      scale_i <- 10^(scale.slope * mean(T_i, T_prev) + scale.int)
      shape_i <- shape.slope * mean(T_i, T_prev) + shape.int
      mean.list <- c(mean.list, ((pi * scale_i) / shape_i) / (sin(pi / shape_i)))
      log.list <- c(log.list, cum.surv[i - 1])
      prob <- 1 - cum.surv[i - 1]
      t_star_i <- scale_i * (prob / (1 - prob))^(1 / shape_i)
      fail.rates[i] <- (1 - pllogis(t_star_i, shape_i, scale_i)) -
                        (1 - pllogis(t_star_i + t_diff, shape_i, scale_i))
      cum.surv[i] <- cum.surv[i - 1] - fail.rates[i]
    }
    med.tdt <- temp.dat$Time_min[which.max(cum.surv < 0.5)]
    cum.surv.post <- na.omit(cum.surv)
    times.post <- temp.dat$Time_min[1:length(cum.surv.post)]
    expected.tdt <- mean(
      -sum((cum.surv.post[2:length(cum.surv.post)] - cum.surv.post[1:(length(cum.surv.post) - 1)]) * times.post[2:length(cum.surv.post)]),
      -sum((cum.surv.post[2:length(cum.surv.post)] - cum.surv.post[1:(length(cum.surv.post) - 1)]) * times.post[1:(length(cum.surv.post) - 1)])
    )
  }

  if (shape == "constant") {
    cum.surv <- rep(NA, nrow(temp.dat))
    cum.surv[1] <- 1
    fail.rates <- rep(0, nrow(temp.dat))
    mean.list <- c(); shape.list <- c(); log.list <- c()

    for (i in 2:nrow(temp.dat)) {
      T_i <- temp.dat$Temperature[i]
      T_prev <- temp.dat$Temperature[i - 1]
      t_i <- temp.dat$Time_min[i]
      t_prev <- temp.dat$Time_min[i - 1]
      t_diff <- t_i - t_prev
      scale_i <- 10^(scale.slope * mean(T_i, T_prev) + scale.int)
      shape_i <- cons.shape.int
      mean.list <- c(mean.list, ((pi * scale_i) / shape_i) / (sin(pi / shape_i)))
      log.list <- c(log.list, cum.surv[i - 1])
      prob <- 1 - cum.surv[i - 1]
      t_star_i <- scale_i * (prob / (1 - prob))^(1 / shape_i)
      fail.rates[i] <- (1 - pllogis(t_star_i, shape_i, scale_i)) -
                        (1 - pllogis(t_star_i + t_diff, shape_i, scale_i))
      cum.surv[i] <- cum.surv[i - 1] - fail.rates[i]
    }
    med.tdt <- temp.dat$Time_min[which.max(cum.surv < 0.5)]
    cum.surv.post <- na.omit(cum.surv)
    times.post <- temp.dat$Time_min[1:length(cum.surv.post)]
    expected.tdt <- mean(
      -sum((cum.surv.post[2:length(cum.surv.post)] - cum.surv.post[1:(length(cum.surv.post) - 1)]) * times.post[2:length(cum.surv.post)]),
      -sum((cum.surv.post[2:length(cum.surv.post)] - cum.surv.post[1:(length(cum.surv.post) - 1)]) * times.post[1:(length(cum.surv.post) - 1)])
    )
  }

  output <- list(med.tdt, expected.tdt, cum.surv.post, shape.pars, logscale.pars)
  names(output) <- c("med.tdt", "expected.tdt", "cum.surv", "shape.pars", "logscale.pars")
  return(output)
}

# ------------------------------------------------------------
# 2. ADAPT YOUR `aphid` DATA to their required TDT.dat format:
#    columns Tref (temperature), time (failure/censor time),
#    event (1 = death observed, 0 = censored), one row per
#    INDIVIDUAL (not counts). This mirrors how the repo expands
#    Jorgensen et al.'s count data into individual-level rows
#    before calling flexsurvreg(Surv(time, event) ~ 1, ...).
# ------------------------------------------------------------
# aphid: spec, age, temp, dur, rep, surv, total, survival_rate, dur_h
# We expand each temp x rep trial into one row per individual:
#   - each newly dead individual gets event = 1 at the dur when it was
#     first recorded dead (matches how flexsurvreg(Surv(time,event)~1,
#     dist="llogis") is used elsewhere in the repo -- exact/right-
#     censored time coding, not interval-censored)
#   - remaining survivors at the final dur get event = 0 (censored)

expand_to_individual <- function(df) {
  df <- df %>% arrange(dur)
  n_prev_dead <- 0
  out <- list()
  for (i in seq_len(nrow(df))) {
    n_dead_now <- df$total[i] - df$surv[i]
    n_new_deaths <- n_dead_now - n_prev_dead
    if (n_new_deaths > 0) {
      out[[length(out) + 1]] <- tibble(
        time = rep(df$dur[i], n_new_deaths),
        event = rep(1, n_new_deaths)
      )
    }
    n_prev_dead <- n_dead_now
  }
  n_survivors <- df$surv[nrow(df)]
  if (n_survivors > 0) {
    out[[length(out) + 1]] <- tibble(
      time = rep(df$dur[nrow(df)], n_survivors),
      event = rep(0, n_survivors)
    )
  }
  bind_rows(out)
}

# Build TDT.dat for a single species (adjust `spec_name` as needed;
# subset first if aphid contains multiple species/ages you want to
# model separately)
spec_name <- "SA"   # <-- set to the species/age combination you want
aphid_sub <- aphid %>% filter(spec == spec_name)

TDT.dat <- aphid_sub %>%
  group_by(temp, rep) %>%
  group_modify(~ expand_to_individual(.x)) %>%
  ungroup() %>%
  rename(Tref = temp) %>%
  dplyr::select(Tref, time, event)

#subset to warm temperatures?
TDT.dat <- TDT.dat[TDT.dat$Tref>30,]

# ------------------------------------------------------------
# 3. ADAPT temperature data to
#    their required temp.dat format: columns Temperature, Time_min
# ------------------------------------------------------------
# temp.dat <- Wuhan_summer_temps %>%
#   arrange(datetime) %>%
#   mutate(Time_min = as.numeric(difftime(datetime, min(datetime), units = "mins"))) %>%
#   rename(Temperature = temp) %>%
#   dplyr::select(Time_min, Temperature)

#start at time 0
temps.all$time <- temps.all$time-1

#convert time in hours to minutes
temps.all$Time_min <- temps.all$time *60
temps.all$Temperature <- temps.all$temp
  
# ------------------------------------------------------------
# 4. RUN the model -- linear (increasing-variance) shape version
#    and constant-shape version for comparison, exactly as the
#    repo's workflow does (e.g., ll.spec vs ll.cons.* objects)
# ------------------------------------------------------------
# ll.linear  <- llogis.mod(TDT.dat, temp.dat, shape = "linear")
# ll.const   <- llogis.mod(TDT.dat, temp.dat, shape = "constant")
# 
# ll.linear$med.tdt        # predicted median time-to-death (minutes) under fluctuating regime
# ll.linear$expected.tdt   # predicted expected time-to-death
# ll.linear$shape.pars     # c(shape.slope, shape.int, cons.shape.int)
# ll.linear$logscale.pars  # c(scale.slope, scale.int) for log10(scale) ~ Tref

temps.all$cumsurv<- NA
temps.all$dr<- NA
temps.all$surv<- NA
temps.all$fec<- NA

#estimate survival by experiment and treatment
for(expt.k in 1:7){
  tdat<- temps.all[temps.all$expt==expt.k,]
  treats<- unique(tdat$treatment)
  
  for(treat.k in 1:length(treats) ){
    ttdat<- tdat[tdat$treatment==treats[treat.k],]
    
    ll.const   <- llogis.mod(TDT.dat, ttdat, shape = "constant")
    
    #add cum surv to dataframe
    inds<- which(temps.all$expt== expt.k & temps.all$treatment== treats[treat.k] )
    temps.all$cumsurv[inds]<- ll.const$cum.surv
    
    #estimate fecundity and development
    temps.all$dr[inds]<- dr(temps.all$Temperature[inds])
    temps.all$surv[inds]<- sur(temps.all$Temperature[inds])
    temps.all$fec[inds]<- fec(temps.all$Temperature[inds])
  
  }} #end loop

#PLOT
#to long format
t.long <- temps.all %>%
  dplyr::select(-X, -Time_min, -temp) %>%                       # drop row-index column, if not needed
  pivot_longer(
    cols = -c(treatment, expt, time),
    names_to = "variable",
    values_to = "value"
  )

#just survival
ggplot(temps.all, aes(x = time, y = cumsurv, color = treatment)) +
  geom_line() +
  theme_bw() + facet_wrap(.~expt, scales = "free_y")

#plot all estimates
ggplot(t.long[t.long$expt==1,], aes(x = time, y = value, color = treatment)) +
  geom_line() +
  facet_wrap(.~variable, scales = "free_y")+
  theme_bw() 

#----------------
#compare to empirical data

#assume max fecundity above topt
temps.all$fecmax<- temps.all$fec
temps.all$fecmax[temps.all$temp>topt]<- fecmax

#estimate fitness as cumulative survival 
temps.all$fit<- temps.all$fec*temps.all$cumsurv
#fit with max fecundity
temps.all$fitmax<- temps.all$fecmax*temps.all$cumsurv

#cumulative fitness
t.cumfit <- temps.all %>%
  arrange(expt, treatment, time) %>%
  group_by(expt, treatment) %>%
  mutate(cum_fit = cumsum(fit),
         cum_fec = cumsum(fec),
         cum_fitmax = cumsum(fitmax)
         ) %>%
  ungroup()

#compare to observed fitness at time 200
fdat<- t.cumfit[t.cumfit$time==200,c("treatment","expt","cumsurv","cum_fit","cum_fec", "cum_fitmax")]

#subset to lab for comparison
PerfDat.sub <- PerfDat[-which(PerfDat$population=="field"),]

#mean across replicates
PerfDat.avg <- PerfDat.sub %>%
  group_by(expt, treatment, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = n(),
    .groups = "drop"
  )

#match observations
fdat$ex_tr<- paste(fdat$expt, fdat$treatment, sep="_")
PerfDat.avg$ex_tr<- paste(PerfDat.avg$expt, PerfDat.avg$treatment, sep="_")

#just fecundity
PerfDat.avg <- PerfDat.avg[-which(PerfDat.avg$metric=="fecundity"),]

match1<- match(fdat$ex_tr, PerfDat.avg$ex_tr) 
fdat$ObsFec <- PerfDat.avg$mean_value[match1]

#plot observed against estimates
ggplot(fdat, aes(x = ObsFec, y = cum_fit, color = treatment, group=treatment)) +
  geom_point() +
  geom_smooth(method="lm")+
  facet_wrap(.~expt, scales = "free")+
  theme_bw() 

ggplot(fdat, aes(x = ObsFec, y = cum_fec, color = treatment, group=treatment)) +
  geom_point() +
  geom_smooth(method="lm")+
  facet_wrap(.~expt, scales = "free")+
  theme_bw() 

#just cumulative survival
ggplot(fdat, aes(x = cumsurv, y = ObsFec, color = treatment, group=treatment)) +
  geom_point() +
  geom_smooth(method="lm")+
  facet_wrap(.~expt, scales = "free")+
  theme_bw() 

#assume fecundity only declines due to thermal stress
ggplot(fdat, aes(x = cum_fitmax, y = ObsFec, color = treatment, group=treatment)) +
  geom_point() +
  geom_smooth(method="lm")+
  facet_wrap(.~expt, scales = "free")+
  theme_bw() 

