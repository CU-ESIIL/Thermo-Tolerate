library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridisLite)
library(patchwork)
library(ggplot2)
library(deSolve)
library(pracma)				# contains sigmoid function
library(TrenchR)
library(zoo)

#+++++++++++++++++++++++++++++
#Set up temperature and fitness data for experiments
#++++++++++++++++++++++++++++
#Analysis for English grain aphid, Sitobion avenae
#Processes environmental and fitness component data for the following experiments:
#Expt 1: vary min
#Expt 2: vary max
#Expt 3: vary variance (expt 3, mild means)
#Expt 4: vary means and variance (expt 3, high means)
#Expt 5: vary duration of heatwave (expt 2)
#Expt 6: vary length of heatwave and timing (expt 5 adult)
#Expt 7: vary length of heatwave and timing (expt 5 nymph)

#=====
#Expt 1. vary min
#Zhao et al. 2014. Night warming on hot days produces novel impacts on development, survival and reproduction in a small arthropod
#Dryad data: http://doi.org/10.5061/dryad.q2070 
#English grain aphid, Sitobion avenae

adat2.dt<- read.csv("3_theory/data/Zhaoetal2014/Zhaoetal2014_devtime.csv")
adat2.p<- read.csv("3_theory/data/Zhaoetal2014/Zhaoetal2014_AdPerf.csv")
adat2.lt<- read.csv("3_theory/data/Zhaoetal2014/Zhaoetal2014_LifeTable.csv")
adat2.sur<- read.csv("3_theory/data/Zhaoetal2014/Zhaoetal2014_SurvNymph.csv")
#could add survival but format unclear
#----
#construct temperatures
#currently sawtooth, make sine wave?

#day: increased from 27 °C to peak at 35 °C (at 09.00–15.00 h) and then decreased to 27 °C (at 15.00–21.00 h)
#temperatures in the six chambers dropped from 27 °C to different minima (13, 16, 19, 21, 23 and 25 °C) at 21.00–03.00 h and then rose to 27 °C (at 03.00–09.00 h). 

tmin= c(13, 16, 19, 21, 23, 25)

t.hrs= c(22:24,1:8)
tmin.fun= function(tmin) c( 27-1*(27-tmin)/6, 27-2*(27-tmin)/6, 27-3*(27-tmin)/6, 27-4*(27-tmin)/6, 27-5*(27-tmin)/6, tmin, 27-5*(27-tmin)/6, 27-4*(27-tmin)/6, 27-3*(27-tmin)/6, 27-2*(27-tmin)/6, 27-1*(27-tmin)/6)

temps<- matrix(NA, nrow=length(tmin), ncol=24)
for(t.k in 1:length(tmin)){
  ts= tmin.fun(tmin[t.k])
  temps[t.k, 1:8]= ts[match(1:8, t.hrs)]
  temps[t.k, 9:21]= c(27, 27+1*8/6, 27+2*8/6, 27+3*8/6, 27+4*8/6, 27+5*8/6, 35, 27+5*8/6, 27+4*8/6, 27+3*8/6, 27+2*8/6, 27+1*8/6, 27)
  temps[t.k, 22:24]= ts[match(22:24, t.hrs)]
}

temps= as.data.frame(cbind(tmin, do.call("cbind", rep(list(temps), 30))))
colnames(temps)[2:ncol(temps)]<- 1:(ncol(temps)-1)

#to long format
temps.l<- melt(temps, id.vars = c("tmin"), variable.name = "time")

#plot temperatures, varying night temperatures
ggplot(data=temps.l, aes(x=time, y =value, color=tmin, group=tmin))+geom_line()

#combine temp data
colnames(temps.l)<- c("treatment","time","temp")
temps.l$expt<- 1

temps.all<- temps.l
temps.all$time<- as.numeric(temps.all$time)

#-----
#Assemble performance
#to long format
adat2.p.l<- melt(adat2.p[,c("NTmin","Fecundtiy")], id.vars = c("NTmin"), variable.name = "metric")

#Survival format unclear

#Development format
#just nymphal development
#subset lines since report dt and rv
adat2.dt.m1 <- na.omit(cbind( adat2.dt$NTminN[1:132], "dev_rate", 1/adat2.dt$Nymph[1:132]))
colnames(adat2.dt.m1)<- c("NTmin","metric","value")

#--------
#combine dataframes
obs<- adat2.p.l[,c("NTmin","metric","value")]
obs<- rbind(obs, adat2.dt.m1)

colnames(obs)[1]<- c("treatment")
obs$expt<- 1

PerfDat<- obs

#plot
PerfDat$value <- as.numeric(PerfDat$value)
ggplot(data=PerfDat, aes(x=treatment, y =value))+geom_point()+facet_wrap(.~metric)

#=========================================
#Expt 2. vary max
#Ma et al. 2015. Daily temperature extremes play an important role in predicting thermal effects. The Journal of Experimental Biology 218 (14), 2289-2296
#https://doi.org/10.1242/jeb.122127, provided data
#Vary daily maximum temperatures, while holding night-time temperatures constant

#read data
adat5.t<- read.csv("data/Maetal2015/Maetal2015_temps.csv")
adat5.p<- read.csv("data/Maetal2015/Maetal2015_perf.csv")

#--------------
# construct temperatures
#to long format
temps.l<- melt(adat5.t, id.vars = c("Hour"), variable.name = "tmax")

#plot temperatures, varying night temperatures
ggplot(data=temps.l, aes(x=Hour, y =value, color=tmax, group=tmax))+geom_line()

#expand to 30 days
temps= do.call("rbind", rep(list(adat5.t), 30))
colnames(temps)[1]<-"time"
temps$time<- 1:720

#to long format
temps.l<- melt(temps, id.vars = c("time"), variable.name = "tmax")
temps.l$time= as.numeric(temps.l$time)
temps.l$value= as.numeric(temps.l$value)

ggplot(data=temps.l, aes(x=time, y =value, color=tmax, group=tmax))+geom_line()

#combine temp data
temps.l<- temps.l[,c("tmax","time","value")]
colnames(temps.l)<- c("treatment","time","temp")
temps.l$expt<- 2
#update treatment name
temps.l$treatment<- gsub("Dmax_", "", temps.l$treatment)

temps.all<- rbind(temps.all, temps.l)

#--------------
# assemble performance metrics

#developmental rate
adat5.p$dr= 1/(adat5.p$DurL/24) #check conversion from hourly
#check low durations at high temperatures
#nymphal data

# #Fecundity, Longevity
# adat5.p.m= adat5.p %>%
#   group_by(Dmax_C) %>%
#   summarise(lon= mean(longevity, na.rm = TRUE), fec=mean(fecundity, na.rm = TRUE), dr=mean(dr, na.rm = TRUE) )

#use all data
adat5.p.m= adat5.p[,c("Dmax_C","longevity","fecundity","dr")]

#to long format
adat5.p.l<- melt(adat5.p.m, id.vars = c("Dmax_C"), variable.name = "metric")
adat5.p.l$treatment<- adat5.p.l$Dmax_C

obs<- adat5.p.l[,c("treatment","metric","value")]
obs$expt<- 2

#drop high outliers
obs<- obs[-which(obs$metric=="dr" & obs$value>0.15),]

PerfDat<- rbind(PerfDat, obs[,c("treatment","metric","value","expt")])

ggplot(data=adat5.p, aes(x=Dmax_C, y =dr))+geom_point()

#=====
#Expt 3: vary variance (expt 3, mild means)
#Expt 4: vary means and variance (expt 3, high means)

#Wang, XJ., Ma, CS. Can laboratory-reared aphid populations reflect the thermal performance of field populations in studies on pest science and climate change biology?. J Pest Sci 96, 509–522 (2023). https://doi.org/10.1007/s10340-022-01565-6
#Sitobion avenae

#read data
#mild means
adat4.var<- read.csv("data/WangMa2023/WangMa2023_temp22mean_diffvar.csv")
#high means
adat4.mean<- read.csv("data/WangMa2023/WangMa2023_diffmeans.csv")
#popgrowth
adat4.r<- read.csv("data/WangMa2023/WangMa2023_popgrowth.csv")

#----------
#set up temperatures
#mean: 22 °C; ±0 °C,±5 °C,±9 °C,±13 °C as fluctuating ranges, and designed four #3 high constant (28 °C, 30 °C, 32 °C) and fluctuating regimes (28 ± 5 °C, 30±5 °C, 32±5 °C)

temps= matrix(NA, nrow= 10, ncol=24)
temps[1,]<- 22
temps[2,]<- diurnal_temp_variation_sine(T_max = 27, T_min = 17, 1:24)
temps[3,]<- diurnal_temp_variation_sine(T_max = 31, T_min = 13, 1:24)
temps[4,]<- diurnal_temp_variation_sine(T_max = 35, T_min = 9, 1:24)
temps[5,]<- 28
temps[6,]<- 30
temps[7,]<- 32
#28 ± 5
temps[8,]<- diurnal_temp_variation_sine(T_max = 33, T_min = 23, 1:24)
#30±5
temps[9,]<- diurnal_temp_variation_sine(T_max = 35, T_min = 25, 1:24)
#32±5
temps[10,]<- diurnal_temp_variation_sine(T_max = 37, T_min = 27, 1:24)

##expand to 30 days
temps= do.call("cbind", rep(list(temps), 30))
colnames(temps)[2:ncol(temps)]<- 1:(ncol(temps)-1)

Tmean= c(22, 22, 22, 22, 28, 30, 32, 28, 30, 32)
Tvar= c(0, 5, 9, 13, 0, 0, 0, 5, 5, 5)
treat= c("mild means", "mild means", "mild means", "mild means", "high means", "high means", "high means", "high means", "high means", "high means")
temps= as.data.frame(cbind(Tmean, Tvar, treat, temps))
temps$Tmean_var= paste(temps$Tmean, temps$Tvar, sep="_")
temps$expt<-3
temps$expt[which(temps$treat=="high means")]<-4

#to long format
temps.l<- melt(temps, id.vars = c("Tmean", "Tvar", "treat", "Tmean_var","expt"), variable.name = "time")
temps.l$time= as.numeric(temps.l$time)
temps.l$value= as.numeric(temps.l$value)

ggplot(data=temps.l, aes(x=time, y =value, color=Tvar, group=Tmean_var))+geom_line()+
  facet_grid(.~treat, scale="free_y")

#combine temp data
temps.l<- temps.l[,c("Tmean_var","time","value","expt")]
colnames(temps.l)<- c("treatment","time","temp","expt")

temps.all<- rbind(temps.all, temps.l)

#----------
#Tmean of 22, difference variance

#drop NA
adat4.var.m<- adat4.var[-which(is.na(adat4.var$value)),]

adat4.var.m$treatment= paste(adat4.var.m$Tmean, adat4.var.m$Tvar, sep="_")

ggplot(data=adat4.var.m, aes(x=Tvar, y =value, color=population))+geom_point()+ geom_line()+
  facet_grid(metric~Tmean, scale="free_y")

#---------------
#combine performance metrics

#add performance data
obs<- adat4.var.m[,c("treatment","metric","value","population")]
obs$expt<- 3
obs$expt[which(obs$treatment %in% c("28_0", "28_5", "30_0", "30_5", "32_0", "32_5"))]<- 4

PerfDat$population<- NA
PerfDat<- rbind(PerfDat, obs[,c("treatment","metric","value","expt","population")])

#=========================================
#Expt 5
#Ma CS, Wang L, Zhang W, Rudolf V 2018. Resolving biological impacts of multiple heat waves: interaction of hot and recovery days. Oikos 127:622–33
#https://doi.org/10.1111/oik.04699
#data: http://dx.doi.org/10.5061/dryad.5qk4s
#Sitobion avenae

#vary number of successive hot days (1–3 days) and normal interval days (1–3 days)
#either exposed aphids first to hot days or normal days

#Finding: Increasing the duration of hot days in heat waves had a negative effect on various demographic rates and life-time fitness of individuals, but magnitude of this effect was typically contingent on the temporal clustering of hot periods.

#load biological data
adat3.dev<- read.csv("data/Maetal2017/Development.csv")
adat3.rep<- read.csv("data/Maetal2017/Reproduction.csv")
adat3.tr<- read.csv("data/Maetal2017/Traits.csv")

#----
#generate temps
#hot: 35°C and 20°C as the daily maximum and minimum temperatures
#normal: 28°C and 13°C as the daily maximum and minimum temperatures
#sine curve to simulate the diurnal fluctuations with a magnitude of 15°C 

hrs<- 0:23
plot(hrs, diurnal_temp_variation_sine(T_max = 35, T_min = 20, hrs))

#make temperature sequences
hdt= c(diurnal_temp_variation_sine(T_max = 35, T_min = 20, 4:23), diurnal_temp_variation_sine(T_max = 35, T_min = 20, 0:3))
ndt= c(diurnal_temp_variation_sine(T_max = 28, T_min = 13, 4:23), diurnal_temp_variation_sine(T_max = 28, T_min = 13, 0:3))
#fix transition

x<- temps[1,1:3]
x[1:2]<- as.numeric(x[1:2])
tz= rep(c(rep(ndt, x[2]), rep(hdt, x[1])), ceiling(40/(x[1]+x[2])))
plot(1:100, tz[1:100], type="l")

#set up matrix of 30 days of hourly temps
trt<- expand.grid(hd=1:3, nd=1:3, first= 1:2) #first: 1 is n, 2 is h
temps<- matrix(NA, nrow=nrow(trt), ncol=24*30)
temps= as.data.frame(cbind(trt, temps))

build.temps<- function(x){
  x[1]<- as.numeric(x[1]); x[2]<- as.numeric(x[2])
if(x[3]==1) ts<- rep(c(rep(ndt, x[2]), rep(hdt, x[1])), ceiling(40/(x[1]+x[2])))
if(x[3]==2) ts<- rep(c(rep(hdt, x[1]), rep(ndt, x[2])), ceiling(40/(x[1]+x[2])))
  return(ts[1:720])
}

#set up temperatures  
tz= t(apply(temps[,1:3], MARGIN=1, FUN=build.temps))
temps[, 4:ncol(temps)]<- tz

#to long format
temps.l<- melt(temps, id.vars = c("hd","nd","first"), variable.name = "time")
temps.l$treat<- paste(temps.l$hd, temps.l$nd, sep="_")
temps.l$treat.nh<- paste(temps.l$hd, temps.l$nd, temps.l$first, sep="_")
temps.l$time<-as.numeric(temps.l$time)

#plot temperatures
plot(1:100, tz[1,1:100], type="l")
plot(1:200, temps[1,4:203], type="l")
ggplot(data=temps.l, aes(x=time, y =value, color=hd, group=treat.nh))+geom_line()+facet_wrap(.~treat)+xlim(0,100)

#combine temp data
temps.l<- temps.l[,c("treat.nh","time","value")]
colnames(temps.l)<- c("treatment","time","temp")
temps.l$expt<- 5

temps.all<- rbind(temps.all, temps.l)

#----
#Traits
#use development time for nymphs
adat3.tr$dt= adat3.tr$NymphDur

#developmental rate
adat3.tr$dr= 1/adat3.tr$dt
#fedundity rate
adat3.tr$fec.rate= adat3.tr$Fecundity / adat3.tr$Longevity

#to long format
adat3.l<- melt(adat3.tr[,c("Treatment","ID","H_C","CycleDays","ContinueNormalDay","ContinueHotday","Longevity","Fecundity","BirthDur","dr","fec.rate")], 
               id.vars = c("Treatment","ID","H_C","CycleDays","ContinueNormalDay","ContinueHotday"), variable.name = "metric")
adat3.l$ContinueNormalDay <- as.factor(adat3.l$ContinueNormalDay)

#NymphDur, Longevity, Fecundity, BirthDur
fig.trait<- ggplot(data=adat3.l, aes(x=ContinueHotday, y =value, color=ContinueNormalDay, group=ContinueNormalDay))+
  geom_smooth(method='lm') +geom_point()+
  facet_grid(H_C~metric, scales="free_y")

#----
#Other data sets

#development:nymphal period
#find first zero 
dt= apply(adat3.dev[,5:ncol(adat3.dev)], MARGIN=1, FUN=function(x)max(which(x==1)) )
dt[is.infinite(dt)]<-NA
adat3.dev$dev<- dt

#mean metric
adat3.dev.m= adat3.dev %>%
  group_by(TRT, Hot_Cold) %>%
  summarise( dev= mean(dev, na.rm=TRUE) )

ggplot(data=adat3.dev, aes(x=TRT, y =dev, color=Hot_Cold))+geom_point()

ggplot(data=adat3.dev.m, aes(x=TRT, y =dev, color=Hot_Cold))+geom_point()

#----
#Reproduction
adat3.rep$rep<- colSums(adat3.rep[,5:ncol(adat3.rep)])

#mean metric
adat3.rep.m= adat3.rep %>%
  group_by(Treatment, Hot_Cold) %>%
  summarise( rep= mean(rep, na.rm=TRUE) )

ggplot(data=adat3.rep, aes(x=Treatment, y =rep, color=Hot_Cold))+geom_point()

ggplot(data=adat3.rep.m, aes(x=Treatment, y =rep, color=Hot_Cold))+geom_point()

#--------------------------
#add performance data
obs<- adat3.l[,c("H_C","ContinueNormalDay","ContinueHotday","metric","value")]

obs$hc<- match(obs$H_C, c("C","H") )
obs$treatment<- paste(obs$ContinueHotday, obs$ContinueNormalDay, obs$hc, sep="_")
obs$expt<- 5
obs$population <- NA

PerfDat<- rbind(PerfDat, obs[,c("treatment","metric","value","expt","population")])

#=========================================
#Expt 6: vary length of heatwave and timing (expt 5 adult)
#Expt 7: vary length of heatwave and timing (expt 5 nymph)

#Zhao et al. 2019. The importance of timing of heat events for predicting the dynamics of aphid pest populations. Pest management science
#https://doi.org/10.1002/ps.5344, provided data
#survival and productivity
#Sitobion avenae 
#using heat stress (20–35°C diurnal cycle) across the nymph and adult stages
#four timings [early nymph (NE), late nymph (NL), early adult (AE) as well as late adult (AL)] and six durations (1, 2, 3, 4, 5 and 6 consecutively hot days)
#Fecundity (nymphs / adult): focus on NE or NL
#Hot day: 20–35C; Normal day: 13–28C

#read data
adat6.t<- read.csv("data/Zhaoetal2019/Zhaoetal2019_temps.csv")
adat6.p<- read.csv("data/Zhaoetal2019/Zhaoetal2019_perf.csv")

#--------------
# construct temperatures
adat6.t$Temp...C<- as.numeric(adat6.t$Temp...C)

min<- as.numeric(format(strptime(adat6.t$Time..GMT.08.00, format="%y/%m/%d %H:%M"),'%M'))
# #restrict to 20 minute interval
# adat6.t<- adat6.t[min %in% c(0,20,40),]
# min<- min[min %in% c(0,20,40)]
#restrict to hourly data
adat6.t<- adat6.t[min %in% c(0),]
min<- min[min %in% c(0)]

hrs<- as.numeric(format(strptime(adat6.t$Time..GMT.08.00, format="%y/%m/%d %H:%M"),'%H'))
adat6.t$doy<- as.numeric(format(strptime(adat6.t$Time..GMT.08.00, format="%y/%m/%d %H:%M"),'%j'))
# #restrict to 20 minute interval
# adat6.t<- adat6.t[min %in% c(0,20,40),]
#restrict to hourly data
adat6.t<- adat6.t[min %in% c(0),]

time<- hrs+min/60

runs<- rle(adat6.t$am_pm)

time[which(adat6.t$am_pm=="pm" & hrs<12)]= time[which(adat6.t$am_pm=="pm" & hrs<12)] +12
time[which(adat6.t$am_pm=="am" & hrs==12)]= time[which(adat6.t$am_pm=="am" & hrs==12)] -12
adat6.t$time<- time

#plot normal and high treatment temperatures
ggplot(data=adat6.t, aes(x=time, y =Temp...C, color=treatment, group=treatment))+geom_line()

#duration normal temps
unique(adat6.t$doy[which(adat6.t$treatment=="normal")])
#experiment run until all dead, but truncate at 30 days

#Focus on AE: adult early and NE: nymphal late
temps.n= adat6.t[which(adat6.t$treatment=="normal"),]
temps.n$day<- temps.n$doy -temps.n$doy[1] + 1
#restrict to 30 days
temps.n= temps.n[which(temps.n$day<=30),]
temps.n<- temps.n[,c("ind","Temp...C","time","day")]
temps.n$dt<- temps.n$day +temps.n$time/24

temps.h= adat6.t[which(adat6.t$treatment=="high"),]
temps.h= temps.h[which(temps.h$doy>=251),]
temps.h$day<- temps.h$doy -temps.h$doy[1] + 1
temps.h<- temps.h[,c("ind","Temp...C","time","day")]
temps.h$dt<- temps.h$day +temps.h$time/24

#fill nas
temps.h$Temp...C<- na.approx(temps.h$Temp...C)
temps.n$Temp...C<- na.approx(temps.n$Temp...C)

#AE1: 8-9
temps.n$AE1 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(8:9))
temps.n$AE1[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#AE2: 8-10
temps.n$AE2 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(8:10))
temps.n$AE2[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#AE3: 8-11
temps.n$AE3 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(8:11))
temps.n$AE3[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#AE4: 8-12
temps.n$AE4 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(8:12))
temps.n$AE4[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#AE5: 8-13
temps.n$AE5 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(8:13))
temps.n$AE5[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#AE6: 8-14
temps.n$AE6 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(8:14))
temps.n$AE6[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#Late nymphal
#NL6: 1-7
temps.n$NL6 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(1:7))
temps.n$NL6[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#NL5: 2-7
temps.n$NL5 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(2:7))
temps.n$NL5[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#NL4: 3-7
temps.n$NL4 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(3:7))
temps.n$NL4[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#NL3: 4-7
temps.n$NL3 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(4:7))
temps.n$NL3[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#NL2: 5-7
temps.n$NL2 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(5:7))
temps.n$NL2[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#NL1: 6-7
temps.n$NL1 <- temps.n$Temp...C
inds<- which(temps.n$day %in% c(6:7))
temps.n$NL1[inds]<- temps.h$Temp...C[match(temps.n$dt[inds], temps.h$dt)]

#to long format
temps.n <- temps.n[,c("ind","AE1","AE2","AE3","AE4","AE5","AE6","NL1","NL2","NL3","NL4","NL5","NL6")]
#fix hour counter
temps.n$ind=1:nrow(temps.n)
temps.l<- melt(temps.n, id.vars = c("ind"), variable.name = "treatment")

#fill 1 nas
temps.l$value <- na.approx(temps.l$value)

#plot temperatures, varying night temperatures
ggplot(data=temps.l, aes(x=ind, y =value, color=treatment, group=treatment))+geom_line()

#----
#combine temp data
colnames(temps.l)<- c("time","treatment","temp")
temps.l$expt<- 6
temps.l$expt[which(temps.l$treatment %in% c("NL1","NL2","NL3","NL4","NL5","NL6"))]<- 7

temps.all<- rbind(temps.all, temps.l)

#--------------
# assemble performance metrics

adat6.p$dr <- 1/adat6.p$Nymph.duration

#Fecundity, Longevity
#Also life period and immediate.death

#reduce to late nymphal and early adult treatments
adat6.p.m<- adat6.p[adat6.p$Treatments %in% c("AE1","AE2","AE3","AE4","AE5","AE6","NL1","NL2","NL3","NL4","NL5","NL6"),]
adat6.p.m$fec <- adat6.p.m$Productivity
adat6.p.m<- adat6.p.m[,c("Treatments","Longevity","fec","dr")]

#to long format
adat6.p.l<- melt(adat6.p.m, id.vars = c("Treatments"), variable.name = "metric")

#drop NAs
adat6.p.l<- adat6.p.l[which(!is.na(adat6.p.l$value)),]

obs<- adat6.p.l[,c("Treatments","metric","value")]
obs$treatment<- obs$Treatments
obs$expt<- 6
obs$expt[which(obs$treatment %in% c("NL1","NL2","NL3","NL4","NL5","NL6"))]<- 7

obs$population<- NA

PerfDat<- rbind(PerfDat, obs[,c("treatment","metric","value","expt","population")])

#====================
#write out data sets
PerfDat$metric= as.character(PerfDat$metric)

#align names
PerfDat$metric[which(PerfDat$metric=="fec")]="fecundity"
PerfDat$metric[which(PerfDat$metric=="Fecundity")]="fecundity"
PerfDat$metric[which(PerfDat$metric=="Fecundtiy")]="fecundity"

PerfDat$metric[which(PerfDat$metric=="lon")]<-"longevity"
PerfDat$metric[which(PerfDat$metric=="Longevity")]="longevity"

PerfDat$metric[which(PerfDat$metric=="dr")]="dev_rate"

#drop NAs
PerfDat<- PerfDat[which(!is.na(PerfDat$value)),]

#write out
write.csv(temps.all, "processed_data/TempTimeSeries.csv")
write.csv(PerfDat, "processed_data/PerformanceData.csv")

temps.all <-read.csv("3_theory/data/buckley_data/TempTimeSeries.csv")
PerfDat <- read.csv("3_theory/data/buckley_data/PerformanceData.csv")

#plot temperatures for seven days
ggplot(data=temps.all[which(temps.all$time<169),], aes(x=time, y =temp, color=treatment))+geom_line()+facet_wrap(.~expt)

#plot fecundity
ggplot(data=PerfDat[which(PerfDat$metric=="fecundity"),], aes(x=treatment, y =value))+geom_point()+facet_wrap(.~expt)

ggplot(data=PerfDat[which(PerfDat$metric=="dev_rate"),], aes(x=treatment, y =value))+geom_point()+facet_wrap(.~expt)


