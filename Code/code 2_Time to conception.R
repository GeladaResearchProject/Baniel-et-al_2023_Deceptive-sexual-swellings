rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(lubridate)
library(dplyr)
library(plyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(plyr)
library(survival)
library(survminer)
library(coxme) 
theme_set (theme_classic(base_size=20))


### Table S3.COX PROPORTIONAL HAZARDS MODELS ###
cox<-read.csv("~/Dataset/dataset 3_Cox model_Time to conception.txt", header = TRUE, sep= ",",na.strings="NA")
head(cox)
dim(cox)#3746
summary(cox)
length(unique(cox$Infant))
length(unique(cox$Female))
length(unique(cox$Unit))
length(unique(cox$Infant[cox$Takeover.fix %in% 1]))#39 infants with takeover before 18 mos
cox$Conception<-as.integer(as.character(cox$Conception))
cox$Sex.infant<-as.factor(as.character(cox$Sex.infant))
cox$Parity<-as.factor(as.character(cox$Parity))
cox$Takeover.fix<-as.factor(as.character(cox$Takeover.fix))
cox$Infant<-as.factor(as.character(cox$Infant))
cox$Female<-as.factor(as.character(cox$Female))
cox$Season<-as.factor(as.character(cox$Season))
cox$Rainfall.cat<-as.factor(as.character(cox$Rainfall.cat))
cox$Temperature.cat<-as.factor(as.character(cox$Temperature.cat))
cox$Unit<-as.factor(as.character(cox$Unit))
cox$log.Cumulative.rainfall<-log10(cox$Cumulative.rainfall+0.1)
cox$log.Avg.temperature<-log10(cox$Avg.temperature)


### MODEL 2a - Takeover.fix as fixed predictor ####
fit1 <- coxph(Surv(Time1, Time2, Conception) ~ Takeover.fix + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)#Takeover.fix violates proportional hazard assumption
plot(cox.zph(fit1)[1],ylab="Beta(t) for Takeover.fix (fixed predictor)",xlab="Time to conception (mos)")
abline(h=0, col="blue",lty=2)

#Final model (w/ time-transformation of Takeover.fix)
fit2 <- coxph(Surv(Time1, Time2, Conception) ~ Takeover.fix + tt(Takeover.fix) + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit2)

#Conservative dataset
table(cox$Dataset)
cox.cons<-droplevels(subset(cox,cox$Dataset %in% "conservative"))
dim(cox.cons)
fit1 <- coxph(Surv(Time1, Time2, Conception) ~ Takeover.fix + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit1)
cox.zph(fit1)
fit2 <- coxph(Surv(Time1, Time2, Conception) ~ Takeover.fix + tt(Takeover.fix) + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit2)


### PLOT SURVIVAL CURVES ###

#Takeover.fix
a <- survfit(Surv(Time1, Time2, Conception) ~ Takeover.fix, data = cox)
ggsurvplot(a,conf.int = TRUE)+ ylab("Proportion of mos without Conception") + xlab("Time to conception (mos)") 

#Parity
b <- survfit(Surv(Time1, Time2, Conception) ~ Parity, data = cox)
ggsurvplot(b,conf.int = TRUE)+ ylab("Proportion of mos without Conception") + xlab("Time to conception (mos)") 

#Sex.infant
c <- survfit(Surv(Time1, Time2, Conception) ~ Sex.infant, data = cox)
ggsurvplot(c,conf.int = TRUE)+ ylab("Proportion of mos without Conception") + xlab("Time to conception (mos)") 

#Seasonal effects
d <- survfit(Surv(Time1, Time2, Conception) ~ Season, data = cox)
ggsurvplot(d,conf.int = TRUE)+ ylab("Proportion of mos without Conception") + xlab("Time to conception (mos)") 

e<- survfit(Surv(Time1, Time2, Conception) ~  Rainfall.cat, data = cox)
ggsurvplot(e,conf.int = TRUE)+ ylab("Proportion of mos without Conception") + xlab("Time to conception (mos)") 

f<- survfit(Surv(Time1, Time2, Conception) ~ Temperature.cat, data = cox)
ggsurvplot(f,conf.int = TRUE)+ ylab("Proportion of mos without Conception") + xlab("Time to conception (mos)") 




### Table S4.LMM ANALYSIS ###
source("~/code/vif.R") 
df<- read.csv("~/Dataset/dataset 4_LMM_Time to conception.txt", header = TRUE, sep= ",",na.strings="NA")
head(df) 
dim(df)#127
summary(df)
df$Takeover.before18mos<-factor(Takeover.before18mos,levels=c("No takeover","Takeover"))
length(unique(df$Unit))
length(unique(df$Female))
length(unique(df$Infant))
table(df$Takeover.before18mos)#88 non-takeover infants and 39 takeover infants

#When infants experienced the takeover before birth, set age at takeover at zero month
sort(df$Age.takeover)
df$Age.takeover[df$Age.takeover<0 & df$Age.takeover %!in% NA]<-0

#Time to conception for non-takeover and takeover infants
summary(df$Time.to.conception[df$Takeover.before18mos %in% "No takeover"])
summary(df$Time.to.conception[df$Takeover.before18mos %in% "Takeover"])
 

### MODEL 2b - ALL INFANTS ###
m1<-lmer(Time.to.conception ~  Takeover.before18mos + Parity + Sex.infant + (1|Female) + (1|Unit), data=df)
summary(m1)
vif.mer(m1)
hist(resid(m1),nclass=100)
plot(resid(m1)~fitted(m1),nclass=100)
qqnorm(resid(m1))
qqline(resid(m1))

#Conservative dataset
table(df$Dataset)
df.cons<-droplevels(subset(df,df$Dataset %in% "conservative"))
dim(df.cons)#122
table(df.cons$Takeover.before18mos)
m1.cons<-lmer(Time.to.conception ~  Takeover.before18mos + Parity + Sex.infant + (1|Female) + (1|Unit), data=df.cons)
summary(m1.cons)

#Plot
a<-ggplot(df, aes(x=Time.to.conception, color=Takeover.before18mos))+ geom_histogram(binwidth=0.5) +xlab("Time.to.conception  (mos)")+ scale_color_manual(values=c("#999999","coral3"))
b<-ggplot(df,aes(y=Time.to.conception, x=Takeover.before18mos,fill=Takeover.before18mos))+ geom_boxplot()+ geom_jitter()+ylab("Time to conception (mos)")+xlab ("")+ scale_fill_manual(values=c("#999999","coral3"))
c<-ggplot(df,aes(y=Time.to.conception, x=Parity))+ geom_boxplot() + geom_jitter()+ ylab("Time to conception (mos)")+xlab ("Maternal parity")
d<-ggplot(df,aes(y=Time.to.conception, x=Sex.infant))+ geom_boxplot()+ geom_jitter()+ ylab("Time to conception (mos)")+xlab ("Sex infant")
ggarrange(a,b,c,d,nrow=2,ncol=2)


### MODEL 2c -  ONLY TAKEOVER INFANTS  ###
df2<-droplevels(subset(df,df$Takeover.before18mos %in% "Takeover"))
dim(df2)#39
m2<-lmer(Time.to.conception ~  scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit), data=df2)
summary(m2)

#Conservative dataset 
df2.cons<-droplevels(subset(df2,df2$Dataset %in% "conservative"))
dim(df2.cons)
m2.cons<-lmer(Time.to.conception ~  scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit), data=df2.cons)
summary(m2.cons)

#Plot age at takeover (continuous)
ggplot(df2,aes(y=Time.to.conception, x=Age.takeover))+ geom_point(size=2) + ylab("Time to conception (mos)")+ xlab ("Age of infant at Takeover.fix (mos)")+ geom_smooth(method='lm', formula= y~x)


