rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(lubridate)
library(dplyr)
library(plyr)
library(lme4)
library(lmerTest)
library(jtools)
library(ggplot2)
library(ggpubr)
library(plyr)
library(survival)
library(survminer)
library(coxme) 
library(ggbeeswarm)
theme_set (theme_classic(base_size=15))


### Table S8-S9 - COX PROPORTIONAL HAZARDS MODELS ###
cox<-read.csv("~/Dataset/dataset 7_Cox model_Post-takeover infant survival.txt", header = TRUE, sep= ",",na.strings="NA")
head(cox)
dim(cox)#975
summary(cox)
length(unique(cox$Infant))
length(unique(cox$Female))
length(unique(cox$Unit))
length(unique(cox$Infant[cox$Death %in% 1]))#33 infant deaths

#When infants experienced the takeover before birth, set age at takeover at zero month
sort(unique(cox$Age.takeover))
cox$Age.takeover[cox$Age.takeover<0]<-0

cox$Death<-as.integer(as.character(cox$Death))
cox$Swelling<-as.factor(as.character(cox$Swelling))
cox$Sex.infant<-as.factor(as.character(cox$Sex.infant))
cox$Parity<-as.factor(as.character(cox$Parity))
cox$Infant<-as.factor(as.character(cox$Infant))
cox$Female<-as.factor(as.character(cox$Female))
cox$Unit<-as.factor(as.character(cox$Unit))

#Conservative dataset
table(cox$Dataset)
cox.cons<-droplevels(subset(cox,cox$Dataset %in% "conservative"))
dim(cox.cons)

#Remove infants dying within 30 days of a takeover
cox$Date.takeover<-as.Date(cox$Date.takeover,"%Y-%m-%d")
cox$Date.death<-as.Date(cox$Date.death,"%Y-%m-%d")
cox$Post.takeover.survival<-as.numeric(cox$Date.death-cox$Date.takeover)
rem<-unique(cox$Infant[cox$Post.takeover.survival %!in% NA & cox$Post.takeover.survival<=30])
length(rem)#9
cox1<-droplevels(subset(cox,cox$Infant %!in% rem))
length(unique(cox1$Infant))
length(unique(cox1$Infant[cox1$Death %in% 1]))


### MODEL 5a - AGE AT TAKEOVER (fix) ### 
fit1 <- coxph(Surv(Time1, Time2, Death) ~ Swelling + Age.takeover + Parity + Sex.infant  + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)

#Conservative dataset
fit2 <- coxph(Surv(Time1, Time2, Death) ~ Swelling + Age.takeover + Parity + Sex.infant  + cluster(Female), data=cox.cons)
summary(fit2)
cox.zph(fit2)


### MODEL 5b - INTERACTION SWELLING * AGE AT TAKEOVER (fix) ### 
fit1 <- coxph(Surv(Time1, Time2, Death) ~ Swelling * Age.takeover + Parity + Sex.infant  + cluster(Female), data=cox)
summary(fit1)

#Conservative dataset
fit2 <- coxph(Surv(Time1, Time2, Death) ~ Swelling * Age.takeover + Parity + Sex.infant  + cluster(Female), data=cox.cons)
summary(fit2)


### MODELS 5c and 5d - Remove infants dying within 40 days of a takeover ### 
fit1 <- coxph(Surv(Time1, Time2, Death) ~ Swelling + Age.takeover + Parity + Sex.infant  + cluster(Female), data=cox1)
summary(fit1)

fit2 <- coxph(Surv(Time1, Time2, Death) ~ Swelling * Age.takeover + Parity + Sex.infant  + cluster(Female), data=cox1)
summary(fit2)


### MODEL 5e - AGE OF INFANT (fix) ### 
fit1 <- coxph(Surv(Time1, Time2, Death) ~ Swelling + Age.infant + Parity + Sex.infant  + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)

#Conservative dataset
fit2 <- coxph(Surv(Time1, Time2, Death) ~ Swelling + Age.infant + Parity + Sex.infant  + cluster(Female), data=cox.cons)
summary(fit2)
cox.zph(fit2)


### MODEL 5f - INTERACTION SWELLING * AGE AT TAKEOVER (fix) ### 
fit1 <- coxph(Surv(Time1, Time2, Death) ~ Swelling * Age.infant + Parity + Sex.infant  + cluster(Female), data=cox)
summary(fit1)

#Conservative dataset
fit2 <- coxph(Surv(Time1, Time2, Death) ~ Swelling * Age.infant + Parity + Sex.infant  + cluster(Female), data=cox.cons)
summary(fit2)


### MODELS 5g and 5h - Remove infants dying within 40 days of a takeover ### 
fit1 <- coxph(Surv(Time1, Time2, Death) ~ Swelling + Age.infant + Parity + Sex.infant  + cluster(Female), data=cox1)
summary(fit1)

fit2 <- coxph(Surv(Time1, Time2, Death) ~ Swelling * Age.infant + Parity + Sex.infant  + cluster(Female), data=cox1)
summary(fit2)



### PLOT SURVIVAL CURVES ###

#Maternal swellings
a <- survfit(Surv(Time1, Time2, Death) ~ Swelling, data = cox)
ggsurvplot(a,conf.int = TRUE)+ ylab("Proportion of surviving infants") + xlab("Time post-takeover (mos)") 

#Age at takeover
cox$Age.takeover.cat<-rep("Med",length(cox$Infant))
cox$Age.takeover.cat[cox$Age.takeover<=6]<-"Young"
cox$Age.takeover.cat[cox$Age.takeover>12]<-"Old"
cox$Age.takeover.cat<-factor(cox$Age.takeover.cat,levels=c("Young","Med","Old"))
tab<-unique(cox[c("Infant","Age.takeover.cat")])
table(tab$Age.takeover.cat)
b <- survfit(Surv(Time1, Time2, Death) ~ Age.takeover.cat, data = cox)
ggsurvplot(b,conf.int = T)+ ylab("Proportion of surviving infants") + xlab("Time post-takeover (mos)") 

#Age infant
cox$Age.infant.cat<-rep("Med",length(cox$Infant))
cox$Age.infant.cat[cox$Age.infant<=6]<-"Young"
cox$Age.infant.cat[cox$Age.infant>=12]<-"Old"
cox$Age.infant.cat<-factor(cox$Age.infant.cat,levels=c("Young","Med","Old"))
tab<-unique(cox[c("Infant","Age.infant.cat")])
table(tab$Age.infant.cat)
c <- survfit(Surv(Time1, Time2, Death) ~ Age.infant.cat, data = cox)
ggsurvplot(c,conf.int = T)+ ylab("Proportion of surviving infants") + xlab("Time post-takeover (mos)") 

#Mum swell * age at takeover
cox$Age.takeover.cat<-rep("Old",length(cox$Infant))
cox$Age.takeover.cat[cox$Age.takeover<=12]<-"Young"
cox$Interaction1<-paste(cox$Swelling,cox$Age.takeover.cat,sep="_")
table(cox$Interaction1)
inter<-droplevels(subset(cox,cox$Interaction1 %!in% c("0_Old")))
d<- survfit(Surv(Time1, Time2, Death) ~ Interaction1, data = inter)
ggsurvplot(d,conf.int = T)+ ylab("Proportion of surviving infants") + xlab("Time post-takeover (mos)") 

#Mum swell * age infant
cox$Age.infant.cat<-rep("Old",length(cox$Infant))
cox$Age.infant.cat[cox$Age.infant<=12]<-"Young"
cox$Interaction2<-paste(cox$Swelling,cox$Age.infant.cat,sep="_")
table(cox$Interaction2)
inter<-droplevels(subset(cox,cox$Interaction2 %!in% c("0_Old")))
e <- survfit(Surv(Time1, Time2, Death) ~ Interaction2, data = inter)
ggsurvplot(e,conf.int = T)+ ylab("Proportion of surviving infants") + xlab("Time post-takeover (mos)") 

#Parity
f <- survfit(Surv(Time1, Time2, Death) ~ Parity, data = cox)
ggsurvplot(f,conf.int = TRUE)+ ylab("Proportion of unweaned infants") + xlab("Age of infant (mos)") 

#Sex infant
g <- survfit(Surv(Time1, Time2, Death) ~ Sex.infant, data = cox)
ggsurvplot(g,conf.int = TRUE)+ ylab("Proportion of unweaned infants") + xlab("Age of infant (mos)") 




### Table S10.LMM ANALYSIS ###
source("~/code/vif.R") 
df<- read.csv("~/Dataset/dataset 8_LMM_Post-takeover infant survival.txt", header = TRUE, sep= ",",na.strings="NA")
head(df) 
dim(df)#105
df$Survive.24mos<-as.factor(as.character(df$Survive.24mos))

#When infants experienced the takeover before birth, set age at takeover at zero month
sort(df$Age.takeover)
df$Age.takeover[df$Age.takeover<0 & df$Age.takeover %!in% NA]<-0

#Conservative dataset
table(df$Dataset)
df.cons<-droplevels(subset(df,df$Dataset %in% "conservative"))
dim(df.cons)

#Plot - infant survival to 2y according to whether mum swell within 3 mos of a takeover
table(df$Swelling.within.3mos)
p<-ggplot(df,aes(y=Age.takeover, x=Survive.24mos,colour=Swelling.within.3mos))+
  geom_quasirandom(size=4)+ ylim(0,20)+
  xlab ("Infant survived to 24 mos?") +ylab("Age of infant at takeover (mos)")+
  theme(legend.position="top",
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18,color="black"))+
  scale_color_manual(values=c("royalblue3","deeppink4"))
p+labs(colour="Mother swells within 3 mos of takeover")

#Plot - Latency to return to swelling following takeovers
sort(df$Latency.takeover.to.swelling)
sw<-droplevels(subset(df,df$Latency.takeover.to.swelling %!in% NA))
ggplot(sw,aes(y=Latency.takeover.to.swelling, x=Age.takeover,colour=Survive.24mos)) + geom_point(size=4,alpha=0.8)+
  xlab ("Age of infant at takeover (mos)") + ylab("Number mos to return to swelling following a takeover")+ 
  theme(legend.position="top")+
  scale_color_manual(values=c("mediumpurple3","lightseagreen"))+labs(colour="Infant survives to 24 mos?")


### MODEL 5i - ALL INFANTS ###
m1<-glmer(Survive.24mos ~ Swelling.within.3mos  + scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit),family=binomial,data=df,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
summary(m1)
vif.mer(m1)

#Conservative dataset
m1<-glmer(Survive.24mos ~ Swelling.within.3mos  + scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit),family=binomial,data=df.cons,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
summary(m1)


### MODEL 5j - ALL INFANTS ###
m2<-glmer(Survive.24mos ~ Swelling.within.3mos  * scale(Age.takeover) + Parity + Sex.infant + (1|Female)+ (1|Unit),family=binomial,data=df,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
anova(m1,m2,test="Chisq")
summary(m2)

#Plot interaction
m2<-glmer(Survive.24mos ~ Swelling.within.3mos * Age.takeover + Parity + Sex.infant + (1|Female)+ (1|Unit),family=binomial,data=df,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
interact_plot(m2, pred = Age.takeover, modx = Swelling.within.3mos)

#Conservative dataset
m2<-glmer(Survive.24mos ~ Swelling.within.3mos * scale(Age.takeover) + Parity + Sex.infant + (1|Female)+ (1|Unit),family=binomial,data=df.cons,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
anova(m1,m2,test="Chisq")
summary(m2)

### MODEL 5k and 5l - Remove infants dying within 30 days of a takeover ### 
df$Date.takeover<-as.Date(df$Date.takeover,"%Y-%m-%d")
df$Date.death<-as.Date(df$Date.death,"%Y-%m-%d")
df$Post.takeover.survival<-as.numeric(df$Date.death-df$Date.takeover)
rem<-unique(df$Infant[df$Post.takeover.survival %!in% NA & df$Post.takeover.survival<=30])
length(rem)#9
df1<-droplevels(subset(df,df$Infant %!in% rem))
length(unique(df1$Infant))
length(unique(df1$Infant[df1$Survive.24mos %in% 1]))

m1<-glmer(Survive.24mos ~ Swelling.within.3mos  + scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit),family=binomial,data=df1,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
summary(m1)

m2<-glmer(Survive.24mos ~ Swelling.within.3mos  * scale(Age.takeover) + Parity + Sex.infant + (1|Female)+ (1|Unit),family=binomial,data=df1,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),na.action=na.fail)
anova(m1,m2,test="Chisq")
summary(m2)

