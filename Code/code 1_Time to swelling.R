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


### Table S1. COX PROPORTIONAL HAZARDS MODELS ###
cox<-read.csv("~/Dataset/dataset 1_Cox model_Time to swelling.txt", header = TRUE, sep= ",",na.strings="NA")
head(cox)
dim(cox)#3746
summary(cox)
length(unique(cox$Infant))
length(unique(cox$Female))
length(unique(cox$Unit))
length(unique(cox$Infant[cox$Takeover.fix %in% 1]))#104 infants with takeover before maternal swelling
cox$Swelling<-as.integer(as.character(cox$Swelling))
cox$Sex<-as.factor(as.character(cox$Sex))
cox$Parity<-as.factor(as.character(cox$Parity))
cox$Takeover.fix<-as.factor(as.character(cox$Takeover.fix))
cox$Takeover.timevarying<-as.factor(as.character(cox$Takeover.timevarying))
cox$Infant<-as.factor(as.character(cox$Infant))
cox$Female<-as.factor(as.character(cox$Female))
cox$Season<-as.factor(as.character(cox$Season))
cox$Rainfall.cat<-as.factor(as.character(cox$Rainfall.cat))
cox$Temperature.cat<-as.factor(as.character(cox$Temperature.cat))
cox$Unit<-as.factor(as.character(cox$Unit))
cox$log.Cumulative.rainfall<-log10(cox$Cumulative.rainfall+0.1)
cox$log.Avg.temperature<-log10(cox$Avg.temperature)


### MODEL 1a - Takeover as fixed predictor ####
fit1 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.fix + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)#Takeover.fix violates proportional hazard assumption
plot(cox.zph(fit1)[1],ylab="Beta(t) for takeover (fixed predictor)",xlab="Time to swelling (mos)")
abline(h=0, col="blue",lty=2)

#Final model (w/ time-transformation of Takeover.fix)
fit2 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.fix + tt(Takeover.fix) + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit2)

#Conservative dataset
#remove cases where we're uncertain whether the swelling occurred a bit
#before or a bit after the takeover or takeovers with uncertain dates. 
table(cox$Dataset)
cox.cons<-droplevels(subset(cox,cox$Dataset %in% "conservative"))
dim(cox.cons)

fit1 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.fix + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit1)
cox.zph(fit1)

fit2 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.fix + tt(Takeover.fix) + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit2)


### MODEL 1b - Takeover as time-varying predictor ###
fit1 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.timevarying + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)
plot(cox.zph(fit1)[1])

#Final model (w/ time-transformation of Takeover.timevarying)
fit2 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.timevarying + tt(Takeover.timevarying) + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit2)

#Conservative dataset
fit1 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.timevarying + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit1)
cox.zph(fit1)
fit2 <- coxph(Surv(Time1, Time2, Swelling) ~ Takeover.timevarying + tt(Takeover.fix) + Parity + Sex + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit2)


### PLOT SURVIVAL CURVES ###
#Takeover.fix
a <- survfit(Surv(Time1, Time2, Swelling) ~ Takeover.fix, data = cox)
ggsurvplot(a,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 

#Takeoveover.timevarying 
b <- survfit(Surv(Time1, Time2, Swelling) ~ Takeover.timevarying, data = cox)
ggsurvplot(b,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 

#Parity
c <- survfit(Surv(Time1, Time2, Swelling) ~ Parity, data = cox)
ggsurvplot(c,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 

#Sex
d <- survfit(Surv(Time1, Time2, Swelling) ~ Sex, data = cox)
ggsurvplot(d,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 

#Seasonal effects
e <- survfit(Surv(Time1, Time2, Swelling) ~ Season, data = cox)
ggsurvplot(e,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 

f<- survfit(Surv(Time1, Time2, Swelling) ~ Rainfall.cat, data = cox)
ggsurvplot(f,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 

g<- survfit(Surv(Time1, Time2, Swelling) ~ Temperature.cat, data = cox)
ggsurvplot(g,conf.int = TRUE)+ ylab("Proportion of mos where mothers had not yet resumed swelling") + xlab("Time to swelling (mos)") 



### Table S2.LMM ANALYSES ###
source("~/code/vif.R") 
df<- read.csv("~/Dataset/dataset 2_LMM_Time to swelling.txt", header = TRUE, sep= ",",na.strings="NA")
head(df) 
dim(df)#218
summary(df)
df$Takeover.before.swelling<-factor(df$Takeover.before.swelling,levels=c("No takeover","Takeover"))
length(unique(df$Unit))
length(unique(df$Female))
length(unique(df$Infant))
table(df$Takeover.before.swelling)#115 non-takeover infants and 103 takeover infants

#When infants experienced the takeover before birth, set age at takeover at zero month
sort(df$Age.takeover)
df$Age.takeover[df$Age.takeover<0 & df$Age.takeover %!in% NA]<-0

#Time to swelling for non-takeover and takeover infants
summary(df$Time.to.swell[df$Takeover.before.swelling %in% "No takeover"])
summary(df$Time.to.swell[df$Takeover.before.swelling %in% "Takeover"])


### MODEL 1c - ALL INFANTS ###
m1<-lmer(Time.to.swell ~  Takeover.before.swelling + Parity + Sex.infant + (1|Female) + (1|Unit), data=df)
summary(m1)
vif.mer(m1)
hist(resid(m1),nclass=100)
plot(resid(m1)~fitted(m1),nclass=100)
qqnorm(resid(m1))
qqline(resid(m1))

#Conservative dataset 
table(df$Dataset)
df.cons<-droplevels(subset(df,df$Dataset %in% "conservative"))
dim(df.cons)#204
table(df.cons$Takeover.before.swelling)
m1.cons<-lmer(Time.to.swell ~  Takeover.before.swelling + Parity + Sex.infant + (1|Female) + (1|Unit), data=df.cons)
summary(m1.cons)

#Plot
a<-ggplot(df, aes(x=Time.to.swell, color=Takeover.before.swelling))+ geom_histogram(binwidth=0.5) +xlab("Time to swelling  (mos)")+ scale_color_manual(values=c("#999999","coral3"))
b<-ggplot(df,aes(y=Time.to.swell, x=Takeover.before.swelling,fill=Takeover.before.swelling))+ geom_boxplot()+ geom_jitter()+ylab("Time to swelling  (mos)")+xlab ("")+ scale_fill_manual(values=c("#999999","coral3"))
c<-ggplot(df,aes(y=Time.to.swell, x=Parity))+ geom_boxplot() + geom_jitter()+ ylab("Time to swelling (mos)")+xlab ("Maternal parity")
d<-ggplot(df,aes(y=Time.to.swell, x=Sex.infant))+ geom_boxplot()+ geom_jitter()+ ylab("Time to swelling  (mos)")+xlab ("Infant sex")
ggarrange(a,b,c,d,nrow=2,ncol=2)


### MODEL 1d - ONLY TAKEOVER INFANTS  ###
df2<-droplevels(subset(df,df$Takeover.before.swelling %in% "Takeover"))
dim(df2)#103

m2<-lmer(Time.to.swell ~  scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit), data=df2)
summary(m2)

#Conservative dataset 
df2.cons<-droplevels(subset(df2,df2$Dataset %in% "conservative"))
dim(df2.cons)#89
m2.cons<-lmer(Time.to.swell ~  scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit), data=df2.cons)
summary(m2.cons)

#Plot 
#Age at takeover (continuous)
ggplot(df2,aes(y=Time.to.swell, x=Age.takeover))+ geom_point(size=2) + ylab("Time to swelling (mos)")+ xlab ("Age of infant at takeover (mos)")+ geom_smooth(method='lm', formula= y~x)

#Age at takeover (categorical)
sort(to$Age.takeover)
no.to<-droplevels(subset(df,df$Takeover.before.swelling %in% "No takeover"));dim(no.to) 
to<-droplevels(subset(df,df$Takeover.before.swelling %in% "Takeover"));dim(to)
to_3<-droplevels(subset(to,to$Age.takeover <= 3));dim(to_3)
to_6<-droplevels(subset(to,to$Age.takeover > 3 & to$Age.takeover <= 6));dim(to_6)
to_9<-droplevels(subset(to,to$Age.takeover > 6 & to$Age.takeover <= 9));dim(to_9)
to_12<-droplevels(subset(to,to$Age.takeover > 9 & to$Age.takeover <= 12));dim(to_12)
to_15<-droplevels(subset(to,to$Age.takeover > 12 & to$Age.takeover <= 15));dim(to_15)
to_18<-droplevels(subset(to,to$Age.takeover > 15 & to$Age.takeover <= 18));dim(to_18)

j1<-cbind(no.to$Time.to.swell,rep("no\ntakeover",length(no.to$Time.to.swell)))
j2<-cbind(to_3$Time.to.swell,rep("takeover\n0-3 mos",length(to_3$Time.to.swell)))
j3<-cbind(to_6$Time.to.swell,rep("takeover\n3-6 mos",length(to_6$Time.to.swell)))
j4<-cbind(to_9$Time.to.swell,rep("takeover\n6-9 mos",length(to_9$Time.to.swell)))
j5<-cbind(to_12$Time.to.swell,rep("takeover\n9-12 mos",length(to_12$Time.to.swell)))
j6<-cbind(to_15$Time.to.swell,rep("takeover\n12-15 mos",length(to_15$Time.to.swell)))
j7<-cbind(to_18$Time.to.swell,rep("takeover\n15-18 mos",length(to_18$Time.to.swell)))
jd<-data.frame(rbind(j1,j2,j3,j4,j5,j6,j7))                                                
colnames(jd)<-c("Time.to.swell","Context")      
jd$Time.to.swell<-as.numeric(as.character(jd$Time.to.swell))
jd$Context<-factor(jd$Context,levels=c("takeover\n0-3 mos","takeover\n3-6 mos","takeover\n6-9 mos","takeover\n9-12 mos","takeover\n12-15 mos","takeover\n15-18 mos","no\ntakeover"))
table(jd$Context)

sum1<-ddply(jd, "Context", summarise, Time.to.swell=mean(Time.to.swell))
sum2<-ddply(jd, "Context",summarise, sd=sd(Time.to.swell, na.rm = TRUE))
summary<-merge(sum1,sum2)

col<-c(rep("coral3",6),"#999999")
ggplot(jd, aes(Context,Time.to.swell,colour=Context)) +
  geom_jitter(position = position_jitter(0.2),size=2) + 
  scale_colour_manual(values=col)+
  geom_pointrange(aes(ymin = Time.to.swell-sd, ymax = Time.to.swell+sd),data = summary,colour="black",size=1.3,alpha=0.8)+
  xlab("")+ ylab("Time to swelling (mos)")+
  theme(legend.position = "none",
        axis.title= element_text(size = 20),
        axis.text.x = element_text(angle=45,size = 18,color="black",vjust=0.7),
        axis.text.y = element_text(size = 20,color="black"))


