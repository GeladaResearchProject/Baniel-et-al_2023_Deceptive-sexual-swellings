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
theme_set (theme_classic(base_size=15))


### Table S5.COX PROPORTIONAL HAZARDS MODELS ###
cox<-read.csv("~/Dataset/dataset 5_Cox model_IBI.txt", header = TRUE, sep= ",",na.strings="NA")
head(cox)
dim(cox)#3466
summary(cox)
length(unique(cox$Infant))
length(unique(cox$Female))
length(unique(cox$Unit))
length(unique(cox$Infant[cox$Takeover.fix %in% 1]))#33 infants with takeover before 18 mos
cox$Birth<-as.integer(as.character(cox$Birth))
cox$Sex.infant<-as.factor(as.character(cox$Sex.infant))
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


### MODEL 3a - Takeover as fixed predictor ####
fit1 <- coxph(Surv(Time1, Time2, Birth) ~ Takeover.fix + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)

#Conservative dataset
table(cox$Dataset)
cox.cons<-droplevels(subset(cox,cox$Dataset %in% "conservative"))
dim(cox.cons)
fit2 <- coxph(Surv(Time1, Time2, Birth) ~ Takeover.fix + Parity + Sex.infant +log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit2)


### MODEL 3b - Takeover as time-varying predictor ####
fit1 <- coxph(Surv(Time1, Time2, Birth) ~ Takeover.timevarying + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox)
summary(fit1)
cox.zph(fit1)

#Conservative dataset
fit2 <- coxph(Surv(Time1, Time2, Birth) ~ Takeover.timevarying + Parity + Sex.infant + log.Avg.temperature + log.Cumulative.rainfall + cluster(Female), data=cox.cons)
summary(fit2)


### Plot Survival curves ###

#Takeover.fix
a <- survfit(Surv(Time1, Time2, Birth) ~ Takeover.fix, data = cox)
ggsurvplot(a,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 

#Takeoveover.timevarying 
b <- survfit(Surv(Time1, Time2, Birth) ~ Takeover.timevarying, data = cox)
ggsurvplot(b,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 

#Parity
c <- survfit(Surv(Time1, Time2, Birth) ~ Parity, data = cox)
ggsurvplot(c,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 

#Sex infant
d <- survfit(Surv(Time1, Time2, Birth) ~ Sex.infant, data = cox)
ggsurvplot(d,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 

#Seasonal effects
e <- survfit(Surv(Time1, Time2, Birth) ~ Season, data = cox)
ggsurvplot(e,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 

f<- survfit(Surv(Time1, Time2, Birth) ~ Rainfall.cat, data = cox)
ggsurvplot(f,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 

g<- survfit(Surv(Time1, Time2, Birth) ~ Temperature.cat, data = cox)
ggsurvplot(g,conf.int = TRUE)+ ylab("Proportion of mos without birth") + xlab("IBI (mos)") 




### Table S6.LMM ANALYSIS ###
source("~/code/vif.R") 
df<- read.csv("~/Dataset/dataset 6_LMM_IBI.txt", header = TRUE, sep= ",",na.strings="NA")
head(df)
dim(df)#119
summary(df)
df$Takeover.before18mos<-factor(df$Takeover.before18mos,levels=c("No takeover","Takeover"))
length(unique(df$Unit))
length(unique(df$Female))
length(unique(df$Infant))
table(df$Takeover.before18mos)#86 non-takeover infants and 33 takeover infants

#When infants experienced the takeover before birth, set age at takeover at zero month
sort(df$Age.takeover)
df$Age.takeover[df$Age.takeover<0 & df$Age.takeover %!in% NA]<-0

#IBI for non-takeover and takeover infants
summary(df$IBI[df$Takeover %in% "No takeover"])
summary(df$IBI[df$Takeover %in% "Takeover"])


### MODEL 3c - ALL INFANTS ###
m1<-lmer(IBI ~  Takeover.before18mos + Parity + Sex.infant + (1|Female) + (1|Unit), data=df)
summary(m1)
vif.mer(m1)
hist(resid(m1),nclass=100)
plot(resid(m1)~fitted(m1),nclass=100)
qqnorm(resid(m1))
qqline(resid(m1))

#Conservative dataset 
table(df$Dataset)
df.cons<-droplevels(subset(df,df$Dataset %in% "conservative"))
dim(df.cons)#114
table(df.cons$Takeover.before18mos)
m1.cons<-lmer(IBI ~  Takeover.before18mos + Parity + Sex.infant + (1|Female) + (1|Unit), data=df.cons)
summary(m1.cons)

#Plot
a<-ggplot(df, aes(x=IBI, color=Takeover.before18mos))+ geom_histogram(binwidth=0.5) +xlab("IBI  (mos)")+ scale_color_manual(values=c("#999999","coral3"))
b<-ggplot(df,aes(y=IBI, x=Takeover.before18mos,fill=Takeover.before18mos))+ geom_boxplot()+ geom_jitter()+ylab("IBI  (mos)")+xlab ("")+ scale_fill_manual(values=c("#999999","coral3"))
c<-ggplot(df,aes(y=IBI, x=Parity))+ geom_boxplot() + geom_jitter()+ ylab("IBI (mos)")+xlab ("Maternal parity")
d<-ggplot(df,aes(y=IBI, x=Sex.infant))+ geom_boxplot()+ geom_jitter()+ ylab("IBI  (mos)")+xlab ("Infant sex")
ggarrange(a,b,c,d,nrow=2,ncol=2)



### MODEL 3d - ONLY TAKEOVER INFANTS  ###
df2<-droplevels(subset(df,df$Takeover.before18mos %in% "Takeover"))
dim(df2)#33
m2<-lmer(IBI ~  scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit), data=df2)
summary(m2)

#Conservative dataset 
df2.cons<-droplevels(subset(df2,df2$Dataset %in% "conservative"))
dim(df2.cons)
m2.cons<-lmer(IBI ~  scale(Age.takeover) + Parity + Sex.infant + (1|Female) + (1|Unit), data=df2.cons)
summary(m2.cons)

#Plot age at takeover (continuous)
ggplot(df2,aes(y=IBI, x=Age.takeover))+ geom_point(size=2) + 
  ylab("IBI (mos)")+ xlab ("Age of infant at takeover (mos)")+
  geom_smooth(method='lm', formula= y~x)e


### MODEL 3e - ALL INFANTS ###
df$Takeover.cat<-rep(NA,length(df$Female))
df$Takeover.cat[df$Takeover.before18mos %in% "No takeover"]<-"No takeover"
df$Takeover.cat[df$Age.takeover <= 6]<-"takeover 0-6mos"
df$Takeover.cat[df$Age.takeover > 6 & df$Age.takeover <= 12]<-"takeover 6-12mos"
df$Takeover.cat[df$Age.takeover > 12]<-"takeover 12-18mos"
table(df$Takeover.cat)
df$Takeover.cat<-as.factor(df$Takeover.cat)

m1<-lmer(IBI ~  Takeover.cat + Parity + Sex.infant + (1|Female) + (1|Unit), data=df)
summary(m1)

#Plot age at takeover (categorical)
no.to<-droplevels(subset(df,df$Takeover.before18mos %in% "No takeover"));dim(no.to) 
to<-droplevels(subset(df,df$Takeover.before18mos %in% "Takeover"));dim(to) 
sort(to$Age.takeover)
to_6<-droplevels(subset(to,to$Age.takeover <= 6));dim(to_6)
to_12<-droplevels(subset(to,to$Age.takeover > 6 & to$Age.takeover <= 12));dim(to_12)
to_18<-droplevels(subset(to,to$Age.takeover > 12 & to$Age.takeover <= 18));dim(to_18)

j1<-cbind(no.to$IBI,rep("no\ntakeover",length(no.to$IBI)))
j2<-cbind(to_6$IBI,rep("takeover\n0-6 mos",length(to_6$IBI)))
j3<-cbind(to_12$IBI,rep("takeover\n6-12 mos",length(to_12$IBI)))
j4<-cbind(to_18$IBI,rep("takeover\n12-18 mos",length(to_18$IBI)))
jd<-data.frame(rbind(j1,j2,j3,j4))                                                
colnames(jd)<-c("IBI","Context")      
jd$IBI<-as.numeric(as.character(jd$IBI))
jd$Context<-factor(jd$Context,levels=c("takeover\n0-6 mos","takeover\n6-12 mos","takeover\n12-18 mos","no\ntakeover"))

sum1<-ddply(jd, "Context", summarise, IBI = mean(IBI))
sum2<-ddply(jd, "Context",summarise, sd=sd(IBI, na.rm = TRUE))
summary<-merge(sum1,sum2)
col<-c(rep("coral3",3),"#999999")

ggplot(jd, aes(Context,IBI,colour=Context)) +
  geom_jitter(position = position_jitter(0.2),size=2) + 
  scale_colour_manual(values=col)+
  geom_pointrange(aes(ymin = IBI-sd, ymax = IBI+sd),data = summary,colour="black",size=1.3,alpha=0.8)+
  xlab("")+ ylab("IBI duration (mos)")+
  theme(legend.position = "none",
        axis.title= element_text(size = 20),
        axis.text.x = element_text(size = 19,color="black"),
        axis.text.y = element_text(size = 20,color="black"))

