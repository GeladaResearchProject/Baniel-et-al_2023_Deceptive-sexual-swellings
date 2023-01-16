rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(lme4)
library(lmerTest)
theme_set (theme_classic(base_size=20))

### Table S7.LMM ANALYSIS ###
df0<- read.csv("~/Dataset/dataset 9_LMM_Estrogen.txt", header = TRUE, sep= ",",na.strings="NA")
head(df0)
dim(df0)#619
summary(df0)
df0$ID<-paste(df0$Female,df0$Infant,df0$ID.takeover,sep="_")
length(unique(df0$ID))#34

#Keep only samples 30 days post-takeover
sort(df0$Nb.days.relative.takeover)
df<-droplevels(subset(df0,df0$Nb.days.relative.takeover<=30))
dim(df)#308

#Only keep females for which we have samples before and after takeover (longitudinal dataset)
t1<-as.data.frame.matrix(table(df$ID,df$Takeover.period));t1
longitudinal<-t1[t1$before %!in% 0 & t1$after %!in% 0,]
df<-droplevels(subset(df,df$ID %in% rownames(longitudinal)))

#Final dataset
head(df)
dim(df)#289
length(unique(df$Female))#22
table(df$Repro.State)
length(unique(df$Infant))#27
df$Takeover.period<-factor(df$Takeover.period,level=c("before","after"))
table(df$Takeover.period)#193 samples before and 96 samples after takeover


### LMM ###
source("~/code/vif.R") 
m1<-lmer(log.Estrogen ~ Takeover.period  + scale(Age.sample) + scale(Avg.temperature) + scale(Cumulative.rainfall) + Wash.step + Laboratory + (1|Female) + (1|Unit),data=df, REML=F, control=lmerControl(optimizer="bobyqa"))
summary(m1)
vif.mer(m1)
hist(resid(m1),nclass=100)
plot(resid(m1)~fitted(m1),nclass=100)
qqnorm(resid(m1))
qqline(resid(m1))


### PLOT ###
#Estrogen & takeover period
ggplot(df,aes(y=Estrogen.ng.g, x=Takeover.period))+geom_boxplot(outlier.shape = NA)+ geom_jitter(size=2, aes(color=Repro.State)) 
ggplot(df,aes(y=log.Estrogen, x=Takeover.period))+geom_boxplot(outlier.shape = NA)+ geom_jitter(size=2, aes(color=Repro.State)) 

#Infant age at sample collection
sort(df$Age.sample)
ggplot(df,aes(y=Estrogen.ng.g, x=Age.sample,colour=Takeover.period))+ geom_point(size=2) + geom_smooth(method='lm', formula= y~x) + xlab("Infant age at sample collection (month)") + scale_colour_brewer(palette="Dark2")
ggplot(df,aes(y=log.Estrogen, x=Age.sample,colour=Takeover.period))+ geom_point(size=2) + geom_smooth(method='lm', formula= y~x) + xlab("Infant age at sample collection (month)") + scale_colour_brewer(palette="Dark2")

table(df$Age.sample.cat)
df$Age.sample.cat<-factor(df$Age.sample.cat,levels = c("0-6 months","6-12 months","12-18 months",">18 months"))
ggplot(df,aes(y=Estrogen.ng.g, x=Takeover.period))+geom_boxplot(outlier.shape = NA)+ geom_jitter(size=2, aes(color=Repro.State))+ facet_wrap(~Age.sample.cat)
ggplot(df,aes(y=log.Estrogen, x=Takeover.period))+geom_boxplot(outlier.shape = NA)+ geom_jitter(size=2, aes(color=Repro.State))+ facet_wrap(~Age.sample.cat)

#Infant age at takeover
sort(df$Age.takeover)
ggplot(df,aes(y=Estrogen.ng.g, x=Age.takeover,colour=Takeover.period))+ geom_point(size=2) + geom_smooth(method='lm', formula= y~x) + xlab("Infant age at takeover (month)") + scale_colour_brewer(palette="Dark2")
ggplot(df,aes(y=log.Estrogen, x=Age.takeover,colour=Takeover.period))+ geom_point(size=2) + geom_smooth(method='lm', formula= y~x) + xlab("Infant age at takeover (month)") + scale_colour_brewer(palette="Dark2")

table(df$Age.takeover.cat)
df$Age.takeover.cat<-factor(df$Age.takeover.cat,levels = c("0-6 months","6-12 months","12-18 months",">18 months"))
ggplot(df,aes(y=Estrogen.ng.g, x=Takeover.period))+geom_boxplot(outlier.shape = NA)+ geom_jitter(size=2, aes(color=Repro.State))+ facet_wrap(~Age.takeover.cat)
ggplot(df,aes(y=log.Estrogen, x=Takeover.period))+geom_boxplot(outlier.shape = NA)+ geom_jitter(size=2, aes(color=Repro.State))+ facet_wrap(~Age.takeover.cat)


##Individual trajectories (max=30 days post-takeover)
ggplot(df,aes(y=Estrogen.ng.g, x=Nb.days.relative.takeover,colour=Repro.State))+xlim(-200,30)+
  geom_point(size=2) + geom_line(color="grey") + facet_wrap(~ID)+ 
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.5)+
  theme(legend.position="top")+ xlab("Days relative to male takeover")+ ylab("Fecal estrogens (ng/g)")


##Individual trajectories
ggplot(df0,aes(y=Estrogen.ng.g, x=Nb.days.relative.takeover,colour=Repro.State))+xlim(-400,400)+
  geom_point(size=2) + geom_line(color="grey") + facet_wrap(~ID)+ 
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.5)+
  theme(legend.position="top")+ xlab("Days relative to male takeover")+ ylab("Fecal estrogens (ng/g)")



