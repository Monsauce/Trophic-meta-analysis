#R script for "Predators have greater impacts on recipient freswater communities" 
library(ggplot2)
library(metafor)
library(RCurl)
library(plyr)

###Coarse Data###
#read in coarse prey categories data scraped from literature from GitHub
coarse.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Meta.Coarse.csv")
coarse<-read.csv(text=coarse.URL)

#calculate effect size using Hedge's D-standardized mean difference (due to negative values)
meta.coarse<-escalc(measure="SMD", m1i=mean, m2i=control, sd1i=mean.SD, sd2i=control.SD,
                        n1i=mean.N, n2i=control.N,data=coarse, var.names=c("SMD","SMD_var"),digits=4)

#set omnivore as the reference level for the models
meta.coarse$trophic.level <- relevel(factor(meta.coarse$trophic.level), ref="omnivore")

#run a random effects model trophic level as a moderator 
Model.1 <- rma(SMD, SMD_var, data = meta.coarse, method="EB")
Model.1

#run a mixed effects model trophic level as a moderator 
Model.2 <- rma(SMD, SMD_var, mods = ~ trophic.level, data = meta.coarse, method="EB")
Model.2

#calculate amount of variance explained by moderator   
anova(Model.1, Model.2)#no interaction 


####Fine Data####
#read in fine prey categories data scraped from literature from GitHub
fine.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Meta.Fine.csv")
fine<-read.csv(text=fine.URL)

#calculate effect size using Hedge's D-standardized mean difference (due to negative values)
meta.fine<-escalc(measure="SMD", m1i=mean, m2i=control, sd1i=mean.SD, sd2i=control.SD,
                    n1i=mean.N, n2i=control.N,data=fine, var.names=c("SMD","SMD_var"),digits=4)

#set omnivore as the reference level for the models
meta.coarse$trophic.level <- relevel(factor(meta.coarse$trophic.level), ref="omnivore")

#run a random effects model 
Model.3 <- rma(SMD, SMD_var, data = meta.fine, method="EB")
Model.3

#run a mixed effects model trophic level as a moderator 
Model.4 <- rma(SMD, SMD_var, mods = ~ trophic.level, data = meta.fine, method="EB")
Model.4

#calculate amount of variance explained by moderator, test if moderator is significant    
anova(Model.3, Model.4)

#run a random effects model with the interaction of trophic level and prey type  
Model.5 <- rma(SMD, SMD_var, mods = ~ trophic.level*prey.type, data = meta.fine, method="REML")

#test if interaction of moderators is significant using likelihood ratio test 
Model.5.additive <- rma(SMD, SMD_var, mods = ~ trophic.level+prey.type, data = meta.fine, method="ML")
Model.5.interaction <- rma(SMD, SMD_var, mods = ~ trophic.level*prey.type, data = meta.fine, method="ML")
anova(Model.5.additive, Model.5.interaction)

####Figure 1####
#make table of effect size values for trophic.level 
Model.4.intercept <- rma(SMD, SMD_var, mods = ~ trophic.level-1, data = meta.fine, method="REML")#run model w/o intercept
SE<-summary(Model.4.intercept)$b
ci_l<-summary(Model.4.intercept)$ci.lb
ci_h<-summary(Model.4.intercept)$ci.ub

Figure1A.table<-data.frame(cbind(SE,ci_l,ci_h))
colnames(Figure1A.table)[1]<-"SE"
colnames(Figure1A.table)[2]<-"ci_l"
colnames(Figure1A.table)[3]<-"ci_h"
Figure1A.table$Trophic.level<-c("Omnivore","Predator")
Figure1A.table$Trophic.level<-as.factor(Figure1A.table$Trophic.level)
Figure1A.table<-ddply(.data=Figure1A.table, .variables=.(SE, Trophic.level), .fun= summarise, CI = abs(SE-ci_h))
Figure1A.table

#plot Fligure 1A
Figure.1A<-ggplot(Figure1A.table, aes(x = Trophic.level, y = SE))+
  geom_point()+geom_errorbar(aes(ymin=SE-CI, ymax=SE+CI, width=0.2))+
  theme_bw()+
  ylab("Effect size")+
  xlab("Trophic level")


#run a random effects model with the interaction of trophic level and species  
Model.6 <- rma(SMD, SMD_var, mods = ~ trophic.level:species - 1, data=meta.fine, method="REML")
Model.6


#make table of effect size values for 
SE<-summary(Model.)$b
ci_l<-summary(Model.)$ci.lb
ci_h<-summary(Model.)$ci.ub

Figure1B.table<-data.frame(cbind(SE,ci_l,ci_h))
colnames(Figure1B.table)[1]<-"SE"
colnames(Figure1B.table)[2]<-"ci_l"
colnames(Figure1B.table)[3]<-"ci_h"
Figure1B.table$Trophic.level<-c("Omnivore:Consumer", "Trophic:Predator", "Prey:Predator", "Prey:Consumer", "Predator:Predator")
Figure1B.table$Trophic.level<-as.factor(Figure1B.table$Trophic.level)
Figure1B.table<-ddply(.data=Figure1B.table, .variables=.(SE, Trophic.level), .fun= summarise, CI = abs(SE-ci_h))
Figure1B.table




