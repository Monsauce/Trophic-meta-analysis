#R script for "Predators have greater impacts on recipient freswater communities" 
library(ggplot2)
library(metafor)
library(RCurl)

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
anova(Model.1, Model.2)

#forest plot of coarse data  
forest(meta.coarse$SMD,meta.coarse$SMD_var,slab=meta.coarse$trophic,pch=19)

####Fine Data####
#read in fine prey categories data scraped from literature from GitHub
fine.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Fine.Coarse.csv")
fine<-read.csv(text=fine.URL)

#calculate effect size using Hedge's D-standardized mean difference (due to negative values)
meta.fine<-escalc(measure="SMD", m1i=mean, m2i=control, sd1i=mean.SD, sd2i=control.SD,
                    n1i=mean.N, n2i=control.N,data=fine, var.names=c("SMD","SMD_var"),digits=4)

#set omnivore as the reference level for the models
meta.coarse$trophic.level <- relevel(factor(meta.coarse$trophic.level), ref="omnivore")

#run a random effects model with the interaction of trophic level and prey type  
Model.3 <- rma(SMD, SMD_var, mods = ~ trophic.level*prey.type, data = meta.fine, method="REML")

#test if interaction of moderators is significant using likelihood ratio test 
Model.4.additive <- rma(SMD, SMD_var, mods = ~ trophic.level+prey.type, data = meta.coarse, method="ML")
Model.5.interaction <- rma(SMD, SMD_var, mods = ~ trophic.level*prey.type, data = meta.coarse, method="ML")
anova(Model.4.additive, Model.5.interaction)




