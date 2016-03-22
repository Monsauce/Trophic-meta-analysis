#R script for "Predators have greater impacts on recipient freswater communities" 
library(ggplot2)
library(metafor)
library(RCurl)
library(plyr)
library(ICC)
library(frair)
library(tidyr)
library(emdbook)
library(gridExtra)
library(nlstools)

####Community Data###
#read in community prey categories data scraped from literature from GitHub
coarse.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Meta.Community.csv")
coarse<-read.csv(text=coarse.URL)

#calculate number of unique studies 
length(unique(coarse$study))

#calculate effect size using Hedge's D-standardized mean difference (due to negative values)
meta.coarse<-escalc(measure="SMD", m1i=mean, m2i=control, sd1i=mean.SD, sd2i=control.SD,
                    n1i=mean.N, n2i=control.N,data=coarse, var.names=c("SMD","SMD_var"),digits=4)

#calculate intraclass corelation coefficient to assess dependence between studies 
ICC.coarse.studies<-ICCbare(study, SMD, data=meta.coarse)

#set omnivore as the reference level for the models
meta.coarse$trophic.level <- relevel(factor(meta.coarse$trophic.level), ref="omnivore")

#run a random effects model 
Model.1 <- rma(SMD, SMD_var, data = meta.coarse, method="REML")
Model.1

#run a mixed effects model trophic level as a moderator 
Model.2 <- rma(SMD, SMD_var, mods = ~ trophic.level, data = meta.coarse, method="REML")
Model.2

#Q-Q plot to test assumptions of data 
Q.plot.Model.2<-qqnorm(Model.2, main="Mixed-Effects Model")


####Population Data####
#read in population prey categories data scraped from literature from GitHub
fine.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Meta.Population.csv")
fine<-read.csv(text=fine.URL)

#calculate number of unique studies 
length(unique(fine$study))

#calculate effect size using Hedge's D-standardized mean difference (due to negative values)
meta.fine<-escalc(measure="SMD", m1i=mean, m2i=control, sd1i=mean.SD, sd2i=control.SD,
                  n1i=mean.N, n2i=control.N,data=fine, var.names=c("SMD","SMD_var"),digits=4)

#calculate intraclass corelation coefficient to assess dependence between studies and species 
ICC.fine.studies<-ICCbare(study, SMD, data=meta.fine)
ICC.fine.species<-ICCbare(species, SMD, data=meta.fine)

#set omnivore as the reference level for the models
meta.coarse$trophic.level <- relevel(factor(meta.coarse$trophic.level), ref="omnivore")

#run a random effects model 
Model.3 <- rma(SMD, SMD_var, data = meta.fine, method="REML")
Model.3

#run a mixed effects model trophic level as a moderator 
Model.4 <- rma(SMD, SMD_var, mods = ~ trophic.level, data = meta.fine, method="REML")
Model.4

#Q-Q plot to test assumptions of data 
Q.plot<-qqnorm(Model.4, main="Mixed-Effects Model")

#test if interaction of moderators is significant using likelihood ratio test 
Model.5.additive <- rma(SMD, SMD_var, mods = ~ trophic.level+prey.type, data = meta.fine, method="ML")
Model.5.interaction <- rma(SMD, SMD_var, mods = ~ trophic.level*prey.type, data = meta.fine, method="ML")
anova(Model.5.additive, Model.5.interaction)

#test if species account for heterogenetity 
Model.6<-rma(SMD, SMD_var, mods = ~ species, data = meta.fine, method="REML")
Model.6

####Figure 1####
Figure1.model <- rma(SMD, SMD_var, mods = ~ trophic.level-1, data = meta.fine, method="REML")#run model w/o intercept with prey.type interaction
ES<-summary(Figure1.model)$b
ci_l<-summary(Figure1.model)$ci.lb
ci_h<-summary(Figure1.model)$ci.ub

Figure1.table<-data.frame(cbind(ES,ci_l,ci_h))
Figure1.table<-data.frame(cbind(ES,ci_l,ci_h))
colnames(Figure1.table)[1]<-"ES"
colnames(Figure1.table)[2]<-"CI_l"
colnames(Figure1.table)[3]<-"CI_h"
Figure1.table$trophic.level<-c("omnivore", "predator")
daply(meta.fine, .(trophic.level), nrow)
Figure1.table$ES.count<-c(81, 75)
Figure1.table$trophic.level<-as.factor(Figure1.table$trophic.level)

#plot Figure 1
Figure.1<-ggplot(data=Figure1.table, aes(x=trophic.level, y=ES))+
  geom_errorbar(aes(ymin=CI_l, ymax=CI_h), width=.200)+
  geom_point(aes(size=ES.count, fill=trophic.level), colour="black",shape=21)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(limits=c(-1,1), name="SMD")+
  scale_x_discrete(name="Trophic level")+
  coord_flip()+
  scale_fill_manual(values=c("black", "white"))+
  theme_bw()

#make table of effect size values for trophic.level and prey type
Figure2.model <- rma(SMD, SMD_var, mods = ~ trophic.level:prey.type-1, data = meta.fine, method="REML")#run model w/o intercept with prey.type interaction
ES<-summary(Figure2.model)$b
ci_l<-summary(Figure2.model)$ci.lb
ci_h<-summary(Figure2.model)$ci.ub

Figure2.table<-data.frame(cbind(ES,ci_l,ci_h))
Figure2.table<-data.frame(cbind(ES,ci_l,ci_h))
colnames(Figure2.table)[1]<-"ES"
colnames(Figure2.table)[2]<-"CI_l"
colnames(Figure2.table)[3]<-"CI_h"
Figure2.table$trophic.level<-c("omnivore", "predator", "omnivore","predator")
Figure2.table$prey.type<-c("omnivore:consumer", "predator:consumer", "omnivore:herbivore","predator:herbivore")
daply(meta.fine, .(trophic.level,prey.type), nrow)
Figure2.table$ES.count<-c(64,46,17,29)
Figure2.table$trophic.level<-as.factor(Figure2.table$trophic.level)
Figure2.table$species<-as.factor(Figure2.table$prey.type)
#reorder factors
Figure2.table<- transform(Figure2.table, prey.type = reorder(prey.type, order(prey.type, decreasing = TRUE)))

#plot Figure 2
Figure.2<-ggplot(data=Figure2.table, aes(x=prey.type, y=ES))+
  geom_errorbar(aes(ymin=CI_l, ymax=CI_h), width=.200)+
  geom_point(aes(size=ES.count, fill=trophic.level), colour="black",shape=21)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(limits=c(-5,5), name="SMD")+
  scale_x_discrete(name="Trophic level")+
  coord_flip()+
  scale_fill_manual(values=c("black", "white"))+
  theme_bw()

####Figure A1####
#make table of effect size values for trophic.level and species 
FigureA1.model <- rma(SMD, SMD_var, mods = ~ species-1, data = meta.fine, method="REML")#run model w/o intercept 
ES<-summary(FigureA1.model)$b
ci_l<-summary(FigureA1.model)$ci.lb
ci_h<-summary(FigureA1.model)$ci.ub

FigureA1.table<-data.frame(cbind(ES,ci_l,ci_h))
colnames(FigureA1.table)[1]<-"ES"
colnames(FigureA1.table)[2]<-"CI_l"
colnames(FigureA1.table)[3]<-"CI_h"
FigureA1.table$trophic.level<-c("predator", "predator", "predator","omnivore","omnivore","omnivore","omnivore","omnivore","predator","omnivore","omnivore","omnivore","omnivore","omnivore","predator")
FigureA1.table$species<-c("Bythotrephes cederstroemi", "Bythotrephes longimanus", "Cercopagis pengoi", "Dikerogammarus villosus", "Gammarus fasciatus", "Gammarus tigrinus", "Hemimysis anomala", "Limnomysis benedeni", "Leptodora kindtii", "Mysis diluviana", "Mysis relicta", "Orconectes rusticus", "Procambarus clarkii", "Pacifastacus leniusculus", "Tortanus dextrilobatus")
FigureA1.table$ES.count<-daply(meta.fine, .(species), nrow)
FigureA1.table$trophic.level<-as.factor(FigureA1.table$trophic.level)
FigureA1.table$species<-as.factor(FigureA1.table$species)
#reorder factors
FigureA1.table<- transform(FigureA1.table, species = reorder(species, order(species, decreasing = TRUE)))

#plot Figure A1
Figure.A1<-ggplot(data=FigureA1.table, aes(x=species, y=ES))+
  geom_errorbar(aes(ymin=CI_l, ymax=CI_h), width=.667)+
  geom_point(aes(size=ES.count), colour="black",shape=21)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(limits=c(-5,5), name="Standardized Mean Difference (SMD)")+
  scale_x_discrete(name="Species")+
  coord_flip()+
  theme_bw()


####Functional response####
#read in functional response data from GitHub
FR.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/FR.csv")
FR<-read.csv(text=FR.URL)

#run fit function using non-linear least for each predator-prey combination 
#subset study 1.1
FR.1.1<-FR[FR$study%in%c("1.1"),]
#check for Type II
FR.1.1.prop<- ddply(FR.1.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.1.1.lm<-lm(prop.consumed~prey.density, data=FR.1.1.prop)
summary(FR.1.1.lm)

FR.rogers.1.1<-frair_fit(prey.consumed~prey.density, data=FR.1.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 1.2
FR.1.2<-FR[FR$study%in%c("1.2"),]
#check for Type II
FR.1.2.prop<- ddply(FR.1.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.1.2.lm<-lm(prop.consumed~prey.density, data=FR.1.2.prop)
summary(FR.1.2.lm)

FR.rogers.1.2<-frair_fit(prey.consumed~prey.density, data=FR.1.2, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.1
FR.2.1<-FR[FR$study%in%c("2.1"),]
#check for Type II
FR.2.1.prop<- ddply(FR.2.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.1.lm<-lm(prop.consumed~prey.density, data=FR.2.1.prop)
summary(FR.2.1.lm)

FR.rogers.2.1<-frair_fit(prey.consumed~prey.density, data=FR.2.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.2
FR.2.2<-FR[FR$study%in%c("2.2"),]
#check for Type II
FR.2.2.prop<- ddply(FR.2.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.2.lm<-lm(prop.consumed~prey.density, data=FR.2.2.prop)
summary(FR.2.2.lm)

FR.rogers.2.2<-frair_fit(prey.consumed~prey.density, data=FR.2.2, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.3
FR.2.3<-FR[FR$study%in%c("2.3"),]
#check for Type II
FR.2.3.prop<- ddply(FR.2.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.3.lm<-lm(prop.consumed~prey.density, data=FR.2.3.prop)
summary(FR.2.3.lm)

FR.rogers.2.3<-frair_fit(prey.consumed~prey.density, data=FR.2.3, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.4
FR.2.4<-FR[FR$study%in%c("2.4"),]
#check for Type II
FR.2.4.prop<- ddply(FR.2.4, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.4.lm<-lm(prop.consumed~prey.density, data=FR.2.4.prop)
summary(FR.2.4.lm)

FR.rogers.2.4<-frair_fit(prey.consumed~prey.density, data=FR.2.4, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.5
FR.2.5<-FR[FR$study%in%c("2.5"),]
#check for Type II
FR.2.5.prop<- ddply(FR.2.5, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.5.lm<-lm(prop.consumed~prey.density, data=FR.2.5.prop)
summary(FR.2.5.lm)

FR.rogers.2.5<-frair_fit(prey.consumed~prey.density, data=FR.2.5, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.6
FR.2.6<-FR[FR$study%in%c("2.6"),]
#check for Type II
FR.2.6.prop<- ddply(FR.2.6, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.6.lm<-lm(prop.consumed~prey.density, data=FR.2.6.prop)
summary(FR.2.6.lm)

FR.rogers.2.6<-frair_fit(prey.consumed~prey.density, data=FR.2.6, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.7
FR.2.7<-FR[FR$study%in%c("2.7"),]
#check for Type II
FR.2.7.prop<- ddply(FR.2.7, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.7.lm<-lm(prop.consumed~prey.density, data=FR.2.7.prop)
summary(FR.2.7.lm)

FR.rogers.2.7<-frair_fit(prey.consumed~prey.density, data=FR.2.7, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 2.8
FR.2.8<-FR[FR$study%in%c("2.8"),]
#check for Type II
FR.2.8.prop<- ddply(FR.2.8, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.2.8.lm<-lm(prop.consumed~prey.density, data=FR.2.8.prop)
summary(FR.2.8.lm)

FR.rogers.2.8<-frair_fit(prey.consumed~prey.density, data=FR.2.8, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 3.1
FR.3.1<-FR[FR$study%in%c("3.1"),]
#check for Type II
FR.3.1.prop<- ddply(FR.3.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.3.1.lm<-lm(prop.consumed~prey.density, data=FR.3.1.prop)
summary(FR.3.1.lm)

FR.rogers.3.1<-frair_fit(prey.consumed~prey.density, data=FR.3.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.67))

#subset study 3.2
FR.3.2<-FR[FR$study%in%c("3.2"),]
#check for Type II
FR.3.2.prop<- ddply(FR.3.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.3.2.lm<-lm(prop.consumed~prey.density, data=FR.3.2.prop)
summary(FR.3.2.lm)

FR.rogers.3.2<-frair_fit(prey.consumed~prey.density, data=FR.3.2, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.67))

#subset study 3.3
FR.3.3<-FR[FR$study%in%c("3.3"),]
#check for Type II
FR.3.3.prop<- ddply(FR.3.3, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.3.3.lm<-lm(prop.consumed~prey.density, data=FR.3.3.prop)
summary(FR.3.3.lm)

rogers.pred = function(prey.density.stand, a, h, P, T) {prey.density.stand - lambertW(a * h * prey.density.stand * exp(-a * (P * T -h * prey.density.stand)))/(a * h)}
params<- list(a=6,h=0.01)

FR.rogers.3.3<-nls(prey.consumed~rogers.pred(prey.density,a,h,P=1,T=0.67),
                   start=params,data=FR.3.3,control=list(maxiter=1000))

#subset study 3.4
FR.3.4<-FR[FR$study%in%c("3.4"),]
#check for Type II
FR.3.4.prop<- ddply(FR.3.4, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.3.4.lm<-lm(prop.consumed~prey.density, data=FR.3.4.prop)
summary(FR.3.4.lm)

FR.rogers.3.4<-nls(prey.consumed~rogers.pred(prey.density,a,h,P=1,T=0.67),
                   start=params,data=FR.3.4,control=list(maxiter=1000))

#subset study 4.1
FR.4.1<-FR[FR$study%in%c("4.1"),]
#check for Type II
FR.4.1.prop<- ddply(FR.4.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.4.1.lm<-lm(prop.consumed~prey.density, data=FR.4.1.prop)
summary(FR.4.1.lm)

FR.rogers.4.1<-frair_fit(prey.consumed~prey.density, data=FR.4.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=1))

#subset study 4.2
FR.4.2<-FR[FR$study%in%c("4.2"),]
#check for Type II
FR.4.2.prop<- ddply(FR.4.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.4.2.lm<-lm(prop.consumed~prey.density, data=FR.4.2.prop)
summary(FR.4.2.lm)

FR.rogers.4.2<-frair_fit(prey.consumed~prey.density, data=FR.4.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=1))

#subset study 4.3
FR.4.3<-FR[FR$study%in%c("4.3"),]
#check for Type II
FR.4.3.prop<- ddply(FR.4.3, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.4.3.lm<-lm(prop.consumed~prey.density, data=FR.4.3.prop)
summary(FR.4.3.lm)

FR.rogers.4.3<-frair_fit(prey.consumed~prey.density, data=FR.4.3, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=1))

#subset study 5.1
FR.5.1<-FR[FR$study%in%c("5.1"),]
#check for Type II
FR.5.1.prop<- ddply(FR.5.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.1.lm<-lm(prop.consumed~prey.density, data=FR.5.1.prop)
summary(FR.5.1.lm)
#non-negative slope

FR.5.1.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.1, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.2
FR.5.2<-FR[FR$study%in%c("5.2"),]
#check for Type II
FR.5.2.prop<- ddply(FR.5.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.2.lm<-lm(prop.consumed~prey.density, data=FR.5.2.prop)
summary(FR.5.2.lm)
#non-negative slope

FR.5.2.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.2, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.3
FR.5.3<-FR[FR$study%in%c("5.3"),]
#check for Type II
FR.5.3.prop<- ddply(FR.5.3, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.3.lm<-lm(prop.consumed~prey.density, data=FR.5.3.prop)
summary(FR.5.3.lm)
#non-negative slope

FR.5.3.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.3, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.4
FR.5.4<-FR[FR$study%in%c("5.4"),]
#check for Type II
FR.5.4.prop<- ddply(FR.5.4, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.4.lm<-lm(prop.consumed~prey.density, data=FR.5.4.prop)
summary(FR.5.4.lm)
#non-negative slope

FR.5.4.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.4, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.5
FR.5.5<-FR[FR$study%in%c("5.5"),]
#check for Type II
FR.5.5.prop<- ddply(FR.5.5, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.5.lm<-lm(prop.consumed~prey.density, data=FR.5.5.prop)
summary(FR.5.5.lm)
#non-negative slope

FR.5.5.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.5, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.6
FR.5.6<-FR[FR$study%in%c("5.6"),]
#check for Type II
FR.5.6.prop<- ddply(FR.5.6, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.6.lm<-lm(prop.consumed~prey.density, data=FR.5.6.prop)
summary(FR.5.6.lm)
#non-negative slope

FR.5.6.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.6, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.7
FR.5.7<-FR[FR$study%in%c("5.7"),]
FR.5.7.prop<- ddply(FR.5.7, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.7.lm<-lm(prop.consumed~prey.density, data=FR.5.7.prop)
summary(FR.5.7.lm)
#non-negative slope

FR.5.7.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.7, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.8
FR.5.8<-FR[FR$study%in%c("5.8"),]
FR.5.8.prop<- ddply(FR.5.8, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.8.lm<-lm(prop.consumed~prey.density, data=FR.5.8.prop)
summary(FR.5.8.lm)
#non-negative slope

FR.5.8.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.8, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.9
FR.5.9<-FR[FR$study%in%c("5.9"),]
FR.5.9.prop<- ddply(FR.5.9, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.9.lm<-lm(prop.consumed~prey.density, data=FR.5.9.prop)
summary(FR.5.9.lm)
#non-negative slope

FR.5.9.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.9, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 5.10
FR.5.10<-FR[FR$study%in%c("5.10."),]
FR.5.10.prop<- ddply(FR.5.10, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.5.10.lm<-lm(prop.consumed~prey.density, data=FR.5.10.prop)
summary(FR.5.10.lm)
#non-negative slope

FR.5.10.Type1<-frair_fit(prey.consumed~prey.density, data=FR.5.10, response="typeI", start=list(a=0.2), fixed=list(T=0.5))

#subset study 6.1
FR.6.1<-FR[FR$study%in%c("6.1"),]
FR.6.1.prop<- ddply(FR.6.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.6.1.lm<-lm(prop.consumed~prey.density, data=FR.6.1.prop)
summary(FR.6.1.lm)

FR.rogers.6.1<-frair_fit(prey.consumed~prey.density, data=FR.6.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=1))

#subset study 6.2
FR.6.2<-FR[FR$study%in%c("6.2"),]
FR.6.2.prop<- ddply(FR.6.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.6.2.lm<-lm(prop.consumed~prey.density, data=FR.6.2.prop)
summary(FR.6.2.lm)
#non-negative slope

FR.6.2.Type1<-frair_fit(prey.consumed~prey.density, data=FR.6.2, response="typeI", start=list(a=0.2), fixed=list(T=1))

#subset study 6.3
FR.6.3<-FR[FR$study%in%c("6.3"),]
FR.6.3.prop<- ddply(FR.6.3, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.6.3.lm<-lm(prop.consumed~prey.density, data=FR.6.3.prop)
summary(FR.6.3.lm)
#non-negative slope

FR.6.3.Type1<-frair_fit(prey.consumed~prey.density, data=FR.6.3, response="typeI", start=list(a=0.2), fixed=list(T=1))

#subset study 7.1
FR.7.1<-FR[FR$study%in%c("7.1"),]
FR.7.1.prop<- ddply(FR.7.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.7.1.lm<-lm(prop.consumed~prey.density, data=FR.7.1.prop)
summary(FR.7.1.lm)

FR.rogers.7.1<-frair_fit(prey.consumed~prey.density, data=FR.7.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=1))

#subset study 8.1
FR.8.1<-FR[FR$study%in%c("8.1"),]
FR.8.1.prop<- ddply(FR.8.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.8.1.lm<-lm(prop.consumed~prey.density, data=FR.8.1.prop)
summary(FR.8.1.lm)

FR.rogers.8.1<-frair_fit(prey.consumed~prey.density, data=FR.8.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 8.2
FR.8.2<-FR[FR$study%in%c("8.2"),]
FR.8.2.prop<- ddply(FR.8.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.8.2.lm<-lm(prop.consumed~prey.density, data=FR.8.2.prop)
summary(FR.8.2.lm)

FR.rogers.8.2<-frair_fit(prey.consumed~prey.density, data=FR.8.2, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 8.3
FR.8.3<-FR[FR$study%in%c("8.3"),]
FR.8.3.prop<- ddply(FR.8.3, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.8.3.lm<-lm(prop.consumed~prey.density, data=FR.8.3.prop)
summary(FR.8.3.lm)

FR.rogers.8.3<-frair_fit(prey.consumed~prey.density, data=FR.8.3, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 9.1
FR.9.1<-FR[FR$study%in%c("9.1"),]
FR.9.1.prop<- ddply(FR.9.1, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.9.1.lm<-lm(prop.consumed~prey.density, data=FR.9.1.prop)
summary(FR.9.1.lm)

FR.rogers.9.1<-frair_fit(prey.consumed~prey.density, data=FR.9.1, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))

#subset study 9.2
FR.9.2<-FR[FR$study%in%c("9.2"),]
FR.9.2.prop<- ddply(FR.9.2, .(species), transform, prop.consumed=prey.consumed/prey.density)
FR.9.2.lm<-lm(prop.consumed~prey.density, data=FR.9.2.prop)
summary(FR.9.2.lm)

FR.rogers.9.2<-frair_fit(prey.consumed~prey.density, data=FR.9.2, response="rogersII", start=list(a=0.2,h=0.05), fixed=list(T=0.5))


####Figure 3####
#make list of models type I and II models  
models.I<-list(FR.5.1.Type1=FR.5.1.Type1,FR.5.2.Type1=FR.5.2.Type1, FR.5.3.Type1=FR.5.3.Type1,
               FR.5.4.Type1=FR.5.4.Type1, FR.5.5.Type1=FR.5.5.Type1, FR.5.6.Type1=FR.5.5.Type1, FR.5.7.Type1=FR.5.7.Type1,
               FR.5.8.Type1=FR.5.8.Type1,FR.5.9.Type1=FR.5.9.Type1, FR.5.10.Type1=FR.5.10.Type1,FR.6.2.Type1=FR.6.2.Type1,
               FR.6.3.Type1=FR.6.3.Type1)

models.II<-list(FR.rogers.1.1=FR.rogers.1.1,FR.rogers.1.2=FR.rogers.1.2, FR.rogers.2.1=FR.rogers.2.1, FR.rogers.2.2=FR.rogers.2.2,
                FR.rogers.2.3=FR.rogers.2.3,FR.rogers.2.4=FR.rogers.2.4,FR.rogers.2.5=FR.rogers.2.5,FR.rogers.2.6=FR.rogers.2.6,
                FR.rogers.2.7=FR.rogers.2.7,FR.rogers.2.8=FR.rogers.2.8, FR.rogers.3.1=FR.rogers.3.1, FR.rogers.3.2=FR.rogers.3.2, 
                FR.rogers.4.1=FR.rogers.4.1,FR.rogers.4.2=FR.rogers.4.2,FR.rogers.4.3=FR.rogers.4.3, FR.rogers.6.1=FR.rogers.6.1, 
                FR.rogers.7.1=FR.rogers.7.1, FR.rogers.8.1=FR.rogers.8.1,FR.rogers.8.2=FR.rogers.8.2,FR.rogers.8.3=FR.rogers.8.3,
                FR.rogers.9.1=FR.rogers.9.1,FR.rogers.9.2=FR.rogers.9.2)

#extract type II coefficients 
Models.I<-as.data.frame(lapply(models.I, coefficients))
#go from wide to long 
FR.models.I<-gather(Models.I[1,], model, a, FR.5.1.Type1:FR.6.3.Type1)
#add model fit 
FR.models.I$fit<-"Type I"
FR.models.I$h<-NA

#extract type II coefficients 
Models.II<-as.data.frame(lapply(models.II, coefficients))
#go from wide to long 
FR.models.II.a<-gather(Models.II[1,], model, a, FR.rogers.1.1:FR.rogers.9.2)
FR.models.II.h<-gather(Models.II[2,], model, h, FR.rogers.1.1:FR.rogers.9.2)
FR.models.II<-merge(FR.models.II.a,FR.models.II.h)
FR.models.II$fit<-"Type II"

#add 3.3 and 3.4
models.II.add<-list(FR.rogers.3.3=FR.rogers.3.3, FR.rogers.3.4=FR.rogers.3.4)
Models.II.add<-as.data.frame(lapply(models.II.add, coefficients))
FR.models.II.add.a<-gather(Models.II.add[1,], model, a, FR.rogers.3.3:FR.rogers.3.4)
FR.models.II.add.h<-gather(Models.II.add[2,], model, h, FR.rogers.3.3:FR.rogers.3.4)
FR.models.II.add<-merge(FR.models.II.add.a,FR.models.II.add.h)
FR.models.II.add$fit<-"Type II"

#merge all DFs
FR.models<-rbind(FR.models.I, FR.models.II, FR.models.II.add)

#add experiment information from GitHub
EXP.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Experiment.Data.csv")
EXP<-read.csv(text=EXP.URL)

#merge experiment info with functional response 
FR.table<- merge(FR.models,EXP,by="model")

#standardize attack rate and handling time 
FR.table<- ddply(FR.table, .(model), transform, a.standard=a/(volume*predator.size), h.standard=h/ratio)

#add effect size from Figure.2.Model
Figure2.table.subset<-Figure2.table[Figure2.table$species%in%c("Tortanus dextrilobatus", "Hemimysis anomala","Mysis diluviana","Dikerogammarus villosus","Leptodora kindtii","Cercopagis pengoi","Pacifastacus leniusculus","Bythotrephes longimanus"),]
Figure3.table<- merge(FR.table,Figure2.table.subset,by="species")

#calculate mean of a across species  
Figure3A.table.summary<- ddply(Figure3.table, .(species, ES,trophic.level), summarize, a.mean=mean(a.standard), a.SE=sd(a.standard)/sqrt(length(a.standard)))

#plot Figure 3A - attack rate 
Figure.3A<-ggplot(Figure3A.table.summary, aes(x=a.mean,y=ES, colour= trophic.level))+geom_point()+
  xlab("Attack rate (a) L-1 mm-1")+ylab("Effect size (SMD)")+
  scale_colour_manual(values=c("black", "dark grey"))+
  geom_errorbarh(aes(xmin=a.mean-a.SE, xmax=a.mean+a.SE),width=.2)+
  stat_smooth(colour="black", method=lm, se=F, aes(group=1))+
  theme_bw()

#plot Figure 3AB - handling time 
#subset Type II data
Figure3B.table.<-Figure3.table[Figure3.table$fit%in%c("Type II"),]

#calculate mean of h across species 
Figure3B.table.summary<- ddply(Figure3B.table, .(species, ES,trophic.level), summarize, h.mean=mean(h.standard), h.SE=sd(h)/sqrt(length(h.standard)))

Figure.3B<-ggplot(Figure3B.table.summary, aes(x=h.mean,y=ES, colour= trophic.level))+geom_point()+
  xlab("Standardized handing time (h)")+ylab("Effect size (SMD)")+
  scale_colour_manual(values=c("black", "dark grey"))+
  geom_errorbarh(aes(xmin=h.mean-h.SE, xmax=h.mean+h.SE),width=.2)+
  stat_smooth(colour="black", method=lm, se=F, aes(group=1))+
  theme_bw()

#knitt graphs together 
Figure.3<-grid.arrange(Figure.3A, Figure.3B,ncol=2)

#test linear models 
#bootstrap coefficents 
Boot.1.1<-frair_boot(FR.rogers.1.1,nboot=20)
Boot.1.2<-frair_boot(FR.rogers.1.2,nboot=20)
Boot.2.1<-frair_boot(FR.rogers.2.1,nboot=20)
Boot.2.2<-frair_boot(FR.rogers.2.2,nboot=20)
Boot.2.3<-frair_boot(FR.rogers.2.3,nboot=20)
Boot.2.4<-frair_boot(FR.rogers.2.4,nboot=20)
Boot.2.5<-frair_boot(FR.rogers.2.5,nboot=20)
Boot.2.6<-frair_boot(FR.rogers.2.6,nboot=20)
Boot.2.7<-frair_boot(FR.rogers.2.7,nboot=20)
Boot.2.8<-frair_boot(FR.rogers.2.8,nboot=20)
Boot.3.1<-frair_boot(FR.rogers.3.1,nboot=20)
Boot.3.2<-frair_boot(FR.rogers.3.2,nboot=20)
Boot.4.1<-frair_boot(FR.rogers.4.1,nboot=20)
Boot.4.2<-frair_boot(FR.rogers.4.2,nboot=20)
Boot.4.3<-frair_boot(FR.rogers.4.3,nboot=20)
Boot.5.1<-frair_boot(FR.5.1.Type1,nboot=20)
Boot.5.2<-frair_boot(FR.5.2.Type1,nboot=20)
Boot.5.3<-frair_boot(FR.5.3.Type1,nboot=20)
Boot.5.4<-frair_boot(FR.5.4.Type1,nboot=20)
Boot.5.5<-frair_boot(FR.5.5.Type1,nboot=20)
Boot.5.6<-frair_boot(FR.5.6.Type1,nboot=20)
Boot.5.7<-frair_boot(FR.5.7.Type1,nboot=20)
Boot.5.8<-frair_boot(FR.5.8.Type1,nboot=20)
Boot.5.9<-frair_boot(FR.5.9.Type1,nboot=20)
Boot.5.10<-frair_boot(FR.5.10.Type1,nboot=20)
Boot.6.1<-frair_boot(FR.rogers.6.1,nboot=20)
Boot.6.2<-frair_boot(FR.6.2.Type1,nboot=20)
Boot.6.3<-frair_boot(FR.6.3.Type1,nboot=20)
Boot.7.1<-frair_boot(FR.rogers.7.1,nboot=20)
Boot.8.1<-frair_boot(FR.rogers.8.1,nboot=20)
Boot.8.2<-frair_boot(FR.rogers.8.2,nboot=20)
Boot.8.3<-frair_boot(FR.rogers.8.3,nboot=20)
Boot.9.1<-frair_boot(FR.rogers.9.1,nboot=20)
Boot.9.2<-frair_boot(FR.rogers.9.2,nboot=20)
Boot.3.3<-nlsBoot(FR.rogers.3.3,niter=20)
Boot.3.4<-nlsBoot(FR.rogers.3.4,niter=20)


Boots<-list(Boot.1.1$bootcoefs[,1:2],Boot.1.2$bootcoefs[,1:2],Boot.1.2$bootcoefs[,1:2], Boot.2.1$bootcoefs[,1:2], Boot.2.3$bootcoefs[,1:2], Boot.2.4$bootcoefs[,1:2],
            Boot.2.5$bootcoefs[,1:2], Boot.2.6$bootcoefs[,1:2], Boot.2.7$bootcoefs[,1:2], Boot.2.8$bootcoefs[,1:2], Boot.3.1$bootcoefs[,1:2], Boot.3.2$bootcoefs[,1:2],
            Boot.3.3$coefboot, Boot.3.4$coefboot, Boot.4.1$bootcoefs[,1:2],Boot.4.2$bootcoefs[,1:2],Boot.4.3$bootcoefs[,1:2],Boot.5.1$bootcoefs, Boot.5.2$bootcoefs, 
            Boot.5.3$bootcoefs,Boot.5.4$bootcoefs,Boot.5.5$bootcoefs,Boot.5.6$bootcoefs,Boot.5.7$bootcoefs,Boot.5.8$bootcoefs,Boot.5.9$bootcoefs,Boot.5.10$bootcoefs,
            Boot.6.1$bootcoefs[,1:2],Boot.6.2$bootcoefs,Boot.6.3$bootcoefs, Boot.7.1$bootcoefs[,1:2],Boot.8.1$bootcoefs[,1:2],Boot.8.2$bootcoefs[,1:2],Boot.8.3$bootcoefs[,1:2],
            Boot.9.1$bootcoefs[,1:2],Boot.9.2$bootcoefs[,1:2])

Boots.DF<-as.data.frame(do.call("rbind", Boots))

#add model names
Boots.names<-(EXP$model)
Boots.DF$model<-rep(Boots.names, each=20)

#add effect size (SMD)
Boots.exp<- merge(Boots.DF,EXP,by="model")
Boot.table<-merge(Boots.exp,Figure2.table.subset,by="species")
Boot.table.standard<-ddply(Boot.table, .(model), transform, a.standard=a/(volume*predator.size), h.standard=h/ratio)

#remove NA
Boot.table.standard.a<-Boot.table.standard[!is.na(Boot.table.standard$a),]

#run linear model for attack rate 
Model.7<-lm(ES~a.standard, data=Boot.table.standard.a)
summary(Model.7)

#run linear model for max feeding rate 
#subset Type II
TypeI<-c("FR.5.1.Type1","FR.5.2.Type1","FR.5.3.Type1","FR.5.4.Type1","FR.5.5.Type1","FR.5.6.Type1","FR.5.7.Type1","FR.5.8.Type1","FR.5.9.Type1","FR.5.10.Type1","FR.6.2.Type1", "FR.6.3.Type1")
Boot.table.standard.h<-Boot.table.standard[ ! Boot.table.standard$model %in% TypeI, ]

Model.8<-lm(ES~h.standard, data=Boot.table.standard.h)
summary(Model.8)

#run t-test on trophic.level
Model.9<-t.test(a.standard~trophic.level, data=Figure3.table)
Model.9

Model.10<-t.test(h.standard~trophic.level, data=Figure3B.table)
Model.10
