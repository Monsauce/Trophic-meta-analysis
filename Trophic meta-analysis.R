#R script for "Predators have greater impacts on recipient freswater communities" 
library(ggplot2)
library(metafor)
library(RCurl)

#read in coarse scraped data from literature from GitHub
coarse.URL <- getURL("https://raw.githubusercontent.com/Monsauce/Trophic-meta-analysis/master/Meta.Coarse.csv")
coarse<-read.csv(text=coarse.URL)

#calculate effect size using Hedge's D-standardized mean difference (due to negative values)
meta.coarse.var<-escalc(measure="SMD", m1i=mean, m2i=control, sd1i=mean.SD, sd2i=control.SD,
                        n1i=mean.N, n2i=control.N,data=meta.coarse, var.names=c("SMD","SMD_var"),digits=4)

#run fixed effects model on trophic.group and field.lab.group 
res <- rma(LRR, LRR_var, mods=cbind(trophic, field.lab.group), data = meta.coarse.var)
