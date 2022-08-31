

# Data cleaning -----------------------------------------------------------


library(readxl)
library(lubridate)
library(dplyr)
library(pammtools)
library(scales)

dt <- read_excel("Cleaned data.xlsx")

dt = dt %>% 
  group_by(infection_id) %>%
  mutate(cumsum=cumsum(clinical),cumsum=ifelse(cumsum>1,1,cumsum),
         cumcause =cumsum(cause))

dt = dt %>% mutate(cause = ifelse(clinical==1 & cumcause==0,2,cause),
                   cause = ifelse(clinical==0 & cumsum==1,0,cause),
                   qpcr = ifelse(clinical==0 & cumsum==1,0,qpcr),
                   Tgam_uL_adj = ifelse(clinical==0 & cumsum==1,0,Tgam_uL_adj))


dtrm = dt %>% filter(!(clinical==0 & cumsum==1))

dtrm= dtrm %>% mutate(male= ifelse(gender=="Male",1,0),
                    hbs.homo = ifelse(HbS=="homo",1,0),
                    hbs.hetero = ifelse(HbS=="hetero",1,0),
                    g6pd376.homo = ifelse(G6PD376=="homo",1,0),
                    g6pd376.hetero = ifelse(G6PD376=="hetero",1,0),
                    g6pd202.homo = ifelse(G6PD202=="homo",1,0),
                    g6pd202.hetero = ifelse(G6PD202=="hetero",1,0),
                    age2 = ifelse(agecat==2,1,0),
                    age3 = ifelse(agecat==3,1,0),
                    moib = ifelse(moi>1,1,0),
                    infection_id=as.factor(infection_id))

saveRDS(as.data.frame(dtrm %>% select(weeks,qpcr,Tgam_uL_adj,infection_id)),"data for figures/Figure5and6data.RDS")

# Analysis ----------------------------------------------------------------

dtrm$logqpcr = log(dtrm$qpcr+0.1)
dtrm$loggam = log(dtrm$Tgam_uL_adj+0.001)

library(mgcv)
library(ggplot2)
library(ggprism)
library(ggpubr)

mod0 = bam(logqpcr ~ s(weeks, bs="bs",k=20)+ s(infection_id, bs="re"), 
           family=gaussian(),
           data=dtrm)
summary(mod0)
gam.vcomp(mod0)

newd = expand.grid(weeks=seq(0,40,0.1), infection_id="1")
preds=predict.gam(mod0,newd,exclude="s(infection_id)",se.fit=TRUE)
newd$logqpcr = preds$fit
newd$logqpcr.lci=preds$fit-1.96*preds$se.fit
newd$logqpcr.uci=preds$fit+1.96*preds$se.fit
newd$qpcr = exp(newd$logqpcr)-0.1
newd$qpcr.lci = exp(newd$logqpcr.lci)-0.1
newd$qpcr.uci = exp(newd$logqpcr.uci)-0.1


plot1=ggplot(data=newd, aes(x=weeks,y=qpcr))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=qpcr.lci,ymax=qpcr.uci),alpha=0.2)+
  theme_prism()+
  ylab("Parasite density/\U00B5L")+
  scale_x_continuous(breaks=seq(0,40,by=4))+
  xlab("Time from first detection of incident malaria (weeks)")
plot1


#ggsave("parasite_dens.pdf",plot=plot1,device="pdf",dpi=300, width=6, height=4)


# Covariate effects -------------------------------------------------------


## Univariable ------------------------------------------------------------
PDoutput=data.frame(matrix(NA, nrow=10, ncol=5))
rownames(PDoutput)=c("age2","age3","male", "hbs.homo","hbs.hetero","g6pd376.homo","g6pd376.hetero", "g6pd202.homo", "g6pd202.hetero","moib")
colnames(PDoutput) <- c("DR","SE","DR.lci","DR.uci","pval")


mod.c = bam(logqpcr ~ age2+age3+s(weeks, bs="bs")+ s(infection_id, bs="re"), 
           family=gaussian(),
           data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, age2=1, age3=0, infection_id="1") %>%
  add_hazard(mod.c, reference=list(age2=0),exclude="s(infection_id)"))
PDoutput["age2",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, age2=0, age3=1, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(age3=0),exclude="s(infection_id)"))
PDoutput["age3",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 


mod.c = bam(logqpcr ~ male+s(weeks, bs="bs")+ s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, male=1, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(male=0),exclude="s(infection_id)"))
PDoutput["male",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:4)] 

mod.c = bam(logqpcr ~ hbs.homo+hbs.hetero+s(weeks, bs="bs")+ s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, hbs.homo=1, hbs.hetero=0, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(hbs.homo=0),exclude="s(infection_id)"))
PDoutput["hbs.homo",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, hbs.homo=0, hbs.hetero=1, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(hbs.hetero=0),exclude="s(infection_id)"))
PDoutput["hbs.hetero",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 

mod.c = bam(logqpcr ~ g6pd376.homo+g6pd376.hetero+s(weeks, bs="bs")+ s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, g6pd376.homo=1, g6pd376.hetero=0, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(g6pd376.homo=0),exclude="s(infection_id)"))
PDoutput["g6pd376.homo",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, g6pd376.homo=0, g6pd376.hetero=1, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(g6pd376.hetero=0),exclude="s(infection_id)"))
PDoutput["g6pd376.hetero",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 

mod.c = bam(logqpcr ~ g6pd202.homo+g6pd202.hetero+
              s(weeks, bs="bs")+ 
              s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, g6pd202.homo=1, g6pd202.hetero=0, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(g6pd202.homo=0),exclude="s(infection_id)"))
PDoutput["g6pd202.homo",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, g6pd202.homo=0, g6pd202.hetero=1, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(g6pd202.hetero=0),exclude="s(infection_id)"))
PDoutput["g6pd202.hetero",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:5)] 


mod.c = bam(logqpcr ~ moib+
              s(weeks, bs="bs")+ 
              s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(expand.grid(tend=1,weeks=4, moib=1, infection_id="1") %>%
                          add_hazard(mod.c, reference=list(moib=0),exclude="s(infection_id)"))
PDoutput["moib",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:4)] 


PDoutput$pval = 2*pnorm(-abs(log(PDoutput$DR)/PDoutput$SE))

PDoutput2 = cbind(round(PDoutput[,1:4],2),round(PDoutput[,5],4))

PDoutput2

# Multivariable results ---------------------------------------------------

PDoutput.x=data.frame(matrix(NA, nrow=10, ncol=5))
rownames(PDoutput.x)=c("age2","age3","male", "hbs.homo","hbs.hetero","g6pd376.homo","g6pd376.hetero", "g6pd202.homo", "g6pd202.hetero","moib")
colnames(PDoutput.x) <- c("DR","SE","DR.lci","DR.uci","pval")

mod.c = bam(logqpcr ~ age2+ age3+ male+ 
              s(weeks, bs="bs")+ 
              s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)

pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=1,age3=0,male=0, cumclinical=0, infection_id="1") %>%
    add_hazard(mod.c, reference=list(age2=0),
               exclude="s(infection_id)"))
PDoutput.x["age2",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:7)] 
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=1,male=0, cumclinical=0, infection_id="1") %>%
    add_hazard(mod.c, reference=list(age3=0),
               exclude="s(infection_id)"))
PDoutput.x["age3",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:7)] 
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=1, cumclinical=0, infection_id="1") %>%
    add_hazard(mod.c, reference=list(male=0),
               exclude="s(infection_id)"))
PDoutput.x["male",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:7)] 


mod.c = bam(logqpcr ~ age2+ age3+ male+ hbs.homo+hbs.hetero+
              s(weeks, bs="bs")+ s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)

pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,hbs.homo=1,hbs.hetero=0, infection_id="1") %>%
    add_hazard(mod.c, reference=list(hbs.homo=0),
               exclude="s(infection_id)"))
PDoutput.x["hbs.homo",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:9)] 
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,hbs.homo=0,hbs.hetero=1, infection_id="1") %>%
    add_hazard(mod.c, reference=list(hbs.hetero=0),
               exclude="s(infection_id)"))
PDoutput.x["hbs.hetero",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:9)] 



mod.c = bam(logqpcr ~ age2+ age3+ male+ g6pd376.homo+g6pd376.hetero+
              s(weeks, bs="bs")+ s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,g6pd376.homo=1,g6pd376.hetero=0, infection_id="1") %>%
    add_hazard(mod.c, reference=list(g6pd376.homo=0),
               exclude="s(infection_id)"))
PDoutput.x["g6pd376.homo",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:9)] 
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,g6pd376.homo=0,g6pd376.hetero=1, infection_id="1") %>%
    add_hazard(mod.c, reference=list(g6pd376.hetero=0),
               exclude="s(infection_id)"))
PDoutput.x["g6pd376.hetero",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:9)] 


mod.c = bam(logqpcr ~ age2+ age3+ male+ g6pd202.homo+g6pd202.hetero+
              s(weeks, bs="bs")+ 
              s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,g6pd202.homo=1,g6pd202.hetero=0, infection_id="1") %>%
    add_hazard(mod.c, reference=list(g6pd202.homo=0),
               exclude="s(infection_id)"))
PDoutput.x["g6pd202.homo",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:9)] 
pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,g6pd202.homo=0,g6pd202.hetero=1, infection_id="1") %>%
    add_hazard(mod.c, reference=list(g6pd202.hetero=0),
               exclude="s(infection_id)"))
PDoutput.x["g6pd202.hetero",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:9)] 


mod.c = bam(logqpcr ~ age2+ age3+ male+ moib+
              s(weeks, bs="bs")+ 
              s(infection_id, bs="re"), 
            family=gaussian(),
            data=dtrm)
summary(mod.c)

pdpred=suppressWarnings(
  expand.grid(tend=1,weeks=4,age2=0,age3=0,male=0, cumclinical=0,moib=1, infection_id="1") %>%
    add_hazard(mod.c, reference=list(moib=0),
               exclude="s(infection_id)"))
PDoutput.x["moib",c("DR","SE","DR.lci","DR.uci")] = pdpred[1,-(1:8)] 


PDoutput.x$pval = 2*pnorm(-abs(log(PDoutput.x$DR)/PDoutput.x$SE))

PDoutput.x2 = cbind(round(PDoutput.x[,1:4],2),round(PDoutput.x[,5],4))

PDoutput.x2



data.frame(DR=paste0(PDoutput2$DR," (",PDoutput2$DR.lci,", ",PDoutput2$DR.uci,")"),
           pval=PDoutput2[,5])

data.frame(DR=paste0(PDoutput.x2$DR," (",PDoutput.x2$DR.lci,", ",PDoutput.x2$DR.uci,")"),
           pval=PDoutput.x2[,5])




resdt = data.frame(var=rownames(PDoutput2),DR=paste0(PDoutput2$DR," (",PDoutput2$DR.lci,", ",PDoutput2$DR.uci,")"), pval=PDoutput2[,5])

xlsx::write.xlsx2(resdt,"results_estimated.xlsx",sheetName = "parasite density",row.names = FALSE)









