library(readxl)
library(dplyr)
library(pammtools)
library(mgcv)

hbsdata=data.frame(coefs=NA,se=NA,pval=NA)
moidata=data.frame(coefs=NA,se=NA,pval=NA)

for(i in 1:100){

dt <- read_excel("Cleaned data.xlsx")
set.seed(i+1986)
dt$random = runif(nrow(dt))

for(j in 1:nrow(dt)){
  dt$weeks[j]=ifelse(dt$start[j]==dt$date[j],dt$random[j]*4,(dt$date[j]-dt$date[j-1])*dt$random[j]+dt$weeks[j-1])
}

dt = dt %>% 
  group_by(infection_id) %>%
  mutate(time=weeks,cumsum=cumsum(clinical),cumsum=ifelse(cumsum>1,1,cumsum),
         cumcause =cumsum(cause))

#resolve malaria at the moment of symptomatic unless gametocytes are observed.

dt = dt %>% mutate(cause = ifelse(clinical==1 & cumcause==0,2,cause),
                   cause = ifelse(clinical==0 & cumsum==1,0,cause),
                   qpcr = ifelse(clinical==0 & cumsum==1,0,qpcr),
                   Tgam_uL_adj = ifelse(clinical==0 & cumsum==1,0,Tgam_uL_adj))


dtrm = dt %>% filter(!(clinical==0 & cumsum==1))

dtf <- dtrm %>% group_by(infection_id)  %>% mutate(status=ifelse(cause==0,0,1),state=ifelse(sum(cause)==1,1,ifelse(sum(cause)==2,2,0)),symptomatic=ifelse(sum(clinical)>0,1,0),clinical=cumsum(clinical),clinical=ifelse(clinical>1,1,clinical),iniclinical=first(clinical)) %>% ungroup() 
table(dtf$cause)

#asymptomatics only
#dtf = dtf %>% filter(iniclinical==0)


dtf2 = dtf %>% group_by(infection_id) %>% mutate(qpcr=first(qpcr),clinical=first(clinical))

faildt = unique(dtf2[dtf2$state==dtf2$cause | dtf2$state==0,] %>% 
                  group_by(infection_id) %>%
                  mutate(time=ifelse(status==0,max(time),time),date=data.table::fifelse(status==0,last(date),date)) %>% dplyr::select("infection_id","cohortid","time","state","start","gender", "agecat", "HbS", "G6PD376", "G6PD202","moi","qpcr","clinical","iniclinical")) %>%
  mutate(time=ifelse(time==0,1,time),id=infection_id) %>% select(-c("qpcr","clinical"))


cuts=sort(unique(c(faildt$time)))

ped <- as_ped(
  data = faildt,
  formula = Surv(time, state) ~ .,
  cut=cuts,
  id = "id", combine = TRUE) %>%
  mutate(cause=as.factor(cause),
         time=tend,
         time=ifelse(tend==1,0,time))

ped = ped


ped = ped %>% select(-c("time")) %>% mutate(date=as.Date(start)+tstart)

ped = ped %>% select(c(colnames(ped)[-21]),"cause")

ped= ped %>% mutate(male= ifelse(gender=="Male",1,0),
                    hbs.homo = ifelse(HbS=="homo",1,0),
                    hbs.hetero = ifelse(HbS=="hetero",1,0),
                    g6pd376.homo = ifelse(G6PD376=="homo",1,0),
                    g6pd376.hetero = ifelse(G6PD376=="hetero",1,0),
                    g6pd202.homo = ifelse(G6PD202=="homo",1,0),
                    g6pd202.hetero = ifelse(G6PD202=="hetero",1,0),
                    age2 = ifelse(agecat==2,1,0),
                    age3 = ifelse(agecat==3,1,0),
                    moib = ifelse(moi>1,1,0))

# analysis ----------------------------------------------------------------
mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
              (hbs.homo+hbs.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
sm.c=summary(mod.c)

hbsdata[i,]=c(mod.c$coefficients[5],sm.c$se[5],sm.c$p.pv[5])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             (moib):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
sm.c=summary(mod.c)

moidata[i,]=c(mod.c$coefficients[3],sm.c$se[3],sm.c$p.pv[3])

print(i)

}


sum(I(hbsdata$pval<0.05)) #81%
sum(I(moidata$pval<0.05)) #10%


mean1=mean(hbsdata$coefs)
mean2=mean(moidata$coefs)

se1 = (mean(sqrt(hbsdata$se))+(1+1/100)*sd(hbsdata$coefs))^2
se2 = (mean(sqrt(moidata$se))+(1+1/100)*sd(moidata$coefs))^2

zscore1=mean1/se1
zscore2=mean2/se2

2*pnorm(q=abs(zscore1), lower.tail=FALSE)
2*pnorm(q=abs(zscore2), lower.tail=FALSE)

hist(hbsdata$pval)
hist(moidata$pval)

# Longitudinal model ------------------------------------------------------

hbsdatapcr=data.frame(coefs=NA,se=NA,pval=NA)
moidatapcr=data.frame(coefs=NA,se=NA,pval=NA)

hbsdatagam=data.frame(coefs=NA,se=NA,pval=NA)
moidatagam=data.frame(coefs=NA,se=NA,pval=NA)

hbsdatagf=data.frame(coefs=NA,se=NA,pval=NA)
moidatagf=data.frame(coefs=NA,se=NA,pval=NA)

for(i in 1:100){
  
  dt <- read_excel("Cleaned data.xlsx")
  
  dtsub=unique(dt %>% select(infection_id))
    set.seed(i+1986)
    dtsub$random = runif(nrow(dtsub))
  
  
  dt = dt %>% left_join(dtsub) %>% 
    group_by(infection_id) %>%
    mutate(weeks=weeks+4*random,cumsum=cumsum(clinical),cumsum=ifelse(cumsum>1,1,cumsum),
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
  
  
  dtrm$logqpcr = log(dtrm$qpcr+0.1)
  dtrm$loggam = log(dtrm$Tgam_uL_adj+0.001)
  
  mod.c = bam(logqpcr ~ hbs.homo+hbs.hetero+s(weeks, bs="bs")+ s(infection_id, bs="re"), 
              family=gaussian(),
              data=dtrm)
  sm.c = summary(mod.c)
  
  hbsdatapcr[i,]=c(mod.c$coefficients[3],sm.c$se[3],sm.c$p.pv[3])
  
  mod.c = bam(logqpcr ~ moib+s(weeks, bs="bs")+ s(infection_id, bs="re"), 
              family=gaussian(),
              data=dtrm)
  sm.c = summary(mod.c)
  
  moidatapcr[i,]=c(mod.c$coefficients[2],sm.c$se[2],sm.c$p.pv[2])
  
  mod.c = bam(loggam ~  hbs.homo+hbs.hetero+
                s(weeks, bs="bs")+ s(logqpcr,bs="bs")+ s(infection_id, bs="re"), 
              family=gaussian(),
              data=dtrm)
  sm.c = summary(mod.c)
  
  hbsdatagam[i,]=c(mod.c$coefficients[3],sm.c$se[3],sm.c$p.pv[3])
  
  mod.c = bam(loggam ~  moib+
                s(weeks, bs="bs")+ s(logqpcr,bs="bs")+ s(infection_id, bs="re"), 
              family=gaussian(),
              data=dtrm)
  sm.c = summary(mod.c)
  
  moidatagam[i,]=c(mod.c$coefficients[2],sm.c$se[2],sm.c$p.pv[2])
  
  
  dtrm$logqpcr = log(dtrm$qpcr+0.1)
  dtrm$loggamfrac = ifelse(dtrm$qpcr==0,NA,log((dtrm$Tgam_uL_adj/dtrm$qpcr)+0.001))
  
  mod.c = bam(loggamfrac ~  hbs.homo+hbs.hetero+
                s(weeks, bs="bs")+ s(logqpcr,bs="bs")+ s(infection_id, bs="re"), 
              family=gaussian(),
              data=dtrm)
  sm.c = summary(mod.c)
  
  hbsdatagf[i,]=c(mod.c$coefficients[3],sm.c$se[3],sm.c$p.pv[3])
  
  mod.c = bam(loggamfrac ~  moib+
                s(weeks, bs="bs")+ s(logqpcr,bs="bs")+ s(infection_id, bs="re"), 
              family=gaussian(),
              data=dtrm)
  sm.c = summary(mod.c)
  
  moidatagf[i,]=c(mod.c$coefficients[2],sm.c$se[2],sm.c$p.pv[2])
  
  
  print(i)
  
  
}




sum(I(hbsdata$pval<0.05)) #81%
sum(I(moidata$pval<0.05)) #10%
sum(I(hbsdatapcr$pval<0.05)) #0%
sum(I(moidatapcr$pval<0.05)) #100%
sum(I(hbsdatagam$pval<0.05)) #100%
sum(I(moidatagam$pval<0.05)) #96%
sum(I(hbsdatagf$pval<0.05)) #99%
sum(I(moidatagf$pval<0.05)) #0%


library(ggplot2)
library(ggprism)
library(ggpubr)

ggarrange(
ggplot()+
  geom_histogram(data=hbsdatapcr,aes(x=pval,fill="HbAS"),alpha=0.5)+
  geom_histogram(data=moidatapcr,aes(x=pval,fill="MOI"),alpha=0.5)+
  scale_fill_viridis_d(name="")+
  xlab("P-values")+
  ylab("Relative frequency (%)")+
  theme_prism()+
  theme(legend.position = "bottom")+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  ggtitle("Parasite density trajectories over time"),

ggplot()+
  geom_histogram(data=hbsdata,aes(x=pval,fill="HbAS"),alpha=0.5)+
  geom_histogram(data=moidata,aes(x=pval,fill="MOI"),alpha=0.5)+
  scale_fill_viridis_d(name="")+
  xlab("P-values")+
  ylab("Relative frequency (%)")+
  theme_prism()+
  theme(legend.position = "bottom")+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  ggtitle("Time to starting gametocyte production"),


ggplot()+
  geom_histogram(data=hbsdatagam,aes(x=pval,fill="HbAS"),alpha=0.5)+
  geom_histogram(data=moidatagam,aes(x=pval,fill="MOI"),alpha=0.5)+
  scale_fill_viridis_d(name="")+
  xlab("P-values")+
  ylab("Relative frequency (%)")+
  theme_prism()+
  theme(legend.position = "bottom")+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  ggtitle("Gametocyte density trajectories over time"),

ggplot()+
  geom_histogram(data=hbsdatagf,aes(x=pval,fill="HbAS"),alpha=0.5)+
  geom_histogram(data=moidatagf,aes(x=pval,fill="MOI"),alpha=0.5)+
  scale_fill_viridis_d(name="")+
  xlab("P-values")+
  ylab("Relative frequency (%)")+
  theme_prism()+
  theme(legend.position = "bottom")+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  ggtitle("Gametocyte fraction trajectories over time"),

nrow=2, ncol=2, common.legend = TRUE)

ggsave("figures/sensitivity_pvals.pdf",device="pdf",dpi=300, width=12, height=9)



ggarrange(
  ggplot()+
    geom_histogram(data=hbsdatapcr,aes(x=exp(coefs),fill="HbAS"),alpha=0.5)+
    geom_histogram(data=moidatapcr,aes(x=exp(coefs),fill="MOI"),alpha=0.5)+
    scale_fill_viridis_d(name="")+
    xlab("Effect ratio")+
    ylab("Relative frequency (%)")+
    theme_prism()+
    scale_x_log10()+
    theme(legend.position = "bottom")+
    geom_vline(xintercept = 1,linetype="dashed")+
    ggtitle("Parasite density trajectories over time"),
  
  ggplot()+
    geom_histogram(data=hbsdata,aes(x=exp(coefs),fill="HbAS"),alpha=0.5)+
    geom_histogram(data=moidata,aes(x=exp(coefs),fill="MOI"),alpha=0.5)+
    scale_fill_viridis_d(name="")+
    xlab("Effect ratio")+
    ylab("Relative frequency (%)")+
    theme_prism()+
    scale_x_log10()+
    theme(legend.position = "bottom")+
    geom_vline(xintercept = 1,linetype="dashed")+
    ggtitle("Time to starting gametocyte production"),
  
  
  ggplot()+
    geom_histogram(data=hbsdatagam,aes(x=exp(coefs),fill="HbAS"),alpha=0.5)+
    geom_histogram(data=moidatagam,aes(x=exp(coefs),fill="MOI"),alpha=0.5)+
    scale_fill_viridis_d(name="")+
    xlab("Effect ratio")+
    ylab("Relative frequency (%)")+
    theme_prism()+    
    scale_x_log10()+
    theme(legend.position = "bottom")+
    geom_vline(xintercept = 1,linetype="dashed")+
    ggtitle("Gametocyte density trajectories over time"),
  
  ggplot()+
    geom_histogram(data=hbsdatagf,aes(x=exp(coefs),fill="HbAS"),alpha=0.5)+
    geom_histogram(data=moidatagf,aes(x=exp(coefs),fill="MOI"),alpha=0.5)+
    scale_fill_viridis_d(name="")+
    xlab("Effect ratio")+
    ylab("Relative frequency (%)")+
    theme_prism()+
    scale_x_log10()+
    theme(legend.position = "bottom")+
    geom_vline(xintercept = 1,linetype="dashed")+
    ggtitle("Gametocyte fraction trajectories over time"),
  
  nrow=2, ncol=2, common.legend = TRUE)

ggsave("figures/sensitivity_coefs.pdf",device="pdf",dpi=300, width=12, height=9)



