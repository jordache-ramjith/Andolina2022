
library(dplyr)

ped <- readRDS("ped.Rds")

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
library(pammtools)
library(mgcv)
library(ggplot2)
library(ggprism)
library(ggpubr)


# Baseline hazard and cumulative incidence including clinical cases -------

mod0 <- pamm(ped_status ~ cause+s(tend, by = cause,k=5,bs="bs"), 
             data = ped)
summary(mod0)

newd =  expand.grid(tend=1:98,cause=levels(ped$cause),intlen=1) %>% 
  group_by(cause) %>%
  add_hazard(mod0) %>% add_cif(mod0)

plot1=ggplot(data=newd, aes(x=tend/7,y=hazard))+
  geom_line(aes(col=cause),size=1.2)+
  geom_ribbon(aes(ymin=ci_lower,ymax=ci_upper,fill=cause),alpha=0.4)+
  theme_prism()+
  theme(legend.position = "top")+
  ylab("Incidence rate")+
   coord_cartesian(ylim=c(0,0.1))+
  scale_color_brewer(palette="Dark2",labels=c("Gametocyte initiation", "Malaria resolution"))+
  scale_fill_brewer(palette="Dark2",labels=c("Gametocyte initiation", "Malaria resolution"))+
  guides(fill="none")+
  scale_x_continuous(breaks=seq(0,14,by=2))+
  xlab("Weeks from first detection of incident malaria")
plot1

ggsave("figures/gametocyte_incidence_rates_over_time.pdf",device="pdf",dpi=300, width=8, height=6)


plot2=ggplot(data=newd, aes(x=tend/7,y=cif))+
  geom_line(aes(col=cause),size=1.2)+
  geom_ribbon(aes(ymin=cif_lower,ymax=cif_upper,fill=cause),alpha=0.4)+
  theme_prism()+
  theme(legend.position = "top")+
  ylab("Cumulative incidence")+
  scale_color_brewer(palette="Dark2",labels=c("Gametocyte initiation", "Malaria resolution"))+
  scale_fill_brewer(palette="Dark2",labels=c("Gametocyte initiation", "Malaria resolution"))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  guides(fill="none")+
  scale_x_continuous(breaks=seq(0,14,by=2))
plot2

plotc=ggarrange(plot1+ rremove("xlab"),
          plot2+ rremove("xlab"),
          nrow=1,
          common.legend=TRUE,
          legend="top")

plotc=annotate_figure(plotc,
                bottom = text_grob("Weeks from first detection of incident malaria",face="bold",size=14))

plotc

#ggsave("figures/gametocyte_incidence_and_cumulative_incidence.pdf",device="pdf",dpi=300, width=12, height=8)



# Baseline hazard and cumulative incidence EXcluding clinical cases -------

ped2 = ped[ped$clinical==0,]

mod0 <- pamm(ped_status ~ cause+s(tend, by = cause,k=5,bs="bs"), 
             data = ped2)
summary(mod0)

newd.ex =  expand.grid(tend=1:98,cause=levels(ped$cause),intlen=1) %>% 
  group_by(cause) %>%
  add_hazard(mod0) %>% add_cif(mod0)


plot2ex=ggplot()+
  geom_line(data=newd[newd$cause=="1",], 
            aes(x=tend/7,y=cif,col="Symptomatics inc."),size=1.2)+
  geom_ribbon(data=newd[newd$cause=="1",], 
              aes(x=tend/7,y=cif,ymin=cif_lower,ymax=cif_upper,fill="Symptomatics inc."),alpha=0.4)+  
  geom_line(data=newd.ex[newd$cause=="1",], 
            aes(x=tend/7,y=cif,col="Symptomatics ex."),size=1.2)+
  geom_ribbon(data=newd.ex[newd$cause=="1",], 
              aes(x=tend/7,y=cif,ymin=cif_lower,ymax=cif_upper,fill="Symptomatics ex."),alpha=0.4)+
  theme_prism()+
  theme(legend.position = "top")+
  xlab("Time from first detection of incident malaria")+
  ylab("Cumulative incidence")+
  scale_colour_hue(l=45)+
  scale_fill_hue(l=45)+
  guides(fill="none")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_x_continuous(breaks=seq(0,14,by=2))
plot2ex

plotc=ggarrange(ggarrange(plot1+ rremove("xlab"),
                plot2+ rremove("xlab"), nrow=1, legend="top", common.legend = TRUE),
                ggarrange(plot2ex+rremove("xlab"),
                          common.legend = TRUE),
                nrow=1,
                legend="top", widths = c(2, 1))

plotc

annotate_figure(plotc,
                bottom = text_grob("Time from first detection of incident malaria (weeks)",face="bold",size=14))

#ggsave("gam_inc.pdf",device="pdf",dpi=300, width=12, height=8)

ped$logqpcr = log(ped$qpcr+0.1)

# Parasite density influence ----------------------------------------------


#saveRDS((ped %>% select(id, tstart,tend,interval, ped_status,offset,cause,logqpcr)),"data for figures/Figure3_PAMdata.RDS")


pam = pamm(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") + 
            s(logqpcr,by=cause,k=5,bs="bs"), 
          data = ped)
summary(pam)

newd2 =  expand.grid(tend=28,
                     logqpcr=seq(-2,5,by=0.1),
                     cause="1",intlen=1) %>% 
  group_by(cause) %>%
  add_hazard(pam) %>% add_cif(pam)
newd2$qpcr = exp(newd2$logqpcr)-0.1

library(scales)
plot3=ggplot(data=newd2, aes(x=qpcr,y=hazard))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=ci_lower,ymax=ci_upper),alpha=0.2)+
  theme_prism()+
  theme(legend.position = "top")+
  xlab("Parasite density/\U00B5L")+
  ylab("Incidence rate (gametocyte initiation)")+
  coord_cartesian(ylim=c(0,0.08))+
  guides(fill="none",col="none")+
  scale_x_continuous(
    trans = log_trans(),
    breaks = trans_breaks("log", function(x) exp(x)),
    labels = trans_format("log", function(x) round(exp(x),1))
  ) 

plot3

newd2b =  expand.grid(tend=1:100,logqpcr=seq(-2,5,by=0.1),cause=levels(ped$cause),intlen=1) %>% 
  group_by(cause) %>%
  add_hazard(pam) %>% group_by(cause,logqpcr) %>% add_cif(pam)

newd2b$cause2=factor(newd2b$cause,labels=c("Gametocyte initiation", "Malaria resolution"))
newd2b$qpcr=exp(newd2b$logqpcr)-0.1



plot4=ggplot(data=newd2b[newd2b$cause=="1",], aes(x=tend/7,y=exp(logqpcr)-0.1,fill=cif))+
  geom_tile(width=1,height=1)+
  theme_prism()+
  theme(legend.position = "bottom")+
  xlab("Time from first detection of incident malaria (weeks)")+
  ylab("Parasite density/\U00B5L")+
  scale_y_continuous(
    trans = log_trans(),
    breaks = trans_breaks("log", function(x) exp(x)),
    labels = trans_format("log", function(x) round(exp(x),1))
  ) +
scale_fill_viridis_c(name="Cumulative\nIncidence",breaks=seq(0,1,0.2),limits=c(0,1))+
  theme(legend.title = element_text())+
  scale_x_continuous(breaks=seq(0,14,by=2))
plot4


plot.pd=ggarrange(plot3,
                plot4,
                nrow=1,
                legend="top",
                labels=c("A","B"))
plot.pd

ggsave("figures/gametocyte_incidence_by_parasite_density.pdf",device="pdf",dpi=300, width=12, height=8)



# Univariable effects -----------------------------------------------------

PDoutput=data.frame(matrix(NA, nrow=20, ncol=5))
rownames(PDoutput)=c("GAM inc: age2","Res: age2","GAM inc: age3","Res: age3","GAM inc: male","Res: male", "GAM inc: hbs.homo", "Res: hbs.homo","GAM inc: hbs.hetero","Res: hbs.hetero","GAM inc: g6pd376.homo","Res: g6pd376.homo","GAM inc: g6pd376.hetero","Res: g6pd376.hetero", "GAM inc: g6pd202.homo", "Res: g6pd202.homo", "GAM inc: g6pd202.hetero", "Res: g6pd202.hetero","GAM inc: moib","Res: moib")
colnames(PDoutput) <- c("cause","HR","HR.lci","HR.uci","pval")

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             (age2+age3):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = ped %>% make_newdata(age2=c(1),age3=c(0),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(age2=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: age2",] = cbind(newd[1,],summary(mod.c)$p.pv[3])
PDoutput["Res: age2",] = cbind(newd[2,],summary(mod.c)$p.pv[4])

newd = ped %>% make_newdata(age2=c(0),age3=c(1),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(age3=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: age3",] = cbind(newd[1,],summary(mod.c)$p.pv[5])
PDoutput["Res: age3",] = cbind(newd[2,],summary(mod.c)$p.pv[6])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             male:cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = ped %>% make_newdata(male=c(1),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(male=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: male",] = cbind(newd[1,],summary(mod.c)$p.pv[3])
PDoutput["Res: male",] = cbind(newd[2,],summary(mod.c)$p.pv[4])


mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             (hbs.homo+hbs.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = ped %>% make_newdata(hbs.homo=c(1),hbs.hetero=c(0),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(hbs.homo=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: hbs.homo",] = cbind(newd[1,],summary(mod.c)$p.pv[3])
PDoutput["Res: hbs.homo",] = cbind(newd[2,],summary(mod.c)$p.pv[4])

newd = ped %>% make_newdata(hbs.homo=c(0),hbs.hetero=c(1),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(hbs.hetero=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: hbs.hetero",] = cbind(newd[1,],summary(mod.c)$p.pv[5])
PDoutput["Res: hbs.hetero",] = cbind(newd[2,],summary(mod.c)$p.pv[6])


mod.c =gam(ped_status ~cause+ s(tend,by=cause,k=5,bs="bs") +
             (g6pd202.homo+g6pd202.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = ped %>% make_newdata(g6pd202.homo=c(1),g6pd202.hetero=c(0),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(g6pd202.homo=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: g6pd202.homo",] = cbind(newd[1,],summary(mod.c)$p.pv[3])
PDoutput["Res: g6pd202.homo",] = cbind(newd[2,],summary(mod.c)$p.pv[4])


newd = ped %>% make_newdata(g6pd202.homo=c(0),g6pd202.hetero=c(1),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(g6pd202.hetero=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: g6pd202.hetero",] = cbind(newd[1,],summary(mod.c)$p.pv[5])
PDoutput["Res: g6pd202.hetero",] = cbind(newd[2,],summary(mod.c)$p.pv[6])


mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             (g6pd376.homo+g6pd376.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = ped %>% make_newdata(g6pd376.homo=c(1),g6pd376.hetero=c(0),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(g6pd376.homo=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: g6pd376.homo",] = cbind(newd[1,],summary(mod.c)$p.pv[3])
PDoutput["Res: g6pd376.homo",] = cbind(newd[2,],summary(mod.c)$p.pv[4])

newd = ped %>% make_newdata(g6pd376.homo=c(0),g6pd376.hetero=c(1),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(g6pd376.hetero=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: g6pd376.hetero",] = cbind(newd[1,],summary(mod.c)$p.pv[5])
PDoutput["Res: g6pd376.hetero",] = cbind(newd[2,],summary(mod.c)$p.pv[6])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             moib:cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = ped %>% make_newdata(moib=c(1),cause=c(1,2)) %>% add_hazard(mod.c,ref=list(moib=c(0))) %>% select(cause, hazard, ci_lower, ci_upper)

PDoutput["GAM inc: moib",] = cbind(newd[1,],summary(mod.c)$p.pv[3])
PDoutput["Res: moib",] = cbind(newd[2,],summary(mod.c)$p.pv[4])

PDoutput$Var = rownames(PDoutput)
rownames(PDoutput)=as.character(1:nrow(PDoutput))

PDoutput[order(PDoutput$cause),] %>% select(Var,HR,HR.lci,HR.uci,pval)


# Adjusted for pcr age sex ----------------------------


PDoutput.x=data.frame(matrix(NA, nrow=20, ncol=5))
rownames(PDoutput.x)=c("GAM inc: age2","Res: age2","GAM inc: age3","Res: age3","GAM inc: male","Res: male", "GAM inc: hbs.homo", "Res: hbs.homo","GAM inc: hbs.hetero","Res: hbs.hetero","GAM inc: g6pd376.homo","Res: g6pd376.homo","GAM inc: g6pd376.hetero","Res: g6pd376.hetero", "GAM inc: g6pd202.homo", "Res: g6pd202.homo", "GAM inc: g6pd202.hetero", "Res: g6pd202.hetero","GAM inc: moib","Res: moib")
colnames(PDoutput.x) <- c("cause","HR","HR.lci","HR.uci","pval")

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
          (age2+age3+male):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(ped %>% make_newdata(age2=c(1),age3=c(0),male=c(0),clinical=c(0),cause=c(1,2)) %>% 
  add_hazard(mod.c,ref=list(age2=c(0))) %>% 
  select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="age2"),
ped %>% make_newdata(age2=c(0),age3=c(1),male=c(0),clinical=c(0),cause=c(1,2)) %>% 
  add_hazard(mod.c,ref=list(age3=c(0))) %>% 
  select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="age3"),
ped %>% make_newdata(age2=c(0),age3=c(0),male=c(1),clinical=c(0),cause=c(1,2)) %>% 
  add_hazard(mod.c,ref=list(male=c(0))) %>% 
  select(cause, hazard, ci_lower, ci_upper)  %>% mutate(var="male")
)

PDoutput.x["GAM inc: age2",] = cbind(newd[1,1:4],summary(mod.c)$p.pv[3])
PDoutput.x["Res: age2",] = cbind(newd[2,1:4],summary(mod.c)$p.pv[4])
PDoutput.x["GAM inc: age3",] = cbind(newd[3,1:4],summary(mod.c)$p.pv[5])
PDoutput.x["Res: age3",] = cbind(newd[4,1:4],summary(mod.c)$p.pv[6])
PDoutput.x["GAM inc: male",] = cbind(newd[5,1:4],summary(mod.c)$p.pv[7])
PDoutput.x["Res: male",] = cbind(newd[6,1:4],summary(mod.c)$p.pv[8])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+hbs.homo+hbs.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),hbs.homo=c(1),hbs.hetero=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(hbs.homo=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="hbs.homo"),
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),hbs.homo=c(0),hbs.hetero=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(hbs.hetero=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="hbs.hetero"))
newd             

PDoutput.x["GAM inc: hbs.homo",] = cbind(newd[7,1:4],summary(mod.c)$p.pv[9])
PDoutput.x["Res: hbs.homo",] = cbind(newd[8,1:4],summary(mod.c)$p.pv[10])
PDoutput.x["GAM inc: hbs.hetero",] = cbind(newd[9,1:4],summary(mod.c)$p.pv[11])
PDoutput.x["Res: hbs.hetero",] = cbind(newd[10,1:4],summary(mod.c)$p.pv[12])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+g6pd376.homo+g6pd376.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd376.homo=c(1),g6pd376.hetero=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd376.homo=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd376.homo"),
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd376.homo=c(0),g6pd376.hetero=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd376.hetero=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd376.hetero"))

PDoutput.x["GAM inc: g6pd376.homo",] = cbind(newd[11,1:4],summary(mod.c)$p.pv[9])
PDoutput.x["Res: g6pd376.homo",] = cbind(newd[12,1:4],summary(mod.c)$p.pv[10])
PDoutput.x["GAM inc: g6pd376.hetero",] = cbind(newd[13,1:4],summary(mod.c)$p.pv[11])
PDoutput.x["Res: g6pd376.hetero",] = cbind(newd[14,1:4],summary(mod.c)$p.pv[12])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+g6pd202.homo+g6pd202.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd202.homo=c(1),g6pd202.hetero=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd202.homo=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd202.homo"),
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd202.homo=c(0),g6pd202.hetero=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd202.hetero=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd202.hetero"))

PDoutput.x["GAM inc: g6pd202.homo",] = cbind(newd[15,1:4],summary(mod.c)$p.pv[9])
PDoutput.x["Res: g6pd202.homo",] = cbind(newd[16,1:4],summary(mod.c)$p.pv[10])
PDoutput.x["GAM inc: g6pd202.hetero",] = cbind(newd[17,1:4],summary(mod.c)$p.pv[11])
PDoutput.x["Res: g6pd202.hetero",] = cbind(newd[18,1:4],summary(mod.c)$p.pv[12])


mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+moib):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),moib=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(moib=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="moi>1"))

PDoutput.x["GAM inc: moib",] = cbind(newd[19,1:4],summary(mod.c)$p.pv[9])
PDoutput.x["Res: moib",] = cbind(newd[20,1:4],summary(mod.c)$p.pv[10])


PDoutput.x$Var = rownames(PDoutput.x)
rownames(PDoutput.x)=as.character(1:nrow(PDoutput.x))

PDoutput.x[order(PDoutput.x$cause),] %>% select(Var,HR,HR.lci,HR.uci,pval)

  

# Adjusted for pcr ----------------------------


PDoutput.y=data.frame(matrix(NA, nrow=20, ncol=5))
rownames(PDoutput.y)=c("GAM inc: age2","Res: age2","GAM inc: age3","Res: age3","GAM inc: male","Res: male", "GAM inc: hbs.homo", "Res: hbs.homo","GAM inc: hbs.hetero","Res: hbs.hetero","GAM inc: g6pd376.homo","Res: g6pd376.homo","GAM inc: g6pd376.hetero","Res: g6pd376.hetero", "GAM inc: g6pd202.homo", "Res: g6pd202.homo", "GAM inc: g6pd202.hetero", "Res: g6pd202.hetero","GAM inc: moib","Res: moib")
colnames(PDoutput.y) <- c("cause","HR","HR.lci","HR.uci","pval")

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(ped %>% make_newdata(age2=c(1),age3=c(0),male=c(0),clinical=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(age2=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="age2"),
             ped %>% make_newdata(age2=c(0),age3=c(1),male=c(0),clinical=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(age3=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="age3"),
             ped %>% make_newdata(age2=c(0),age3=c(0),male=c(1),clinical=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(male=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper)  %>% mutate(var="male")
)

PDoutput.y["GAM inc: age2",] = cbind(newd[1,1:4],summary(mod.c)$p.pv[3])
PDoutput.y["Res: age2",] = cbind(newd[2,1:4],summary(mod.c)$p.pv[4])
PDoutput.y["GAM inc: age3",] = cbind(newd[3,1:4],summary(mod.c)$p.pv[5])
PDoutput.y["Res: age3",] = cbind(newd[4,1:4],summary(mod.c)$p.pv[6])
PDoutput.y["GAM inc: male",] = cbind(newd[5,1:4],summary(mod.c)$p.pv[7])
PDoutput.y["Res: male",] = cbind(newd[6,1:4],summary(mod.c)$p.pv[8])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+hbs.homo+hbs.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),hbs.homo=c(1),hbs.hetero=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(hbs.homo=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="hbs.homo"),
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),hbs.homo=c(0),hbs.hetero=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(hbs.hetero=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="hbs.hetero"))
newd             

PDoutput.y["GAM inc: hbs.homo",] = cbind(newd[7,1:4],summary(mod.c)$p.pv[9])
PDoutput.y["Res: hbs.homo",] = cbind(newd[8,1:4],summary(mod.c)$p.pv[10])
PDoutput.y["GAM inc: hbs.hetero",] = cbind(newd[9,1:4],summary(mod.c)$p.pv[11])
PDoutput.y["Res: hbs.hetero",] = cbind(newd[10,1:4],summary(mod.c)$p.pv[12])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+g6pd376.homo+g6pd376.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd376.homo=c(1),g6pd376.hetero=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd376.homo=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd376.homo"),
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd376.homo=c(0),g6pd376.hetero=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd376.hetero=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd376.hetero"))

PDoutput.y["GAM inc: g6pd376.homo",] = cbind(newd[11,1:4],summary(mod.c)$p.pv[9])
PDoutput.y["Res: g6pd376.homo",] = cbind(newd[12,1:4],summary(mod.c)$p.pv[10])
PDoutput.y["GAM inc: g6pd376.hetero",] = cbind(newd[13,1:4],summary(mod.c)$p.pv[11])
PDoutput.y["Res: g6pd376.hetero",] = cbind(newd[14,1:4],summary(mod.c)$p.pv[12])

mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+g6pd202.homo+g6pd202.hetero):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd202.homo=c(1),g6pd202.hetero=c(0),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd202.homo=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd202.homo"),
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),g6pd202.homo=c(0),g6pd202.hetero=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(g6pd202.hetero=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="g6pd202.hetero"))

PDoutput.y["GAM inc: g6pd202.homo",] = cbind(newd[15,1:4],summary(mod.c)$p.pv[9])
PDoutput.y["Res: g6pd202.homo",] = cbind(newd[16,1:4],summary(mod.c)$p.pv[10])
PDoutput.y["GAM inc: g6pd202.hetero",] = cbind(newd[17,1:4],summary(mod.c)$p.pv[11])
PDoutput.y["Res: g6pd202.hetero",] = cbind(newd[18,1:4],summary(mod.c)$p.pv[12])


mod.c =gam(ped_status ~ cause+s(tend,by=cause,k=5,bs="bs") +
             s(logqpcr,by=cause,k=5,bs="bs")+
             (age2+age3+male+moib):cause, 
           data = ped,
           family = poisson(), 
           offset = offset)
summary(mod.c)
newd = rbind(newd,
             ped %>% 
               make_newdata(age2=c(0),age3=c(0),male=c(0),clinical=c(0),moib=c(1),cause=c(1,2)) %>% 
               add_hazard(mod.c,ref=list(moib=c(0))) %>% 
               select(cause, hazard, ci_lower, ci_upper) %>% mutate(var="moi>1"))

PDoutput.y["GAM inc: moib",] = cbind(newd[19,1:4],summary(mod.c)$p.pv[9])
PDoutput.y["Res: moib",] = cbind(newd[20,1:4],summary(mod.c)$p.pv[10])


PDoutput.y$Var = rownames(PDoutput.y)
rownames(PDoutput.y)=as.character(1:nrow(PDoutput.y))

PDoutput.y[order(PDoutput.y$cause),] %>% select(Var,HR,HR.lci,HR.uci,pval)


PDoutput.y = PDoutput.y[order(PDoutput.y$cause),] %>% select(Var,HR,HR.lci,HR.uci,pval)

PDoutput.y
PDoutput.y1 = PDoutput.y[1:10,]
PDoutput.y2 = PDoutput.y[11:20,]


resdt = data.frame(var=c("age2", "age3","male" ,"hbs.homo","hbs.hetero", "g6pd376.homo", "g6pd376.hetero", "g6pd202.homo","g6pd202.hetero", "moib"),DR=paste0(round(PDoutput.y1$HR,2)," (",round(PDoutput.y1$HR.lci,2),", ",round(PDoutput.y1$HR.uci,2),")"),pval=round(PDoutput.y1$pval,4))

xlsx::write.xlsx2(resdt,"results_estimated.xlsx",sheetName = "gametocyte incidence",row.names = FALSE,append=TRUE)

resdt = data.frame(var=c("age2", "age3","male" ,"hbs.homo","hbs.hetero", "g6pd376.homo", "g6pd376.hetero", "g6pd202.homo","g6pd202.hetero", "moib"),DR=paste0(round(PDoutput.y2$HR,2)," (",round(PDoutput.y2$HR.lci,2),", ",round(PDoutput.y2$HR.uci,2),")"),pval=round(PDoutput.y2$pval,4))

xlsx::write.xlsx2(resdt,"results_estimated.xlsx",sheetName = "Resolution priot to gams",row.names = FALSE,append=TRUE)

