

# Data cleaning -----------------------------------------------------------



library(readxl)
library(lubridate)
library(dplyr)
library(survival)
library(survminer)
library(cowplot)
library(ggsci)
library(ggprism)
library(ggpubr)
library(mgcv)

source("Pamms competing risks.R")

dt <- read_excel("Cleaned data.xlsx")

dt = dt %>% 
  group_by(infection_id) %>%
  mutate(cumsum=cumsum(clinical),cumsum=ifelse(cumsum>1,1,cumsum),
         cumcause =cumsum(cause))

dtasym = dt

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
                      moib = ifelse(moi>1,1,0))


faildt= faildt %>% mutate(male= ifelse(gender=="Male",1,0),
                      hbs.homo = ifelse(HbS=="homo",1,0),
                      hbs.hetero = ifelse(HbS=="hetero",1,0),
                      g6pd376.homo = ifelse(G6PD376=="homo",1,0),
                      g6pd376.hetero = ifelse(G6PD376=="hetero",1,0),
                      g6pd202.homo = ifelse(G6PD202=="homo",1,0),
                      g6pd202.hetero = ifelse(G6PD202=="hetero",1,0),
                      age2 = ifelse(agecat==2,1,0),
                      age3 = ifelse(agecat==3,1,0),
                      moib = ifelse(moi>1,1,0))

df <- dtrm[dtrm$qpcr>0 & is.na(dtrm$qpcr)==F,]

df <- df %>% mutate(totaltime = last(weeks))

df$weekscat <- factor(cut(df$weeks,
                          breaks=c(-0.00001,0,4,8,12,96)),
                      labels = c("0","(0,4]","(4,8]","(8,12]","(12,96]"))

ggplot(data=df, aes(x=weeks,y=qpcr,col=as.factor(infection_id)),alpha=0.2)+
  geom_point()+
  geom_line()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0.1)+
  geom_hline(yintercept = 0.01,linetype="dashed")+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
#ggsave("figures/pcr_long.pdf",device="pdf",dpi=300,width=8,height=6)


df = df %>% group_by(infection_id) %>% 
  mutate(initial=first(clinical),
         initial=ifelse(initial==0,"Asymptomatic","Symptomatic"),
         id2 = paste0(initial,": ",infection_id),
         gam = ifelse(Tgam_uL_adj>0 & is.na(Tgam_uL_adj)==F,"Positive gametocytes","Negative gametocytes"))

ggplot(data=df, aes(x=weeks,y=qpcr),alpha=0.2)+
  geom_point()+
  geom_line()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0.1)+
  geom_hline(yintercept = 0.01,linetype="dashed")+
  ylab("Parasite density/\U00B5L")+
  xlab("Time since detected incident infection (weeks)")+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  facet_wrap(. ~id2)+
  geom_point(data=df[df$gam=="Positive gametocytes",],aes(x=weeks,y=qpcr), col="red" )+
  theme(axis.title=element_text(size=18,face="bold"))
#ggsave("figures/pcr_long_individual.pdf",device="pdf",dpi=300,width=18,height=12)


df12 = df[df$totaltime>=12,]

df12 = df12 %>% group_by(infection_id) %>% mutate(eversymptomatic=ifelse(any(clinical==1),1,0), initial=ifelse(first(clinical)==1,"Symptomatic","Asymptomatic"),sumgams=sprintf("%.2f",sum(I(Tgam_uL_adj>0)/100,na.rm=T)),gamfrac=Tgam_uL_adj/qpcr,peakgamfrac=max(gamfrac,na.rm=TRUE),order=paste0(sumgams,sprintf("%.5f", round(peakgamfrac/1000,5)),sprintf("%.3f", round(totaltime/100,3)),sprintf("%.3f", round(infection_id/1000,3))))

ggplot()+
  geom_point(data=na.omit(df12%>%select(weeks,order,qpcr,infection_id)),aes(x=weeks,y=qpcr,col="Parasites"),alpha=0.7)+
  geom_line(data=na.omit(df12%>%select(weeks,order,qpcr,infection_id)),aes(x=weeks,y=qpcr,col="Parasites"),alpha=0.7)+
  geom_point(data=na.omit(df12%>%select(weeks,order,Tgam_uL_adj,infection_id)),aes(x=weeks,y=Tgam_uL_adj,col="Gametocytes"),alpha=0.7)+
  geom_line(data=na.omit(df12%>%select(weeks,order,Tgam_uL_adj,infection_id)),aes(x=weeks,y=Tgam_uL_adj,col="Gametocytes"),alpha=0.7)+
  theme(legend.position = "top")+
  geom_hline(yintercept = 0.1,alpha=0.6)+
  geom_hline(yintercept = 0.01,linetype="dashed",alpha=0.6)+
  ylab("Density/\U00B5L")+
  xlab("Time since detected incident infection (weeks)")+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  facet_wrap(. ~order)+
  theme(axis.title=element_text(size=18,face="bold"))+
  scale_color_viridis_d(name="",begin=0.3, end=0.9)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())
#ggsave("figures/pcr_and_gam_long_individual.pdf",device="pdf",dpi=300,width=12,height=8)


ggplot()+
  geom_point(data=na.omit(df12%>%select(weeks,order,Tgam_uL_adj,infection_id,qpcr)),aes(x=weeks,y=Tgam_uL_adj/qpcr,col="Gametocytes"),alpha=0.7)+
  geom_line(data=na.omit(df12%>%select(weeks,order,Tgam_uL_adj,infection_id,qpcr)),aes(x=weeks,y=Tgam_uL_adj/qpcr,col="Gametocytes"),alpha=0.7)+
  theme(legend.position = "none")+
  geom_hline(yintercept = 1,alpha=0.6)+
  #geom_hline(yintercept = 0.01,linetype="dashed")+
  ylab("Gametocyte fraction")+
  xlab("Time since detected incident infection (weeks)")+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  facet_wrap(. ~order)+
  theme(axis.title=element_text(size=18,face="bold"))+
  scale_color_viridis_d(name="",begin=0, end=0)+ 
  theme(strip.background = element_blank(),strip.text.x = element_blank())
#ggsave("figures/gamfrac_long_individual.pdf",device="pdf",dpi=300,width=12,height=8)



# Analysis ----------------------------------------------------------------


# What proportion of the data were asymptomatic/symptomatic initially?

df0 <- df[df$weeks==0,]
table(df0$clinical)
round(prop.table(table(df0$clinical)),3)*100

prop.test(c(88,16), c(104,104), p = NULL, alternative = "two.sided",
          correct = TRUE)


table(df0$moib,df0$clinical,useNA="ifany")
prop.table(table(df0$moib,df0$clinical),2)
fisher.test(table(df0$moib,df0$clinical))

table(I(df0$Tgam_uL_adj>0),df0$clinical,useNA="ifany")
prop.table(table(I(df0$Tgam_uL_adj>0),df0$clinical),2)
fisher.test(table(I(df0$Tgam_uL_adj>0),df0$clinical))


table(df0$HbS,df0$clinical,useNA="ifany")
prop.table(table(df0$HbS,df0$clinical),2)
fisher.test(table(df0$HbS,df0$clinical))


# From those that were initially asymptomatic, 
## what proportion became symptomatic through the duration of the infection?
dtasymsub = dtasym[dtasym$cumsum==1,] %>% mutate(symtime=first(weeks))
dtasymsub2 =dtasymsub[dtasymsub$infection_id %in% df0$infection_id[df0$clinical==0],]
length(unique(dtasymsub2$infection_id))
(12/88)*100

#Gmean at the moment of clinical symptoms
test1=t.test(I(log(dtasymsub2$qpcr[dtasymsub2$clinical==1])))
exp(test1$estimate)
exp(test1$conf.int)

#Gmean at initial visit
pd = df[df$infection_id %in% unique(dtasymsub2$infection_id) & df$weeks==0,]$qpcr
pg = df[df$infection_id %in% unique(dtasymsub2$infection_id) & df$weeks==0,]$Tgam_uL_adj
test1=t.test(I(log(pd)))
exp(test1$estimate)
exp(test1$conf.int)


## what proportion of the 16+12 = 28 clinical infections had subsequent visits with positive parasites?
dtasymsub = dtasymsub %>% mutate(qpcr = ifelse(is.na(qpcr),0,qpcr))
dtasymsub3 = dtasymsub[dtasymsub$qpcr>0,] 
t1 = table(dtasymsub3$infection_id,I(dtasymsub3$symtime-dtasymsub3$weeks==0))
sum(t1[,1]>0)
20/28

## What proportion of the 16 initially symptomatic had gametocytes and how long?
sum(I(dtasymsub3$Tgam_uL_adj[dtasymsub3$weeks==0]>0),na.rm=T)
1

dtasymsub3$weeks[dtasymsub3$infection_id %in% dtasymsub3$infection_id[dtasymsub3$weeks==0 & dtasymsub3$Tgam_uL_adj>0]]

# 1 and only at the moment of detection.

## What proportion of the 12 later symptomatic had gametocytes?
t2 = table(dtasymsub3[dtasymsub3$weeks!=0,]$infection_id,I(dtasymsub3[dtasymsub3$weeks!=0,]$Tgam_uL_adj>0))
sum(t2[,2]>0)
t2[,2]

# 5 had gametocytes, 1 only at the clinical visit, 3 had gametocytes in the visit following clinical inf and 1 had gametocytes in the next two visits after clinical infection
temp = dtasymsub3[dtasymsub3$infection_id %in% c(5,25,65,73) & dtasymsub3$Tgam_uL_adj>0,] 
temp$weeks-temp$symtime
# ids 5 65 73 has two (at the moment of clinical infection and then 1.86, 0.86 and 1.29 weeks later)
#25 had three (at the moment of clinical infection and 1.14 weeks and 5.14 weeks later)
#94 had one at the moment of clinical infection

#Histograms

ggarrange(
ggplot(data=df0)+
  geom_histogram(aes(x=totaltime,fill="All infections"),bins=23.92857)+
  theme_prism()+
  scale_x_continuous(breaks=seq(0,96,4))+
  xlab("")+
  scale_y_continuous(limits=c(0,80))+
  scale_fill_viridis_d()+
  theme(legend.position = "none")+
  ylab("Frequency\n(All infections)")+
  ggtitle(" "),

ggplot(data=df0[df0$clinical==0,])+
  geom_histogram(aes(x=totaltime,fill="Initially\nasymptomatic"),bins=23.92857)+
  theme_prism()+
  scale_x_continuous(breaks=seq(0,96,4))+
  xlab("Duration since detection of incident malaria (weeks)")+
  scale_y_continuous(limits=c(0,80))+
  scale_fill_viridis_d(begin=0.5)+
  theme(legend.position = "none")+
  ylab("Frequency\n(Initially asymptomatic)"),
nrow=2,
ncol = 1,
align = "hv",
labels=c("A","B"),
vjust=1)

#ggsave("figures/histogram_duration.pdf",device="pdf",width=12,heigh=6,dpi=600)

#xlsx::write.xlsx2(as.data.frame(df0 %>% mutate(initially_asymptomatic=clinical) %>% select(totaltime,initially_asymptomatic)),"data for figures/Figure1_histogram.xlsx",sheetName = "ALL",row.names = FALSE)


# Duration of infection for initially asymptomatic infections----------------------------

## Here we define the duration as the time from detection of incident malaria until either the last parasite positive visit or until symptomatic infections which are interrupted due to treatment. n=88.

df0$time_cat <- cut(df0$totaltime,
                       breaks=c(-0.00001,0,4,8,12,96))

table(df0$time_cat,df0$clinical)
round(prop.table(table(df0$time_cat,df0$clinical),2),3)*100


df0$agecat = factor(ifelse(df0$age2==1,2,ifelse(df0$age3==1,3,1)),labels=c("<5 years","5-15 years","16+ years"))
df0 = df0 %>% mutate(hbs=as.factor(ifelse(hbs.homo==1,"SS",ifelse(hbs.hetero==1,"AS","AA"))))
df0$hbs = relevel(df0$hbs,ref="AA")
df0$moib = as.factor(ifelse(df0$moi==1,"1",">1"))

#by age cat and hbs 
table(df0[df0$clinical==0,]$time_cat,df0[df0$clinical==0,]$agecat)
round(prop.table(table(df0[df0$clinical==0,]$time_cat,df0[df0$clinical==0,]$agecat),2),3)*100


table(df0[df0$clinical==0,]$time_cat,df0[df0$clinical==0,]$hbs)
round(prop.table(table(df0[df0$clinical==0,]$time_cat,df0[df0$clinical==0,]$hbs),2),3)*100

table(df0[df0$clinical==0,]$moib,df0[df0$clinical==0,]$time_cat)
prop.table(table(df0[df0$clinical==0,]$moib,df0[df0$clinical==0,]$time_cat),1)

median(df0$moi[df0$clinical==0],na.rm=T)
median(df0$moi[df0$clinical==1],na.rm=T)
range(df0$moi[df0$clinical==0],na.rm=T)
range(df0$moi[df0$clinical==1],na.rm=T)
wilcox.test(df0$moi~df0$clinical)

#t - tests
df0$totaltime2 = ifelse(df0$totaltime==0,0.1,df0$totaltime)
df0$logtottime = log(df0$totaltime2)
library(psych)
anova(lm(df0[df0$clinical==0,]$logtottime ~ df0[df0$clinical==0,]$agecat))

anova(lm(df0[df0$clinical==0,]$logtottime ~ df0[df0$clinical==0,]$hbs))

df0$clin = factor(df0$clinical,labels=c("Asymptomatic","Symptomatic"))
df0$totaltime2 = ifelse(df0$totaltime==0,0.1,df0$totaltime)

options(scipen=10000)
ggarrange(
ggplot(data=df0[df0$clinical==0,]) +
  geom_violin(aes(x="Initally asymptomatic", y=totaltime, fill="Overall"),width=1, alpha=0.7) +
  geom_jitter(aes(x="Initally asymptomatic", y=totaltime, fill="Overall"),height = 0,alpha=0.2,width=0.2)+
  scale_fill_viridis_d(begin=0, end=0.2) +
  #scale_y_continuous(breaks=seq(0,96,4))+
  theme_prism()+
  theme(
    legend.position="none",
  ) +
  ylab("Duration of infection (weeks)")+
  theme(axis.title.x = element_blank())+
  ggtitle("Overall"),

ggplot(data=df0[df0$clinical==0,]) +
  geom_violin(aes(x=agecat, y=totaltime, fill=agecat),width=1, alpha=0.7) +
  geom_jitter(aes(x=agecat, y=totaltime, fill=agecat),height = 0,alpha=0.2,width=0.2)+
  scale_fill_viridis_d(begin=0.4, end=0.6) +
  theme_prism()+
  theme(
    legend.position="none",
  ) +
  ylab("")+
  theme(axis.title.x = element_blank())+
  ggtitle("Age")+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()),

ggplot(data=na.omit(df0[df0$clinical==0,] %>% select(totaltime,hbs))) +
  geom_violin(aes(x=hbs, y=totaltime, fill=hbs),width=1, alpha=0.7) +
  geom_jitter(aes(x=hbs, y=totaltime, fill=hbs),height = 0,alpha=0.2,width=0.2)+
  scale_fill_viridis_d(begin=0.8, end=1) +
  theme_prism()+
  theme(
    legend.position="none",
  ) +
  ylab("")+
  theme(axis.title.x = element_blank())+
  ggtitle("HbS")+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()),

nrow=1)
#ggsave("figures/Descriptive_duration.pdf",device="pdf",width=12,heigh=6,dpi=600)

# durations of infection for those who remained asymptomatic
df0a <- df %>% group_by(infection_id) %>%
  mutate(sumclin=sum(clinical),sumclin=ifelse(sumclin>=1,1,0)) %>%
  filter(weeks==0 & sumclin==0)

table(I(df0a$totaltime==0))
prop.table(table(I(df0a$totaltime==0)))


# Parasite densities: immediately symptomatic and those that were  --------


  ggplot() +
    geom_violin(data=df %>% filter(infection_id %in% df0$infection_id[df0$clinical==0]),
                aes(x="Initially asymptomatic ", y=qpcr, fill="Initially asymptomatic "),width=1, alpha=0.7) +
    geom_jitter(data=df %>% filter(infection_id %in% df0$infection_id[df0$clinical==0]),
                aes(x="Initially asymptomatic ", y=qpcr, fill="Initially asymptomatic "),height = 0,alpha=0.2,width=0.2)+
    geom_violin(data=df %>% filter(infection_id %in% df0$infection_id[df0$clinical==1]),
                aes(x="Initially\nsymptomatic", y=qpcr, fill="Initially\nsymptomatic"),width=1, alpha=0.7) +
    geom_jitter(data=df %>% filter(infection_id %in% df0$infection_id[df0$clinical==1]),
                aes(x="Initially\nsymptomatic", y=qpcr, fill="Initially\nsymptomatic"),height = 0,alpha=0.2,width=0.2)+
  geom_violin(data=df0 %>% filter(clinical==0 & totaltime==0),
              aes(x="Initially asymptomatic \nonly detected once", y=qpcr, fill="Initially asymptomatic \nonly detected once"),width=1, alpha=0.7) +
  geom_jitter(data=df0 %>% filter(clinical==0 & totaltime==0),
              aes(x="Initially asymptomatic \nonly detected once", y=qpcr, fill="Initially asymptomatic \nonly detected once"),height = 0,alpha=0.2,width=0.2)+
    scale_fill_viridis_d() +
    theme_prism()+
    theme(
      legend.position="none",
    ) +
    ylab("Parasite density/\U00B5L")+
    theme(axis.title.x = element_blank())+
    scale_y_log10(breaks=c(0,0.01,0.1,1,10,100,1000,10000))
#ggsave("figures/Descriptive_parasitedens.pdf",device="pdf",width=10,height=6,dpi=600)


# What is the Gmean parasite density amongst at the first visit by clinical status?
table(df0$clinical)
m1=lm(data=df0, I(log(qpcr)) ~ clinical)
newd = data.frame(clinical=c(0,1),infection_id="1")
pred = predict(m1,newdata = newd,se.fit=TRUE)
newd$pred = round(exp(pred$fit),2)
newd$lci = round(exp(pred$fit-1.96*pred$se.fit),2)
newd$uci = round(exp(pred$fit+1.96*pred$se.fit),2)
newd
summary(m1)

test1=t.test(I(log(df0$qpcr[df0$clinical==0])))
test2=t.test(I(log(df0$qpcr[df0$clinical==1])))
exp(test1$estimate)
exp(test2$estimate)
exp(test1$conf.int)
exp(test2$conf.int)




# What is the Gmean parasite density amongst those that had only visit only by clinical status?

df.pcr2 = df0 %>% filter(totaltime==0)
table(df.pcr2$clinical)
m1=lm(data=df.pcr2, I(log(qpcr)) ~ clinical)
newd = data.frame(clinical=c(0,1),infection_id="1")
pred = predict(m1,newdata = newd,se.fit=TRUE)
newd$pred = round(exp(pred$fit),2)
newd$lci = round(exp(pred$fit-1.96*pred$se.fit),2)
newd$uci = round(exp(pred$fit+1.96*pred$se.fit),2)
newd
summary(m1)

test1=t.test(I(log(df.pcr2$qpcr[df.pcr2$clinical==0])))
test2=t.test(I(log(df.pcr2$qpcr[df.pcr2$clinical==1])))
exp(test1$estimate)
exp(test2$estimate)
exp(test1$conf.int)
exp(test2$conf.int)


# Gmean parasite density at initial visit for missing and non-missing moi.
exp(mean(log(df0[is.na(df0$moi)==T,]$qpcr)))
exp(mean(log(df0[is.na(df0$moi)==F,]$qpcr)))


m1=lm(data=df0, I(log(qpcr)) ~ I(is.na(moi)))
summary(m1)

test1=t.test(I(log(df0[is.na(df0$moi)==T,]$qpcr)))
test2=t.test(I(log(df0[is.na(df0$moi)==F,]$qpcr)))
exp(test1$estimate)
exp(test2$estimate)
exp(test1$conf.int)
exp(test2$conf.int)


# Peak parasite density

peakpar = df %>% 
  group_by(infection_id) %>% 
  mutate(parpeak = max(qpcr,na.rm=TRUE)) %>%
  filter(weeks==0)
m1=lm(data=peakpar, I(log(parpeak)) ~ clinical)
newd = data.frame(clinical=c(0,1),infection_id="1")
pred = predict(m1,newdata = newd,se.fit=TRUE)
newd$pred = round(exp(pred$fit),2)
newd$lci = round(exp(pred$fit-1.96*pred$se.fit),2)
newd$uci = round(exp(pred$fit+1.96*pred$se.fit),2)
newd
summary(m1)

test1=t.test(I(log(peakpar$parpeak[peakpar$clinical==0])))
test2=t.test(I(log(peakpar$parpeak[peakpar$clinical==1])))
exp(test1$estimate)
exp(test2$estimate)
exp(test1$conf.int)
exp(test2$conf.int)


#Parasite density of initial visit for those who carried infection for more than one visit
dfnew = df %>% group_by(infection_id) %>% mutate(weekscat2=last(weekscat),weekscat2=ifelse(weekscat2==0,"Once","More")) %>% filter(weekscat=="0" & initial=="Asymptomatic") 

table(dfnew$weekscat2)
m1=lm(data=dfnew, I(log(qpcr)) ~ weekscat2)
summary(m1)

# Gametocyte initiation ---------------------------------------------------

# What % of initially asymptomatic incident infections had gametocytes at the moment of first detection. 

df = df %>% mutate(gambin = ifelse(Tgam_uL_adj==0 | is.na(Tgam_uL_adj),0,1), gambin=cumsum(gambin),gambin=ifelse(gambin>=1,"pos","neg"))

table(df$gambin[df$weeks==0], df$clinical[df$weeks==0], useNA = "ifany")
round(prop.table(table(df$gambin[df$weeks==0], df$clinical[df$weeks==0], useNA = "ifany"),2),3)*100

#Amongst the initially asymptomatic infections that only had one detection time point (n=49), what % had gametocytes?

df0 = df0 %>% mutate(gambin = ifelse(Tgam_uL_adj>0,"pos","neg"))
table(df0$gambin[df0$totaltime==0],df0$clinical[df0$totaltime==0], useNA = "ifany")
prop.table(table(df0$gambin[df0$totaltime==0],df0$clinical[df0$totaltime==0], useNA = "ifany"),2)



# Long infections ---------------------------------------------------------

#How many long infections were there? And what proportion had gametocytes at the initial visit?

dflong = df[df$totaltime>=12,]
table(dflong$gambin[dflong$weeks==0], useNA = "ifany")
prop.table(table(dflong$gambin[dflong$weeks==0], useNA = "ifany"))

table(dflong$gambin[dflong$weeks==0], dflong$clinical[dflong$weeks==0], useNA = "ifany")
prop.table(table(dflong$gambin[dflong$weeks==0], dflong$clinical[dflong$weeks==0], useNA = "ifany"))

#When did these long infections cumulatively get gametocytes?
table(dflong$gambin[dflong$weeks==0])
table(dflong$gambin[dflong$weeks<=4])
table(dflong$gambin[dflong$weeks<=8])
table(dflong$gambin[dflong$weeks<=12])

#who had long infections?
subid=unique(dflong$infection_id)

#Who had initially asymptomatic infections of any length?
subdf3 =  df[df$infection_id %in% df0$infection_id[df0$clinical==0],]
subid3=df0$infection_id[df0$clinical==0]


# Non-parametric analysis -------------------------------------------------
##Asymptomatic at start and long infections
pedsub1 = faildt[faildt$infection_id %in% subid,] 
pedsub1$state = factor(pedsub1$state,levels=c(0,1,2),labels=c("None","gam","res"))
fit <- survfit(Surv(time/7, state, type="mstate") ~ 1, data = pedsub1)
summary(fit)
summary(fit)$table
fit$states = c("Event\nfree","Gametocyte\ninitiation","Malaria\nResolution")
plot1 = ggcompetingrisks(fit,cumevents=TRUE)+ theme_cowplot() + scale_fill_jco() + facet_grid(. ~ "Long infections")+labs(title=NULL)+theme_prism()+theme(legend.position = "bottom")+scale_x_continuous(breaks=seq(0,20,4))
fitdt=summary(fit)
plot1

##Asymptomatic at start and any duration infections
pedsub3 = faildt[faildt$infection_id %in% as.character(subid3),] 
pedsub3$state = factor(pedsub3$state,levels=c(0,1,2),labels=c("None","gam","res"))
table(pedsub3$state)
fit3 <- survfit(Surv(time/7, state, type="mstate") ~ 1, data = pedsub3)
summary(fit3)
summary(fit3)$table
fit3$states = c("Event\nfree","Gametocyte\ninitiation","Malaria\nResolution")
plot3 = ggcompetingrisks(fit3)+ theme_cowplot() + scale_fill_jco() + facet_grid(. ~ "Initially asymptomatic & any duration")+labs(title=NULL)+theme_prism()+theme(legend.position = "bottom")+scale_x_continuous(breaks=seq(0,20,4))
plot3

##ALL infections including symptomatic at start
faildt$state = factor(faildt$state,levels=c(0,1,2),labels=c("None","gam","res"))
fit4 <- survfit(Surv(time/7, state, type="mstate") ~ 1, data = faildt)
summary(fit4)
summary(fit4)$table
fit4$states = c("Event\nfree","Gametocyte\ninitiation","Malaria\nResolution")
plot4 = ggcompetingrisks(fit4)+ theme_cowplot() + scale_fill_jco() + facet_grid(. ~ "All infections")+labs(title=NULL)+theme_prism()+theme(legend.position = "bottom")+scale_x_continuous(breaks=seq(0,20,4))

annotate_figure(ggarrange(
  plot4+rremove("ylab")+rremove("xlab"),
  plot3+rremove("ylab")+rremove("xlab"),
  plot1+rremove("ylab")+rremove("xlab"),
  common.legend = TRUE,
  legend = "top",
  nrow=1,
  labels=c("A","B","C"),
  align = "hv"
), 
bottom=text_grob("Time since detected incident infection (weeks)",face = "bold", size = 14),
left=text_grob("Proportion of infections",face = "bold", size = 14, rot = 90)
)

#ggsave("figures/gametocyte_incidence.pdf",device="pdf",width=12,heigh=6,dpi=600)

plotdt1 = data.frame(
time = fitdt$time,
n.risk = fitdt$n.risk[,1],
n.gam = fitdt$n.event[,2],
n.resolve = fitdt$n.event[,3],
p.eventfree = fitdt$pstate[,1],
p.gam = fitdt$pstate[,2],
p.resolve = fitdt$pstate[,3]
)

plotdt2 = data.frame(
  time = summary(fit3)$time,
  n.risk = summary(fit3)$n.risk[,1],
  n.gam = summary(fit3)$n.event[,2],
  n.resolve = summary(fit3)$n.event[,3],
  p.eventfree = summary(fit3)$pstate[,1],
  p.gam = summary(fit3)$pstate[,2],
  p.resolve = summary(fit3)$pstate[,3]
)


plotdt3 = data.frame(
  time = summary(fit4)$time,
  n.risk = summary(fit4)$n.risk[,1],
  n.gam = summary(fit4)$n.event[,2],
  n.resolve = summary(fit4)$n.event[,3],
  p.eventfree = summary(fit4)$pstate[,1],
  p.gam = summary(fit4)$pstate[,2],
  p.resolve = summary(fit4)$pstate[,3]
)



#xlsx::write.xlsx2(plotdt3,"data for figures/Figure2_comprisks.xlsx",sheetName = "ALL infections",row.names = FALSE)
#xlsx::write.xlsx2(plotdt2,"data for figures/Figure2_comprisks.xlsx",sheetName = "Asymptomatic at start and any duration infections",row.names = FALSE,append=TRUE)
#xlsx::write.xlsx2(plotdt1,"data for figures/Figure2_comprisks.xlsx",sheetName = "Long infections",row.names = FALSE,append=TRUE)


#xlsx::write.xlsx2(as.data.frame(faildt %>% select(time,state)),"data for figures/Figure2_comprisksdata.xlsx",sheetName = "ALL infections",row.names = FALSE)
#xlsx::write.xlsx2(as.data.frame(pedsub3 %>% select(time,state)),"data for figures/Figure2_comprisksdata.xlsx",sheetName = "Asymptomatic at start and any duration infections",row.names = FALSE,append=TRUE)
#xlsx::write.xlsx2(as.data.frame(pedsub1 %>% select(time,state)),"data for figures/Figure2_comprisksdata.xlsx",sheetName = "Long infections",row.names = FALSE,append=TRUE)


# now stratified by HbS hetero --------------------------------------------

##Asymptomatic at start and long infections

pedsub1 = pedsub1 %>% mutate(hbs=factor(ifelse(hbs.homo==1,"SS",ifelse(hbs.hetero==1,"AS","AA")),levels=c("AA","AS","SS")))
pedsub1[(nrow(pedsub1)+1),]=NA
pedsub1$time[(nrow(pedsub1))]=0
pedsub1$state[(nrow(pedsub1))]="res"
pedsub1$hbs[(nrow(pedsub1))]="SS"
fit <- survfit(Surv(time/7, state, type="mstate") ~ hbs, data = pedsub1)
summary(fit)
summary(fit)$table
fit$states = c("Event\nfree","Gametocyte\ninitiation","Malaria\nResolution")
plot1b = ggcompetingrisks(fit,cumevents=TRUE)+ theme_cowplot() + scale_fill_jco()+ggtitle("Initially asymptomatic & long infections")+theme_prism()+theme(legend.position = "bottom")+scale_x_continuous(breaks=seq(0,20,4))

plot1b



##Asymptomatic at start and any duration infections
pedsub3 = pedsub3 %>% mutate(hbs=factor(ifelse(hbs.homo==1,"SS",ifelse(hbs.hetero==1,"AS","AA")),levels=c("AA","AS","SS")))
fit3 <- survfit(Surv(time/7, state, type="mstate") ~ hbs, data = pedsub3)
summary(fit3)
summary(fit3)$table
fit3$states = c("Event\nfree","Gametocyte\ninitiation","Malaria\nResolution")
plot3b = ggcompetingrisks(fit3,cumevents=TRUE)+ theme_cowplot() + scale_fill_jco()+ggtitle("Initially asymptomatic & any duration")+theme_prism()+theme(legend.position = "bottom")+scale_x_continuous(breaks=seq(0,20,4))
plot3b

##ALL infections
faildt = faildt %>% mutate(hbs=factor(ifelse(hbs.homo==1,"SS",ifelse(hbs.hetero==1,"AS","AA")),levels=c("AA","AS","SS")))
fit4 <- survfit(Surv(time/7, state, type="mstate") ~ hbs, data = faildt)
summary(fit4)
summary(fit4)$table
fit4$states = c("Event\nfree","Gametocyte\ninitiation","Malaria\nResolution")
plot4b = ggcompetingrisks(fit4,cumevents=TRUE)+ theme_cowplot() + scale_fill_jco()+ggtitle("All infections")+theme_prism()+theme(legend.position = "bottom")+scale_x_continuous(breaks=seq(0,20,4))
plot4b


annotate_figure(ggarrange(
  plot4b+rremove("ylab")+rremove("xlab"),
  plot3b+rremove("ylab")+rremove("xlab"),
  plot1b+rremove("ylab")+rremove("xlab"),
  common.legend = TRUE,
  legend = "top",
  nrow=3,
  align = "hv"
), 
bottom=text_grob("Time since detected incident infection (weeks)",face = "bold", size = 14),
left=text_grob("Proportion of infections",face = "bold", size = 14, rot = 90)
)

#ggsave("figures/gametocyte_incidence_by_HBS.pdf",device="pdf",width=12,heigh=12,dpi=600)


# Describe gametocyte densities -------------------------------------------


# Gametocyte initiation ---------------------------------------------------

# What proportion of initially asymptomatic ever had gametocytes
peakgam = df %>% 
  group_by(infection_id) %>% 
  mutate(gampeak = max(Tgam_uL_adj,na.rm=TRUE)) %>%
  filter(weeks==0)
gams = peakgam$gampeak[peakgam$gampeak>0 & peakgam$clinical==0]
length(gams)

# What was the mean peak Geometric density amongst initially asymptomatic that ever had gametocytes
quantile(gams)
DescTools::Gmean(gams,conf.level = .95)

## BY HbS
peakgam$HbS = relevel(factor(peakgam$HbS),ref="WT")

m1=lm(data=peakgam[peakgam$gampeak>0,], I(log(gampeak)) ~ HbS)
summary(m1)

test1=t.test(I(log(peakgam[peakgam$gampeak>0,]$gampeak[peakgam[peakgam$gampeak>0,]$HbS=="hetero"])))
test2=t.test(I(log(peakgam[peakgam$gampeak>0,]$gampeak[peakgam[peakgam$gampeak>0,]$HbS=="WT"])))
exp(test1$estimate)
exp(test2$estimate)
exp(test1$conf.int)
exp(test2$conf.int)

quantile(peakgam[peakgam$gampeak>0,]$gampeak[peakgam[peakgam$gampeak>0,]$HbS=="WT"])
quantile(peakgam[peakgam$gampeak>0,]$gampeak[peakgam[peakgam$gampeak>0,]$HbS=="hetero"])

## When was this peak?
peakgam = df %>% 
  group_by(infection_id) %>% mutate(iniclinical=first(clinical)) %>% 
  mutate(gampeak = max(Tgam_uL_adj,na.rm=TRUE)) %>%
  filter(Tgam_uL_adj==gampeak & gampeak>0)
table(peakgam$weekscat[peakgam$iniclinical==0])
prop.table(table(peakgam$weekscat[peakgam$iniclinical==0]))

# What was the mean peak Geometric density amongst initially asymptomatic excluding all symptomatic infections
# two others from a later symptomatic visit
table(peakgam$weekscat[peakgam$iniclinical==0 & peakgam$clinical==0])
prop.table(table(peakgam$weekscat[peakgam$iniclinical==0 & peakgam$clinical==0]))

peakgam$gampeak[peakgam$iniclinical==0 & peakgam$clinical==0]
quantile(peakgam$gampeak[peakgam$iniclinical==0 & peakgam$clinical==0])
DescTools::Gmean(peakgam$gampeak[peakgam$iniclinical==0 & peakgam$clinical==0],conf.level = .95)

## When was this peak?

table(peakgam$weekscat[peakgam$iniclinical==0 & peakgam$clinical==0])
prop.table(table(peakgam$weekscat[peakgam$iniclinical==0 & peakgam$clinical==0]))


### One peak density was at an initially symptomatic visit (week=0) with density=0.80, the other two peak densities came from a later symptomatic visit (between week=0 and week=4) with gametocyte densities 195 and 26 respectively.

# What were the gametocyte densities for initially asymptomatic incident infections that had gametocytes at the moment of first detection. 

df$Tgam_uL_adj[df$weeks==0 & df$clinical==0 & df$gambin=="pos"]

## quantiles and geometric mean
quantile(df$Tgam_uL_adj[df$weeks==0 & df$clinical==0 & df$gambin=="pos"])
DescTools::Gmean(df$Tgam_uL_adj[df$weeks==0 & df$clinical==0 & df$gambin=="pos"],conf.level = .95)

# What was the gametocyte densities for initially symptomatic incident infections that had gametocytes at the moment of first detection. 
df$Tgam_uL_adj[df$weeks==0 & df$clinical==1 & df$gambin=="pos"]


#What were the gametocyte densities for initially asymptomatic incident infections that had gametocytes and were only observed at the moment of infection?
df0$Tgam_uL_adj[df0$totaltime==0 & df0$clinical==0 & df0$gambin=="pos"]


# gametocyte density distributions across duration of infections for initially asymptomatic infections.

dfgam1 = df %>% mutate(iniclinical=first(clinical)) %>% filter(Tgam_uL_adj>0 & iniclinical==0)
dfgam1$loggam = log(dfgam1$Tgam_uL_adj)
hist(dfgam1$loggam)

table(dfgam1$weekscat)

DescTools::Gmean(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="0"],conf.level = .95)
quantile(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="0"])
DescTools::Gmean(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(0,4]"],conf.level = .95)
quantile(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(0,4]"])
DescTools::Gmean(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(4,8]"],conf.level = .95)
quantile(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(4,8]"])
DescTools::Gmean(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(8,12]"],conf.level = .95)
quantile(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(8,12]"])
DescTools::Gmean(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(12,96]"],conf.level = .95)
quantile(dfgam1$Tgam_uL_adj[dfgam1$weekscat=="(12,96]"])


#plot 

ggarrange(
ggplot() +
  geom_violin(data=dfgam1,
              aes(x=weekscat, y=Tgam_uL_adj, fill=weekscat),width=1, alpha=0.7) +
  geom_jitter(data=dfgam1,
              aes(x=weekscat, y=Tgam_uL_adj, fill=weekscat),height = 0,alpha=0.2,width=0.2)+
  scale_fill_viridis_d(end=0.9) +
  theme_prism()+
  theme(
    legend.position="none",
  ) +
  ylab("Gametocyte density/\U00B5L")+
  theme(axis.title.x = element_blank())+
  scale_y_log10(breaks=c(0,0.01,0.1,1,10,100,1000,10000)),
ggplot() +
  geom_violin(data=dfgam1,
              aes(x="Overall", y=Tgam_uL_adj, fill="Overall"),width=1, alpha=0.7) +
  geom_jitter(data=dfgam1,
              aes(x="Overall", y=Tgam_uL_adj, fill="Overall"),height = 0,alpha=0.2,width=0.2)+
  scale_fill_viridis_d(begin=1, end=1) +
  theme_prism()+
  theme(
    legend.position="none",
  ) +
  ylab("Gametocyte density/\U00B5L")+
  theme(axis.title.x = element_blank())+
  scale_y_log10(breaks=c(0,0.01,0.1,1,10,100,1000,10000))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank()),
nrow=1
, widths = c(1, 0.25))



#ggsave("figures/gametocyte_density_initially_asymptomatic.pdf",device="pdf",width=10,height=6,dpi=600)


#xlsx::write.xlsx2(as.data.frame(dfgam1 %>% select(weekscat,Tgam_uL_adj)),"data for figures/Figure4_gamdensity.xlsx",sheetName = "ALL",row.names = FALSE)





