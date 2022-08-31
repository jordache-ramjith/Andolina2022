library(readxl)
dt <- read_excel("Cleaned data.xlsx")

dt = dt %>% 
  group_by(infection_id) %>%
  mutate(cumcause=cumsum(cause),cumsum=cumsum(clinical),cumsum=ifelse(cumsum>1,1,cumsum))%>% 
  filter(!(clinical==0 & cumsum==1)) %>%
  mutate(revorder = n():1,
         initial = as.factor(ifelse(first(clinical)==1,"Symptomatic","Asymptomatic"))) %>% filter(revorder==1) %>% mutate(cleared = ifelse(qpcr==0 | clinical==1,1,0), cleared=ifelse(is.na(cleared)==T,0,cleared))


table(dt$cleared)
prop.table(table(dt$cleared))

table(dt$cleared,dt$initial)
prop.table(table(dt$cleared,dt$initial),2)


dts1 = dt %>% mutate(whichdata = "ALL")
dts2 = dt[dt$initial=="Asymptomatic",] %>% mutate(whichdata = "Initially \nAsymptomatic")

dtcomb = rbind(dts1,dts2)

require("survival")
fit <- survfit(Surv(weeks, cleared) ~ whichdata, data = dtcomb)
summary(fit)


library(ggplot2)
library(ggprism)
library(survminer)

ggsurvplot(fit, data = dtcomb,
           #surv.median.line = "hv", # Add medians survival
           break.time.by=4,
           surv.scale="percent",
           xlim=c(0,20),
           # Change legends: title & labels
           legend.title = "",
           legend.labs = c("All", "Initially \nAsymptomatic"),
           # Add p-value and tervals
           #pval = TRUE,
           
           conf.int = TRUE,
           # Add risk table
           risk.table = FALSE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_prism(), # Change ggplot2 theme,
           xlab="Time to clearance (weeks)",
           ylab="Proportion remaining infected"
)

ggsave("figures/KM_duration.pdf",device="pdf",dpi=300,width=8,height=6)


require("survival")
fit <- survfit(Surv(weeks, cleared) ~ HbS, data = dts2)
summary(fit)

ggsurvplot(fit, data = dts2,
           #surv.median.line = "hv", # Add medians survival
           break.time.by=4,
           surv.scale="percent",
           xlim=c(0,20),
           # Change legends: title & labels
           legend.title = "",
           #legend.labs = c("All", "Initially \nAsymptomatic"),
           # Add p-value and tervals
           pval = TRUE,
           
           conf.int = TRUE,
           # Add risk table
           risk.table = FALSE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           #palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_prism(), # Change ggplot2 theme,
           xlab="Time to clearance (weeks)",
           ylab="Proportion remaining infected"
)
