# Code for time-varying effects of maternal age on time to miscarriage using PAMM

library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


library(pammtools)
library(mgcv)

d= fread(snakemake@input[[1]])

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age),
           !is.na(misc))

set.seed(1234)
df = df %>% group_by(lopnr)  %>% sample_n(1) %>% ungroup()

df1= filter(df, misc== 1) %>% mutate(maternal_tertiles= as.factor(ntile(maternal_age, 3)))

ped = as_ped(df1, Surv(gest_duration, misc) ~ maternal_tertiles, id = "id", cut=c(16,seq(20, 154, by=7)))

mod.tv.age = bam(ped_status ~ ti(tend,bs='cr',k=11)  + 
                   maternal_tertiles + ti(tend, by=as.ordered(maternal_tertiles),
bs='cr'), data=ped, offset=offset, family=poisson())

summary(mod.tv.age)

mod.tv.age.ptable= summary(mod.tv.age)$p.table
mod.tv.age.ptable$outcome= 'Maternal tertiles'

mod.tv.age.stable= summary(mod.tv.age)$s.table
mod.tv.age.stable$outcome= 'Maternal tertiles'

pred_df = make_newdata(ped, tend=unique(tend), maternal_tertiles=factor(maternal_tertiles)) %>%
  add_term(mod.tv.age, term="maternal_tertiles", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

ylimits= if (snakemake@wildcards[['oocytes']]== 'own-oocyte') c(-1.1, 0.9) else c(-3.1, 2.9)

scaleFUN= function(x) sprintf("%.1f", x)

p1= ggplot(pred_df, aes(x=(tmid)/7, y=fit, col=maternal_tertiles, fill=maternal_tertiles)) +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_fill_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_x_continuous(breaks=seq(0, 42, by=2), expand=c(0,0)) +
  ylim(ylimits) + 
  theme_classic(base_size= 10) + 
  xlab("Gestational age, weeks") + 
  ylab("Log hazard ratio") +
  theme(legend.position= "bottom", 
        legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), 
        legend.title = element_text(size=10)) +
scale_y_continuous(labels=scaleFUN)

ggsave(snakemake@output[[1]], p1, width= 120, height= 80, units= 'mm')

#ped = as_ped(df1, Surv(gest_duration, misc) ~ maternal_age, id = "id", cut=c(16,seq(20, 154, by=7)))

#mod.tv.age = bam(ped_status ~ ti(tend,bs='cr',k=11)  +
#                   maternal_age + ti(tend, by=as.ordered(maternal_age),
#bs='cr'), data=ped, offset=offset, family=poisson())

#summary(mod.tv.age)

#mod.tv.age.ptable2= summary(mod.tv.age)$p.table
#mod.tv.age.ptable2$outcome= 'Maternal age'

#mod.tv.age.stable2= summary(mod.tv.age)$s.table
#mod.tv.age.stable2$outcome= 'Maternal age'

#tv.age= rbind(mod.tv.age.ptable2, mod.tv.age.ptable)
#tv.age.s= rbind(mod.tv.age.stable2, mod.tv.age.stable)

fwrite(mod.tv.age.ptable, snakemake@output[[2]], sep= '\t')

fwrite(mod.tv.age.stable, snakemake@output[[3]], sep= '\t')

fwrite(df1, snakemake@output[[4]], sep= '\t')
