# Code for time-varying effects of recurrent and multiple miscarriage on time to miscarriage using PAMM

library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


library(pammtools)
library(mgcv)


#d= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/main_data/QIVF-QC-own-oocyte.txt')
d= fread(snakemake@input[[1]])
d$cat_prev_misc= factor(with(d, ifelse(is.na(prev_misc), NA, ifelse(prev_misc> 3, 3, prev_misc))))

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age),
           !is.na(misc))

set.seed(1234)
df = df %>% group_by(lopnr)  %>% sample_n(1) %>% ungroup()

#df1= filter(df, misc== 1)
df1= df


ped = as_ped(df1, Surv(gest_duration, misc) ~ cat_prev_misc, id = "id", cut=c(16, seq(20, 154, by=7)))

mod.tv.age = bam(ped_status ~ ti(tend,bs='cr',k=11)  + 
                   cat_prev_misc + ti(tend, by=as.ordered(cat_prev_misc),
bs='cr'), data=ped, offset=offset, family=poisson())

summary(mod.tv.age)

mod.tv.age.ptable= summary(mod.tv.age)$p.table
mod.tv.age.ptable$outcome= 'Previous miscarriages'

mod.tv.age.stable= summary(mod.tv.age)$s.table
mod.tv.age.stable$outcome= 'Previous miscarriages'

pred_df = make_newdata(ped, tend=unique(tend), cat_prev_misc=cat_prev_misc) %>%
  add_term(mod.tv.age, term="cat_prev_misc", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

ylimits= if (snakemake@wildcards[['oocytes']]== 'own-oocyte') c(-1.1, 0.9) else c(-3.1, 2.9)

p1= ggplot(pred_df, aes(x=(tmid)/7, y=fit, col=cat_prev_misc, fill=cat_prev_misc)) +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_manual(values=colorBlindBlack8[c(1,6,2, 8)], name="Previous miscarriages") +
  scale_fill_manual(values=colorBlindBlack8[c(1,6,2, 8)], name="Previous miscarriages") +
  scale_x_continuous(breaks=seq(0, 42, by=2), expand=c(0,0)) +
  ylim(ylimits) + 
  theme_classic(base_size= 10) + 
  xlab("Gestational age, weeks") + 
  ylab("Log hazard ratio") +
  theme(legend.position= "bottom", 
        legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), 
        legend.title = element_text(size=10))


ggsave(snakemake@output[[1]], p1, width= 88, height= 80, units= 'mm')


fwrite(data.frame(mod.tv.age.ptable), snakemake@output[[2]], sep= '\t')

fwrite(data.frame(mod.tv.age.stable), snakemake@output[[3]], sep= '\t')

fwrite(df1, snakemake@output[[4]], sep= '\t')
