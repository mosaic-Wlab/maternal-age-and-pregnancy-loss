# Code for time-varying effects of maternal age on time to miscarriage using PAMM in male infertility

library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(cowplot)
library(showtext)
options(warn=-1)

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ICD_COL_NAMES=c("DIA1","DIA2","DIA3","DIA4","DIA5","DIA6","DIA7","DIA8","DIA9","DIA10","DIA11","DIA12","DIA13","DIA14","DIA15","DIA16","DIA17","DIA18","DIA19","DIA20","DIA21","DIA22","DIA23","DIA24","DIA25","DIA26","DIA27","DIA28","DIA29","DIA30")

library(pammtools)
library(mgcv)

d= fread(snakemake@input[[1]])

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age),
           !is.na(misc))

outpatient= fread(snakemake@input[[2]])
outpatient= filter(outpatient, AR> 2006)
outpatient= select(outpatient, lopnr, AR, all_of(ICD_COL_NAMES))
INFERTILITY_CODES=c('N97','N970','N971','N972','N973','N978','N979')

MALE_INFERTILITY_CODES=c('N974')

outpatient= mutate(outpatient, male_infertility= case_when((if_any(ICD_COL_NAMES, ~ . %in% MALE_INFERTILITY_CODES)) ~ TRUE, TRUE ~ FALSE), female_infertility=  case_when((if_any(ICD_COL_NAMES, ~ . %in% INFERTILITY_CODES)) ~ TRUE, TRUE ~ FALSE))
outpatient= filter(outpatient, male_infertility, !female_infertility)

male_inf_lopnr= unique(pull(outpatient, lopnr))

df= filter(df, lopnr %in% male_inf_lopnr)

set.seed(1234)
df = df %>% group_by(lopnr)  %>% sample_n(1) %>% ungroup()

df1= filter(df, misc== 1) %>% mutate(maternal_tertiles= as.factor(ntile(maternal_age, 3)))

if (nrow(df1)< 100) {

file.create(snakemake@output[[1]])
file.create(snakemake@output[[2]])
file.create(snakemake@output[[3]])
file.create(snakemake@output[[4]])
file.create(snakemake@output[[5]])
file.create(snakemake@output[[6]])

} else {

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

p1= ggplot(pred_df, aes(x=(tmid)/7, y=fit, col=maternal_tertiles, fill=maternal_tertiles)) +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_fill_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_x_continuous(breaks=seq(0, 42, by=2), expand=c(0,0)) +
  ylim(c(-1.1, 0.9)) + 
  theme_classic(base_size= 10) + 
  xlab("Gestational age, weeks") + 
  ylab("Log hazard ratio") +
  theme(legend.position= "bottom", 
        legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), 
        legend.title = element_text(size=10))

ggsave(snakemake@output[[1]], p1, width= 120, height= 80, units= 'mm')

df$early_loss= with(df, ifelse(misc == 1 & gest_duration< 7*6, 'Early loss', ifelse(misc==1 & gest_duration>= 7*6, 'Late loss', ifelse(gest_duration>154 & all_loss== 0, 'Live born', NA))))

p2= ggplot(filter(df, !is.na(early_loss)), aes(early_loss, maternal_age, colour= early_loss)) +
geom_boxplot() +
theme_cowplot(font_size= 10) +
theme(legend= element_blank(),
        legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(),
        legend.title = element_text(size=10))

ggsave(snakemake@output[[5]], p2, width= 88, height= 80, units= 'mm')

print(summary(lm(maternal_age ~ early_loss, df)))
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

df1= mutate(df1, maternal_tertiles= ntile(x= maternal_age, 3))

df1$maternal_tertiles= as.character(df1$maternal_tertiles)

m_tertiles= survfit2(Surv(gest_duration, misc)~ maternal_tertiles, filter(df1, misc==1))

p1= ggsurvfit(m_tertiles) +
theme_cowplot(font_size= 10) +
# add_confidence_interval() +
scale_color_manual(values = colorBlindBlack8[c(2,6,1)]) +
  scale_fill_manual(values = colorBlindBlack8[c(2,6,1)]) +
#  add_risktable() +
  add_quantile(color = "gray50", linewidth = 0.2) +
#  scale_ggsurvfit() +
xlab("Duration of gestation, days") +
ylab("Percentage survival") +
theme(legend.title = element_blank(),
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        axis.ticks.x= element_line(color= 'black', size= 0.2),
        legend.position="bottom")

ggsave(snakemake@output[[6]], p1, width= 88, height= 60, units= 'mm')


}
