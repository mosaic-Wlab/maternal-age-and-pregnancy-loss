# Figures for maternal age and time-to-miscarriage

library(data.table)
library("dplyr")
library("ggplot2")
library(survival)
library(ggsurvfit)
library(cowplot)
library(showtext)
options(warn=-1)
library(gtsummary)
library(broom)

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#d= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/main_data/QIVF-QC-own-oocyte.txt')

d= fread(snakemake@input[[1]])
d$implantation_failure= as.numeric(d$Resultfetus1=='')
d$early_miscarriage= ifelse(is.na(d$misc) | is.na(d$gest_duration), NA, ifelse(d$misc== 1 & d$gest_duration< 7*10, 1, 0))
d$cat_prev_misc= factor(with(d, ifelse(is.na(prev_misc), NA, ifelse(prev_misc> 3, 3, prev_misc))))


df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age),
           !is.na(misc))

set.seed(1234)
df = df %>% group_by(lopnr)  %>% sample_n(1) %>% ungroup()

#df1= filter(df, misc== 1)
df1= df

m_tertiles= survfit2(Surv(gest_duration, misc)~ cat_prev_misc, filter(df1, misc==1))

p1= ggsurvfit(m_tertiles) +
theme_cowplot(font_size= 10) +
# add_confidence_interval() +
scale_color_manual(values = colorBlindBlack8[c(2,6,1,8)]) +
  scale_fill_manual(values = colorBlindBlack8[c(2,6,1,8)]) +
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

ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')

fwrite(data.frame(t(summary(m_tertiles)$table)), snakemake@output[[2]], sep= '\t')

m1.zph= cox.zph(coxph(Surv(gest_duration, misc)~ cat_prev_misc, filter(df1, misc== 1)))
x1= data.frame(m1.zph$table)
x1$variable= row.names(x1)

fwrite(x1, snakemake@output[[3]], sep= '\t')

x1= df1 %>%
  tbl_summary(
    include = c(maternal_age, year_transfer, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean} ± {sd}"),
digits = list(everything() ~ c(2))
) %>%
  add_ci(method = list(all_continuous() ~ "t.test")) %>%
  modify_header(ci_stat_0 ~ "**95% CI**")


x2= df1 %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(everything() ~ c(0, 2))
  ) %>%
  add_ci(method = list(all_categorical() ~ "wilson")) %>%
  modify_header(ci_stat_0 ~ "**95% CI**")

x= bind_rows(as.data.frame(x1), as.data.frame(x2))
fwrite(x, snakemake@output[[4]], sep= '\t')
