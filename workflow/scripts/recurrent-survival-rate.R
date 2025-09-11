# Figures for maternal age and time-to-miscarriage

library(data.table)
library("dplyr")
library("ggplot2")
library(survival)
library(ggsurvfit)
library(cowplot)
library(showtext)
options(warn=-1)

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#d= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/main_data/QIVF-QC-own-oocyte.txt')

d= fread(snakemake@input[[1]])

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age),
           !is.na(misc))

set.seed(1234)
df = df %>% group_by(lopnr)  %>% sample_n(1) %>% ungroup()

#df1= filter(df, misc== 1)
df1= df

#recurrent= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/main_data/recurrent-multiple-loss.txt')
recurrent= fread(snakemake@input[[2]])


# Setting any recurrent pregnancy loss and any miscarriage

anymisc= arrange(recurrent, desc(prev_misc)) %>%
  group_by(lopnr) %>% filter(row_number()== 1)

anymisc= anymisc[, c('lopnr', 'prev_misc')]
anymisc$prev_misc= ifelse(is.na(anymisc$prev_misc), 0, anymisc$prev_misc)
names(anymisc)= c('lopnr', 'nmisc')

df1= left_join(df1, anymisc, by= c('lopnr'))

df1$nmisc= df1$nmisc - df1$misc

df1= mutate(df1, anyrecurrent= ifelse(lopnr %in% (filter(recurrent, recurrent== 1) %>% pull(lopnr)), 'Recurrent',
                                   ifelse(nmisc> 0, 'Multiple', 'None')))

df1$anyrecurrent= ifelse(is.na(df1$anyrecurrent), 'None', df1$anyrecurrent)

df1$anyrecurrent= factor(df1$anyrecurrent, levels= c('None', 'Multiple', 'Recurrent'))

x= full_join(df1, recurrent, by= 'lopnr')
x= filter(x, Etdate> date)
x$prev_misc= ifelse(is.na(x$prev_misc), 0, x$prev_misc)
x= group_by(x, lopnr) %>% filter(prev_misc == max(prev_misc))
x= filter(x, !duplicated(lopnr))

x$cat_prev_misc= ifelse(is.na(x$prev_misc), NA, 
		ifelse(x$prev_misc> 3, 3, x$prev_misc))

m_tertiles= survfit2(Surv(gest_duration, misc)~ cat_prev_misc, filter(x, misc==1))

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

m1.zph= cox.zph(coxph(Surv(gest_duration, misc)~ cat_prev_misc, filter(x, misc== 1)))
x= data.frame(m1.zph$table)
x$variable= row.names(x)

fwrite(x, snakemake@output[[3]], sep= '\t')


