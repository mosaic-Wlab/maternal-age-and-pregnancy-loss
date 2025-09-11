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

d= fread(snakemake@input[[1]])

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age), 
           !is.na(misc))

set.seed(1234)
df = df %>% group_by(lopnr) %>% sample_n(1) %>% ungroup()

print(paste('Total sample size for ', snakemake@wildcards[['oocytes']], ':', nrow(df)))
print(paste('Cases Miscarriage sample size for ', snakemake@wildcards[['oocytes']], ':', nrow(filter(df, misc==1))))


if (grepl('raw', snakemake@output[[1]])) {
raw_m= survfit2(Surv(gest_duration, misc)~ 1, df)

p1= ggsurvfit(raw_m) +
theme_cowplot(font_size= 10) +
 add_confidence_interval() +
#  add_risktable() +
# add_quantile(color = "gray50", linewidth = 0.2) +
#  scale_ggsurvfit() +
xlab("Duration of gestation, days") +
ylab("Percentage survival") +
theme(legend.title = element_blank(),
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position="bottom")

raw_m2= survfit2(Surv(gest_duration, misc)~ 1, filter(df, misc==1))

p2= ggsurvfit(raw_m2) +
theme_cowplot(font_size= 10) +
 add_confidence_interval() +
#  add_risktable() +
  add_quantile(color = "gray50", linewidth = 0.2) +
#  scale_ggsurvfit() +
xlab("Duration of gestation, days") +
ylab("Percentage survival") +
theme(legend.title = element_blank(),
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position="bottom")

print(snakemake@output[[1]])
ggsave(snakemake@output[[1]], p1, width= 88, height= 60,  units= 'mm')
ggsave(snakemake@output[[2]], p2, width= 88, height= 60, units= 'mm')

x= data.frame(t(summary(raw_m2)$table))

fwrite(x, snakemake@output[[3]], sep= '\t')
}


if (grepl('maternal-age', snakemake@output[[1]])) {

df= mutate(df, maternal_tertiles= ntile(x= maternal_age, 3))

df$maternal_tertiles= as.character(df$maternal_tertiles)

m_tertiles= survfit2(Surv(gest_duration, misc)~ maternal_tertiles, filter(df, misc==1))

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

ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')

fwrite(summary(m_tertiles)$table, snakemake@output[[2]], sep= '\t')

m1.zph= cox.zph(coxph(Surv(gest_duration, misc)~ maternal_age, filter(df, misc== 1)))
x= data.frame(m1.zph$table)
x$variable= row.names(x)

fwrite(x, snakemake@output[[3]], sep= '\t')

}


