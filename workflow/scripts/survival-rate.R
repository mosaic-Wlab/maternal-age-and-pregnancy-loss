# Figures for maternal age and time-to-miscarriage

library(data.table)
library("dplyr")
library("ggplot2")
library(survival)
library(ggsurvfit)
library(cowplot)
library(showtext)
library(gtsummary)
library(broom)
options(warn=-1)

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sPBIYlGn= c("#FAE9A0FF", "#DBD797FF", "#BCC68DFF", "#9CB484FF", "#7DA37BFF", "#5E9171FF", "#3F7F68FF", "#1F6E5EFF", "#005C55FF")

d= fread(snakemake@input[[1]])
d$implantation_failure= as.numeric(d$Resultfetus1=='')
d$early_miscarriage= ifelse(is.na(d$misc) | is.na(d$gest_duration), NA, ifelse(d$misc== 1 & d$gest_duration< 7*10, 1, 0))


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

x= df %>%
  tbl_summary(
    include = c(maternal_age, year_transfer, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean} ± {sd}"),
digits = list(everything() ~ c(2))
) %>%
  add_ci(method = list(all_continuous() ~ "t.test")) %>%
  modify_header(ci_stat_0 ~ "**95% CI**")


x1= df %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, misc, early_miscarriage),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
digits = list(everything() ~ c(0, 2))
  ) %>%
  add_ci(method = list(all_categorical() ~ "wilson")) %>%
  modify_header(ci_stat_0 ~ "**95% CI**")

x= bind_rows(as.data.frame(x), as.data.frame(x1))


fwrite(x, snakemake@output[[4]], sep= '\t')

}


if (grepl('maternal-age', snakemake@output[[1]])) {

df= mutate(df, maternal_tertiles= ntile(x= maternal_age, 3))

df$maternal_tertiles= as.character(df$maternal_tertiles)

m_tertiles= survfit2(Surv(gest_duration, misc)~ maternal_tertiles, filter(df, misc==1))

p1= ggsurvfit(m_tertiles) +
theme_cowplot(font_size= 10) +
# add_confidence_interval() +
scale_color_manual(values = sPBIYlGn[c(1,5,9)]) +
  scale_fill_manual(values = sPBIYlGn[c(1,5,9)]) +
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


