# Figures for implantation failure and miscsarriage rate per maternal age   

library(data.table)
library("dplyr")
library("ggplot2")
library(cowplot)
library(showtext)
options(warn=-1)
library(gtsummary)
library(broom)
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


sPBIYlGn= c("#FAE9A0FF", "#DBD797FF", "#BCC68DFF", "#9CB484FF", "#7DA37BFF", "#5E9171FF", "#3F7F68FF", "#1F6E5EFF", "#005C55FF")

d= fread(snakemake@input[[1]])

d$implantation_failure= as.numeric(d$Resultfetus1=='')
d$early_miscarriage= ifelse(is.na(d$misc) | is.na(d$gest_duration), NA, ifelse(d$misc== 1 & d$gest_duration< 7*10, 1, 0))


df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age), !is.na(misc))

#Implantation rate

print(paste('Implantation rate sample size for ', snakemake@wildcards[['oocytes']], ':', nrow(d)))

imp_by_age= group_by(d, ifelse(trunc(maternal_age)<= 20, 20,
                                 ifelse(trunc(maternal_age)> 43, 44, trunc(maternal_age)))) %>%
  summarize(age= mean(maternal_age), p_hat= mean(Resultfetus1==''), cases= sum(Resultfetus1==''),
            total= n()) %>%
  mutate(loCI= p_hat - (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total),
         upCI= p_hat + (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total))

names(imp_by_age)[1]= 'maternal_age'
imp_by_age$outcome= 'Implantation failure'

df$Resultfetus1= with(df, ifelse(Resultfetus1== 'Biokemisk graviditet', 'Biochemical loss',
                               ifelse(Resultfetus1== 'Dodfott barn vecka 22+0 \x96 27+6', 'Fetal loss 22-28w',
                                      ifelse(Resultfetus1== 'Dodfott barn vecka 28+0 eller mer', 'Fetal loss >28w',
                                             ifelse(Resultfetus1== 'Levande fott barn', 'Liveborn',
                                                    ifelse(Resultfetus1== 'Spontan abort fore 13 veckor', 'Fetal loss <13w',
                                                           ifelse(Resultfetus1== '', 'Implantation failure', 'Fetal loss 13-22w')))))))

df$Resultfetus1= factor(df$Resultfetus1, levels= c('Implantation failure', 'Biochemical loss',
                                                 'Fetal loss <13w', 'Fetal loss 13-22w',
                                                 'Fetal loss 22-28w', 'Fetal loss >28w',
                                                 'Liveborn'))



### Descriptive characteristics

x= d %>%
  tbl_summary(
    include = c(maternal_age, year_transfer, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean} ± {sd}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test")) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


x1= d %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, misc, early_miscarriage),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
digits = list(everything() ~ c(0, 2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson")) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")



x= bind_rows(as.data.frame(x), as.data.frame(x1))
fwrite(x, snakemake@output[[3]], sep= '\t')

# Miscarriage rate

print(paste('Implantation rate sample size for ', snakemake@wildcards[['oocytes']], ':', nrow(df)))
print(paste('Implantation rate sample size for ', snakemake@wildcards[['oocytes']], ':', nrow(filter(df, misc==1))))


misc_by_age= group_by(df, ifelse(trunc(maternal_age)<= 20, 20,
ifelse(trunc(maternal_age)> 43, 44, trunc(maternal_age)))) %>%
  summarize(age= mean(maternal_age), p_hat= mean(misc), cases= sum(misc), total= n()) %>%
  mutate(loCI= p_hat - (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total),
           upCI= p_hat + (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total))

names(misc_by_age)[1]= 'maternal_age'

misc_by_age$outcome= 'Miscarriage'
misc_by_age= rbind(misc_by_age, imp_by_age)

misc_by_age$maternal_age= ifelse(misc_by_age$maternal_age== 20, '<=20', ifelse(misc_by_age$maternal_age== 44,
                                                                               '>=44', misc_by_age$maternal_age))

misc_by_age$maternal_age= factor(misc_by_age$maternal_age, levels= c('<=20', seq(21,43, 1), '>=44'),
                                 labels= c('<=20', seq(21,43, 1), '>=44'))

misc_by_age$loCI= ifelse(misc_by_age$loCI< 0, 0, misc_by_age$loCI)
misc_by_age$upCI= ifelse(misc_by_age$upCI> 1, 1, misc_by_age$upCI)

misc_by_age= filter(misc_by_age, !is.na(maternal_age))

p1= ggplot(misc_by_age, aes(x= factor(maternal_age), y= p_hat, colour= outcome, fill= outcome)) +
  geom_pointrange(aes(ymin= loCI, ymax=upCI), size = .1, fatten= 0.2,
                position=position_dodge(0.05), shape= 21 ) +
  scale_color_manual(values= sPBIYlGn[c(2, 8)]) +
  scale_fill_manual(values= sPBIYlGn[c(2, 8)]) +
  theme_cowplot(font_size= 10) +
  ylab('Probability (95% CI)') +
  xlab('Maternal age, years') +
  scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]}) +
  theme(legend.title = element_blank(), 
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position="bottom")

ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')

fwrite(misc_by_age, snakemake@output[[2]], sep= '\t')

x= df %>%
  tbl_summary(include = c(maternal_age, year_transfer, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean} ± {sd}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test")) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


x1= df %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, misc, early_miscarriage),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
digits = list(everything() ~ c(0, 2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson")) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")



x= bind_rows(as.data.frame(x), as.data.frame(x1))
fwrite(x, snakemake@output[[4]], sep= '\t')

