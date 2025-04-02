# Figures for implantation failure and miscsarriage rate per maternal age   

library(data.table)
library("dplyr")
library("ggplot2")

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Implantation rate

imp_by_age= group_by(d, ifelse(trunc(maternal_age)<= 20, 20,
                                 ifelse(trunc(maternal_age)> 43, 44, trunc(maternal_age)))) %>%
  summarize(age= mean(maternal_age), p_hat= mean(Resultfetus1==''), cases= sum(Resultfetus1==''),
            total= n()) %>%
  mutate(loCI= p_hat - (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total),
         upCI= p_hat + (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total))

names(imp_by_age)[1]= 'maternal_age'
imp_by_age$outcome= 'Implantation failure'

df$Resultfetus1= factor(df$Resultfetus1, levels= c('Biokemisk graviditet', 'Spontan abort fore 13 veckor',
                                                   'Spontan abort vecka 13+0 \x96 21+6',
                                                   'Dodfott barn vecka 22+0 \x96 27+6',
                                                   'Dodfott barn vecka 28+0 eller mer', 'Levande fott barn'))

# Miscarriage rate

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

p1= ggplot(misc_by_age, aes(x= factor(maternal_age), y= p_hat, colour= outcome)) +
  geom_point() +
  geom_errorbar(aes(ymin= loCI, ymax=upCI), width=.1,
                position=position_dodge(0.05)) +
  theme_classic() +
  ylab('Probability (95% CI)') +
  xlab('Maternal age, years') +
  scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]})

ggsave(snakemake@output[[1]], p1)

fwrite(misc_by_age, snakemake@output[[2]], sep= '\t')
