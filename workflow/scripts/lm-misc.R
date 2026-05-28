library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(cowplot)
library(lmerTest)
library(lme4)


sPBIYlGn= c("#FAE9A0FF", "#DBD797FF", "#BCC68DFF", "#9CB484FF", "#7DA37BFF", "#5E9171FF",
            "#3F7F68FF", "#1F6E5EFF", "#005C55FF")


#### This script has two sections. One where we estimate the effects of maternal age and amternal age tertiles on biochemical pregnancy loss (among all pregnancy losses), and the second where we estimate the effect of maternal age (and tertiles) on the duration of pregnancy losses.


d= fread(snakemake@input[[1]])

d$Resultfetus1= with(d, ifelse(Resultfetus1== 'Biokemisk graviditet', 'Biochemical loss',
  ifelse(Resultfetus1== 'Dodfott barn vecka 22+0 \x96 27+6', 'Fetal loss 22-28w',
    ifelse(Resultfetus1== 'Dodfott barn vecka 28+0 eller mer', 'Fetal loss >28w',
      ifelse(Resultfetus1== 'Levande fott barn', 'Liveborn',
        ifelse(Resultfetus1== 'Spontan abort fore 13 veckor', 'Fetal loss <13w',
          ifelse(Resultfetus1== '', 'Implantation failure', 'Fetal loss 13-22w')))))))

d$Resultfetus1= factor(d$Resultfetus1, levels= c('Implantation failure', 'Biochemical loss',
                                                 'Fetal loss <13w', 'Fetal loss 13-22w', 
                                                 'Fetal loss 22-28w', 'Fetal loss >28w', 
                                                 'Liveborn'))

df= filter(d, Resultfetus1 != 'Implantation failure', !is.na(maternal_age))

if (grepl('male-infertility', snakemake@output[[1]]) ){

df= filter(df, male_infer== 1 & female_infer== 0)

} else if (grepl('sperm-recipient', snakemake@output[[1]]) ){
df = filter(df, Spermowndonated== 'Donerade')

}

#set.seed(1234)
#df = df  %>% group_by(lopnr) %>% sample_n(1) %>% ungroup()

df= mutate(df, maternal_tertiles= as.factor(ntile(maternal_age, 3)))

table(df$Resultfetus1)


#### we first plot the estimated probability of biochemical pregnancy loss among all losses


df$biok= ifelse(df$Resultfetus1== 'Biochemical loss', 1, 0)

df$abort= ifelse(df$Resultfetus1 %in% c('Fetal loss <13w', 
                                        'Fetal loss 13-22w', 'Fetal loss 22-28w', 'Fetal loss >28w'), 1, 0)

df$biok_loss= ifelse(df$biok== 1, 1, ifelse(df$abort== 1, 0, NA))


z= group_by(df, maternal_tertiles) %>% summarize(total= n(), 
                                                 abort= sum(abort),
                                                 biok= sum(biok),
						p_hat= sum(biok) / sum(abort)) %>%
mutate(loCI= p_hat - (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total),
upCI= p_hat + (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total))

p1= ggplot(z, aes(x= factor(maternal_tertiles), y= p_hat, colour=maternal_tertiles , fill= maternal_tertiles)) +
  geom_pointrange(aes(ymin= loCI, ymax=upCI), size = .1, fatten= 0.2,
                position=position_dodge(0.05), shape= 21 ) +
  scale_color_manual(values= sPBIYlGn[c(1, 5, 9)]) +
  scale_fill_manual(values= sPBIYlGn[c(1, 5, 9)]) +
  theme_cowplot(font_size= 10) +
  ylab('Probability (95% CI)') +
  xlab('Maternal age tertiles') +
  theme(
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position="none")

#ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')


#### Use estimated marginal means to get the estimated probability for each tertile
m1= (glmer(biok_loss ~ maternal_tertiles + prev_misc + year_transfer + Incubationdays + (1|lopnr), data=df, family= binomial(link = "logit")))

library(emmeans)

emm_df <- as.data.frame(emmeans(m1, ~ maternal_tertiles, type = "response"))

p1= ggplot(emm_df, aes(x = maternal_tertiles, y = prob, fill= maternal_tertiles, colour= maternal_tertiles)) +
geom_pointrange(aes(ymin= asymp.LCL,, ymax= asymp.UCL), size = .1, fatten= 0.2,
                position=position_dodge(0.05), shape= 21 ) +   
  scale_fill_manual(values = sPBIYlGn[c(1, 5, 9)]) +
  scale_colour_manual(values = sPBIYlGn[c(1, 5, 9)]) +
  labs(x = "Maternal age tertiles",
    y = "Predicted probability\n(95% CI)") +
  theme_cowplot(10) +
  theme(legend.position="none")

ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')

# Estimate the effect of maternal age on biochemical pregnancy loss

m1= (glmer(biok_loss ~ maternal_age + prev_misc + year_transfer + Incubationdays + (1|lopnr), data=df, family= binomial(link = "logit")))

glm_confs= data.frame(confint.merMod(m1, parm= c("maternal_age", "prev_misc"), method= 'Wald'))
names(glm_confs)= c('lo95', 'up95')
glm_confs$exposure= row.names(glm_confs)
glm_confs$lo95= exp(glm_confs$lo95)
glm_confs$up95= exp(glm_confs$up95) 
glm_sums= data.frame(summary(m1)$coefficients)
glm_sums$exposure= row.names(glm_sums)
glm_sums$Estimate= exp(glm_sums$Estimate)
glm_sums_cont= inner_join(glm_sums, glm_confs, by= 'exposure')
print('We are here')
glm_sums_cont$cases = sum(model.frame(m1)$biok_loss == 1, na.rm = TRUE)
glm_sums_cont$sample_size= nobs(m1)
glm_sums_cont$individuals= length(unique(model.frame(m1)$lopnr))

# Estimate the effect of maternal age tertiles on biochemical pregnancy loss

m1= (glmer(biok_loss ~ maternal_tertiles + prev_misc + year_transfer + Incubationdays + (1|lopnr), data=df, family= binomial(link = "logit")))

glm_confs= data.frame(confint(m1, parm= c('maternal_tertiles2', 'maternal_tertiles3', 'prev_misc'), method = 'Wald'))

names(glm_confs)= c('lo95', 'up95')
glm_confs$lo95= exp(glm_confs$lo95)
glm_confs$up95= exp(glm_confs$up95)

glm_confs$exposure= row.names(glm_confs) 
glm_sums= data.frame(summary(m1)$coefficients)
glm_sums$Estimate= exp(glm_sums$Estimate)
glm_sums$exposure= row.names(glm_sums)

mat_tert= group_by(df, maternal_tertiles) %>% filter(!is.na(maternal_age), !is.na(prev_misc), !is.na(biok_loss), !is.na(year_transfer), !is.na(Incubationdays)) %>% summarize(cases= sum(biok_loss), sample_size= n())
mat_tert$individuals= length(unique(model.frame(m1)$lopnr))
glm_sums_bin= inner_join(glm_sums, glm_confs, by= 'exposure')

print('We are here 2')
names(mat_tert)= c('exposure', 'cases', 'sample_size', 'individuals')
print(mat_tert)
mat_tert$exposure= paste0('maternal_tertiles', mat_tert$exposure)

glm_sums_bin= left_join(glm_sums_bin, mat_tert, by= 'exposure' )

glm_df= bind_rows(glm_sums_bin, glm_sums_cont)

fwrite(glm_df, snakemake@output[[2]], sep= '\t')

### Estimate maternal age effects on duration of pregnancy loss

m1= (lmer(gest_duration~ maternal_age + prev_misc + year_transfer + Incubationdays + (1|lopnr), filter(df, misc == 1)))
lm_confs= data.frame(confint(m1, parm= c('maternal_age', 'prev_misc'), method = 'Wald'))
names(lm_confs)= c('lo95', 'up95')
lm_confs$exposure= row.names(lm_confs)
 
lm_sums= data.frame(summary(m1)$coefficients)
lm_sums$exposure= row.names(lm_sums)

lm_sums_cont= inner_join(lm_sums, lm_confs, by= 'exposure')
lm_sums_cont$sample_size= nobs(m1)
lm_sums_cont$individuals= length(unique(model.frame(m1)$lopnr))
print('we are here 3')
### Estimate maternal age tertiles effects on duration of pregnancy loss using linear mixed models

m1= (lmer(gest_duration~ maternal_tertiles + prev_misc + year_transfer + Incubationdays + (1|lopnr), filter(df, !is.na(biok_loss))))
lm_confs= data.frame(confint(m1, parm= c('maternal_tertiles2', 'maternal_tertiles3', 'prev_misc'), method = 'Wald'))
lm_confs$exposure= row.names(lm_confs)

lm_sums= data.frame(summary(m1)$coefficients)
lm_sums$exposure= row.names(lm_sums)

lm_sums_bin= inner_join(lm_sums, lm_confs, by= 'exposure')

mat_tert= filter(df, misc== 1, !is.na(maternal_age), !is.na(gest_duration), !is.na(prev_misc), !is.na(year_transfer), !is.na(Incubationdays)) %>% group_by(maternal_tertiles) %>% summarize(sample_size= n())
names(mat_tert)= c('exposure', 'sample_size')
mat_tert$exposure= paste0('maternal_tertiles', mat_tert$exposure)
mat_tert$individuals= length(unique(model.frame(m1)$lopnr))


lm_sums_bin= left_join(lm_sums_bin, mat_tert, by= 'exposure' )

lm= bind_rows(lm_sums_cont, lm_sums_bin)
fwrite(lm, snakemake@output[[3]], sep= '\t')

p1= ggplot(filter(df, misc == 1), aes(gest_duration, group= maternal_tertiles, fill = maternal_tertiles)) +
geom_density(alpha = 0.6) +
  scale_color_manual(values= sPBIYlGn[c(1, 5, 9)], name= 'Maternal age tertiles') +
  scale_fill_manual(values= sPBIYlGn[c(1, 5, 9)], name= 'Maternal age tertiles') +
  theme_cowplot(font_size= 10) +
  ylab('Density') +
  xlab('Gestational duration, days') + 
  theme(#legend.title = element_blank(),
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position="bottom")

ggsave(snakemake@output[[4]], p1, width= 88, height= 60, units= 'mm')

