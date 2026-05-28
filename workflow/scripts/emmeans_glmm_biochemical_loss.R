# Packages

library(dplyr)
library(ggplot2)
library(data.table)
library(lme4)
library(lmerTest)
library(cowplot)

sPBIYlGn= c("#FAE9A0FF", "#DBD797FF", "#BCC68DFF", "#9CB484FF", "#7DA37BFF", "#5E9171FF",
            "#3F7F68FF", "#1F6E5EFF", "#005C55FF")


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

df= filter(d, Resultfetus1 != 'Liveborn',
           Resultfetus1 != 'Implantation failure', 
	!is.na(maternal_age),
	!is.na(prev_misc))

table(df$Resultfetus1)

df$biok= ifelse(df$Resultfetus1== 'Biochemical loss', 1, 0)

df$abort= ifelse(df$Resultfetus1 %in% c('Fetal loss <13w',
'Fetal loss 22-28w', 'Fetal loss >28w', 'Fetal loss 13-22w'), 1, 0)

df$biok_loss= ifelse(df$biok== 1, 1, ifelse(df$abort== 1, 0, NA))

library(lme4)

fit <- glmer(
  biok_loss ~ maternal_age + prev_misc + Incubationdays + year_transfer + 
    (1 | lopnr),
  data = df,
  family = binomial()
)

fit_own= fit

coefs <- summary(fit)$coefficients
coefs= data.frame(coefs)
coefs$exposure= rownames(coefs)

ci <- confint(fit, parm = c("maternal_age" , 'prev_misc'), method = "Wald")
ci= data.frame(ci)
names(ci)= c('lo95', 'up95')
ci$exposure= rownames(ci)
tab= inner_join(coefs, ci, by= 'exposure')

tab$oocyte= 'Own'

age_grid <- seq(
  min(df$maternal_age, na.rm = TRUE),
  max(df$maternal_age, na.rm = TRUE),
  length.out = 200
)

library(emmeans)

emm <- emmeans(
  fit,
  ~ maternal_age,
  at = list(maternal_age = age_grid),
  type = "response"
)

emm_df <- as.data.frame(emm)
emm_df$oocyte= 'Own'

d= fread(snakemake@input[[2]])

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

df= filter(d, Resultfetus1 != 'Liveborn',
           Resultfetus1 != 'Implantation failure', !is.na(prev_misc), !is.na(maternal_age))


table(df$Resultfetus1)

df$biok= ifelse(df$Resultfetus1== 'Biochemical loss', 1, 0)

df$abort= ifelse(df$Resultfetus1 %in% c('Fetal loss <13w', 'Fetal loss 22-28w', 'Fetal loss >28w', 'Fetal loss 13-22w'), 1, 0)

df$biok_loss= ifelse(df$biok== 1, 1, ifelse(df$abort== 1, 0, NA))

library(lme4)

fit <- glmer(
  biok_loss ~ maternal_age + prev_misc + Incubationdays + year_transfer +
    (1 | lopnr),
  data = df,
  family = binomial()
)

fit_recipient= fit

age_grid <- seq(
  min(df$maternal_age, na.rm = TRUE),
  max(df$maternal_age, na.rm = TRUE),
  length.out = 200
)

library(emmeans)

emm <- emmeans(
  fit,
  ~ maternal_age,
  at = list(maternal_age = age_grid),
  type = "response"
)

emm_df_recipient <- as.data.frame(emm)
emm_df_recipient$oocyte= 'Recipient'

emm_df= bind_rows(emm_df, emm_df_recipient)

p1= ggplot(emm_df, aes(x = maternal_age, y = prob, fill= oocyte, colour= oocyte)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.2) +
  scale_fill_manual(values = sPBIYlGn[c(2, 8)]) + 
  scale_colour_manual(values = sPBIYlGn[c(2,8)]) +
  geom_line(linewidth = 0.5) +
  labs(x = "Maternal age (years)",
    y = "Predicted probability\n(95% CI)") +
  theme_cowplot(10) +
  theme(legend.position="bottom")

ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')


ci <- confint(fit, parm = c("maternal_age", "prev_misc"), method = "Wald")
coefs <- data.frame(summary(fit)$coefficients)

coefs$exposure = rownames(coefs)
ci= data.frame(ci)
names(ci)= c('lo95', 'up95')

ci$exposure= rownames(ci)

tab2 = inner_join(coefs, ci, by= 'exposure')

tab2$oocyte= 'recipient'

tab= bind_rows(tab, tab2)
tab$Estimate= exp(tab$Estimate)
tab$lo95= exp(tab$lo95)
tab$up95= exp(tab$up95)

fwrite(tab, snakemake@output[[2]], sep= '\t')

