library(data.table)
library(lmerTest)
library(lme4)
library(dplyr)
library(ggplot2)
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


d= fread(snakemake@input[[1]], h=T)
d$implantation_failure= as.numeric(d$Resultfetus1=='')
d= arrange(d, Etdate) %>% group_by(lopnr) %>% mutate(transfer= row_number()) %>% ungroup()
df= filter(d, !is.na(maternal_age), !is.na(prev_misc), !is.na(transfer) ,!is.na(liveborn))

m1= (glmer(implantation_failure~ maternal_age + prev_misc + transfer + liveborn + (1|lopnr), df, binomial(link = "logit")))

ci= data.frame(confint.merMod(m1, parm = c('maternal_age', 'prev_misc'), method= 'Wald'))
names(ci)= c('lo95', 'up95')

ci$dependent= row.names(ci)
betas= data.frame(summary(m1)$coefficients)
names(betas)= c('beta', 'se', 'z', 'pvalue')
betas$dependent= row.names(betas)

betas= inner_join(betas, ci, by= 'dependent')
betas$cases = sum(model.frame(m1)$implantation_failure == 1, na.rm = TRUE)
betas$sample_size= nobs(m1)
betas$individuals= length(unique(model.frame(m1)$lopnr))
print(str(betas))
betas$beta= exp(betas$beta)
betas$lo95= exp(betas$lo95)
betas$up95= exp(betas$up95)
betas$outcome= "Implantation failure"

implantation= betas

set.seed(1234)
df= sample_frac(df, size= 1) %>% filter(!duplicated(lopnr))

m1= (glm(implantation_failure~ maternal_age + prev_misc + transfer + liveborn, data= df, family= binomial(link = "logit")))

ci= data.frame(confint(m1, parm= c('maternal_age', 'prev_misc'), method= 'Wald'))
names(ci)= c('lo95', 'up95')

ci$dependent= row.names(ci)
betas= data.frame(summary(m1)$coefficients)
names(betas)= c('beta', 'se', 'z', 'pvalue')
betas$dependent= row.names(betas)

betas= inner_join(betas, ci, by= 'dependent')
betas$cases = sum(model.frame(m1)$implantation_failure == 1, na.rm = TRUE)
betas$sample_size= nobs(m1)
betas$individuals= length(unique(model.frame(m1)$lopnr))

print(str(betas))
betas$beta= exp(betas$beta)
betas$lo95= exp(betas$lo95)
betas$up95= exp(betas$up95)
betas$outcome= "(standard)Implantation failure"

implantation_st= betas

df= filter(d, Resultfetus1 != '', !is.na(maternal_age), !is.na(prev_misc), !is.na(misc), !is.na(transfer) ,!is.na(liveborn))

m1= (glmer(misc~ maternal_age + prev_misc + transfer + liveborn + (1|lopnr), df, binomial(link = "logit")))

ci= data.frame(confint.merMod(m1, parm= c('maternal_age', 'prev_misc'), method= 'Wald'))
names(ci)= c('lo95', 'up95')

ci$dependent= row.names(ci)
betas= data.frame(summary(m1)$coefficients)
names(betas)= c('beta', 'se', 'z', 'pvalue')
betas$cases = sum(model.frame(m1)$misc == 1, na.rm = TRUE)
betas$sample_size= nobs(m1)
betas$individuals= length(unique(model.frame(m1)$lopnr))

betas$dependent= row.names(betas)

betas= inner_join(betas, ci, by= 'dependent')

betas$beta= exp(betas$beta)
betas$lo95= exp(betas$lo95)
betas$up95= exp(betas$up95)
betas$outcome= "Miscarriage"

d= bind_rows(implantation, betas)
d=bind_rows(d, implantation_st)

set.seed(1234)
df= sample_frac(df, size= 1) %>% filter(!duplicated(lopnr))

m1= (glm(misc~ maternal_age + prev_misc + transfer + liveborn, data= df, family= binomial(link = "logit")))

ci= data.frame(confint(m1, parm= c('maternal_age', 'prev_misc'), method= 'Wald'))
names(ci)= c('lo95', 'up95')

ci$dependent= row.names(ci)
betas_st= data.frame(summary(m1)$coefficients)
names(betas_st)= c('beta', 'se', 'z', 'pvalue')
betas_st$cases = sum(model.frame(m1)$misc == 1, na.rm = TRUE)
betas_st$sample_size= nobs(m1)
betas_st$individuals= length(unique(model.frame(m1)$lopnr))

betas_st$dependent= row.names(betas_st)

betas_st= inner_join(betas_st, ci, by= 'dependent')

betas_st$beta= exp(betas_st$beta)
betas_st$lo95= exp(betas_st$lo95)
betas_st$up95= exp(betas_st$up95)
betas_st$outcome= "(standard)Miscarriage"

d= bind_rows(d, betas_st)

d= filter(d, dependent %in% c('maternal_age', 'prev_misc'))

d$dependent= with(d, ifelse(dependent== 'maternal_age', 'Maternal age', 'Previous miscarriages'))

fwrite(d, snakemake@output[[2]], sep= '\t')
d= filter(d, outcome %in% c('Miscarriage', 'Implantation failure'))
p1= ggplot(d, aes(x= dependent, y= beta, colour= outcome, fill= outcome)) +
  geom_pointrange(aes(ymin= lo95, ymax=up95), size = .1, fatten= 0.2,
                position=position_dodge(0.1), shape= 21 ) +
  scale_color_manual(values= sPBIYlGn[c(2, 8)]) +
  scale_fill_manual(values= sPBIYlGn[c(2, 8)]) +
  theme_cowplot(font_size= 10) +
geom_hline(yintercept= 1, color= 'grey', linetype = 'dashed', size= 0.6) + 
  ylab('OR (95% CI)') +
  xlab('Outcome, years') +
  theme(legend.title = element_blank(),
        axis.text= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position="bottom")

ggsave(snakemake@output[[1]], p1, width= 88, height= 60, units= 'mm')

