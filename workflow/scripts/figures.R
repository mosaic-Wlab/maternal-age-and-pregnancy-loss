# Figures for maternal age and time-to-miscarriage

library(data.table)
library("dplyr")
library("ggplot2")
library(survival)
library(ggsurvfit)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Read data from women using own oocytes

d= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/maternal_age_miscarriage.txt')

# Filter out implantation failures, biochemical pregnancies with gestational duration< 15 days

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age), 
           !is.na(misc))

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
ggplot(misc_by_age, aes(x= factor(maternal_age), y= p_hat, colour= outcome)) +
  geom_point() +
  geom_errorbar(aes(ymin= loCI, ymax=upCI), width=.1,
                position=position_dodge(0.05)) +
  theme_bw() +
  ylab('Probability (95% CI)') + 
  xlab('Maternal age, years') +
  scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]})

###################### Survival time ### This section needs cleaning

df = df %>% group_by(lopnr) %>% sample_n(1) %>% ungroup()

raw_m= survfit2(Surv(gest_duration, misc)~ 1, df)

ggsurvfit(raw_m)

raw_m= survfit2(Surv(ifelse(gest_duration>150, 150, gest_duration), misc)~ 1, df)

ggsurvfit(raw_m)

raw_m= survfit2(Surv(gest_duration, misc)~ 1, filter(df, misc==1))

ggsurvfit(raw_m)

coxph(Surv(gest_duration, misc)~ maternal_age, df)
coxph(Surv(ifelse(gest_duration>150, 150, gest_duration), misc)~ maternal_age, df)
zph= cox.zph(coxph(Surv(gest_duration, misc)~ maternal_age, filter(df, misc==1)))

plot(zph, resid=F)
abline(coef(fit3)[1:2], lwd= 2, lty= 3, col=2)


fit3 <- coxph(Surv(gest_duration, misc) ~ maternal_age + tt(maternal_age),
               data=filter(df, misc==1),
               tt = function(x, t, ...) x * log(t+3))

zp3= cox.zph(fit3)


x= data.frame(summary(raw_m)[1:10])

x$cum_freq= x$cumhaz / max(x$cumhaz)

raw_misc= survfit2(Surv(gest_duration-15, misc)~ 1, filter(df, misc==1))
raw_misc2= survfit2(Surv(ifelse(gest_duration>150, 150, gest_duration), misc)~ 1, df)

ggsurvfit(raw_misc2)

x= data.frame(summary(raw_misc)[1:10])


df= filter(df, misc==1) %>% mutate(maternal_tertiles= ntile(x= maternal_age, 3))
raw_misc= survfit2(Surv(ifelse(gest_duration>150, 150, gest_duration), misc)~ maternal_tertiles, df)
ggsurvfit(raw_misc)


###################### Time varying effects using cox.zph

df= filter(df, misc== 1) # df= filter(df, !is.na(dat$Resultmiscarriagedate))

m1.zph= cox.zph(coxph(Surv(gest_duration, misc)~ maternal_age, df))

plot(m1.zph, resid= F)

###################### Time varying effects using PAM models
library(pammtools)
library(mgcv)

df1= mutate(df, maternal_tertiles= as.factor(ntile(maternal_age, 3)))

ped = as_ped(df1, Surv(gest_duration, all_loss) ~ maternal_age + maternal_tertiles + Embryofreezetype ,
             id = "id", cut=c(16,seq(20, 308, by=7)))
nrow(ped)  # 389543

# run PAMM
mod.tv.age = bam(ped_status ~ ti(tend,bs='cr',k=11)  + 
                   maternal_tertiles + ti(tend, by=as.ordered(maternal_tertiles),
bs='cr'), data=ped, offset=offset, family=poisson())

summary(mod.tv.age)

# Plot partial effects
pred_df = make_newdata(ped, tend=unique(tend), maternal_tertiles=factor(maternal_tertiles)) %>%
  add_term(mod.tv.age, term="maternal_tertiles", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

ggplot(pred_df, aes(x=(tmid)/7, y=fit, col=maternal_tertiles, fill=maternal_tertiles)) +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_fill_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_x_continuous(breaks=seq(0, 42, by=2), expand=c(0,0)) +
  ylim(c(-1.1, 0.9)) + 
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.75, 0.88), legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), legend.title = element_text(size=10))

################ Same analyses for oocyte recipients
################

egg_rec_all= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/maternal_age_miscarriage_donated_eggs.txt')

# Remove implantation failures 
egg_rec= filter(egg_rec_all, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age), 
                   !is.na(misc))

nrow(egg_rec) # 3170

# Implantation rate
imp_by_age= group_by(egg_rec_all, ifelse(trunc(maternal_age)<= 20, 20, 
                                      ifelse(trunc(maternal_age)> 43, 44, trunc(maternal_age)))) %>% 
  summarize(age= mean(maternal_age), p_hat= mean(Resultfetus1==''), cases= sum(Resultfetus1==''), 
            total= n()) %>% 
  mutate(loCI= p_hat - (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total), 
         upCI= p_hat + (qnorm(1-0.05/2)) * sqrt(p_hat * (1-p_hat)/total))

names(imp_by_age)[1]= 'maternal_age'
imp_by_age$outcome= 'Implantation failure'

#  Miscarriage rate

egg_rec_all$Resultfetus1= factor(egg_rec_all$Resultfetus1, levels= c('Biokemisk graviditet', 'Spontan abort fore 13 veckor', 
                                                                   'Spontan abort vecka 13+0 \x96 21+6', 
                                                                   'Dodfott barn vecka 22+0 \x96 27+6', 
                                                                   'Dodfott barn vecka 28+0 eller mer', 'Levande fott barn'))

misc_by_age= group_by(egg_rec_all, ifelse(trunc(maternal_age)<= 20, 20, 
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

ggplot(misc_by_age, aes(x= factor(maternal_age), y= p_hat, colour= outcome)) +
  geom_point() +
  geom_errorbar(aes(ymin= loCI, ymax=upCI), width=.1,
                position=position_dodge(0.05)) +
  theme_bw() +
  ylab('Probability (95% CI)') + 
  xlab('Maternal age, years') +
  scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]})

###### Survival analysis

egg_rec= group_by(egg_rec, lopnr) %>% sample_n(1)
nrow(egg_rec) # 2636

egg_rec_raw= survfit2(Surv(gest_duration, misc)~ 1, egg_rec)

ggsurvfit(egg_rec_raw)

egg_rec_raw2= survfit2(Surv(gest_duration, misc)~ 1, filter(egg_rec, misc==1))

ggsurvfit(egg_rec_raw2)

zph= cox.zph(coxph(Surv(gest_duration, misc)~ maternal_age, filter(egg_rec, misc==1)))

fit3 <- coxph(Surv(gest_duration, misc) ~ maternal_age + tt(maternal_age),
              data=filter(df, misc==1),
              tt = function(x, t, ...) x * log(t+3))


x= data.frame(summary(raw_m)[1:10])

x$cum_freq= x$cumhaz / max(x$cumhaz)

raw_misc= survfit2(Surv(gest_duration-15, misc)~ 1, filter(df, misc==1))
raw_misc2= survfit2(Surv(ifelse(gest_duration>150, 150, gest_duration), misc)~ 1, df)

ggsurvfit(raw_misc2)

x= data.frame(summary(raw_misc)[1:10])


df= filter(df, misc==1) %>% mutate(maternal_tertiles= ntile(x= maternal_age, 3))
raw_misc= survfit2(Surv(ifelse(gest_duration>150, 150, gest_duration), misc)~ maternal_tertiles, df)
ggsurvfit(raw_misc)

## Test for time-varying effects in oocyte recipients

egg_rec= filter(egg_rec, misc== 1)

m1.oocyte.zph= cox.zph(coxph(Surv(gest_duration, misc)~ maternal_age, egg_rec))

plot(m1.oocyte.zph, resid= F)

# run PAMM for oocyte recipients

egg_rec$maternal_tertiles= as.factor(ntile(egg_rec$maternal_age, 3))

ped_oocyte = as_ped(egg_rec, Surv(gest_duration, misc) ~ maternal_age + maternal_tertiles,
                    id = "id", cut=c(16,seq(20, 150, by=7)))
nrow(ped_oocyte)  # 3039

# run PAMM
mod.tv.age.oocyte = bam(ped_status ~ ti(tend,bs='cr',k=11) + maternal_tertiles + 
                          ti(tend, by=as.ordered(maternal_tertiles), bs='cr'), data=ped_oocyte, offset=offset, 
                        family=poisson())

summary(mod.tv.age.oocyte)

# Plot partial effects
pred_df_oocyte = make_newdata(ped_oocyte, tend=unique(tend), maternal_tertiles=factor(maternal_tertiles)) %>%
  add_term(mod.tv.age.oocyte, term="maternal_tertiles", se_mult=1.96) %>%
  mutate(tmid = 0.5*tstart+0.5*tend)

ggplot(pred_df_oocyte, aes(x=(tmid)/7, y=fit, col=maternal_tertiles, fill=maternal_tertiles)) +
  geom_line(lwd=0.6) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),alpha=0.02,lwd=0.2,lty="dashed") +
  scale_color_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_fill_manual(values=colorBlindBlack8[c(1,6,2)], name="Maternal age tertiles") +
  scale_x_continuous(breaks=seq(0, 23, by=2), expand=c(0,0)) +
  ylim(c(-3.1, 2.9)) + 
  theme_bw() + xlab("gestational age, weeks") + ylab("log hazard ratio") +
  theme(legend.position=c(0.75, 0.88), legend.box.background= element_rect(colour="black"),
        panel.grid.minor.x=element_blank(), legend.title = element_text(size=10))


####### Logistic regression

d= fread('/mnt/hdd/common/pol/maternal_age_miscarriage/results/maternal_age_miscarriage.txt')

# Filter out implantation failures, biochemical pregnancies with gestational duration< 15 days

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age), 
           !is.na(misc))


df$early_loss= ifelse(df$gest_duration< 6* 7, 'early', 
                      ifelse(df$gest_duration < 154, 'late', 
                             ifelse(is.na(df$gest_duration), NA, 'alive')))

df$maternal_tertiles= ntile(df$maternal_age, 3)

df= filter(df, misc==1)
x= group_by(df, maternal_tertiles)  %>% summarize(tot_n= n(), m1= sum(as.numeric(early_loss== 'early')), m2= 'early',
                                                  m3= sum(as.numeric(early_loss=='alive')))
x1= group_by(df, maternal_tertiles) %>% summarize(tot_n= n(), m1= sum(as.numeric(early_loss== 'late')), m2= 'late',
                                                  m3= sum(as.numeric(early_loss== 'alive')))

x= rbind(x, x1)
x$rate= x$m1 / x$tot_n  

x$loCI= x$rate - qnorm(1-0.05/2)*sqrt(x$rate*(1-x$rate)/x$tot_n)
x$upCI= x$rate + qnorm(1-0.05/2)*sqrt(x$rate*(1-x$rate)/x$tot_n)
ggplot(x, aes(as.factor(maternal_tertiles), rate, ymin= loCI, ymax= upCI, colour= m2)) + 
  geom_pointrange() + 
  theme_classic()

