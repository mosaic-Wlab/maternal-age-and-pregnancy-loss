library(data.table)
library("dplyr")


d= fread(snakemake@input[[1]]) #d= fread('/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt')

# Reformatting data

d$Resultmiscarriagedate= as.Date(d$Resultmiscarriagedate)
d$Resultdeliverydate= as.Date(d$Resultdeliverydate)
d$Etdate= as.Date(d$Etdate)
d$Cyclestartdate= as.Date(d$Cyclestartdate)

d$year_transfer= format(as.Date(d$Etdate, format="%Y-%m-%d"),"%Y")

# Formatting data from MBR

mbr= fread(snakemake@input[[2]]) #mbr= fread('/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt')

mbr= select(mbr, mlopnr, blopnr, MFODDAT, BFODDAT, PARITET, TIDSPOAB, TIDDODF) # Think about what to do with parity, etc

mbr$MFODDAT= as.Date(as.character(mbr$MFODDAT), format= "%Y%m%d")
mbr$BFODDAT= as.Date(as.character(mbr$BFODDAT), format= "%Y%m%d")

# For each individual we want to get their birth date

# We can do much more to clean MFODDAT/ BFODDAT

yob_m= select(mbr, mlopnr, MFODDAT) %>% filter(!is.na(MFODDAT)) %>% filter(!duplicated(mlopnr))

yob_b= select(mbr, blopnr, BFODDAT) %>% filter(!is.na(BFODDAT)) %>% filter(!duplicated(blopnr))

names(yob_m)= c('lopnr', 'FODDAT')
names(yob_b)= c('lopnr', 'FODDAT')

yob= bind_rows(yob_m, yob_b)
yob= filter(yob, !duplicated(lopnr))

d= left_join(d, yob, by= c('lopnr'))

################################################## Defining events, censoring and time-to-miscarriage
# Exclusion of implantation failure and biochemical pregnancy should be done downstream

# Clinical miscarriage (this includes some miscarriages occuring < gestational 6 weeks)
# All other losses are censored (implantation failure, biochemical pregnancies and abortions after gw 22)

d$clin_misc= with(d, ifelse(Resultfetus1 %in% c('Spontan abort vecka 13+0 \x96 21+6', 'Spontan abort fore 13 veckor'), 1, 0))

# Miscarriage (includes both miscarriages < 6 gw and biochemical pregnancies)
# All implantation failure and losses > 22gw are censored

d$misc= with(d, ifelse(Resultfetus1 %in% c('Spontan abort vecka 13+0 \x96 21+6', 'Spontan abort fore 13 veckor',
                                           'Biokemisk graviditet'), 1, 0))

# clinical pregnancy loss (as Clinical miscarriage, but includes all losses during pregnancy)
# All implantation failure and biochemical pregnancies are censored

d$clin_loss= with(d, ifelse(Resultfetus1 %in% c('Spontan abort vecka 13+0 \x96 21+6', 'Spontan abort fore 13 veckor', 
                                                'Dodfott barn vecka 22+0 \x96 27+6', 'Dodfott barn vecka 28+0 eller mer'), 
                            1, 0))
# all pregnancy losses (as Miscarriage, but includes all losses during pregnancy)
# Only implantation failures are censored

d$all_loss= with(d, ifelse(Resultfetus1 %in% c('Biokemisk graviditet', 
                                               'Spontan abort vecka 13+0 \x96 21+6', 'Spontan abort fore 13 veckor', 
                                               'Dodfott barn vecka 22+0 \x96 27+6', 'Dodfott barn vecka 28+0 eller mer'),
                           1, 0))

# Gestational duration --- Beware that it includes clinical pregnancies < 6 gestational weeks

d$gest_duration= ifelse(d$all_loss== 0, d$Resultdeliverydate - d$Etdate + d$Incubationdays,
                        ifelse(d$all_loss== 1 & !is.na(d$Resultmiscarriagedate), d$Resultmiscarriagedate - d$Etdate + 
                                 d$Incubationdays,
                               ifelse(!is.na(d$Incubationdays), d$Incubationdays, NA)))

d$misc= d$all_loss
# Remove gestational duration values that do not correspond to the result of the fetus

d= mutate(d, gest_duration= ifelse(gest_duration < 16 | gest_duration > 308 | (Resultfetus1 == 'Levande fott barn' & gest_duration < 154) | (Resultfetus1 == '' & gest_duration > 7) | (Resultfetus1 == 'Biokemisk graviditet' & gest_duration < 10) | ((!is.na(Resultmiscarriagedate) & !is.na(Resultdeliverydate)) & (Resultmiscarriagedate != Resultdeliverydate)), NA, gest_duration))

d$gest_duration2= with(d, ifelse(is.na(gest_duration), NA, 
                                 ifelse(gest_duration < 154, gest_duration, 154)))

########################################## Add male and female infertility

ICD_COL_NAMES=c("DIA1","DIA2","DIA3","DIA4","DIA5","DIA6","DIA7","DIA8","DIA9","DIA10","DIA11","DIA12","DIA13","DIA14","DIA15","DIA16","DIA17","DIA18","DIA19","DIA20","DIA21","DIA22","DIA23","DIA24","DIA25","DIA26","DIA27","DIA28","DIA29","DIA30")

outpatient= fread(snakemake@input[[3]])

outpatient= select(outpatient, lopnr, AR, all_of(ICD_COL_NAMES))
INFERTILITY_CODES=c('N97','N970','N971','N972','N973','N978','N979')

MALE_INFERTILITY_CODES=c('N974')
outpatient= mutate(outpatient, male_infertility= case_when((if_any(ICD_COL_NAMES, ~ . %in% MALE_INFERTILITY_CODES)) ~ TRUE, TRUE ~ FALSE), female_infertility=  case_when((if_any(ICD_COL_NAMES, ~ . %in% INFERTILITY_CODES)) ~ TRUE, TRUE ~ FALSE))

female_infer= unique(filter(outpatient, female_infertility) %>% pull(lopnr))
male_infer= unique(filter(outpatient, male_infertility) %>% pull(lopnr))

d$male_infer= with(d, ifelse(lopnr %in% male_infer, 1, 0))
d$female_infer= with(d, ifelse(lopnr %in% female_infer, 1, 0))

############################## Inclusions / exclusions #########################
# keep singletons, exclude legal abortions, ongoing pregnancy, unkown status of fetus, ectopic pregnancies and missing transfer date
# Exclusion of implantation failure and biochemical pregnancy should be done downstream
ninitial= nrow(d)

d= filter(d, !is.na(Etdate))

n_missing= ninitial - nrow(d)
d= filter(d, Resultfetus2== '', Resultfetus3== '')
n_twins= ninitial - n_missing - nrow(d)

d= filter(d, Resultfetus1!= 'Legal abort')

n_abort= ninitial - n_missing - n_twins - nrow(d)

d= filter(d, Resultfetus1 != 'Ektopisk graviditet')

n_ectopic= ninitial - n_missing - n_twins - n_abort - nrow(d)

d= filter(d, Resultfetus1 != 'Pagaende viabel graviditet', Resultfetus1 != 'Okand')

n_other= ninitial - n_missing - n_twins - n_abort - n_ectopic - nrow(d)

d= filter(d, !is.na(year_transfer))

n_year= ninitial - n_missing - n_twins - n_abort - n_ectopic - n_other - nrow(d)
d = filter(d,  !is.na(Incubationdays))

n_incubation= ninitial - n_missing - n_twins - n_abort - n_ectopic - n_other - n_year - nrow(d)
print(paste('Number of missing Embryo transfer data', n_missing))

print(paste('Number of twin pregnancies excluded', n_twins))
print(paste('Number of elective abortions excluded', n_abort))
print(paste('Number of ectopic pregnancies excluded', n_ectopic))
print(paste('Number of unknown pregnancy status', n_other))
print(paste('Number of missing year of transfer', n_year))
print(paste('Number of missing embryo incubation days', n_incubation))

########################### Add previous miscarriages

recurrent= fread(snakemake@input[[4]])#fread('results/main_data/recurrent-multiple-loss.txt')

recurrent= filter(recurrent, !is.na(recurrent))

x= full_join(d, recurrent, by= 'lopnr')

x= filter(x, (Etdate> date) | is.na(date) )
x$prev_misc= ifelse(is.na(x$prev_misc), 0, x$prev_misc)

x= group_by(x, lopnr, Etdate) %>% filter(date == max(date))


x$cat_prev_misc= ifelse(is.na(x$prev_misc), NA,
                        ifelse(x$prev_misc> 3, 3, x$prev_misc))

d$lopnr=  as.integer(d$lopnr)
d$Etdate= as.Date(d$Etdate)

x$lopnr= as.integer(x$lopnr)
x$Etdate= as.Date(x$Etdate)

d= left_join(d, x[, c('lopnr', 'Etdate', 'date', 'prev_misc')], by = c('lopnr', 'Etdate'))

d$prev_misc= ifelse(is.na(d$prev_misc), 0, d$prev_misc)


### Add number of previous live births

mbr= fread(snakemake@input[[2]])

mbr$BFODDAT = as.Date(mbr$BFODDAT, format='%Y%m%d')
mbr = mbr %>% filter(is.na(DODFOD), !is.na(blopnr), !is.na(mlopnr)) %>%  arrange(mlopnr, BFODDAT) %>%  select(mlopnr, BFODDAT)

d = left_join(d, mbr, by = c('lopnr' = 'mlopnr'), relationship = "many-to-many") %>% mutate(prior_birth = if_else(is.na(BFODDAT) | Etdate <= BFODDAT, FALSE, TRUE)) %>%   group_by(lopnr, Etdate) %>% summarise(liveborn = sum(prior_birth, na.rm = TRUE), .groups = 'drop') %>%  right_join(d, by = c('lopnr', 'Etdate')) %>% mutate(liveborn = if_else(is.na(liveborn), 0, liveborn))

########################### Define maternal age

d$maternal_age= as.numeric(difftime(d$Cyclestartdate, d$FODDAT, units= 'weeks'))/52.25

d= group_by(d, lopnr, Etdate) %>% filter(row_number()== 1)

df= filter(d, Oocytowndonated=='Egna')
d_donated= filter(d, Oocytowndonated!='Egna')

df$maternal_tertiles= ntile(df$maternal_age, 3)
d_donated$maternal_tertiles= ntile(d_donated$maternal_age, 3)

fwrite(df, snakemake@output[[1]], sep= '\t')
fwrite(d_donated, snakemake@output[[2]], sep= '\t')
fwrite(d, snakemake@output[[3]], sep ='\t')
