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

d$all_loss= with(d, ifelse(Resultfetus1 %in% c('', 'Biokemisk graviditet', 
                                               'Spontan abort vecka 13+0 \x96 21+6', 'Spontan abort fore 13 veckor', 
                                               'Dodfott barn vecka 22+0 \x96 27+6', 'Dodfott barn vecka 28+0 eller mer'),
                           1, 0))

# Gestational duration --- Beware that it includes clinical pregnancies < 6 gestational weeks

d$gest_duration= ifelse(d$all_loss== 0, d$Resultdeliverydate - d$Etdate + d$Incubationdays,
                        ifelse(d$all_loss== 1 & !is.na(d$Resultmiscarriagedate), d$Resultmiscarriagedate - d$Etdate + 
                                 d$Incubationdays,
                               ifelse(!is.na(d$Incubationdays), d$Incubationdays, NA)))


# Remove gestational duration values that do not correspond to the result of the fetus

d= mutate(d, gest_duration= ifelse((misc == 1 & gest_duration > 154) | gest_duration < 0 | gest_duration > 308 | 
                                     (Resultfetus1 == 'Levande fott barn' & gest_duration < 154) |
                                     (Resultfetus1 == '' & gest_duration > 7) | (Resultfetus1 == 'Biokemisk graviditet' & 
                                                                                   gest_duration < 10) | 
                                     ((Resultfetus1 == 'Dodfott barn vecka 22+0 \x96 27+6' | 
                                        Resultfetus1 == 'Dodfott barn vecka 28+0 eller mer') & gest_duration< 154) |
                                     ((!is.na(Resultmiscarriagedate) & !is.na(Resultdeliverydate)) & 
                                        (Resultmiscarriagedate != Resultdeliverydate)), NA, gest_duration))

d$gest_duration2= with(d, ifelse(is.na(gest_duration), NA, 
                                 ifelse(gest_duration < 154, gest_duration, 154)))

############################## Inclusions / exclusions #########################
# keep singletons, exclude legal abortions, ongoing pregnancy, unkown status of fetus, ectopic pregnancies and missing transfer date
# Exclusion of implantation failure and biochemical pregnancy should be done downstream

d= filter(d, Resultfetus2== '', Resultfetus3== '', Resultfetus1!= 'Legal abort', Resultfetus1 != 'Pagaende viabel graviditet',
          Resultfetus1 != 'Okand', Resultfetus1 != 'Ektopisk graviditet', !is.na(Etdate))


########################### Add previous miscarriages

recurrent= fread('results/main_data/recurrent-multiple-loss.txt')

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
