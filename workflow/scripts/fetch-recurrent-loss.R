
library(data.table)
library(dplyr)
library(tidyverse)
library(lubridate)

qivf= fread(snakemake@input[[1]])#fread('/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt')
mfr= fread(snakemake@input[[2]])#fread('/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt')
patreg= fread(snakemake@input[[3]])#fread('/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Patientregistret/UT_R_PAR_OV_M_14517_2022.txt')

patreg <- patreg[, c("lopnr","INDATUM","hdia","AR","DIA1","DIA2","DIA3","DIA4","DIA5","DIA6","DIA7","DIA8","DIA9","DIA10","DIA11","DIA12","DIA13","DIA14","DIA15","DIA16","DIA17","DIA18","DIA19","DIA20","DIA21","DIA22","DIA23","DIA24","DIA25","DIA26","DIA27","DIA28","DIA29","DIA30")]

mfr <- mfr[, c("mlopnr","AR","BFODDAT","TIDSPOAB","MDIAG1","MDIAG2","MDIAG3","MDIAG4","MDIAG5","MDIAG6","MDIAG7","MDIAG8","MDIAG9","MDIAG10","MDIAG11","MDIAG12","MDIAGNOS")]
colnames(mfr)[which(names(mfr)=="mlopnr")]<-"lopnr"

qivf <- qivf[, c("lopnr","Resultmiscarriagedate")]
qivf = filter(qivf, !is.na(Resultmiscarriagedate))
names(qivf)= c('lopnr', 'INDATUM')

ICD_COL_NAMES_MFR=c("MDIAG1","MDIAG2","MDIAG3","MDIAG4","MDIAG5","MDIAG6","MDIAG7","MDIAG8","MDIAG9","MDIAG10","MDIAG11","MDIAG12","MDIAGNOS")
ICD_COL_NAMES_PATREG=c("hdia","DIA1","DIA2","DIA3","DIA4","DIA5","DIA6","DIA7","DIA8","DIA9","DIA10","DIA11","DIA12","DIA13","DIA14","DIA15","DIA16","DIA17","DIA18","DIA19","DIA20","DIA21","DIA22","DIA23","DIA24","DIA25","DIA26","DIA27","DIA28","DIA29","DIA30")
RECURRENT_MISCARRIAGE_CODES=c('N96','N969','646D','63490','O262')
MISCARRIAGE_CODES=c('O03','O030','O031','O032','O033','O034','O035','O036','O037','O038','O039','O03','O020','O021','O028','O029','632','634','63460','64300','64310','64320','64399')

temp <- mfr %>% 
mutate(RPL=case_when((if_any(ICD_COL_NAMES_MFR, ~ . %in% RECURRENT_MISCARRIAGE_CODES)) ~ TRUE, TRUE ~ FALSE))

temp1=subset(temp, RPL==TRUE)         #cases from MFR with ICD code in MFR for RPL (3 consecutive)

temp2=subset(temp, TIDSPOAB>=2)       #cases from MFR marked in MFR as 2 or more previous miscarriages

temp3=subset(temp, TIDSPOAB>=3)       #cases from MFR marked in MFR as 3 or more previous miscarriages

x=subset(temp1, TIDSPOAB<3)           #cases from MFR with incoherent coding in MFR, ie there is an ICD code for RPL (3 consecutive), but at the same time marked as "less than two previous miscarriages"

temp1 <- temp1[, c("lopnr")]

temp2 <- temp2[, c("lopnr")]

x <- x[, c("lopnr")]

mfr_multiplemisc= rbind(temp1, temp2, x) 

mfr_multiplemisc= mfr_multiplemisc %>% distinct() #all mothers from MFR where multiplemisc==TRUE, including all where RPL==TRUE, no duplications. No regard given to the order or when she qualifies for multiplemisc or RPL. I'd say it is not reasonable to try to extract information about the order of previous miscarriages from MFR, so this is as much as we can extract here.

mfr_RPL= temp1
mfr_RPL= mfr_RPL %>% distinct() #all mothers from MFR where RPL==TRUE, no duplications

patreg_RPL <- patreg %>% 
mutate(RPL=case_when((if_any(ICD_COL_NAMES_PATREG, ~ . %in% RECURRENT_MISCARRIAGE_CODES)) ~ TRUE, TRUE ~ FALSE))

patreg_RPL= subset(patreg_RPL, RPL==TRUE)

patreg_RPL <- patreg_RPL[, c("lopnr")]

patreg_RPLcoded= patreg_RPL %>% distinct()          #all mothers from PATREG where there is a code registered for RPL (3 consecutive), no duplications

rm(temp,temp1,temp2,temp3,x,patreg_RPL)


##identifying individual miscarriages from PATREG, defined as at least 42 days between different miscarriages (own definition, but 6 weeks should be a reasonable amount of time for treating one miscarriage)

temp <- patreg %>% 
mutate(misc= case_when((if_any(ICD_COL_NAMES_PATREG, ~ . %in% MISCARRIAGE_CODES)) ~ TRUE, TRUE ~ FALSE))

temp= subset(temp, misc== TRUE)

temp <- temp[, c("lopnr","INDATUM")]

temp= rbind(qivf, temp)

temp= temp %>% 
  arrange(lopnr, INDATUM)
  
temp= temp2temp= temp %>% 
  group_by(lopnr) %>%
  mutate(lagged_date = lag(INDATUM)) %>%
  ungroup()
  

temp$INDATUM2 = temp$INDATUM

temp$lagged_date2 = temp$lagged_date

temp$INDATUM2= as.numeric(temp$INDATUM2)
temp$lagged_date2=as.numeric(temp$lagged_date2)

temp= temp %>% 
  mutate(date_diff = INDATUM2 - lagged_date2)
  
temp=subset(temp, is.na(date_diff) | date_diff>41)

temp <- temp[, c("lopnr","INDATUM")]    #mothers from PATREG with each individual miscarriage represented with one row, and the date it was first documented (defined as a new miscarriage if a miscarriage-code was given the mother at least 42 days after a previous miscarriage-code was given to her)


#define multiple miscarriage group (any 2 miscarriages regardless of any live births) from temp-file above from PATREG

temp <- temp %>%
  group_by(lopnr) %>%
  mutate(index = 1:n()) %>%
  ungroup()

total_misc= group_by(temp, lopnr) %>% summarize(n_misc= n())

patreg_multiplemisc= filter(total_misc, n_misc>= 2) #mothers from PATREG with at least two separate miscarriages at any time (and in the columns when they happened)
 

#add the live births to the miscarriage info to be able to construct a manually created RPL variable (3 consecutive miscarriages with no live births in-between)

mfr <- mfr[, c("lopnr","BFODDAT")]

# Previous miscarriages
prev_misc= inner_join(temp, mfr, by= 'lopnr')

prev_misc$BFODDAT= as.Date(prev_misc$BFODDAT, format= '%Y%m%d')

prev_misc= filter(prev_misc, INDATUM< BFODDAT)
prev_misc= group_by(prev_misc, lopnr, BFODDAT) %>% filter(index== max(index))

# Get recurrent miscarriages

names(temp)= c('lopnr', 'date', 'index')
temp$live_misc= 0
temp$date= as.Date(temp$date, format= '%Y-%m-%d')
names(mfr)= c('lopnr', 'date')
mfr$date= as.Date(mfr$date, format= '%Y%m%d')
mfr$live_misc= 1

rec_misc= bind_rows(temp, mfr)

rec_misc= group_by(rec_misc, lopnr) %>% 
	arrange(date) %>%
	mutate(recurrent= cumsum(ifelse(!is.na(lag(lag(index))) & !is.na(lag((index))) & !is.na(((index))), 1, 0)))

names(rec_misc)[names(rec_misc) == 'index'] <- 'prev_misc'

rec_misc$recurrent= ifelse(rec_misc$recurrent > 1, 1, 0)

fwrite(rec_misc, snakemake@output[[1]], sep= '\t')
