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

style_number_2digits <- purrr::partial(gtsummary::style_number, digits = 2)
style_percent_2digits <- purrr::partial(gtsummary::style_percent, digits = 2)

d= fread(snakemake@input[[1]])

d$implantation_failure= as.numeric(d$Resultfetus1=='')
d$early_miscarriage= ifelse(is.na(d$misc) | is.na(d$gest_duration), NA, ifelse(d$misc== 1 & d$gest_duration< 7*10, 1, 0))
d$cat_prev_misc= with(d, ifelse(is.na(prev_misc), NA, ifelse(prev_misc> 3, 3, prev_misc)))
d$Oocytowndonated= ifelse(d$Oocytowndonated== '', NA, d$Oocytowndonated)
d$Spermowndonated= ifelse(d$Spermowndonated== '', NA, d$Spermowndonated)
d$Performedtreatmenttype= ifelse(d$Performedtreatmenttype== '', NA, d$Performedtreatmenttype)

d$gest_duration= ifelse(d$gest_duration<= 15, NA, d$gest_duration)

df= filter(d, Resultfetus1 != '', gest_duration > 15, !is.na(gest_duration), !is.na(maternal_age),
           !is.na(misc))


if (!grepl('tertiles|multiple', snakemake@output[[1]])) {


### Descriptive characteristics

desc_all= d %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
	style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_all= d %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
	style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

# Only in implantation success


desc_succ= df %>%
  tbl_summary(
    include = c(maternal_age,  gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
	style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_succ= df %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
	style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

fwrite(as.data.frame(desc_all), snakemake@output[[1]], sep = '\t', col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_all), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_succ), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_succ), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)


} else if (grepl('tertiles', snakemake@output[[1]])) {

d$maternal_tertiles= ntile(d$maternal_age, 3)

df$maternal_tertiles= ntile(df$maternal_age, 3)

desc_all1= d %>%
  filter(maternal_tertiles== 1) %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_all1= d %>%
  filter(maternal_tertiles== 1) %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

desc_all2= d %>%
  filter(maternal_tertiles== 2) %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_all2= d %>%
  filter(maternal_tertiles== 2) %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

desc_all3= d %>%
  filter(maternal_tertiles== 3) %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_all3= d %>%
  filter(maternal_tertiles== 3) %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


# Only in implantation success


desc_succ1= df %>%
  filter(maternal_tertiles== 1) %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_succ1= df %>%
  filter(maternal_tertiles== 1) %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

desc_succ2= df %>%
  filter(maternal_tertiles== 2) %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_succ2= df %>%
  filter(maternal_tertiles== 2) %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

desc_succ3= df %>%
  filter(maternal_tertiles==3) %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_succ3= df %>%
  filter(maternal_tertiles== 3) %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

write("########## All embryo transfers#########", file = snakemake@output[[1]], append = TRUE)
write("#### First maternal age tertile\n\n\n", file = snakemake@output[[1]], append = TRUE)
fwrite(as.data.frame(desc_all1), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_all1), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)


write("\n\n\n#### Second maternal age tertile\n\n\n", file = snakemake@output[[1]], append = TRUE)
fwrite(as.data.frame(desc_all2), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_all2), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)

write("\n\n\n#### Third maternal age tertile\n\n\n", file = snakemake@output[[1]], append = TRUE)
fwrite(as.data.frame(desc_all3), snakemake@output[[1]], sep = '\t', append=TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_all3), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)

write("\n\n\n########## Only after successfull embryo implantation #########", file = snakemake@output[[1]], append = TRUE)
write("#### First maternal age tertile\n\n\n", file = snakemake@output[[1]], append = TRUE)
fwrite(as.data.frame(desc_succ1), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_succ1), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)

write("\n\n\n#### Second maternal age tertile\n\n\n", file = snakemake@output[[1]], append = TRUE)
fwrite(as.data.frame(desc_succ2), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_succ2), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)

write("\n\n\n#### Third maternal age tertile\n\n\n", file = snakemake@output[[1]], append = TRUE)
fwrite(as.data.frame(desc_succ3), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_succ3), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)

} else {
d$cat_prev_misc= factor(with(d, ifelse(is.na(prev_misc), NA, ifelse(prev_misc> 3, 3, prev_misc))))
df$cat_prev_misc= factor(with(df, ifelse(is.na(prev_misc), NA, ifelse(prev_misc> 3, 3, prev_misc))))

write("########## All embryo transfers#########", file = snakemake@output[[1]], append = TRUE)

for (i in 0:3) {
tertile= ifelse(i == 0, 'No previous', ifelse(i == 1, 'one previous', ifelse(i == 2, 'two previous', 'three or more previous')))

write(paste("\n\n\n####", tertile, "pregnancy losses\n\n\n"), file = snakemake@output[[1]], append = TRUE)

temp_d= d[d$cat_prev_misc== i, ]
temp_df= df[df$cat_prev_misc== i, ]

desc_all= temp_d %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_all= d %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

fwrite(as.data.frame(desc_all), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_all), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)

}
# Only in implantation success

write("\n\n\n########## Only after successfull embryo implantation #########", file = snakemake@output[[1]], append = TRUE)

for (i in 0:3) {
tertile= ifelse(i == 0, 'No previous', ifelse(i == 1, 'one previous', ifelse(i == 2, 'two previous', 'three or more previous')))

write(paste("\n\n\n####", tertile, "pregnancy losses\n\n\n"), file = snakemake@output[[1]], append = TRUE)

temp_d= d[d$cat_prev_misc== i, ]
temp_df= df[df$cat_prev_misc== i, ]

desc_succ= df %>%
  tbl_summary(
    include = c(maternal_age, gest_duration), # your continuous variables
statistic = list(all_continuous() ~ "{mean}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_continuous() ~ "t.test"),
        style_fun= list(gtsummary::all_continuous() ~ style_number_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")


desc_cat_succ= df %>%
  tbl_summary(
    include = c(Incubationdays, implantation_failure, Oocytowndonated, Spermowndonated, Performedtreatmenttype, Embryotreatmenttype, misc, early_miscarriage, cat_prev_misc),
    statistic = list(all_categorical() ~ "{p}"),
digits = list(everything() ~ c(2))) %>%
  add_ci(method = list(all_categorical() ~ "wilson"),
        style_fun= list(gtsummary::all_categorical() ~ style_percent_2digits)) %>%
  gtsummary::modify_header(ci_stat_0 ~ "**95% CI**")

fwrite(as.data.frame(desc_succ), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
fwrite(as.data.frame(desc_cat_succ), snakemake@output[[1]], sep = '\t', append= TRUE, col.names = TRUE, row.names = FALSE)
} 

} 


