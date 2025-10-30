rule fetch_pregn_loss:
        'Generate data for recurrent and multiple miscarriage.'
        input:
                '/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt',
                '/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt',
                '/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Patientregistret/UT_R_PAR_OV_M_14517_2022.txt'
        output:
                'results/main_data/recurrent-multiple-loss.txt'
        script:
                '../scripts/fetch-recurrent-loss.R'

rule format_main_data:
	'Format main data set.'
	input:
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt',
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt'
	output:
		'results/main_data/QIVF-QC-own-oocyte.txt',
		'results/main_data/QIVF-QC-oocyte-recipient.txt',
		'results/main_data/all-QIVF-QC.txt',
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/format-data.R'

rule plot_rates_by_age:
	'Plot of implantation failure and miscarriage rate per maternal age.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt'
	output:
		'results/figures/rates-by-maternal-age-{oocytes}.pdf',
		'results/figures/data/rates-by-maternal-age-{oocytes}-miscarriages.txt',
		'results/figures/data/descriptive-total-{oocytes}-miscarriages.txt',
		'results/figures/data/descriptive-implantation-success-{oocytes}-miscarriages.txt',
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/rates-per-age.R'

rule plot_survival_time:
	'Plot of survival rate for time-to-miscarriage.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt'
	output:
		'results/figures/raw-survival-rate-{oocytes}.pdf',
		'results/figures/raw-survival-rate-{oocytes}-only-miscarriages.pdf',
		'results/figures/data/raw-survival-rate-{oocytes}.txt',
		'results/figures/data/descriptive-{oocytes}.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/survival-rate.R'


rule plot_survival_age:
	'Survival rates by maternal age tertiles.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt'
	output:
		'results/figures/survival-rate-maternal-age-{oocytes}.pdf',
		'results/figures/data/survival-rate-maternal-age-{oocytes}.txt',
		'results/figures/data/time-varying-zph-maternal-age-{oocytes}.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/survival-rate.R'

rule pamm_models:
	'Analyse time-to-event data using PAMM'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt'
	output:
		'results/figures/time-varying-PAMM-{oocytes}.pdf',
		'results/figures/data/time-varying-PAMM-{oocytes}-parametric.txt',
		'results/figures/data/time-varying-PAMM-{oocytes}-smooth.txt',
		'results/figures/data/PAMM-{oocytes}-data.txt'
	conda:
                '../envs/plots.yml'
	script:
		'../scripts/time-varying-models.R'

rule pamm_models_male_infertility:
	'Analyse timeto-event data using PAMM in cases with male infertility.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt',
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Patientregistret/UT_R_PAR_OV_M_14517_2022.txt'
	output:
		'results/figures/male-infertility-time-varying-PAMM-{oocytes}.pdf',
		'results/figures/data/male-infertility-time-varying-PAMM-{oocytes}-parametric.txt',
		'results/figures/data/male-infertility-time-varying-PAMM-{oocytes}-smooth.txt',
		'results/figures/data/male-infertility-PAMM-{oocytes}-data.txt',
		'results/figures/male-infertility-boxplot-{oocytes}.pdf',
		'results/figures/male-infertility-survival-rate-{oocytes}.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/male-infertility.R'

rule pamm_models_recurrent_loss:
	'Analyse time-to-event data using PAMM in cases with recurrent pregnancy loss or multiple miscarriage.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt',
		'results/main_data/recurrent-multiple-loss.txt',
	output:
		'results/figures/time-varying-PAMM-recurrent-{oocytes}.pdf',
		'results/figures/data/time-varying-PAMM-recurrent-{oocytes}-parametric.txt',
		'results/figures/data/time-varying-PAMM-recurrent-{oocytes}-smooth.txt',
		'results/figures/data/PAMM-recurrent-{oocytes}-data.txt',
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/time-varying-recurrent-loss.R'

rule plot_survival_recurrent:
	'Survival rates by maternal age tertiles.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt',
		'results/main_data/recurrent-multiple-loss.txt'
	output:
		'results/figures/survival-rate-recurrent-{oocytes}.pdf',
		'results/figures/data/survival-rate-recurrent-{oocytes}.txt',
		'results/figures/data/time-varying-zph-recurrent-{oocytes}.txt',
		'results/figures/data/recurrent-descriptive-{oocytes}-data.txt',
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/recurrent-survival-rate.R'

rule describe_all:
	''
	input:
		'results/main_data/all-QIVF-QC.txt',
	output:
		'results/description/descriptive-stats-all.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/description.R'

rule describe_data_by_oocytes:
	''
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt',
	output:
		'results/description/descriptive-stats-{oocytes}.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/description.R'

rule describe_by_exposure:
	''
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt',
	output:
		'results/description/descriptive-stats-{exposure}-{oocytes}.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/description.R'
