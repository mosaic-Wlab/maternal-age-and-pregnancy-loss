rule fetch_pregn_loss:
        'Generate data for recurrent and multiple miscarriage.'
        input:
                '/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt',
                '/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt',
                '/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Patientregistret/UT_R_PAR_OV_M_14517_2022.txt'
        output:
                'results/main_data/recurrent-multiple-loss.txt'
        conda:
                 '../envs/plots.yml'
        script:
                '../scripts/fetch-recurrent-loss.R'

rule format_main_data:
	'Format main data set.'
	input:
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt',
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt',
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Patientregistret/UT_R_PAR_OV_M_14517_2022.txt',
		'results/main_data/recurrent-multiple-loss.txt'
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

rule glm_mixed_model:
	'Implement mixed glm for associations between maternal age and history of miscarriage on implantation and miscarriage.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt',
	output:
                'results/figures/mixed-glm-{oocytes}.pdf',
                'results/figures/data/mixed-glm-{oocytes}.txt',
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/mixed-glm-implantation-preg-loss.R'

rule glm_mixed_model_plots:
	'GLM mixed models for associations between maternal age and biochemical loss among all losses.'
	input:
		'results/main_data/QIVF-QC-own-oocyte.txt',
		'results/main_data/QIVF-QC-oocyte-recipient.txt'
	output:
                'results/figures/mixed-glm-maternal_age.pdf',
                'results/figures/data/mixed-glm-maternal_age_miscarriages.txt',	
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/emmeans_glmm_biochemical_loss.R'

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

rule plot_biochemical_vs_clinical_loss:
        'Plot of biochemical to clinical loss by maternal age.'
        input:
                'results/main_data/QIVF-QC-{oocytes}.txt'
        output:
                'results/figures/biochemical-to-clinical-loss-by-maternal-age-{oocytes}.pdf',
                'results/figures/data/glm-biochemical-to-clinical-loss-by-maternal-age-{oocytes}.txt',
                'results/figures/data/lm-gest_duration-by-maternal-age-{oocytes}.txt',
                'results/figures/lm-gest_duration-by-maternal-age-{oocytes}.pdf'
        conda:
                '../envs/plots.yml'
        script:
                '../scripts/lm-misc.R'

rule plot_biochemical_vs_clinical_loss_standard:
        'Plot of biochemical to clinical loss by maternal age.'
        input:
                'results/main_data/QIVF-QC-{oocytes}.txt'
        output:
                'results/figures/standard-biochemical-to-clinical-loss-by-maternal-age-{oocytes}.pdf',
                'results/figures/data/standard-glm-biochemical-to-clinical-loss-by-maternal-age-{oocytes}.txt',
                'results/figures/data/standard-lm-gest_duration-by-maternal-age-{oocytes}.txt',
                'results/figures/standard-lm-gest_duration-by-maternal-age-{oocytes}.pdf'
        conda:
                '../envs/plots.yml'
        script:
                '../scripts/lm-misc.R'


rule plot_biochemical_vs_clinical_loss_male_infer:
        'Plot of biochemical to clinical loss by maternal age.'
        input:
                'results/main_data/QIVF-QC-{oocytes}.txt'
        output:
                'results/figures/{male_infer}-biochemical-to-clinical-loss-by-maternal-age-{oocytes}.pdf',
                'results/figures/data/glm-biochemical-to-clinical-loss-by-maternal-age-{oocytes}-{male_infer}.txt',
                'results/figures/data/lm-gest_duration-by-maternal-age-{oocytes}-{male_infer}.txt',
                'results/figures/lm-gest_duration-by-maternal-age-{male_infer}-{oocytes}.pdf'
        conda:
                '../envs/plots.yml'
        script:
                '../scripts/lm-misc.R'

