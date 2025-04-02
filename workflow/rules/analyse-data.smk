rule format_main_data:
	'Format main data set.'
	input:
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/QIVF/SOS_IMPORTFIL_FINAL_230613.txt',
		'/mnt/hdd/data/swed/Graviditetsrelaterade_infektioner_BoJ/SoS/Medicinska_Födelseregistret/UT_R_MFR_14517_2022.txt'
	output:
		'results/main_data/QIVF-QC-own-oocytes.txt',
		'results/main_data/QIVF-QC-oocyte-recipients.txt'
	script:
		'../workflow/scripts/format-data.R'

rule plot_rates_by_age:
	'Plot of implantation failure and miscarriage rate per maternal age.'
	input:
		'results/main_data/QIVF-QC-{oocytes}.txt'
	output:
		'results/figures/rates-by-maternal-age-{oocytes}.png'
	script:
		'../workflow/scripts/rates-per-age.R'
	
	



