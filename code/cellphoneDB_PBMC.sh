#!/bin/sh

source /hdd/yuan/venv/cpdb/bin/activate
nohup cellphonedb method statistical_analysis Normal_female_cpdb_meta.txt Normal_female_cpdb_count.txt --counts-data=gene_name --output-path=./out/Normal_female > my_out_Normal_female.log 2>&1 &
nohup cellphonedb method statistical_analysis Normal_male_cpdb_meta.txt Normal_male_cpdb_count.txt --counts-data=gene_name --output-path=./out/Normal_male > my_out_Normal_male.log 2>&1 &
nohup cellphonedb method statistical_analysis COV_Severe_female_cpdb_meta.txt COV_Severe_female_cpdb_count.txt --counts-data=gene_name --output-path=./out/COV_Severe_female > my_out_COV_Severe_female.log 2>&1 &
nohup cellphonedb method statistical_analysis COV_Severe_male_cpdb_meta.txt COV_Severe_male_cpdb_count.txt --counts-data=gene_name --output-path=./out/COV_Severe_male > my_out_COV_Severe_male.log 2>&1 &
nohup cellphonedb method statistical_analysis COV_Mild_female_cpdb_meta.txt COV_Mild_female_cpdb_count.txt --counts-data=gene_name --output-path=./out/COV_Mild_female > my_out_COV_Mild_female.log 2>&1 &
nohup cellphonedb method statistical_analysis COV_Mild_male_cpdb_meta.txt COV_Mild_male_cpdb_count.txt --counts-data=gene_name --output-path=./out/COV_Mild_male > my_out_COV_Mild_male.log 2>&1 &

