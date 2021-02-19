# spite-SB-internship

These scripts have been written by Pedro Batista Tan, for his major internship of the master program Systems biology and Bioinformatics, from the Vrije University of Amsterdam (VU) and Universiteit van Amsterdam (UvA). The code for the models has been adapted from a previous version written by Silvia Espada Burriel. It has been written with the supervision of Rutger Hermsen, from the Theoretical biology group of the University of Utrecht.

"stepwise_spite_model.py" and "dichotomous_spite_model.py" are python 3.8 scripts which perform simulations on a minimal models of spite evolution and selection inspired by allelopathic bacteriocin production. For details on the models, please check the internship report. Running these scripts will generate folders with results for each parameter combination, determined with dictionaries within the scripts. 

Simulations results are parsed with the "stepwise_spite_model_analysis_get_stats.py" and "dichotomous_spite_model_analysis_get_stats.py" scripts. These take a results folder as input and output a csv file containing parsed results for all experiments within.

The R script "Spite_plot_results.R" was used to generate figures for the report. It uses csv files obtained with the "get_stats" python scripts. The csv files containing results used for the main figures are available at the "Results_csv" folder.

Fig 1
Spite_results_df_uu.csv

Fig 2
Spite_results_df_uu_diff.csv
Spite_results_uu_selfspite.csv
Spite_results_selfspite_res.csv

Fig 3
Results_dichotomous_l.csv

Fig 4
Results_dichotomous_l_freq.csv

Fig 5
Results_dichotomous_l_kd.csv

Fig 6
Results_dichotomous_l_kc_kd.csv

Fig 7
Results_dichotomous_l_kc_freq.csv
