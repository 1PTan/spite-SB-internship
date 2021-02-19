"""This code has been written by Pedro Batista Tan , for his major internship of the master program Systems biology and
Bioinformatics, from the Vrije University of Amsterdam (VU) and Universiteit van Amsterdam (UvA). The code has been
written with the supervision of Rutger Hermsen, from the Theoretical biology group of the University of Utrecht.

The objective is to parse simulation results from multiple runs in an input folder and output a csv file containing the
final carrier frequency for each experiment. A naming convention for folders was followed to parse simulation parameters
and include them in the csv.

input: main folder containing simulation results of multiple experiments
output: csv file containing parsed simulation results and parameters
"""
# IMPORTS ---------------------------------------------------------------------
import sys
import numpy as np
from scipy.stats import norm
from time import time
from datetime import timedelta
import os
import pathlib
import matplotlib.pyplot as plt
import pandas as pd


results_folder = "Results_uu"  # Input folder containing simulation results
export_file_name = "Spite_results_df_uu.csv"  # Output name for csv file
process_all_folders = True  # Choose whether to include all simulation results from the input folder

# If process_all_folders is set to false, a list of the folder names to parse can be provided instead
folder_list = ["d4_c8_s8_cost0.005"]

if process_all_folders:
    directory = pathlib.Path(__file__).parent.absolute()
    directory = str(directory) + "/" + results_folder
    folder_list = [f.name for f in os.scandir(directory) if f.is_dir()]

# Create a dictionary structure that will contain the parsed results and convert it into a pandas data frame
stats_dict = {'Experiment': [], "spite_sigma": [], "comp_sigma": [], "diff_sigma": [], "cost": [],
              'average_n_individuals': [],
              'average_mean_spite': [], "average_max_spite": [], "average_sd_spite": [],
              "average_mean_density": [], "final_spite_10_quantile": [], "final_spite_25_quantile": [],
              "final_spite_50_quantile": [], "final_spite_75_quantile": [], "final_spite_90_quantile": []}
df = pd.DataFrame(stats_dict)

for folder in folder_list:
    print("Processing results for folder: ", folder)
    # Parse the relevant simulation parameters from the file name and print them below
    folder_name = folder.split("\\")[-1]
    folder_split = folder_name.split("_")
    cost = float(folder_split[3][4:])
    comp_sigma = int(folder_split[1][1:])
    spite_sigma = int(folder_split[2][1:])
    try:
        diff_sigma = int(folder_split[0][1:])
    except:
        diff_sigma = folder_split[0][1:]

    print(f"Parsed spite sigma: {spite_sigma}")
    print(f"Parsed comp sigma: {comp_sigma}")
    print(f"Parsed spite sigma: {diff_sigma}")
    print(f"Parsed cost: {cost}")

    # load simulation results from the folder
    path = f"{results_folder}/{folder_name}/"
    stats = np.load(f"{path}stats.npz")
    max_t = len(stats['t']) - 1
    final_field = np.load(f"{path}time{max_t}.npz")

    # Parse results, averaging over the last time steps
    average_over_last = 150000

    average_mean_spite = np.mean(stats['mp'][max_t - average_over_last:max_t])
    average_max_spite = np.mean(stats['maxp'][max_t - average_over_last:max_t])
    average_sd_spite = np.mean(stats['std'][max_t - average_over_last:max_t])
    average_n_individuals = np.mean(stats['n'][max_t - average_over_last:max_t])

    average_mean_density = np.mean(stats['d'][max_t - average_over_last:max_t])
    final_spite_quantiles = np.percentile(final_field['p'], [10, 25, 50, 75, 90])
    final_spite_10_quantile = final_spite_quantiles[0]
    final_spite_25_quantile = final_spite_quantiles[1]
    final_spite_50_quantile = final_spite_quantiles[2]
    final_spite_75_quantile = final_spite_quantiles[3]
    final_spite_90_quantile = final_spite_quantiles[4]

    # Format folder results in the same dictionary structure
    folder_stats = {'Experiment': folder,
                    "spite_sigma": spite_sigma, "comp_sigma": comp_sigma, "diff_sigma": diff_sigma,
                    "cost": cost, "average_n_individuals": average_n_individuals,
                    "average_mean_spite": average_mean_spite, "average_max_spite": average_max_spite,
                    "average_sd_spite": average_sd_spite, "average_mean_density": average_mean_density,
                    "final_spite_10_quantile": final_spite_10_quantile,
                    "final_spite_25_quantile": final_spite_25_quantile,
                    "final_spite_50_quantile": final_spite_50_quantile,
                    "final_spite_75_quantile": final_spite_75_quantile,
                    "final_spite_90_quantile": final_spite_90_quantile}

    # Add results to the pandas data frame containing all results
    folder_df = pd.DataFrame(folder_stats, index=[0])
    df = pd.concat([df, folder_df], axis=0)

# Export data frame to csv
df.to_csv(export_file_name, index=False)




