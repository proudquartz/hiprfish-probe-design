"""
Collect HiPRFISH probe design results
Hao Shi 2019
De Vlaminck Lab
Cornell University
"""

import argparse
import os
import re
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO

###############################################################################################################
# helper functions
###############################################################################################################

def collect_taxon_best_probes(design_directory, sim_input_filename, taxon_best_probes_filename, taxon_best_probes_filtered_filename, output_probes_summary_filename, bot):
    simulation_directory, design_id = os.path.split(design_directory)
    data_dir = os.path.split(simulation_directory)[0]
    sim_tab = pd.read_csv(sim_input_filename)
    sample = sim_tab[sim_tab.DESIGN_ID == design_id].SAMPLE.values[0]
    target_rank = sim_tab[sim_tab.DESIGN_ID == design_id].TARGET_RANK.values[0]
    similarity = sim_tab[sim_tab.DESIGN_ID == design_id].SIMILARITY.values[0]
    taxon_evaluation_directory = data_dir + '/%s/%s/s_%s/evaluation/' % (sample, target_rank, str(similarity))
    cluster_lookup_filename = data_dir + '/%s/%s/s_%s/consensus/cluster_lookup.tab' % (sample, target_rank, str(similarity))
    taxon_evaluation_filename_list = glob.glob(design_directory + '/*_probe_selection.csv')
    best_probes_list = [pd.read_csv(filename) for filename in taxon_evaluation_filename_list]
    best_probes_quality_sorted_list = [df.sort_values(by = ['quality'], ascending = True) for df in best_probes_list]
    best_probes_df = pd.concat(best_probes_quality_sorted_list)
    best_probes_summary = best_probes_df['target_taxon_full'].value_counts().reset_index()
    best_probes_summary.columns = ['strain', 'probe_counts']
    cluster_lookup = pd.read_csv(cluster_lookup_filename, sep = '\t')
    best_probes_cluster_summary = best_probes_summary.merge(cluster_lookup, on = 'strain')
    best_probes_cluster_summary = best_probes_cluster_summary[['molecule_id', 'strain', 'probe_counts']]
    best_probes_filtered = best_probes_df[best_probes_df['blast_on_target_rate'] > bot]
    best_probes_filtered_summary = best_probes_filtered['target_taxon_full'].value_counts()
    best_probes_df.to_csv(taxon_best_probes_filename, index = False)
    best_probes_filtered.to_csv(taxon_best_probes_filtered_filename, index = False)
    best_probes_filtered_summary.to_csv(output_probes_summary_filename)
    return

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Collect summary statistics of HiPRFISH probes for a complex microbial community')

    parser.add_argument('design_directory', type = str, help = 'Directory to the folder for each design')

    parser.add_argument('sim_input_filename', type = str, help = 'Input csv table containing simulation information')

    parser.add_argument('taxon_best_probes_filename', type = str, help = 'CSV file that include all the best probes for each taxon')

    parser.add_argument('taxon_best_probes_filtered_filename', type = str, help = 'CSV file that include all the best probes for each taxon filtered by blast on target rate')

    parser.add_argument('output_probes_summary_filename', type = str, help = 'Output csv table containing simulation results')

    parser.add_argument('bot', type = float, help = 'Minimum blast on target rate')

    args = parser.parse_args()

    collect_taxon_best_probes(args.design_directory, args.sim_input_filename, args.taxon_best_probes_filename, args.taxon_best_probes_filtered_filename, args.output_probes_summary_filename, args.bot)

if __name__ == '__main__':
    main()
