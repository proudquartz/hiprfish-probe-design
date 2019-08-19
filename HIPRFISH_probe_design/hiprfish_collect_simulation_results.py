"""
Collect HiPRFISH probe design results
Hao Shi 2019
De Vlaminck Lab
Cornell University
"""

import os
import argparse
import pandas as pd
import numpy as np

###############################################################################################################
# helper functions
###############################################################################################################

def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    array = array.flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1,array.shape[0]+1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array)))

def get_cumulative_coverage(taxon_best_probes, taxon_best_probes_filtered):
    n_total = taxon_best_probes.taxon_abundance.drop_duplicates().sum()
    taxon_best_probes['taxon_coverage_absolute'] = taxon_best_probes.taxon_coverage.values*taxon_best_probes.taxon_abundance.values
    taxon_best_probes_sorted = taxon_best_probes.sort_values(['taxon_coverage_absolute'], ascending = False)
    taxon_best_probes_sorted_filtered = pd.merge(taxon_best_probes_filtered.ix[:,['target_taxon_full']], taxon_best_probes, how = 'left', on = 'target_taxon_full')
    taxon_best_probes_sorted_filtered.sort_values(['taxon_coverage_absolute'], ascending = False, inplace = True)
    taxon_cumsum = taxon_best_probes_sorted_filtered.taxon_coverage_absolute.drop_duplicates().cumsum()
    cumcov = taxon_cumsum.values/n_total
    return(cumcov)

def collect_probe_coverage_results(data_dir, simulation_table, output_filename):
    print('Loading samples table: %s' % (simulation_table))
    sim_tab = pd.read_csv(simulation_table)
    sim_tab['Gini_Index'] = np.nan
    sim_tab['Inverse_Simpson_Index'] = np.nan
    sim_tab['Shannon_Index'] = np.nan
    sim_tab['COVERED_TAXA_RICHNESS'] = np.nan
    sim_tab['COVERED_TAXA_RICHNESS_FRACTION'] = np.nan
    sim_tab['CUMCOV_all'] = np.nan
    sim_tab['CUMCOV_TOP10'] = np.nan
    sim_tab['CUMCOV_TOP20'] = np.nan
    sim_tab['OFF_TARGET_PHYLUM'] = np.nan
    sim_tab['OFF_TARGET_CLASS'] = np.nan
    sim_tab['OFF_TARGET_ORDER'] = np.nan
    sim_tab['OFF_TARGET_FAMILY'] = np.nan
    sim_tab['OFF_TARGET_GENUS'] = np.nan
    sim_tab['OFF_TARGET_SPECIES'] = np.nan
    print('Loading result files:')
    for i in range(0, sim_tab.shape[0]):
        sample = sim_tab.SAMPLE[i]
        target_rank = sim_tab.TARGET_RANK[i]
        similarity = sim_tab.SIMILARITY[i]
        taxon_best_probes_filename = data_dir + '/simulation/%s/taxon_best_probes.csv' % (sim_tab.DESIGN_ID[i])
        taxon_best_probes_filtered_filename = data_dir + '/simulation/%s/taxon_best_probes_filtered.csv' % (sim_tab.DESIGN_ID[i])
        if os.path.exists(taxon_best_probes_filename):
            taxon_best_probes = pd.read_csv(taxon_best_probes_filename)
            taxon_best_probes_filtered = pd.read_csv(taxon_best_probes_filtered_filename)
            sim_tab.COVERED_TAXA_RICHNESS.loc[i] = taxon_best_probes_filtered.target_taxon.drop_duplicates().shape[0]
            sim_tab.COVERED_TAXA_RICHNESS_FRACTION.loc[i] = taxon_best_probes_filtered.target_taxon.drop_duplicates().shape[0]/taxon_best_probes.target_taxon.drop_duplicates().shape[0]
            cumcov = get_cumulative_coverage(taxon_best_probes, taxon_best_probes_filtered)
            sim_tab.CUMCOV_all.loc[i] = cumcov[-1]
            if cumcov.shape[0] >= 10:
                sim_tab.CUMCOV_TOP10.loc[i] = cumcov[9]
            if cumcov.shape[0] >= 20:
                sim_tab.CUMCOV_TOP20.loc[i] = cumcov[19]
            print('Saving collected results to %s...' % (output_filename))
        else:
            print('Sample result file %s does not exist' % taxon_best_probes_filename)
        taxon_abundance_filename = '{}/{}/{}/s_{}/consensus/taxon_abundance.csv'.format(data_dir, sample, target_rank, similarity)
        taxon_abundance = pd.read_csv(taxon_abundance_filename)
        taxon_abundance['rel_freq'] = taxon_abundance.counts.values/taxon_abundance.counts.sum()
        sim_tab.loc[i, 'Gini_Index'] = gini(taxon_abundance.rel_freq.sort_values(ascending = True).values)
        sim_tab.loc[i, 'Inverse_Simpson_Index'] = 1/(np.sum(taxon_abundance.rel_freq**2))
        sim_tab.loc[i, 'Shannon_Index'] = -np.sum(taxon_abundance.rel_freq*np.log(taxon_abundance.rel_freq))
        sim_tab.to_csv(output_filename, index = False, header = True)
    return(sim_tab)

###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Collect summary statistics of HiPRFISH probes for a complex microbial community')

    # data directory
    parser.add_argument('data_dir', type = str, help = 'Directory of the data files')

    # input simulation table
    parser.add_argument('simulation_table', type = str, help = 'Input csv table containing simulation information')

    # output simulation results table
    parser.add_argument('simulation_results', type = str, help = 'Output csv table containing simulation results')

    args = parser.parse_args()

    sim_tab = collect_probe_coverage_results(args.data_dir, args.simulation_table, args.simulation_results)

if __name__ == '__main__':
    main()
