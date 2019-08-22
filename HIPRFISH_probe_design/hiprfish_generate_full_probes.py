"""
Collect HiPRFISH probe design results
Hao Shi 2019
De Vlaminck Lab
Cornell University
"""

import argparse
import pandas as pd
import os
import glob
import re
import numpy as np
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from matplotlib import pyplot as plt

###############################################################################################################
# helper functions
###############################################################################################################

def add_spacer(taxon_best_probes, consensus_directory, output_filename):
    probes = pd.read_csv(taxon_best_probes)
    probes['seqrcsa'] = ''
    for i in range(0, probes.shape[0]):
        probe_seq = Seq(probes['seq'][i], generic_dna)
        rrna_file = consensus_directory + '/' + str(probes['target_taxon'][i]) + '.consensus.fasta'
        rrna_file_length = sum(1 for record in SeqIO.parse(rrna_file, 'fasta'))
        if rrna_file_length > 1:
            cluster = re.sub('.*_', '', probes['target_taxon_full'][i])
            rrna_seq = SeqIO.to_dict(SeqIO.parse(rrna_file, 'fasta'))[cluster].seq
        else:
            rrna_seq = SeqIO.read(rrna_file, 'fasta').seq
        sstart = probes['p_start'][i]
        probe_length = probes['ln'][i]
        send = sstart + probe_length - 1
        right_spacer = rrna_seq[sstart - 2] + rrna_seq[sstart - 3] + rrna_seq[sstart - 4]
        left_spacer = rrna_seq[send + 2] + rrna_seq[send + 1] + rrna_seq[send]
        probe_seq_rcsa = str(left_spacer).upper() + str(probe_seq.reverse_complement()) + str(right_spacer).upper()
        probes.ix[i, 'seqrcsa'] = probe_seq_rcsa
    probes.loc[:,'quadg'] = (probes['seqrcsa'].str.upper().str.count('GGGG')) + (probes['seqrcsa'].str.upper().str.count('GGGGG'))
    probes = probes[probes['quadg'] == 0]
    probes.to_csv(output_filename, index = False)
    return

def convert_numeric_barcode(num, nbit):
    code = re.sub('0b', '', format(num, '#0' + str(nbit+2) + 'b'))
    return(code)

def count_number_of_bits(binary_barcode):
    bin_list = list(binary_barcode)
    bin_int_list = [int(index) for index in bin_list]
    return(np.sum(bin_int_list))

def generate_full_probes(design_dir, bot, plf = 'T', primer = 'T', primerset = 'B', barcode_selection = 'MostSimple'):
    design_id = os.path.basename(design_dir)
    RP = ['GATGATGTAGTAGTAAGGGT',
          'AGGTTAGGTTGAGAATAGGA',
          'AGGGTGTGTTTGTAAAGGGT',
          'TTGGAGGTGTAGGGAGTAAA',
          'AGAGTGAGTAGTAGTGGAGT',
          'ATAGGAAATGGTGGTAGTGT',
          'TGTGGAGGGATTGAAGGATA']
    nbit = len(RP)
    if plf == 'T':
        probe_length_filler = 'GTCTATTTTCTTATCCGACG'
    else:
        probe_length_filler = ''
    if primer == 'T':
        if primerset == 'A':
            forward_primer = 'CGCGGGCTATATGCGAACCG'
            reverse_primer = 'GCGTTGTATGCCCTCCACGC'
            # TAATACGACTCACTATAGGGCGTGGAGGGCATACAACGC
        elif primerset == 'B':
            forward_primer = 'CGATGCGCCAATTCCGGTTC'
            reverse_primer = 'CAACCCGCGAGCGATGATCA'
            # TAATACGACTCACTATAGGGTGATCATCGCTCGCGGGTTG
        elif primerset == 'C':
            forward_primer = 'GTTGGTCGGCACTTGGGTGC'
            reverse_primer = 'CCACCGGATGAACCGGCTTT'
            # TAATACGACTCACTATAGGGAAAGCCGGTTCATCCGGTGG
    else:
        forward_primer = ''
        reverse_primer = ''
    taxon_best_probes_sa_filenames = glob.glob('{}/*_probe_selection_sa.csv'.format(design_dir))
    if '{}/0_probe_selection_sa.csv'.format(design_dir) in taxon_best_probes_sa_filenames:
        taxon_best_probes_sa_filenames.remove('{}/0_probe_selection_sa.csv'.format(design_dir))
    taxon_best_probes_sa_filenames.sort()
    taxon_best_probes_list = [pd.read_csv(filename) for filename in taxon_best_probes_sa_filenames]
    taxon_best_probes_filtered_list = [x for x in taxon_best_probes_list if x.blast_on_target_rate.max() > bot]
    oligo_list = []
    blocking_probe_list = []
    barcodes = pd.DataFrame(np.arange(1,2**nbit))
    barcodes.columns = ['NumericBarcode']
    barcodes['BinaryBarcode'] = barcodes.NumericBarcode.apply(convert_numeric_barcode, args = (7,))
    barcodes['NumberBits'] = barcodes.BinaryBarcode.apply(count_number_of_bits)
    if barcode_selection == 'MostSimple':
        barcodes_sorted = barcodes.sort_values(by = ['NumberBits', 'NumericBarcode'])
    elif barcode_selection == 'MostComplex':
        barcodes_sorted = barcodes.sort_values(by = ['NumberBits', 'NumericBarcode'], ascending = [False, False])
    elif barcode_selection == 'Random':
        barcodes_sorted = barcodes.reindex(np.random.permutation(barcodes.index))
    elif barcode_selection == 'Sequential':
        barcodes_sorted = barcodes.copy()
    else:
        barcode_selection = 'Sequential'
        barcodes_sorted == barcodes.copy()
    for i in range(len(taxon_best_probes_filtered_list)):
        probes = taxon_best_probes_filtered_list[i]
        probes = probes[(probes['blast_on_target_rate'] > bot) & (probes['off_target_full_qcovhsp_fraction'] < 0.01)]
        assigned = barcodes_sorted.NumericBarcode.values[i]
        if barcodes_sorted.NumberBits.values[i] > 2:
            barcode_repeat = np.round(15/barcodes_sorted.NumberBits.values[i]).astype(int)
        else:
            barcode_repeat = 15
        if probes.shape[0] > 0:
            for k in range(barcode_repeat):
                probes = probes.sort_values(by = 'quality', ascending = True).reset_index().drop(columns = ['index'])
                for i in range(0, probes.shape[0]):
                    tarseq = probes['seqrcsa'][i]
                    if 'N' not in list(tarseq):
                        code = np.asarray(list(re.sub('0b', '', format(assigned, '#0' + str(7+2) + 'b'))), dtype = np.int64)
                        indices = np.where(code == 1)[0]
                        if len(indices) > 2:
                            indices = np.append(indices, indices[0])
                            subcode = np.zeros((len(indices)-1, nbit), dtype = np.int64)
                            for j in range(0, len(indices) - 1):
                                subcode[j, indices[j]] = 1
                                subcode[j, indices[j+1]] = 1
                            oligo = [[str(probes['target_taxon_full'][i]), probes['p_start'][i], probes['ln'][i], probes['taxon_abundance'][i], tarseq, probes['Tm'][i], probes['GC'][i], re.sub('\[|\]| ','',str(code)), assigned, probes['probe_id'][i], str(probes['target_taxon_full'][i]) + '_' + str(assigned) + '_' + re.sub('\[|\]| ','',str(subcode[k])) + '_' + str(probes['probe_id'][i]), forward_primer + RP[np.where(subcode[k] == 1)[0][0]] + tarseq + RP[np.where(subcode[k] == 1)[0][1]] + reverse_primer] for k in range(0, len(subcode))]
                        elif len(indices) == 2:
                            oligo = [[str(probes['target_taxon_full'][i]), probes['p_start'][i], probes['ln'][i], probes['taxon_abundance'][i], tarseq, probes['Tm'][i], probes['GC'][i], re.sub('\[|\]| ','',str(code)), assigned, probes['probe_id'][i], str(probes['target_taxon_full'][i]) + '_' + str(assigned) + '_' + re.sub('\[|\]| ','', str(code)) + '_' + str(probes['probe_id'][i]), forward_primer + RP[np.where(code == 1)[0][0]] + tarseq + RP[np.where(code == 1)[0][1]] + reverse_primer]]
                        else:
                            oligo = [[str(probes['target_taxon_full'][i]), probes['p_start'][i], probes['ln'][i], probes['taxon_abundance'][i], tarseq, probes['Tm'][i], probes['GC'][i], re.sub('\[|\]| ','',str(code)), assigned, probes['probe_id'][i], str(probes['target_taxon_full'][i]) + '_' + str(assigned) + '_' + re.sub('\[|\]| ','',str(code)) + '_' + str(probes['probe_id'][i]), forward_primer + RP[np.where(code == 1)[0][0]] + tarseq + probe_length_filler + reverse_primer]]
                        oligo_list = oligo_list + oligo
    oligo_df = pd.DataFrame(oligo_list)
    print(oligo_df.shape)
    oligo_df.columns = ['target_taxon', 'p_start', 'ln', 'abundance', 'rna_seq', 'Tm', 'GC', 'code', 'numeric_code', 'probe_numeric_id', 'probe_id', 'probe_full_seq']
    return(oligo_df)

def generate_blocking_probes(design_dir, bplc, target_rank, plf = 'T', primer = 'T', primerset = 'B', barcode_selection = 'MostSimple'):
    design_dir_folder = os.path.split(design_dir)[1]
    blocking_probes_filenames = glob.glob('{}/*_off_target_summary_info.csv'.format(design_dir))
    probes_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    probes = pd.read_csv(probes_filename)
    if '{}/0_off_target_summary_info.csv'.format(design_dir) in blocking_probes_filenames:
        blocking_probes_filenames.remove('{}/0_off_target_summary_info.csv'.format(design_dir))
    blocking_probes_list = []
    for f in blocking_probes_filenames:
        target_taxon = re.sub('_off_target_summary_info.csv', '', os.path.basename(f))
        if not probes.loc[probes.target_taxon.values.astype(str) == target_taxon, :].empty:
            off_target_summary_full = pd.read_csv(f)
            if not off_target_summary_full.empty:
                off_target_summary_full.loc[:,'max_encoding_interference_fraction'] = 0.0
                current_taxon_probes = probes.loc[probes.target_taxon.values.astype(str) == target_taxon, :]
                for bp in current_taxon_probes.probe_numeric_id.drop_duplicates().values:
                    off_target_summary = off_target_summary_full.loc[off_target_summary_full.probe_id.values == bp, :]
                    blocked_taxa = off_target_summary[target_rank].drop_duplicates()
                    for tt in blocked_taxa.values:
                        taxon_probes = probes.loc[probes.target_taxon.values == tt, :]
                        off_target_summary_blocked_taxon = off_target_summary.loc[off_target_summary[target_rank].values == tt, :]
                        off_target_summary_p_start = off_target_summary_blocked_taxon.loc[:, ['molecule_start', 'molecule_end']].drop_duplicates()
                        encoding_p_start = taxon_probes.loc[:,['p_start', 'ln']].drop_duplicates()
                        if not encoding_p_start.empty:
                            significant_overlap_fraction = np.zeros(encoding_p_start.shape[0])
                            for i in range(encoding_p_start.shape[0]):
                                max_start = off_target_summary_p_start.molecule_end.apply(max, args = (encoding_p_start.p_start.values[i] - encoding_p_start.ln.values[i], ))
                                min_end = off_target_summary_p_start.molecule_end.apply(min, args = (encoding_p_start.p_start.values[i], ))
                                overlap = min_end - max_start
                                significant_overlap_fraction[i] = np.sum(overlap > 0)/overlap.shape[0]
                            off_target_summary.loc[off_target_summary[target_rank].values == tt, 'max_encoding_interference_fraction'] = np.max(significant_overlap_fraction)
                        else:
                            off_target_summary.loc[off_target_summary[target_rank].values == tt, 'max_encoding_interference_fraction'] = 0
                    off_target_summary_full.loc[off_target_summary_full.probe_id.values == bp, 'max_encoding_interference_fraction'] = off_target_summary.max_encoding_interference_fraction.values
            blocking_probes_list.append(off_target_summary_full)
    blocking_probes = pd.concat(blocking_probes_list)
    blocking_probes.loc[:,'quadg'] = (blocking_probes['sseq'].str.upper().str.count('GGGG')) + (blocking_probes['sseq'].str.upper().str.count('GGGGG'))
    print(blocking_probes.shape)
    blocking_probes = blocking_probes.loc[blocking_probes.quadg.values == 0,:]
    print(blocking_probes.shape)
    blocking_probes_order_format = blocking_probes[['sseq', 'length']]
    blocking_probes_order_format = blocking_probes_order_format.drop_duplicates()
    blocking_probes_order_format = blocking_probes_order_format.sort_values(by = 'length', ascending = False)
    blocking_probes_order_format = blocking_probes_order_format.loc[blocking_probes_order_format.length.values >= bplc]
    blocking_probes_order_format.insert(0, 'blocking_probe_name', ['blocking_probe_{}'.format(i) for i in range(blocking_probes_order_format.shape[0])])
    blocking_probes_order_format = blocking_probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    blocking_probes_order_format.loc[blocking_probes_order_format.length.values < 15, 'Amount'] = '100nm'
    blocking_probes_order_format = blocking_probes_order_format[['blocking_probe_name', 'sseq', 'Amount', 'Purification', 'length']]
    blocking_probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    probe_length_filler = 'GTCTATTTTCTTATCCGACG'
    # 'GTCTATTTTCTTATCCGACGTGTTG'
    if primerset == 'A':
        forward_primer = 'CGCGGGCTATATGCGAACCG'
        reverse_primer = 'GCGTTGTATGCCCTCCACGC'
        # TAATACGACTCACTATAGGGCGTGGAGGGCATACAACGC
    elif primerset == 'B':
        forward_primer = 'CGATGCGCCAATTCCGGTTC'
        reverse_primer = 'CAACCCGCGAGCGATGATCA'
        # TAATACGACTCACTATAGGGTGATCATCGCTCGCGGGTTG
    elif primerset == 'C':
        forward_primer = 'GTTGGTCGGCACTTGGGTGC'
        reverse_primer = 'CCACCGGATGAACCGGCTTT'
        # TAATACGACTCACTATAGGGAAAGCCGGTTCATCCGGTGG
    else:
        forward_primer = ''
        reverse_primer = ''
    blocking_probe_seq_list = []
    for i in range(blocking_probes_order_format.shape[0]):
        bp_seq = blocking_probes_order_format.sseq.values[i]
        if len(bp_seq) == 15:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTGTTG'
        elif len(bp_seq) == 16:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTGTT'
        elif len(bp_seq) == 17:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTGT'
        elif len(bp_seq) == 18:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGTG'
        elif len(bp_seq) == 19:
            probe_length_filler = 'GTCTATTTTCTTATCCGACGT'
        elif len(bp_seq) >= 20:
            probe_length_filler = 'GTCTATTTTCTTATCCGACG'
        blocking_probe = [[forward_primer + probe_length_filler + bp_seq + probe_length_filler+ reverse_primer]]
        blocking_probe_seq_list = blocking_probe_seq_list + blocking_probe
    blocking_probe_seq_list = pd.DataFrame(blocking_probe_seq_list)
    blocking_probe_seq_list.columns = ['blocking_sequence']
    print(blocking_probe_seq_list.shape[0])
    blocking_probe_seq_list['blocking_sequence'].str.upper().to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probe_sequences.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = False, index = False)
    return

def generate_blocking_probes_plate_format(design_dir, plf = 'T', primer = 'T', primerset = 'B', barcode_selection = 'MostSimple'):
    design_dir_folder = os.path.split(design_dir)[1]
    blocking_probes_filenames = glob.glob('{}/*_off_target_summary_info.csv'.format(design_dir))
    if '{}/0_off_target_summary_info.csv'.format(design_dir) in blocking_probes_filenames:
        blocking_probes_filenames.remove('{}/0_off_target_summary_info.csv'.format(design_dir))
    blocking_probes_list = [pd.read_csv(filename) for filename in blocking_probes_filenames]
    blocking_probes = pd.concat(blocking_probes_list)
    blocking_probes_order_format = blocking_probes[['sseq', 'length']]
    blocking_probes_order_format.columns = ['Sequence', 'length']
    blocking_probes_order_format = blocking_probes_order_format.drop_duplicates()
    blocking_probes_order_format = blocking_probes_order_format.sort_values(by = 'length', ascending = False)
    blocking_probes_order_format.insert(0, 'Name', ['blocking_probe_{}'.format(i) for i in range(blocking_probes_order_format.shape[0])])
    blocking_probes_order_format = blocking_probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    blocking_probes_order_format.loc[blocking_probes_order_format.length.values < 15, 'Amount'] = '100nm'
    plate_row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    blocking_probes_order_format['Well'] = ''
    blocking_probes_order_format = blocking_probes_order_format.reset_index().drop(columns = ['index'])
    if blocking_probes_order_format.shape[0] > 384:
        for i in range(384):
            plate_row_number, plate_column_number = np.divmod(i, 24)
            if plate_row_number >= 16:
                plate_row_number = plate_row_number - 16
            blocking_probes_order_format.loc[i, 'Well'] = plate_row[plate_row_number] + str(plate_column_number+1)
        blocking_probes_order_format = blocking_probes_order_format[['Well', 'Name', 'Sequence', 'Amount', 'Purification', 'length']]
        blocking_probes_order_format.to_excel('{}/{}_blocking_probes_plate_1_order_format.xlsx'.format(design_dir, design_dir_folder), index = False)
        for i in range(384, blocking_probes_order_format.shape[0]):
            plate_row_number, plate_column_number = np.divmod(i-384, 12)
            if plate_row_number >= 8:
                plate_row_number = plate_row_number - 8
            blocking_probes_order_format.loc[i, 'Well'] = plate_row[plate_row_number] + str(plate_column_number+1)
        blocking_probes_order_format = blocking_probes_order_format[['Well', 'Name', 'Sequence', 'Amount', 'Purification', 'length']]
        blocking_probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes_plate_2_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection), index = False)
    else:
        for i in range(blocking_probes_order_format.shape[0]):
            plate_row_number, plate_column_number = np.divmod(i, 24)
            if plate_row_number >= 16:
                plate_row_number = plate_row_number - 16
            blocking_probes_order_format.loc[i, 'Well'] = plate_row[plate_row_number] + str(plate_column_number+1)
        blocking_probes_order_format = blocking_probes_order_format[['Well', 'Name', 'Sequence', 'Amount', 'Purification', 'length']]
        blocking_probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_blocking_probes_plate_1_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection), index = False)
    return

def write_final_probes_fasta(probes, design_dir, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    probes_fasta_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.fasta'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    probes_list = [SeqRecord(Seq(probes['probe_full_seq'][i]), id = str(probes['probe_id'][i]), description = '') for i in range (0, probes.shape[0])]
    SeqIO.write(probes_list, probes_fasta_filename, 'fasta')
    return

def write_final_unique_probes_fasta(probes, design_dir, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    probes = probes.drop_duplicates().reset_index().drop(['index'], axis = 1)
    probes_fasta_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_unique.fasta'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    probes_list = [SeqRecord(Seq(probes['probe_full_seq'][i]), id = str(probes['probe_id'][i]), description = '') for i in range (0, probes.shape[0])]
    SeqIO.write(probes_list, probes_fasta_filename, 'fasta')
    return

def blast_final_probes(design_dir, primerset, barcode_selection, blast_database):
    design_dir_folder = os.path.split(design_dir)[1]
    infile = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_unique.fasta'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    blast_output_filename = '{}/{}_full_length_probes_unique.blast.out'.format(design_dir, design_dir_folder)
    out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq'
    return_code = subprocess.call(['blastn', '-db', blast_database, '-query', infile, '-out', blast_output_filename, '-outfmt', out_format, '-task', 'blastn-short', '-max_hsps', '1', '-max_target_seqs', '100000', '-strand', 'minus', '-evalue', '100', '-num_threads', '1'])
    return(return_code)

def mch_test_v4(df):
    qseq_array = df.qseq.values
    sseq_array = df.sseq.values
    mch_array = np.zeros(len(qseq_array))
    for i in range(len(qseq_array)):
        qseq = qseq_array[i]
        sseq = sseq_array[i]
        if qseq != sseq:
            snp_indices = np.where(np.array(list(qseq)) != np.array(list(sseq)))[0]
            diffs = np.diff(snp_indices)
            mch_array[i] = np.max(np.append(diffs,[snp_indices[0], len(qseq) - 1 - snp_indices[-1]]))
        else:
            mch_array[i] = len(qseq)
    return(mch_array)

def summarize_final_probe_blast(df, mch, target_rank = None):
    df['mch'] = mch_test_v4(df)
    df_filtered = df.loc[df.mch >= mch]
    if df_filtered.shape[0] > 0:
        probe_id = np.unique(df_filtered['probe_id'])[0]
        probe_length = np.unique(df_filtered['probe_length'])[0]
        check = np.sum(((df_filtered['probe_start'] >= 44) & (df_filtered['probe_end'] <= 44 + probe_length - 1) & (df_filtered['mch'] <= probe_length)) | (df_filtered['target_taxon'] == df_filtered[target_rank]))/df_filtered.shape[0]
    else:
        print('{} has no blast hits...'.format(df.probe_id.values[0]))
        check = 0
    return(check)

def check_final_probes_blast(design_dir, probes, mch, bot, target_rank, blast_lineage_filename, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    probes_blast_filename = '{}/{}_full_length_probes_unique.blast.out'.format(design_dir, design_dir_folder)
    probes_blast = pd.read_table(probes_blast_filename, header = None)
    probes_blast.columns = ['probe_id', 'molecule_id', 'pid', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'probe_start', 'probe_end', 'molecule_start', 'molecule_end', 'evalue', 'bitscore', 'staxids', 'qseq', 'sseq']
    probes['probe_length'] = probes['rna_seq'].apply(len) - 6
    print(probes.shape)
    probes_length_df = probes[['probe_id', 'target_taxon', 'probe_length']]
    probes_blast = probes_blast.merge(probes_length_df, on = 'probe_id', how = 'left')
    blast_lineage = pd.read_table(blast_lineage_filename)
    blast_lineage_slim = blast_lineage.loc[:,['molecule_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']]
    probes_blast = probes_blast.merge(blast_lineage_slim, on = 'molecule_id', how = 'left')
    probes_blast_summary = probes_blast.groupby('probe_id', axis = 0).apply(summarize_final_probe_blast, mch = mch, target_rank = target_rank)
    probes_blast_summary = probes_blast_summary.reset_index()
    probes_blast_summary.columns = ['probe_id', 'final_check']
    probes_final_check = probes.merge(probes_blast_summary, on = 'probe_id', how = 'left')
    print(probes_final_check.shape)
    problematic_probes = probes_final_check.loc[(probes_final_check.final_check < bot)][['target_taxon', 'probe_numeric_id', 'probe_id', 'final_check']].drop_duplicates()
    print(problematic_probes.shape)
    probes_final_filter = probes_final_check.loc[~(probes_final_check.target_taxon.isin(problematic_probes.target_taxon.values) & probes_final_check.probe_numeric_id.isin(problematic_probes.probe_numeric_id.values))]
    probes_final_filter.to_csv('{}/{}_primerset_{}_barcode_selection_{}_probes_final_check_filter_N.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    print(probes_final_filter.shape)
    probes_final_filter.loc[:,'probe_full_seq_length'] = probes_final_filter.probe_full_seq.apply(len)
    probes_final_filter.to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_probes.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = True, index = False)
    probes_final_filter_summary = probes_final_filter['target_taxon'].value_counts()
    probes_final_filter_summary.to_csv('{}/{}_primerset_{}_barcode_selection_{}_probes_final_filter_summary.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    probes_final_filter['probe_full_seq'].str.upper().to_csv('{}/{}_primerset_{}_barcode_selection_{}_full_length_probe_sequences.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection), header = False, index = False)
    probes_order_format = probes_final_filter[['abundance', 'probe_id', 'probe_full_seq']]
    probes_order_format = probes_order_format.assign(Amount = '25nm', Purification = 'STD')
    probes_order_format.loc[probes_order_format['probe_full_seq'].str.len() > 60]['Amount'] = '100nm'
    probes_order_format.loc[probes_order_format['probe_full_seq'].str.len() > 60]['Purification'] = 'PAGE'
    probes_order_format = probes_order_format.sort_values(by = ['abundance', 'probe_id'], ascending = [False,True]).reset_index().drop(columns = ['index'])
    probes_order_format.to_excel('{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_order_format.xlsx'.format(design_dir, design_dir_folder, primerset, barcode_selection))
    return(probes_final_filter)

def generate_probe_statistics_plots(design_dir, probes_final_filter, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    fig = plt.figure()
    fig.set_size_inches(5,3)
    plt.hist(probes_final_filter.Tm.values, bins = 20)
    plt.xlabel(r'Melting Temperature [$^\circ$C]')
    plt.ylabel('Frequency')
    plt.tight_layout()
    probes_tm_histogram_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_tm_histogram.png'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(probes_tm_histogram_filename, dpi = 300)
    plt.close()
    fig = plt.figure()
    fig.set_size_inches(5,3)
    plt.hist(probes_final_filter.GC.values, bins = 20)
    plt.xlabel('GC Content [-]')
    plt.ylabel('Frequency')
    plt.tight_layout()
    probes_gc_histogram_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_gc_histogram.png'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(probes_gc_histogram_filename, dpi = 300)
    plt.close()
    fig = plt.figure()
    fig.set_size_inches(5,5)
    abundance_probe_multiplexity = probes_final_filter.groupby('target_taxon').agg({'abundance': np.unique, 'probe_numeric_id': lambda x: np.count_nonzero(np.unique(x))}).reset_index()
    abundance_probe_multiplexity = abundance_probe_multiplexity.rename(columns = {'probe_numeric_id': 'probe_multiplexity'})
    abundance_probe_multiplexity['relative_abundance'] = abundance_probe_multiplexity.abundance.values/abundance_probe_multiplexity.abundance.sum()
    plt.plot(abundance_probe_multiplexity.relative_abundance, abundance_probe_multiplexity.probe_multiplexity, 'o', color = 'orange', alpha = 0.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Relative Abundance')
    plt.ylabel('Probe Plurality')
    plt.tight_layout()
    abundance_plurality_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_abundance_plurality.png'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(abundance_plurality_filename, dpi = 300)
    plt.close()
    abundance_plurality_table_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_abundance_plurality.csv'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    abundance_probe_multiplexity.to_csv(abundance_plurality_table_filename)
    fig = plt.figure()
    fig.set_size_inches(5,5)
    plt.hist(probes_final_filter.probe_full_seq_length.values)
    plt.xlabel('Probe Length [bp]')
    plt.ylabel('Frequency')
    plt.tight_layout()
    probe_length_histogram_filename = '{}/{}_primerset_{}_barcode_selection_{}_full_length_probes_length_histogram.png'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    fig.savefig(probe_length_histogram_filename, dpi = 300)
    plt.close()
    return

def generate_probe_summary_file(design_dir, probes_final_filter, primerset, barcode_selection):
    design_dir_folder = os.path.split(design_dir)[1]
    probes_summary_filename = '{}/{}_full_length_probes_summary.txt'.format(design_dir, design_dir_folder, primerset, barcode_selection)
    with open(probes_summary_filename, 'w') as tf:
        tf.write('Probe design complete.')
    return
###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    # input blast filename
    parser.add_argument('design_dir', type = str, help = 'Directory to store all probe information')

    parser.add_argument('consensus_directory', type = str, help = 'Directory to the consensus folder created during probe design')

    parser.add_argument('blast_database', type = str, help = 'Blast database to check probes against')

    parser.add_argument('bot', type = float, help = 'Minimum blast on target rate ')

    parser.add_argument('mch', type = int, help = 'Maximum continout homology cut off')

    parser.add_argument('bplc', type = int, help = 'Blocking probe length cut off')

    parser.add_argument('-ps', '--primerset', dest = 'primerset', default = 'B', type = str, help = 'Primer sets to use')

    parser.add_argument('-plf', '--plf', dest = 'plf', default = 'T', type = str, help = 'Boolean to indicate whether to include probe length filler')

    parser.add_argument('-p', '--primer', dest = 'primer', default = 'T', type = str, help = 'Boolean to indicate whether to add primers to the full length sequences')

    parser.add_argument('-t', '--target_rank', dest = 'target_rank', default = '', type = str, help = 'Target taxonomic rank at which the probes were designed')

    parser.add_argument('-bs', '--barcode_selection', dest = 'barcode_selection', default = 'MostSimple', type = str, help = 'Method for taxa barcode assignment')

    args = parser.parse_args()

    print('Generating full probes for design {}...'.format(os.path.basename(args.design_dir)))
    taxon_best_probes = glob.glob('{}/*_probe_selection.csv'.format(args.design_dir))
    blast_lineage_filename = '{}/blast_lineage_strain.tab'.format(args.consensus_directory)
    for file in taxon_best_probes:
        taxon_best_probes_sa_filename = re.sub('_probe_selection.csv', '_probe_selection_sa.csv', file)
        add_spacer(file, args.consensus_directory, taxon_best_probes_sa_filename)
    oligo_df = generate_full_probes(args.design_dir, args.bot, plf = args.plf, primer = args.primer, primerset = args.primerset, barcode_selection = args.barcode_selection)
    write_final_probes_fasta(oligo_df, args.design_dir, args.primerset, args.barcode_selection)
    write_final_unique_probes_fasta(oligo_df, args.design_dir, args.primerset, args.barcode_selection)
    blast_final_probes(args.design_dir, args.primerset, args.barcode_selection, args.blast_database)
    probes_final_filter = check_final_probes_blast(args.design_dir, oligo_df, args.mch, args.bot, args.target_rank, blast_lineage_filename, args.primerset, args.barcode_selection)
    generate_blocking_probes(args.design_dir, args.bplc, args.target_rank, plf = args.plf, primer = args.primer, primerset = args.primerset, barcode_selection = args.barcode_selection)
    generate_probe_statistics_plots(args.design_dir, probes_final_filter, args.primerset, args.barcode_selection)
    generate_probe_summary_file(args.design_dir, probes_final_filter, args.primerset, args.barcode_selection)
    return

if __name__ == '__main__':
    main()
