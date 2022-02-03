#!/usr/bin/env python
# -*- coding: UTF-8

# imports
import os
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
import logomaker as lm

from Bio import SeqIO
from inputs.fastq_loader import result_directory


# start


def data_to_df_trans(data_matrix):  # covert matrix to a transposed pandas dataframe for heatmap building
    """Generating a transposed pandas dataframe"""
    data_df = pd.DataFrame(data_matrix)
    data_df_transposed = data_df.T
    data_df_transposed.columns = data_df_transposed.columns + 1
    return data_df_transposed


def data_to_df(data_matrix):  # covert matrix to a pandas dataframe
    """Generating a pandas dataframe from a matrix"""
    data_df = pd.DataFrame(data_matrix)
    data_df.columns = data_df.columns + 1
    return data_df


def convert_to_fasta(translated_list):
    """convert list of AA sequences back to a .fasta file for clustering"""
    save_path_fasta = result_directory
    result_fasta_in_file = 'results.fasta'
    complete_file = os.path.join(save_path_fasta, result_fasta_in_file)  # save path to results directory
    n = 0
    with open(complete_file, 'w') as out:
        for line in translated_list:
            n += 1
            out.write(">%d\n%s\n" % (n, line))
    return out


def get_seqs_back_to_list(FASTA):
    new_seqList = []
    with open(FASTA, 'r') as inFile:
        for record in SeqIO.parse(inFile, 'fasta'):
            new_seqList.append(str(record.seq))
    return new_seqList


# heatmaps
try:  # heatmaps from NNK pssm dataframe
    from isolate_translate import pssm_nnk_lib
except ImportError as i_error:
    print("\n", i_error)
else:
    save_path = result_directory
    save_file_plt1 = input("\nGive NNK heatmap a name and safe as .png, .jpg or .pdf, "
                           "for example: 'example_heatmap.pdf'.")

    complete_heatmap1 = os.path.join(save_path, save_file_plt1)

    pssm_nnk = data_to_df_trans(pssm_nnk_lib)
    matplotlib.use('TkAgg')
    # start
    fig1 = plt.figure(figsize=(12, 12))
    sb.set(font_scale=1.1)
    r1 = sb.heatmap(pssm_nnk, cmap='OrRd', annot=True, cbar_kws={'label': 'Color scale'})
    r1.set_yticklabels(r1.get_yticklabels(), rotation=0)
    r1.set_title("Heatmap of an amino acid position-specific scoring matrix"
                 "\n"
                 "generated with pseudo-count and NNK codon usage",
                 fontsize=18)

    plt.xlabel("Amino acid position in library", fontsize=18)
    plt.ylabel("Amino acid residue", fontsize=18)
    plt.savefig(complete_heatmap1, transparent=True)

try:  # heatmaps from NNN pssm dataframe
    from isolate_translate import pssm_nnn_lib
except ImportError as i_error:
    print("\n", i_error)
else:
    save_path = result_directory
    save_file_plt2 = input("\nGive NNN heatmap a name and safe as .png, .jpg or .pdf, "
                           "for example: 'example_heatmap.pdf'.")

    complete_heatmap2 = os.path.join(save_path, save_file_plt2)

    pssm_nnn = data_to_df_trans(pssm_nnn_lib)
    matplotlib.use('TkAgg')
    # start
    fig1 = plt.figure(figsize=(12, 12))
    sb.set(font_scale=1.1)
    r1 = sb.heatmap(pssm_nnn, cmap='coolwarm', annot=True, cbar_kws={'label': 'Color scale'})
    r1.set_yticklabels(r1.get_yticklabels(), rotation=0)
    r1.set_title("Heatmap of an amino acid position-specific scoring matrix"
                 "\n"
                 " generated with pseudo-count and NNN codon usage",
                 fontsize=18)
    plt.xlabel("Amino acid position in library", fontsize=18)
    plt.ylabel("Amino acid residue", fontsize=18)
    plt.savefig(complete_heatmap2, transparent=True)

try:
    from isolate_translate import translate_nnk_lib
except ImportError as i_error:
    print("\n", i_error)
else:
    try:
        save_path = result_directory
        save_file_logo1 = input("\nGive the NNK sequence logo a name and save as a .pdf. E.g example_logo.pdf")
        complete_path_logo = os.path.join(save_path, save_file_logo1)
    except FileExistsError as exists:
        print("\n", exists)
    else:
        # convert list of translated sequences back to fasta file for clustering
        new_fasta = convert_to_fasta(translate_nnk_lib)

        os.system('/home/christoph/anaconda3/envs/statistical_analysis_test/bin/cd-hit '
                  '-i results/results.fasta -o clustered_nnk -c 0.9 -n 5 -T 7 -M 16000 -l 6')

        try:
            with open('clustered_nnk', 'r') as f:
                raw_seqs = f.readlines()
                seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and '>' not in seq]
            print('\nThere are %d sequences, all of length %d' % (len(seqs), len(seqs[0])))
        except FileNotFoundError as not_found:
            print("\n", not_found)
        else:
            background_dic_nnk = {"A": 0.1, "C": 0.05, "D": 0.05, "E": 0.05, "F": 0.05, "G": 0.15, "H": 0.05, "I": 0.05,
                                  "K": 0.05, "L": 0.15, "M": 0.05, "N": 0.05, "P": 0.1, "Q": 0.05, "R": 0.15, "S": 0.15,
                                  "T": 0.1, "V": 0.1, "W": 0.05, "Y": 0.05,
                                  }
            background_tmp_nnk = list(background_dic_nnk.values())
            background_array_nnk = np.array(background_tmp_nnk)

            counts_mat = lm.alignment_to_matrix(seqs)
            print(counts_mat)
            prob_mat = lm.transform_matrix(counts_mat, normalize_values=True, pseudocount=0)
            print(prob_mat)
            info_mat = lm.transform_matrix(counts_mat,
                                           from_type='counts',
                                           to_type='information',
                                           background=background_array_nnk,
                                           pseudocount=0)
            info_mat.index = np.arange(1, len(info_mat) + 1)
            print(info_mat)
            info_logo = lm.Logo(info_mat, color_scheme='chemistry',
                                baseline_width=1.0,
                                figsize=(12, 5),
                                vpad=.1,
                                width=1)
            # set axes label
            info_logo.ax.set_xlabel('Positions', fontsize=16)
            info_logo.style_xticks(anchor=0, spacing=1)
            info_logo.ax.set_xlim([0, len(info_mat) + 1])
            info_logo.ax.set_ylabel('Information (bits)', labelpad=-1, fontsize=16)
            plt.title('Sequence motif information logo with NNK codon usage')
            plt.savefig(complete_path_logo)

            try:
                pd.set_option('display.max_columns', None)
                pd.set_option('display.width', None)
                pd.set_option('display.max_colwidth', None)

                save_path = result_directory
                save_file_mat = "nnk_matrix_report.txt"
                complete_path_mat = os.path.join(save_path, save_file_mat)

                with open(complete_path_mat, 'w') as f:
                    print('\nThis is the count matrix of the clustered sequences:\n', counts_mat, file=f)
                    print('\nThis is the probability matrix of the clustered sequences:\n', prob_mat, file=f)
                    print('\nThis is the information content matrix of the clustered sequences:\n', info_mat, file=f)
            except FileExistsError as exists:
                print("\n", exists)

try:
    from isolate_translate import translate_nnn_lib
except ImportError as i_error:
    print("\n", i_error)
else:
    try:
        save_path = result_directory
        save_file_logo1 = input("\nGive the NNN sequence logo a name and save as a .pdf. E.g example_logo.pdf")
        complete_path_logo = os.path.join(save_path, save_file_logo1)
    except FileExistsError as exists:
        print("\n", exists)
    else:
        # convert list of translated sequences back to fasta file for clustering
        new_fasta = convert_to_fasta(translate_nnn_lib)

        os.system('/home/christoph/anaconda3/envs/statistical_analysis_test/bin/cd-hit '
                  '-i results/results.fasta -o clustered_nnn -c 0.9 -n 5 -T 7 -M 16000 -l 6')

        try:
            with open('clustered_nnn', 'r') as f:
                raw_seqs = f.readlines()
                seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and '>' not in seq]
            print('\nThere are %d sequences, all of length %d' % (len(seqs), len(seqs[0])))
        except FileNotFoundError as not_found:
            print("\n", not_found)
        else:
            background_dic_nnn = {"A": 0.1, "C": 0.05, "D": 0.05, "E": 0.05, "F": 0.05, "G": 0.15, "H": 0.05, "I": 0.05,
                                  "K": 0.05, "L": 0.15, "M": 0.05, "N": 0.05, "P": 0.1, "Q": 0.05, "R": 0.15, "S": 0.15,
                                  "T": 0.1, "V": 0.1, "W": 0.05, "Y": 0.05
                                  }
            background_tmp_nnn = list(background_dic_nnn.values())
            background_array_nnn = np.array(background_tmp_nnn)

            counts_mat = lm.alignment_to_matrix(seqs)
            print(counts_mat)
            prob_mat = lm.transform_matrix(counts_mat, normalize_values=True, pseudocount=0)
            print(prob_mat)
            info_mat = lm.transform_matrix(counts_mat,
                                           from_type='counts',
                                           to_type='information',
                                           background=background_array_nnn,
                                           pseudocount=0)
            info_mat.index = np.arange(1, len(info_mat) + 1)
            print(info_mat)
            info_logo = lm.Logo(info_mat, color_scheme='chemistry',
                                baseline_width=1.0,
                                figsize=(12, 5),
                                vpad=.1,
                                width=1)
            # set axes label
            info_logo.ax.set_xlabel('Positions', fontsize=16)
            info_logo.style_xticks(anchor=0, spacing=1)
            info_logo.ax.set_xlim([0, len(info_mat) + 1])
            info_logo.ax.set_ylabel('Information (bits)', labelpad=-1, fontsize=16)
            plt.title('Sequence motif information logo with NNN codon usage')
            plt.savefig(complete_path_logo)

            try:
                pd.set_option('display.max_columns', None)
                pd.set_option('display.width', None)
                pd.set_option('display.max_colwidth', None)

                save_path = result_directory
                save_file_mat = "nnn_matrix_report.txt"
                complete_path_mat = os.path.join(save_path, save_file_mat)

                with open(complete_path_mat, 'w') as f:
                    print('\nThis is the count matrix of the clustered sequences:\n', counts_mat, file=f)
                    print('\nThis is the probability matrix of the clustered sequences:\n', prob_mat, file=f)
                    print('\nThis is the information content matrix of the clustered sequences:\n', info_mat, file=f)
            except FileExistsError as exists:
                print("\n", exists)
