#!/usr/bin/env python
# -*- coding: UTF-8

# imports ##############################################################################################################
import os.path
import re
import regex

import pandas as pd

from collections import Counter
from Bio import motifs

from inputs.fastq_loader import seq_List
from inputs.fastq_loader import result_directory

# isolate random peptide insert ########################################################################################

# add option: if insert.txt already exists --> skip
try:
    text_name = input("\nGive the .txt (text) file where the inserts are going to be saved a name:")
    file_name = '{}.txt'.format(text_name)
except FileExistsError as exists:
    print("\n", exists)
else:
    complete_path = os.path.join(result_directory, file_name)

    while True:
        ask_operation = input('\nDo you want to start the pipeline and '
                              'isolate the random sequence portion of experiment, yes or no?')
        if ask_operation == 'yes':
            ask_task = input("\nDo you use a Ph.D.7 library, a Ph.D.12 library or a custom library, "
                             "enter 'A for Ph.D.7', 'B for Ph.D.12' or ' C for custom':")

            seqList_string = str(seq_List)
            rev_seqList_string = seqList_string[::-1]

            if ask_task == "A":

                phd_up_flank = 'TCT'
                print('\nUpstream flanking region is:', phd_up_flank)
                phd_down_flank = 'GGTGGAGGT'
                print('\nDownstream flanking region is:', phd_down_flank)
                standard_outs = open(complete_path, 'w')
                capture = regex.compile("TCT(.{21})GGT")
                for match in capture.finditer(seqList_string):
                    if match:
                        insert_phd = match.group(1)
                        standard_outs.write(insert_phd + "\n")
                    else:
                        print('\nNo forward matches found:')

                print('\nForward inserts found:', len([*capture.finditer(seqList_string)]))

                for match in capture.finditer(rev_seqList_string):
                    if match:
                        rev_insert_phd = match.group(1)
                        standard_outs.write(rev_insert_phd + "\n")
                    else:
                        print('\nNo reverse matches found')
                standard_outs.close()
                print('\nReverse inserts found:', len([*capture.finditer(rev_seqList_string)]))

                print('\nTogether:', len([*capture.finditer(seqList_string)])
                      + len([*capture.finditer(rev_seqList_string)]))

                break

            if ask_task == "B":

                phd_up_flank = 'TCT'
                print('\nUpstream flanking region is:', phd_up_flank)
                phd_down_flank = 'GGTGGAGGT'
                print('\nDownstream flanking region is:', phd_down_flank)
                standard_outs = open(complete_path, 'w')
                capture = regex.compile("TCT(.{36})GGT")
                for match in capture.finditer(seqList_string):
                    if match:
                        insert_phd = match.group(1)
                        standard_outs.write(insert_phd + "\n")
                    else:
                        print('\nNo forward matches found:')

                print('\nForward inserts found:', len([*capture.finditer(seqList_string)]))

                for match in capture.finditer(rev_seqList_string):
                    if match:
                        rev_insert_phd = match.group(1)
                        standard_outs.write(rev_insert_phd + "\n")
                    else:
                        print('\nNo reverse matches found')
                standard_outs.close()
                print('\nReverse inserts found:', len([*capture.finditer(rev_seqList_string)]))

                print('\nTogether',
                      len([*capture.finditer(seqList_string)]) + len([*capture.finditer(rev_seqList_string)]))

                break

            if ask_task == "C":

                print('\nYou will now be asked to enter your custom up- and downstream markers')
                seqList_string = str(seq_List)
                upstream_marker = input('\nPlease enter the upstream marker sequence that you want to use '
                                        'in this run. Please only use uppercase letters:')
                downstream_marker = input('\nPlease enter the downstream marker sequence that you want to use '
                                          'in this run. Please only use uppercase letters:')
                insert_length = input('\nPlease enter the custom length of the insert that you are searching for, '
                                      'not more than three digits allowed:')

                print("\nUser selected:'{}'as upstream marker."
                      " Input type: {}.".format(upstream_marker, type(upstream_marker)))
                print("\nUser selected:'{}' as downstream marker."
                      " Input type: {}.".format(downstream_marker, type(downstream_marker)))
                print("\nUser selected: '{}' as random insert length".format(insert_length))


                def validation_check_marker(input_string):
                    regex_check = re.compile("[ATGC]+")
                    hit = regex_check.match(str(input_string))
                    return bool(hit)


                def validation_check_length(input_string):
                    regex_check = re.compile("[0-9]{1,3}")
                    hit = regex_check.match(str(input_string))
                    return bool(hit)


                # add function for calculation of Hamming/Levenshtein distance

                validate_1 = validation_check_marker(upstream_marker)
                if validate_1:
                    print('\nThe input {} as a upstream marker is valid and will be used!'.format(upstream_marker))
                else:
                    print('\nThe input {} is not valid'.format(upstream_marker))

                validate_1 = validation_check_marker(downstream_marker)
                if validate_1:
                    print('\nThe input {} as a downstream marker is valid and will be used!'.format(downstream_marker))
                else:
                    print('\nThe input {} is not valid'.format(downstream_marker))

                validate_2 = validation_check_length(insert_length)
                if validate_2:
                    print('\nThe input {} will be used as the length of the '
                          'random insert portion of the sequence!'.format(insert_length))

                capture_rex = "."
                generic_re = regex.compile("%s{s<=1}(%s{%s})%s{s<=2}" %
                                           (upstream_marker, capture_rex, insert_length, downstream_marker),
                                           flags=regex.V1)
                custom_outs = open(complete_path, 'w')
                for match in generic_re.finditer(seqList_string):
                    if match:
                        insert_custom = match.group(1)
                        custom_outs.write(insert_custom + "\n")
                    else:
                        print('\nNo forward matches found!')

                print('\nForward inserts found:', len([*generic_re.finditer(seqList_string)]))

                for match in generic_re.finditer(rev_seqList_string):
                    if match:
                        rev_insert_custom = match.group(1)
                        custom_outs.write(rev_insert_custom + "\n")
                    else:
                        print("\nNo reverse matches found!")
                custom_outs.close()

                print("\nReverse inserts found:", len([*generic_re.finditer(rev_seqList_string)]))

                print("\nTogether:", len([*generic_re.finditer(seqList_string)]) +
                      len([*generic_re.finditer(rev_seqList_string)]))
        else:
            print("\nOkay, I will ask no questions and stop the process!")

            break
        break

    # convert DNA to amino acid sequence ##############################################################################

    input_file = complete_path
    with open(input_file, 'r') as f:
        seq = f.readlines()
    seq = [s.replace(" ", "").replace(",", "").replace("'", "").replace("\n", "") for s in seq]


    def translate_nnn(seq):
        """Return NNN codon usage translation for every DNA sequence in list of DNA sequences """
        nnn_table = {
            'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F',
            'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S',
            'TAA': '*', 'TGA': '*', 'TTG': 'L', 'TCG': 'S', 'TAG': '*',
            'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
            'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L',
            'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CTG': 'L', 'CCG': 'P',
            'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N',
            'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
            'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'ATG': 'M',
            'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A',
            'GAT': 'D', 'GGT': 'G', 'GTC': 'V', 'GCC': 'A', 'GAC': 'D',
            'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
            'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
        }
        nnn_aa_seq = []
        nnn_aa_seq_not_length = []
        print("\nStarting to translate:")
        for dna in seq:
            protein_seq = ""
            if len(dna) % 3 != 0:
                nnn_aa_seq_not_length.append(protein_seq)
                continue  # go to next element of seq
            for i in range(0, len(dna), 3):
                nnn_codon = nnn_table[dna[i:i + 3]]
                protein_seq += nnn_codon
            nnn_aa_seq.append(protein_seq)
        return nnn_aa_seq


    def translate_nnk(seq):
        """Return NNK codon usage translation for every DNA sequence in list of DNA sequences """
        nnk_table = {
            'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTG': 'L', 'TCG': 'S', 'TAG': 'G',
            'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 'CTG': 'L', 'CCG': 'P',
            'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 'ATG': 'M',
            'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
            'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
        }
        nnk_aa_seq = []
        nnk_aa_seq_not_length = []
        print("\nStarting to translate:")
        for dna in seq:
            protein_seq = ""
            if len(dna) % 3 != 0:
                nnk_aa_seq_not_length.append(protein_seq)
                continue
            for i in range(0, len(dna), 3):
                nnk_codon = nnk_table.get(dna[i:i + 3], 'X')
                protein_seq += nnk_codon
            if protein_seq.count('X') == 0:
                nnk_aa_seq.append(protein_seq)
        return nnk_aa_seq


    def pos_specific_score_matrix_nnn(translated_list):
        """Build position specific scoring matrix for AA residues with pseudo-count generated with NNN codon usage"""
        m = motifs.create(translated_list, alphabet="ACDEFGHILKMNPQRSTVWY*")
        background = {"F": 0.1, "S": 0.3, "Y": 0.1, "C": 0.1, "L": 0.3, "G": 0.2,
                      "W": 0.05, "P": 0.2, "H": 0.1, "R": 0.3, "Q": 0.1, "I": 0.15,
                      "T": 0.2, "N": 0.1, "M": 0.05, "K": 0.1, "V": 0.2, "A": 0.2,
                      "D": 0.1, "E": 0.1, "*": 0.15
                      }
        pwm_nnn = m.counts.normalize(pseudocounts=background)
        print(pwm_nnn)
        pssm_nnn = pwm_nnn.log_odds()
        print(pssm_nnn)
        return pssm_nnn


    def pos_specific_score_matrix_nnk(translated_list):
        """Build position specific scoring matrix for AA residues with pseudo-count generated with NNK codon usage"""
        m = motifs.create(translated_list, alphabet="ACDEFGHILKMNPQRSTVWY")
        background = {"F": 0.05, "S": 0.15, "Y": 0.05, "C": 0.05, "L": 0.15, "G": 0.15,
                      "W": 0.05, "P": 0.1, "H": 0.05, "R": 0.15, "Q": 0.05, "I": 0.05,
                      "T": 0.1, "N": 0.05, "M": 0.05, "K": 0.05, "V": 0.1, "A": 0.1,
                      "D": 0.05, "E": 0.05
                      }
        pwm_nnk = m.counts.normalize(pseudocounts=background)
        print(pwm_nnk)
        pssm_nnk = pwm_nnk.log_odds()
        print(pssm_nnk)
        return pssm_nnk


    # create results

    while True:
        ask_translation = input("\nYou will now be asked if you want to "
                                "\ntranslate with NNN codon usage, "
                                "NNK codon usage or both, type 'NNN', 'NNK' or 'both'")
        if ask_translation == 'NNN':
            translate_nnn_lib = translate_nnn(seq)
            count_nnn = Counter(translate_nnn_lib)
            ordered_count_nnn = dict(count_nnn.most_common())
            print("\nCount of amino acid sequences with NNN codon usage:", ordered_count_nnn)
            ask_excel_nnn = input("\nDo you want to save your results in an excel table, yes or no?")
            if ask_excel_nnn == "yes":
                name_excel = input("\nGive the .xlsx (excel) file a name")
                ordered_count_nnn_file = '{}.xlsx'.format(name_excel)
                complete_path_count_nnn = os.path.join(result_directory, ordered_count_nnn_file)
                with pd.ExcelWriter(complete_path_count_nnn) as writer:
                    df = pd.DataFrame.from_dict(ordered_count_nnn, orient='index', columns=['count'])
                    df.to_excel(writer, sheet_name='sheet_1')
                print("\nSaving results to {} file".format('ordered_count_nnn.xlsx'))
            else:
                print("\nResults won't be saved in an excel table")

            print("\nCreating the positional weight matrix for the options you have chosen")

            pssm_nnn_lib = pos_specific_score_matrix_nnn(translate_nnn_lib)

            break

        if ask_translation == 'NNK':
            translate_nnk_lib = translate_nnk(seq)
            count_nnk = Counter(translate_nnk_lib)
            ordered_count_nnk = dict(count_nnk.most_common())
            print("\nCount of amino acid sequences translated with NNK codon usage:", ordered_count_nnk)
            ask_excel_nnk = input("\nDo you want to save your results in an excel table, yes or no?")
            if ask_excel_nnk == "yes":
                name_excel = input("\nGive the .xlsx (excel) file a name")
                ordered_count_nnk_file = '{}.xlsx'.format(name_excel)
                complete_path_count_nnk = os.path.join(result_directory, ordered_count_nnk_file)
                with pd.ExcelWriter(complete_path_count_nnk) as writer:
                    df = pd.DataFrame.from_dict(ordered_count_nnk, orient='index', columns=['count'])
                    df.to_excel(writer, sheet_name='sheet_1')

                print("\nSaving results to {} file".format('ordered_count_nnk.xlsx'))
            else:
                print("\nResults won't be saved in an excel table")

            print("\nCreating the positional weight matrix for the options you have chosen:")

            pssm_nnk_lib = pos_specific_score_matrix_nnk(translate_nnk_lib)

            break

        if ask_translation == 'both':
            translate_nnn_lib = translate_nnn(seq)
            count_nnn = Counter(translate_nnn_lib)
            ordered_count_nnn = dict(count_nnn.most_common())
            print("\nCount of amino acid sequences translated with NNN codon usage:", ordered_count_nnn)

            translate_nnk_lib = translate_nnk(seq)
            count_nnk = Counter(translate_nnk_lib)
            ordered_count_nnk = dict(count_nnk.most_common())
            print("\nCount of amino acid sequences translated with NNK codon usage:", ordered_count_nnk)

            ask_excel_both = input("\nDo you want to save your results in an excel table, yes or no?")
            if ask_excel_both == "yes":
                name_excel = input("\nGive the .xlsx (excel) file a name")
                ordered_count_translation = '{}.xlsx'.format(name_excel)
                complete_path_count_translation = os.path.join(result_directory, ordered_count_translation)
                with pd.ExcelWriter(complete_path_count_translation) as writer:
                    df1 = pd.DataFrame.from_dict(ordered_count_nnn, orient='index', columns=['count'])
                    df2 = pd.DataFrame.from_dict(ordered_count_nnk, orient='index', columns=['count'])
                    df1.to_excel(writer, sheet_name='sheet_nnn')
                    df2.to_excel(writer, sheet_name='sheet_nnk')

                    print("\nSaving results to {} file".format(ordered_count_translation))
            else:
                print("\nResults won't be saved in an excel table")

            print("\nCreating the positional weight matrix for the options you have chosen:")

            pssm_nnk_lib = pos_specific_score_matrix_nnk(translate_nnk_lib)
            pssm_nnn_lib = pos_specific_score_matrix_nnn(translate_nnn_lib)

        else:
            print("\nInvalid input! Please write which codon usage you want to use.")

        break
