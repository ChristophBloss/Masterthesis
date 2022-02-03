#!/usr/bin/env python
# -*- coding: UTF-8

# imports
import os.path
import sys
import shutil

from Bio import SeqIO
from tqdm import tqdm

from inputs.user_input import fast_file

current_directory = os.getcwd()
result_directory = os.path.join(current_directory, r'results')
if not os.path.exists(result_directory):
    os.makedirs(result_directory)

file_name = 'copy_file.fastq'
complete_path_file = os.path.join(result_directory, file_name)


# set user working directory. aka directory of the script
os.chdir(os.path.dirname(sys.argv[0]))
print("\nCurrent working directory: {0}".format(os.getcwd()))

# load joined fastq file into cwd
scr = fast_file
dest = os.getcwd()
try:
    shutil.copy(scr, complete_path_file)
except shutil.SameFileError:
    pass
# check if file is in cwd
with os.scandir((os.getcwd())) as entries:
    for entry in entries:
        if entry.is_file():
            print('\nFiles in the current working directory', entry.name)


# parse fastq files and create list of sequences
seq_List = []
with open(complete_path_file, 'r') as inFile:
    sequences = SeqIO.parse(inFile, 'fastq')
    for record in tqdm(sequences):
        if record.seq.count('N') == 0:
            seq_List.append(str(record.seq))
    print('\nNumber of sequences found in the input file', len(seq_List))
