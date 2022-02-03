import os

found = False
dest_folder = None
fastq_file = None


# set destination folder
found = False
while not found:
    dest_folder = input("\nEnter full path of destination folder: ")
    if not os.path.isdir(dest_folder):
        print(dest_folder, '\nThis folder has not been found. Enter correct path. ')
    else:
        print("\nFolder path was found!", dest_folder)
        found = True

# set input file or files
found = False
while not found:
    fast_file = input("\nEnter full path your fasta, fastq ,fastq.gz files (best paired-end joined): ")
    if not os.path.isfile(fast_file):
        print(fastq_file, '\nThis file has not been found. Enter correct path. ')
    else:
        print("\nFile path was found!", fast_file)
        found = True



