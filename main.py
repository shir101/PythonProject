# Import parts of Biopython
from Bio import SeqIO
# import anytree
# import collections
import os
import numpy as np


# the function takes a fasta file and extracts the genome as a string
def extractFasta(path_to_file):
    sequence = ''
    with open(path_to_file, mode='r') as handle:
        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        for record in SeqIO.parse(handle, 'fasta'):
            # Extract individual parts of the FASTA record
            identifier = record.id
            description = record.description
            sequence = record.seq

            # Example: adapt to extract features you are interested in
            print('----------------------------------------------------------')
            print('Processing the record {}:'.format(identifier))
            print('Its description is: \n{}'.format(description))
            amount_of_nucleotides = len(sequence)
            print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides))

    sequence = str(sequence)
    print(type(sequence))
    # print(sequence)
    return sequence[0:49]


# the function cuts the sequence to kmers of length k and step
# returns a list with all the sequences
def kmers(k, step, sequence):
    # string_name[start:end:step]
    list1 = []
    for index in range(0, len(sequence), step):
        # print(sequence[index:index+k:step],'\n')
        list1.append(sequence[index:index + k:1])
    return list1


# !!!extend function to work on comparing the new genome to deeper lineages
# (not ust the basic genome, but the one that is the closest to him on the phylogenetic tree)
# the function compares the new genome to the base covid19 and saves the differences
def find_differntial(kmers_list, basic, families):
    for lineage in range(0, len(families)):  # iterating over he lineages (pages) - rows of the "families" array
        for position in range(0, len(families[lineage])):  # iterating over the kmers in the lineage - cols of the "families" array
            if (families[lineage][position] != kmers_list[position]):
                # FIX INSERTING DIFFERENT KMER INTO FAMILIES !!!!!!!!!! 
                np.insert(families,position,families[lineage][position])  # !!!!check if inserts into the next page end of kmers_list


if __name__ == '__main__':
    # input_file =
    # output_file =
    # fasta_sequenses = SeqIO.parse(open(input_file),'fasta')
    # with open(output_file) as out_file:
    #     for fasta in fasta_sequences:
    #         name, sequence = fasta.id, str(fasta.seq)
    #         new_sequence = some_function(sequence)
    #         write_fasta(out_file)

    # Install the biopython library (if not already installed)

    families = []
    path_to_file = 'C:/Users/User/PycharmProjects/DatabaseBuild/fasta/EPI_ISL_402124.fasta'  # <--- substitute by your local path
    sequence = extractFasta(path_to_file)
    kmers_list = kmers(4, 4, sequence)
    print(kmers_list)
    # print(type(kmers_list[7]))
    families.append(kmers_list)
    # families=np.array(families)

    directory = 'fasta'
    for filename in os.scandir(directory):
        if filename.is_file() and filename.path != 'fasta\EPI_ISL_402124.fasta':
            print(filename.path)
            # File path to your FASTA file
            path_to_file = 'C:/Users/User/PycharmProjects/DatabaseBuild/' + filename.path  # <--- substitute by your local path
            sequence = extractFasta(path_to_file)
            kmers_list = kmers(4, 4, sequence)
            # saving the differences between the new genome and the basic genome - !!compare new genome to the closest one in the lineage that is already in the database
            find_differntial(kmers_list, families[0], families)
            print(kmers_list)
            # print(type(kmers_list[7]))
            #families.append(kmers_list)

    print(families[0])  # the first genome is the reference basic covid19
    print(len(families))  # num of rows
    print(len(families[0]))  # num of columns

    print(families)
