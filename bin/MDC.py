#!/usr/bin/env python3
"""
Version: 1.0
Author: María Halldórsdóttir

List of functions:
This program has 3 functions:
    populate_lists_dictionary > Populates 2 lists and 1 dictionary used for the clustering and 1 dictionary for dereplication.
    percentage_similarity > calculates % similarity between 2 sequences.
    remove_empty_from_dict > removes empty keys from dictionary.
"""
import argparse
import itertools

usage = "This program dereplicates sequences, or clusters together similar sequences"
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-i',metavar='InFile',help='A fasta input file',type=argparse.FileType('r'),required=1)
parser.add_argument('-o',metavar='OutFile_1',help='A fasta output file with dereplicated sequences',type=argparse.FileType('w'),required=1)
parser.add_argument('-u',metavar='OutFile_2',help='A log file with ID´s of replicate sequences',type=argparse.FileType('w'),required=1)
parser.add_argument('-p',metavar='Percentage',help='optional; choose % similarity of sequences if clustering',type=float,default=1.0)
args=parser.parse_args()


def populate_lists_dictionary(input):
    IDS = []
    Sequences = []
    dict_replication = {}
    dict_file = {}
    with input as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                line = line.split()
                ID = line[0]
                IDS.append(ID)
            else:
                Sequences.append(line)
                if line not in dict_file:
                    dict_file[ID] = line
                if line in dict_replication:
                    dict_replication[line].append(ID)
                else:
                    dict_replication[line] = [ID]
        return(IDS, Sequences, dict_replication, dict_file)

def percentage_similarity(seq1, seq2):
    if len(seq1) == len(seq2):
        count = 0.00
        for bp in range(len(seq1)):
            if seq1[bp] == seq2[bp]:
                count += 1
        return(count/len(seq1))
    if len(seq1) > len(seq2):
        count = 0.00
        for bp in range(len(seq2)):
            if seq1[bp] == seq2[bp]:
                count +=1
        return(count/len(seq1))
    if len(seq1) < len(seq2):
        count = 0.00
        for bp in range(len(seq1)):
            if seq2[bp] == seq1[bp]:
                count +=1
        return(count/len(seq2))

def remove_empty_from_dict(d):
    if type(d) is dict:
        return dict((k, remove_empty_from_dict(v)) for k, v in d.items() if v and remove_empty_from_dict(v))
    elif type(d) is list:
        return [remove_empty_from_dict(v) for v in d if v and remove_empty_from_dict(v)]
    else:
        return d

IDS, Sequences, dict_replicates, dict_file = populate_lists_dictionary(args.i)  # open input file and populate all the lists and dictionaries.

if args.p == 1: # For dereplication of the fasta file

    Replicates = dict_replicates    # get the dictionary that has sequences as key and list of IDs as values.
    with args.o as fout, args.u as fout1:   # print the dereplicated sequences.
        print("This is a log file, which has the information of which ID´s have identical/similar sequences :)", file=fout1)
        for k in sorted (Replicates, key=lambda k: len(Replicates[k]), reverse=True):   # sort the dictionary with so that keys with more replicates print first.
            ID = Replicates[k][0]
            ID_string = Replicates[k]
            ID_string = ",".join(ID_string)
            replicates = len(Replicates[k])
            print("{} ;Size={}\n{}".format(ID, replicates, k), file=fout)
            print("*Sequence:{}\nIDs:{}".format(k, ID_string), file=fout1)

else: # For clustering of fasta file

    dict_clusters = {} # dictionary for clustering
    list_with_IDS = [] # list with IDS
    IDS_in_dict = [] # IDS that are already in dict_clusters
    seq_in_dict = [] # sequences that are in dict_clusters
    for i in range(len(Sequences)): # access the indexes in the list with all the sequences
        for k in range(len(Sequences)): # then compare every sequence against each other
            if i != k:  # dont compare the same sequence aginast it-self.
                seq1 = Sequences[i] # sequence 1 to compare
                seq2 = Sequences[k] # agianst sequence 2
                ID_1 = IDS[i]   # ID for sequence 1
                ID_2 = IDS[k]   # ID for sequence 2
                score = percentage_similarity(seq1, seq2)   # compare sequence1 and sequence2 for similarity score.
                if score >= args.p: # if score bigger or equal to the input percentage
                    if seq1 not in dict_clusters:   # if sequence1 is not in dict
                        dict_clusters[seq1] = []    # add to dictionary key
                        list_with_IDS.append(ID_1)  # add the ID to a list
                        list_with_IDS.append(ID_2)
                        seq_in_dict.append(seq1)    # add the sequence to a list
                        if ID_1 not in IDS_in_dict: # if ID not in a list with IDS that have been added to the clusters dictionary -
                            dict_clusters[seq1].append(ID_1)    # then add the ID to the dictionary value. This is done to prevent duplication of the IDs.
                            IDS_in_dict.append(ID_1)    # add ID to a list
                        if ID_2 not in IDS_in_dict:
                            dict_clusters[seq1].append(ID_2)
                            IDS_in_dict.append(ID_2)
                    else:
                        list_with_IDS.append(ID_1)  # if sequence1 in dict // same procedure as if seq1 not in dict
                        if ID_1 not in IDS_in_dict:
                            dict_clusters[seq1].append(ID_1)
                            IDS_in_dict.append(ID_1)
                        if ID_2 not in IDS_in_dict:
                            dict_clusters[seq1].append(ID_2)
                            IDS_in_dict.append(ID_2)

    dictionary_of_file= dict_file

    for key in dictionary_of_file:  # populate the main clusters dictionary with those sequences that did not have
        sequence = dictionary_of_file[key] # any match in the percentage_similarity function.
        ID = key
        if ID not in list_with_IDS:
            dict_clusters[sequence] = [ID]

    for k, v in dict_clusters.items(): # remove duplicate values in dictionary for each key.
        v.sort()
        dict_clusters[k] = [item for item, _ in itertools.groupby(v)]


    dict_clusters1 = {tuple(v): k for k, v in dict_clusters.items()} # remove keys from dictionary with the same value.
    dict_clusters = {v: list(k) for k, v in dict_clusters1.items()}

    dict_clusters = remove_empty_from_dict(dict_clusters) # remove empty keys from dictionary.


    with args.o as fout, args.u as fout1:
        print("This is a log file, which has the information of which ID´s have identical/similar sequences :)", file=fout1)
        for k in sorted (dict_clusters, key=lambda k: len(dict_clusters[k]), reverse=True): # sort the dictionary with so that keys with more replicates print first.
            ID = dict_clusters[k][0] # select the first ID in the value list
            ID_string = dict_clusters[k]
            ID_string = ",".join(ID_string)
            Replicates = len(dict_clusters[k]) # count how many IDs
            print("{} ;Size={}\n{}".format(ID, Replicates, k), file=fout)
            print("*Sequence;{}\nIDs;{}".format(k, ID_string), file=fout1)
