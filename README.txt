
Title: MDC.py
Version: 1.0
Date: 16-03-2020
Author: María Halldórsdóttir

Description:
  MDC program performs dereplication and clustering of sequences.

Methodology:
  List of functions:
  This program has 3 functions:
    - populate_lists_dictionary > Populates 2 lists and 1 dictionary used for the
      clustering and 1 dictionary for dereplication.
    - percentage_similarity > calculates % similarity between 2 sequences.
    - remove_empty_from_dict > removes empty keys from dictionary.

Installation:
  No installation needed, just download file and it is ready for use.
  Requirements: Python version 3.7

Usage:
  ./MDC.py -i "input file" -o "output Fasta file" -u "output log file" -p "sequence percentage similarity"
  If no -p is selected then the program performs dereplication.

Support:
  ma6364ha-s@student.lu.se

Future releases:
  Updated clustering algorithm, which will perform clustering more
  accurately.
