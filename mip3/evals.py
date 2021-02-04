# Copyright (c) 2014 - Niek de Klein, Enrico Magnani, Seung Y. Rhee
#
#     This file is part of microProtein Prediction Program (miP3) version 2.
#
#     microProtein Prediction Program (miP3) version 2 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     microProtein Prediction Program (miP3) version 2 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with miP3.py.  If not, see <http://www.gnu.org/licenses/>.")


import csv
import os
import argparse
import newblast
from Bio import SeqIO
import re


def fix_file_names(value):
    """
    Normalizes string, converts to lowercase, removes non-alpha characters,
    and converts spaces to hyphens.
    """
    filename, file_extension = os.path.splitext(value)
    value = unicode(re.sub('[^\w\s-]', '', filename).strip().lower())
    value = re.sub('[-\s]+', '-', value)
    value += file_extension
    return(value)

def validate(seq):
    alphabets = re.compile('^[acdefghiklmnpqrstvwyx]*$', re.I)
    
    if alphabets.search(seq) is not None:
        return True
    else:
        return False


parser = argparse.ArgumentParser(description='evals', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--all_protein_file', help='Fasta file containing all proteins', required=True)                 # all_proteins_path
parser.add_argument('-i','--interest_file', help='Fasta file containing proteins of interest', required=True)    # interest_protein_path
parser.add_argument('-o','--output_file', help='Output file name', required=True)
parser.add_argument('-b','--blast_folder', help='Folder that contains makeblastd, blastp and rpsblast')
parser.add_argument('-r','--results_file', help='Results file name', required=True)
args = vars(parser.parse_args())
all_proteins_path = args['all_protein_file']
interest_proteins_path = args['interest_file']
output_path = args['output_file']
blast_folder = args['blast_folder'].strip('"')
results_path = args['results_file']
if len(blast_folder) > 0 and not blast_folder.endswith(os.sep):
        blast_folder+=os.sep

script_path = os.path.realpath(__file__).rstrip(__file__.split('/')[-1])      # find path where miP3 is, so that files can be saved relative to this script
        

# start program
''' <original description of function. code is taken out of function, description
     is relevant till output.write_fasta(small_proteins..)>
    Takes the fasta files given by the user and makes 2 dictionaries, one containing interesting proteins (query proteins) and the other containing
    all proteins.
    Both dictionaries have as key the name of the protein and value the protein sequence.

    Result:
    2 dictionaries, one containing interesting proteins (query proteins) and the other containing all proteins
    with as key the name of the protein and value the protein sequence.
    Also writes out a fasta file with all proteins that are 200a.a. or smaller.
    '''
#blastdata = newblast.blastp(all_proteins_path, interest_proteins_path, blast_folder, output_path)
#blast_eval = newblast.getSubjectInfo(blastdata)
#print blast_eval

#count the number of sequences in fasta file
countall = 0
for record in SeqIO.parse(open(all_proteins_path),'fasta'):
    countall = countall + 1
countinterest = 0
for record in SeqIO.parse(open(interest_proteins_path),'fasta'):
    countinterest = countinterest + 1
matrix = [[1 for i in range(countall)] for j in range(countinterest)]

countall = 0
fasta_sequences = SeqIO.parse(open(all_proteins_path),'fasta')
for fasta in fasta_sequences:
    output_handle = open("example.fasta", "w")
    SeqIO.write(fasta, output_handle, "fasta")
    output_handle.close()
    fasta_sequences2 = SeqIO.parse(open(interest_proteins_path),'fasta')
    countinterest = 0
    for fasta2 in fasta_sequences2:
        output_handle2 = open("example2.fasta", "w")
        SeqIO.write(fasta2, output_handle2, "fasta")
        output_handle2.close()
        blastdata = newblast.blastp("example.fasta", "example2.fasta", blast_folder, output_path)
        blast_eval = newblast.getSubjectInfo(blastdata)
        print "fbox protein number"
        print countinterest
        print "decoy protein number"
        print countall
        print blast_eval
        matrix[countinterest][countall] =  blast_eval
        countinterest = countinterest + 1
    countall = countall + 1

#write results to csv file
with open(results_path, 'wb') as csvfile:
    mywriter = csv.writer(csvfile)
    for i in range(len(matrix)):
        mywriter.writerow(matrix[i])
