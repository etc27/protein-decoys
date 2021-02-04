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




import timing
import os
# imports from predictmiP
import argparse
import blast
import interpro_scanner
import output
from Bio import SeqIO
import pickle
import sys
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


#not explained:
#'AT4G25690.1', 'AT4G25690.2','AT3G60890.1']
not_found = ['AT1G74500', 'AT5G07260']
# AT2G35940 can't be found because it's an alternative splice mip, but the alternative splice form is not found tair database
# AT4G25670 and AT5G52550 cannot be found with this method because they are too dissimalar to the transcription factors they diverged from to find with BLAST
# not in TAIR all proteins file (no protein data at tair, only genomic)
# AT2G45450 don't know why not found, homolog of AT1G52150
# AT2G35940 can't be found because it's an alternative splice mip, but the alternative splice form is not found tair database
# Parsing commandline
parser = argparse.ArgumentParser(description='microProtein prediction program. See the README.txt for a more detailed description.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--all_protein_file', help='Fasta file containing all proteins', required=True)                 # all_proteins_path
parser.add_argument('-i','--interest_file', help='Fasta file containing proteins of interest', required=True)    # interest_protein_path
# changed pfam domain names to domain names of any database in InterproScan, but don't have time now to change the variable names
parser.add_argument('-f','--pfam_file', help='Text file with unwanted domain names', required=True)
parser.add_argument('-o','--output_file', help='Output file name', required=True)
parser.add_argument('-b','--blast_folder', help='Folder that contains makeblastd, blastp and rpsblast')
parser.add_argument('-a','--all_size', help='Max size of all proteins to search for',default=550)
parser.add_argument('-e','--eval_all', help='E-value limit to use with blast against all proteins',default=0.0000001)
parser.add_argument('-d','--developing', help='When developing flag is set, extra info is printed out and it uses the blast result from previous run.', action='store_true' )
parser.add_argument('-m','--email', help='Email is necessary for InterproScan webservice', required=True)
args = vars(parser.parse_args())
all_proteins_path = args['all_protein_file']
interest_proteins_path = args['interest_file']
pfam_domains_file_path = args['pfam_file']
outfile_path = args['output_file']
blast_folder = args['blast_folder'].strip('"')
print 'thresholds are:'
print 'all size (-s):',args['all_size']
print 'evalue all (-e):', args['eval_all']
print 'email (-m):', args['email']
print 'developing:', args['developing']
developing = args['developing']
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
proteins_of_interest = {}
interest_list = []
for seq_record in SeqIO.parse(interest_proteins_path, "fasta"):
    proteins_of_interest[' '.join(seq_record.description.split())] = seq_record.seq
    interest_list.append(seq_record.seq)
all_proteins = {}
proteins_500aa = {}
for seq_record in SeqIO.parse(all_proteins_path, 'fasta'):
    seq = str(seq_record.seq).replace('*','').strip()
    if not validate(seq):
        print 'Not a valid protein sequence'
        print seq_record.description
        print seq
    ##    sys.exit()
    else:
        all_proteins[' '.join(seq_record.description.split())] = seq
        if len(seq_record.seq) <= int(args['all_size']):
            proteins_500aa[' '.join(seq_record.description.split())] = seq
# also add the transcription factors because they can be from two different sources, so difference in title
proteins_of_interest = {}
for seq_record in SeqIO.parse(interest_proteins_path, "fasta"):
    seq = str(seq_record.seq).replace('*','')
    if not validate(seq):
        print 'Not a valid protein sequence'
        print seq_record.description
        print seq
    ##    sys.exit()
    else:
        all_proteins[' '.join(seq_record.description.split())] = seq
        proteins_of_interest[' '.join(seq_record.description.split())] = seq_record.seq
        if len(seq_record.seq)  <= int(args['all_size']):
            proteins_500aa[' '.join(seq_record.description.split())] = seq_record.seq

organism_name = all_proteins_path.split(os.sep)[-1].split('.')[0]
if not os.path.exists(script_path+'fasta_files'):
    os.makedirs(script_path+'fasta_files')
# proteins_500aa and small_proteins given twice, first to get ids, second to get
# the proteins (see function)
output.write_fasta(proteins_500aa, proteins_500aa, script_path+'fasta_files'+os.sep+organism_name + '_all_proteins_smaller_than_500aa.fasta')
output.write_fasta(proteins_of_interest, proteins_of_interest, script_path+'fasta_files'+os.sep+organism_name+'_proteins_of_interest.fasta')

if not os.path.exists(script_path+'databases'):
    os.makedirs(script_path+'databases')
if not os.path.exists(script_path+'blast_results'):
    os.makedirs(script_path+'blast_results')
blast_all_pickle = script_path+'blast_results'+os.sep+fix_file_names(organism_name+'_blast_all_'+str(args['eval_all'])+'.p')
if developing:
    if os.path.isfile(blast_all_pickle):
        print 'BLAST records with this e-value against ALL database already exists, loading: '+blast_all_pickle
        try:
            subject_info_allDB = pickle.load( open(blast_all_pickle, 'rb' ))
        except (EOFError, KeyError) as e:
            repr(e)
            print('Removing pickle file, redoing blast')
            os.remove(blast_all_pickle)
###### use if not instead of else incase the error above got thrown, and file got removed
if not developing or not os.path.isfile(blast_all_pickle):
    print 'blasting against all'
    blast.makeBLASTdb(script_path+'fasta_files'+os.sep+organism_name + '_all_proteins_smaller_than_500aa.fasta', script_path+'databases'+os.sep+'allDB_'+organism_name, blast_folder)       # make all proteins database
    blast_records = blast.blastp(interest_proteins_path, script_path+'databases'+os.sep+'allDB_'+organism_name, args['eval_all'], blast_folder, script_path+'blast_results/'+organism_name+'_blastpAllOutput.xml')
    subject_info_allDB = blast.getSubjectInfo(blast_records, proteins_of_interest, args['eval_all'])
    if developing:
        print 'saving in '+blast_all_pickle
        if not os.path.exists(script_path+'blast_results'):
            os.makedirs(script_path+'blast_results')
        f = open(blast_all_pickle, 'wb' )
        pickle.dump( subject_info_allDB, f )
print('len subject_info_allDB after < '+str(args['all_size'])+'a.a. Protein blast: '+str(len(subject_info_allDB))+'\n')



subject_info_total = dict(subject_info_allDB.items())
if len(subject_info_total) == 0:
    print('No putative mip found, exiting program')                              # tell about the error
    sys.exit()        

print('len subject_info_total: '+str(len(subject_info_total))+'\n') 
#print('len subject_info_total after combining smallDB and allDB results: '+str(len(subject_info_total))+'\n')               # just a print >>

# delete the old dicts
del(subject_info_allDB)
##del(subject_info_smallDB)

# Remove TFs that are exact match with the miP
subject_info_filtered_on_compared_length = {}
mip_longer = False
for subject in subject_info_total:
    for query in subject_info_total[' '.join(subject.split())]['query_title']:
        ###adding step to check against all protein subjects instead of just the one it matched with
        if all_proteins[' '.join(subject.split())] in interest_list:
            continue
        #    continue if the subject is in the interest proteins list
        else:
            subject_info_filtered_on_compared_length[' '.join(subject.split())] = subject_info_total[' '.join(subject.split())]
if developing:
    for subject in subject_info_filtered_on_compared_length:
        for prot in not_found:
            if prot.lower() in subject.lower():
                print('found in subject_info_filtered_on_compared_length')
                print(subject)

del(subject_info_total)

print('len subject_info_filtered_on_compared_length after filtering on comparative length: '+str(len(subject_info_filtered_on_compared_length)))
###### Pfam search is removed because it could be difficult for users to get it working. Instead, uses the interpro scan server. If this takes to long
###### try uncommenting this line and adding 'import pfam' at the top of the file (check out pfam.py before)
#subject_info_filter_on_pfam = pfam.pfamSearch(subject_info_filtered_on_compared_length, all_proteins)
########################################################################################
subject_info_filtered_on_domains = interpro_scanner.interproScan(subject_info_filtered_on_compared_length, all_proteins,pfam_domains_file_path, args['email'], developing, script_path)  # filter on domain information gained from interproscan
if developing:
    for subject in subject_info_filtered_on_domains:
        for prot in not_found:
            if prot.lower() in subject.lower():
                print('found in subject_info_filtered_on_domains')
                print(subject)

output.write_fasta(subject_info_filtered_on_domains, all_proteins, outfile_path)

    



