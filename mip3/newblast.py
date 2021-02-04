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


from Bio.Blast.Applications import NcbiblastpCommandline
import ntpath
from Bio.Blast import NCBIXML


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def blastp(blast_file, blast_db, blast_path, blastoutput_xml):
    """
    Run a BLASTP search with blast_file against blast_db
    
    Args:
        blast_file (str):	  Path to fasta file used to BLAST against blast_db
        blast_db (str):		  Name of the database to BLAST against 
        blast_path (str):	  Path to the blastp program
        blastoutput_csv (str): location of the blast output csv file
    
    Returns:
        An iterable of blast records as returned by NCBIXML.parse
    """
    def cline(blastp_cline, blastoutput_xml):
        blastp_cline()
        result_handle = open(blastoutput_xml)
        blast_records = NCBIXML.parse(result_handle)
        return blast_records
    blastp_cline = NcbiblastpCommandline(blast_path+'/blastp', query=blast_file, subject=blast_db,
                                                                      outfmt=5, out=blastoutput_xml)

    return cline(blastp_cline, blastoutput_xml)


def getSubjectInfo(blast_records):
    '''Run through an iterator of blast records and save the results to a dictionary
    
       Args:
           blast_records (iterable):	An iterable of blast records
           prots_of_interest (dict):	A dictionary with proteins of interest. 
                                                                            This contains as key protein name (has to be same as protein names found in blast_records)
                                                                            and as value inform
    
       Returns:
           Dictionary with as key name of the subject (found) protein of a BLAST records and as value a dict with as key query (searched) protein from BLAST records
    '''
    
    #see if target is in TF list
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                return hsp.expect
