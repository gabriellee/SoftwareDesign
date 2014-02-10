# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: YOUR NAME HERE
"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    output = ""
    for s in L:
        output = output + s
    return output


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    # YOUR IMPLEMENTATION HERE
    cdn = [None]*int(len(dna)/3)
    aminoAcids = list([None]*int(len(dna)/3))
    for i in range(len(dna)/3):
    	cdn[i] = dna[i:i+3]
    	if cdn[i] == 'ATT' or cdn[i] == 'ATC' or cdn[i] == 'ATA':
    		aminoAcids[i] = 'I'
    	elif dna[i:i+2] == 'CT' or cdn[i] == 'TTA' or cdn[i] == 'TTG':
    		aminoAcids[i] = 'L'
    	elif dna[i:i+2] == 'GT':
    		aminoAcids[i]  = 'V'
    	elif dna[i:i+2] == 'TT':
    		aminoAcids[i] = 'F'
    	elif cdn[i] == 'ATG':
    		aminoAcids[i] = 'M'
    	elif cdn[i] == 'TGT' or cdn[i] == 'TGC':
    		aminoAcids[i] = 'C'
    	elif dna[i:i+2] == 'GC':
    		aminoAcids[i] = 'A'
    	elif dna[i:i+2] == 'GG':
    		aminoAcids[i] = 'G'
    	elif dna[i:i+2] == 'CC':
    		aminoAcids[i] = 'P'
    	elif dna[i:i+2] == 'AC':
    		aminoAcids[i] = 'T'
    	elif dna[i:i+2] == 'TC' or cdn[i] == 'AGT' or cdn[i] == 'AGC':
    		aminoAcids[i] = 'S'
    	elif dna[i:i+2] == 'TA':
    		aminoAcids[i] = 'Y'
    	elif cdn[i] == 'TGG':
    		aminoAcids[i] = 'W'
    	elif cdn[i] == 'CAA' or cdn[i] == 'CAG':
    		aminoAcids[i] = 'Q'
    	elif cdn[i] == 'AAT' or cdn[i] == 'AAC':
    		aminoAcids[i] = 'N'
    	elif dna[i:i+2] == 'CA':
    		aminoAcids[i] = 'H'
    	elif cdn[i] == 'GAA' or cdn[i] == 'GAG':
    		aminoAcids[i] = 'E'
    	elif cdn[i] == 'GAT' or cdn[i] == 'GAC':
    		aminoAcids[i] = 'D'
    	elif dna[i:i+2] == 'AA':
    		aminoAcids[i] = 'K'
    	else:
    		aminoAcids[i] = 'R'
    return ''.join(aminoAcids)

def coding_strand_to_AA_unit_tests():
    """ Unit tests for the coding_strand_to_AA function """
        
    # YOUR IMPLEMENTATION HERE
    dnacheck1 = "AAAAAAGAAAAGGACTCCTGTATG"
    out_hyp1 = "KKEKDSCM"
    out_act1 = coding_strand_to_AA(dnacheck1)
    print "input: " + dnacheck1 +", " +"expected output: "
    print out_hyp1 + ", actual output: "
    print out_act1
    
    dnacheck2 = "CAAATTCGT"
    out_hyp2 = "QIR"
    out_act2 = coding_strand_to_AA(dnacheck2)
    print "input: " + dnacheck2 +", " +"expected output: "
    print out_hyp2 + ", actual output: "
    print out_act2
    
    dnacheck3 = "CTTGTTCCTTAT"
    out_hyp3 = "LVPY"
    out_act3 = coding_strand_to_AA(dnacheck3)
    print "input: " + dnacheck3 +", " +"expected output: "
    print out_hyp3 + ", actual output: "
    print out_act3
    
    

def get_complement(dna):
    """ Computes the complementary sequence of DNA for the specfied DNA sequence
    
        dna: a DNA sequence represented as a string
        returns: the complementary DNA sequence represented as a string
    """
    
    # YOUR IMPLEMENTATION HERE
    rvs_dna = dna[::-1]

    rep = {'A':'T','T':'A','G':'C','C':'G'}
    rep_dict = dict(rep)
    repfun = lambda match: rep_dict[match.group(0)]
    pattern = rvs_dna.compile("|".join([rvs_dna.escape(i) for i, j in rep]), rvs_dna.M)
    rvs_cmpl = lambda string: pattern.sub(repfun, string)
    
    
    
    
    for i,j in dic.iteritems():
        rvs_cmpl = rvs_dna.replace(i,j)
def get_complement_unit_tests():
    """ Unit tests for the get_complement function """
        
    # YOUR IMPLEMENTATION HERE
    

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    
    # YOUR IMPLEMENTATION HERE

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """
        
    # YOUR IMPLEMENTATION HERE
        
def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  By non-nested we mean that an ORF that occurs entirely within
        another ORF will not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
        
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """

    # YOUR IMPLEMENTATION HERE

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""

    # YOUR IMPLEMENTATION HERE

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    # YOUR IMPLEMENTATION HERE

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """

    # YOUR IMPLEMENTATION HERE
