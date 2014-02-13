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
    #import pdb
    for s in L:
        #pdb.set_trace()
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
    #import re
    rvs_dna = dna[::-1]
    rvs_cmpl = list(dna)

    #rep = {'A':'T','T':'A','G':'C','C':'G'}
    for i in range(len(dna)):
        if rvs_dna[i] == 'A':
            rvs_cmpl[i] = 'T'
        elif rvs_dna[i] == 'T':
            rvs_cmpl[i] = 'A'
        elif rvs_dna[i] == 'G':
            rvs_cmpl[i] = 'C'
        else:
            rvs_cmpl[i] = 'G'
    outcmp = ''.join(rvs_cmpl)
    return outcmp

    #rep_dict = dict(rep)
    #repfun = lambda match: rep_dict[match.group(0)]
    #pattern = re.compile("|".join([re.escape(i) for i, j in rep]), re.M)
    #rvs_cmpl = lambda rvs_dna: pattern.sub(repfun, rvs_dna)
    
    
    
    
    #for i,j in dic.iteritems():
        #rvs_cmpl = rvs_dna.replace(i,j)
def get_complement_unit_tests():
    """ Unit tests for the get_complement function """
        
    # YOUR IMPLEMENTATION HERE
    rvscheck1 = "CAAATTCGT"
    out_hyp1 = "ACGAATTTG"
    out_act1 = get_complement(rvscheck1)
    print "input: " + rvscheck1 +", " +"expected output: "
    print out_hyp1 + ", actual output: "
    print out_act1
    
    rvscheck2 = "CTTGTTCCTTAT"
    out_hyp2 = "ATAAGGAACAAG"
    out_act2 = get_complement(rvscheck2)
    print "input: " + rvscheck2 +", " +"expected output: "
    print out_hyp2 + ", actual output: "
    print out_act2

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    
    # YOUR IMPLEMENTATION HERE
    #import pdb
    
    cdn = [None]*int(len(dna)/3)
    ORF = list()
    for n in range(len(dna)/3):
        cdn[n] = dna[n*3:n*3+3]
        #pdb.set_trace()
        if cdn[n] == 'TAG' or cdn[n] == 'TAA' or cdn[n] == 'TGA':
            break
        else:
            ORF.append(cdn[n])

    outORF = ''.join(ORF)
    return outORF
    

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """
        
    # YOUR IMPLEMENTATION HERE
    check3 = "CTTAGTGTTCCTTATAAATAG"
    out_hyp3 = "CTTAGTGTTCCTTATAAA"
    out_act3 = rest_of_ORF(check3)
    print "input: " + check3 +", " +"expected output: "
    print out_hyp3 + ", actual output: "
    print out_act3
    
    check4 = "CTTGTTCCTTATTGAAAATAGGTTTAA"
    out_hyp4 = "CTTGTTCCTTAT"
    out_act4 = rest_of_ORF(check4)
    print "input: " + check4 +", " +"expected output: "
    print out_hyp4 + ", actual output: "
    print out_act4
def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  By non-nested we mean that an ORF that occurs entirely within
        another ORF will not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE
    import pdb
    ind = 0
    
    # truncate the length of the DNA sequence to be a multiple of 3
    dna = dna[:len(dna)-len(dna)%3]
    ORFlist = list()
    #ind is the index of the value in dna, it goes by 3s
    while ind < len(dna):
            cdn = [None]*int(len(dna)/3)
            #pdb.set_trace()
            for n in range(ind/3,len(dna)/3):#  look for a start codon until you get to the last codon, then restart the loop at the next codon after the reading frame.  If you get to the last codon and do not find a start codon, end the while loop. n is the index in cdn.
                cdn[n] = dna[n*3:n*3+3]
                #pdb.set_trace()
                if cdn[n] == 'ATG':
                    ORF = rest_of_ORF(dna[3*n:len(dna)])
                    ind = len(ORF)+3*n
                    ORFlist.append(ORF)
                    break
                if n == len(dna)/3 - 1:
                    ind = len(dna)
    #pdb.set_trace()
            
    return ORFlist
        

def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
        
    # YOUR IMPLEMENTATION HERE
    check4 = "GCAATGAAATGAAATATGGGGTGA"
    out_hyp4 = "'ATGAAA', 'ATGGGG'"
    out_act4 = find_all_ORFs(check4)
    print "input: " + check4 +", " +"expected output: "
    print out_hyp4 + ", actual output: "
    print out_act4
    
    check3 = "CAGTTTATGAGTGTTTAGCCTATGTATAAA"
    out_hyp3 = "ATGAGTGTT, ATGTATAAA"
    out_act3 = find_all_ORFs(check3)
    print "input: " + check3 +", " +"expected output: "
    print out_hyp3 + ", actual output: "
    print out_act3
    
    check5 = "ATATATATATATATATATATATATA"
    out_hyp5 = ""

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    # YOUR IMPLEMENTATION HERE
    #import pdb
    cmp = get_complement(dna)
    combolist = list()
    dnaORFlist = find_all_ORFs(dna)
    cmpORFlist = find_all_ORFs(cmp)
    #pdb.set_trace()
    combolist = dnaORFlist + cmpORFlist
    return combolist
    

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """
    

    # YOUR IMPLEMENTATION HERE
    check3 = "GGCATGCTAAAATTACTATAGCAT"
    out_hyp3 = "ATGCTAAAATTACTA ATGCTA"
    out_act3 = find_all_ORFs_both_strands(check3)
    print "input: " + check3 +", " +"expected output: "
    print out_hyp3 + ", actual output: "
    print out_act3
    
    check4 = "'ATGGCCCATTAGCTAATG'"
    out_hyp4 = "'ATGGCCCAT','ATG', 'ATGGCCAT'"
    out_act4 = find_all_ORFs_both_strands(check4)
    print "input: " + check4 +", " +"expected output: "
    print out_hyp4 + ", actual output: "
    print out_act4

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""

    # YOUR IMPLEMENTATION HERE
    import pdb
    allORFs = find_all_ORFs_both_strands(dna)
    mymax = 0
    #pdb.set_trace()
    for n in range(len(allORFs)):
        if len(allORFs[n]) > mymax:
            mymax = len(allORFs[n])
    return mymax
    

def longest_ORF_unit_tests():
    
    """ tests"""

    # YOUR IMPLEMENTATION HERE
    check3 = "GGCATGCTAAAATTACTATAGCAT"
    out_hyp3 = "ATGCTAAAATTACTA"
    out_act3 = longest_ORF(check3)
    print "input: " + check3 +", " +"expected output: "
    print out_hyp3 + ", actual output: "
    print out_act3
    
    check4 = "ATGGCCCATTAGCTAATG"
    out_hyp4 = "ATGGCCCAT"
    out_act4 = longest_ORF(check4)
    print "input: " + check4 +", " +"expected output: "
    print out_hyp4 + ", actual output: "
    print out_act4

import random


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    # YOUR IMPLEMENTATION HERE
    dna_list = list(dna)
    maxORF = 0
    mew=''
    import pdb
    #pdb.set_trace()
    for i in range(num_trials):
        if i%10 == 0:
            print i
        random.shuffle(dna_list)
        #print dna_list
        mew=''.join(dna_list)
        #pdb.set_trace()   
        #print me
        ORF = longest_ORF(mew)
        #print ORF
        #print maxORF
        if (ORF > maxORF):
            maxORF = ORF   
    return maxORF
    
#    collapse(dna_list)
#    g=collapse(dna_list)
#    hi=list(dna)
#    ji=random.shuffle(dna_list)
#    mew=collapse(dna_list)
#    find_all_ORFs(mew)
    

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
    ORFlist = find_all_ORFs_both_strands(dna)
    AAlist = list()
    for n in range(len(ORFlist)):
        if len(ORFlist[n]) > threshold:
            AAlist.append(coding_strand_to_AA(ORFlist[n]))
    return AAlist
    
from load import load_seq
dna = load_seq("./data/X73525.fa")
    
h = longest_ORF_noncoding(dna,1500)
print h
#find_all_ORFs_unit_tests()