# -*- coding: utf-8 -*-
# ptw -c --runner='python gene_finder.py'
"""
Gene Finder - SoftDes Mini Project 1

@author: Lauren Gulland

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide): #week 1
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
        All four cases are simple to test, so I decided to test all 4 cases to 
        make sure there were no simple fallacies
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """
    if nucleotide=='C':
        return 'G'
    elif nucleotide=='G':
        return 'C'
    elif nucleotide=='A':
        return 'T'
    elif nucleotide=='T':
        return 'A'
    else:
        return None


def get_reverse_complement(dna): #week 1
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement('ATGCGAATGTAGCATCAAA')
    'TTTGATGCTACATTCGCAT'
    """
    reverse=''

    for i in reversed(range(0, len(dna))):
        reverse+=get_complement(dna[i])

    return reverse
    


def rest_of_ORF(dna): #week 1
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    for i in range(0, len(dna), 3):
        if len(dna[i:])>=3 and (dna[i:i+3]=='TAG' or dna[i:i+3]=='TGA' or dna[i:i+3]=='TAA'):
                return dna[:i]
    return dna



def find_all_ORFs_oneframe(dna): #week 1
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    orfList=[]
    i=0
    while i < len(dna)-2:
        if dna[i:i+3]=='ATG':
            orf=rest_of_ORF(dna[i:])
            orfList.append(orf)
            i+=len(orf)
        else:
            i+=3
    return orfList
    


def find_all_ORFs(dna): #week 1
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    orfList=[]
    for i in range(0, 3):
        orfList.extend(find_all_ORFs_oneframe(dna[i:]))
    return orfList
    


def find_all_ORFs_both_strands(dna): #week 1
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    orfList=[]
    orfList.extend(find_all_ORFs(dna))
    orfList.extend(find_all_ORFs(get_reverse_complement(dna)))
    return orfList
    


def longest_ORF(dna): #week 2
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfList = find_all_ORFs_both_strands(dna)
    longest = orfList[0]
    for i in orfList:
        if len(i)>len(longest):
            longest = i
    return longest




def longest_ORF_noncoding(dna, num_trials): #week 2
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    maxLength = 0
    longestORF = ''
    for i in range(0, num_trials):
        orfTest = longest_ORF(shuffle_string(dna))
        if len(orfTest)>maxLength:
            longestORF = orfTest
            maxLength= len(orfTest)
    return longestORF


def coding_strand_to_AA(dna): #week 2
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    aaStrand = ''
    for i in range(0, len(dna)-2, 3):
        aaStrand += aa_table[dna[i:i+3]]   
    return aaStrand



def gene_finder(dna): #week 2
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = len(longest_ORF_noncoding(dna, 1500))
    orfList = find_all_ORFs_both_strands(dna)
    aminoList=[]
    for i in orfList:
        if len(i)>= threshold:
            aminoList.append(coding_strand_to_AA(i))
    return aminoList

from load import load_seq
dna = load_seq("./data/X73525.fa")
salmonellaGenes = gene_finder(dna) #produces list of candidate genes for the provided salmonella gene
print salmonellaGenes


if __name__ == "__main__":
    import doctest
    doctest.testmod()
