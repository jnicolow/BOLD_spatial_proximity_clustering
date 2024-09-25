"""
This module contains functions used for converting DNA sequences to amino acid sequences and then extracting protine sequences from the aminoacid sequences

Author: Joel Nicolow, Department of Information and Computer Sciences, University of Hawaii at Manoa (October 20, 2023)
"""

from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
def DNA_to_aminoacid(sequence, table="Invertebrate Mitochondrial"):
    """
    convert DNA to amino acid sequences
    :param sequence: Biopython Seq
    :param table: String which conversion table fi be used (default Invertebrate Mitochondrial)
    """
    aminoAcidSeq = sequence.translate(table)# choose which translation (these should be invertabrits)
    return(aminoAcidSeq)

def longest_p_sequence(sequence):
    sequences = sequence.split('*')
    sequences = [seq[seq.find('M'):] if 'M' in seq else '' for seq in sequences]  # trim until the first 'M'
    if len(sequences) == 0: return None
    longest = max(sequences, key=len)  # find the longest sequence
    return longest


def get_aminoacid_seq(DNAseq, table="Invertebrate Mitochondrial", trim=False):
    """
    try the open reading frames and the reverse complement and which ever returns the longest sequence take that as the protein sequence

    :param DNAseq: String DNA sequence (ATG kine)
    :param table: String name of conversion table used in Biopython Bio.Seq.translate() (default "Invertebrate Mitochondrial")
    :param trim: Boolean trim DNAseq to be length multiple of 3 (default False)

    """
    if trim:
        DNAseq = DNAseq[:len(DNAseq)//3*3]
    sequence = Seq(DNAseq)
    sequenceR = sequence.reverse_complement() # this is the other side of the DNA strand
    protineSeqs = []
    # itterate through the three reading frames (start from begining, drop first letter, or drop first two)
    for i in range(0, 3, 1):
        F = DNA_to_aminoacid(sequence[i:], table) # forword reading frame
        F_longest_sequence = longest_p_sequence(F)
        protineSeqs.append(F_longest_sequence)
        R = DNA_to_aminoacid(sequenceR[i:], table) # reverse reading frome
        R_longest_sequence = longest_p_sequence(R) # reverse reading frame
        protineSeqs.append(R_longest_sequence)

    return(str(max(protineSeqs, key=len))) # convert it back to a string
