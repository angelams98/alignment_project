#!/usr/bin/env python3
from functions_align import *
import re
import sys


#PSEUDOCODE

accepted_filetypes = ['.fsa', '.fasta', '.fna']

#Read the file
if len(sys.argv) == 2:
    infile = sys.argv[1]
else:
    infile = input("Give the name of the infile: ")



# Check that the specified file is a fasta file
#O(m), m is the number of elements in accepted_filetypes
if re.search(r'\.\w+$', infile).group(0) in accepted_filetypes:

    #Try to open the file and generates an error message if it fails
    try:
        #O(1), open a file
        infile = open(infile, 'r')

    except IOError as err:
        print("Can't open file, reason:", str(err))
        sys.exit(1)

else:
    print('The given file was not a fasta-file.')
    sys.exit(1)


#Initialize the variables
seq_list = []
title = []
sequences = ""
nucleotides = ["A", "T", "G", "C"]
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
amino_acids_checker = ["D", "E", "F", "H", "I", "K", "L", "M", "N", 
                        "P", "Q", "R", "S", "V", "W", "Y"]
blosum_files = ["30", "35", "40", "45", "50", "55", "62", "65", "70", 
                "75", "80", "85", "90", "100"]
is_protein = False
pos = 0


#Saves the header lines and the sequences in lists
#O(n), n is the number of lines in the file
for line in infile:
    if line.startswith('>'):
        #O(1), it is appending an element
        title.append(line[:-1])

        #If in the next line after ">", sequences is not empty, it is added to seq_list
        if sequences != "":
            #O(1), it is appending an element
            seq_list.append(sequences)
        sequences = ""

    else:
        sequences += line[:-1] 

if sequences != "":
    #O(1), it is appending an element
    seq_list.append(sequences)


#Create an only string with the sequences
seq = "".join(seq_list)


#Check if sequence from infile is either a DNA or a protein sequence
#O(o), o is the number of sequences (theoretically only 2)
while pos < len(seq):
    
    #O(p), p is the number of elements in amino_acid_checker list
    if seq[pos] in amino_acids_checker:
        is_protein = True

    #O(q), q is the number of elements in amino_acids list
    elif seq[pos] not in amino_acids:
        print("The sequence found in the file", infile.name, "contained an impure sequence.")
        sys.exit(1)
        
    pos = pos + 1


#Once we know if it is DNA or protein, we can decide which alignment methods should the program use

#They are DNA sequences
if is_protein == False:
    
    dna_prot = "dna"
    #This variable is not used in DNA
    blosum = 0

    print("What you have given me is a DNA sequence \n")
    alignment = input("So do you want to do local or global alignment?\n")
    alignment = alignment.upper()
    print("Please tell me which parameters you want to use\n")

    match = float(input("Give me the match value: "))
    mismatch = float(input("Give me the mismatch value: "))
    opening = float(input("Give me the indel value: "))
    extension = float(input("Give me the extension value: "))
    print("\n")


    #O(r), r is the strings we are comparing in the conditional
    if alignment == "GLOBAL":

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("And the parameters are match:{}, mismatch:{}, opening score:{} and extension score:{}\n".format(match, mismatch, opening, extension))
        
        (matrix, alignment) = alignment_nw(seq_list[0], seq_list[1], dna_prot, match, mismatch, opening, extension)
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)

    #O(q), q is the strings we are comparing in the conditional
    elif alignment == "LOCAL":

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("And the parameters are match:{}, mismatch:{}, opening score:{} and extension score:{}\n".format(match, mismatch, opening, extension))

        (matrix, alignment, score) = alignment_sw(seq_list[0], seq_list[1], dna_prot, match, mismatch, opening, extension)

        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)
        print("The score of the best alignment is", score[0])
    
    else:
        print("I cannot recognize this alignment")
        sys.exit(1)


#They are protein sequences
elif is_protein == True:

    dna_prot = "protein"
    #Those variables are not used in proteins
    match = 0 
    mismatch = 0
    
    print("What you have given me is a protein sequence \n")
    alignment = input("So do you want to do local or global alignment?\n")
    alignment = alignment.upper()

    print("Please tell me which parameters you want to use\n")
    opening = float(input("Give me the indel value: "))
    extension = float(input("Give me the extension value: "))
    blosum_int = input("Give me which blosum matrix you want to use (introduce the number only): ")
    print("\n")

    #O(s), s is the number of the elements in the blosum_files list
    if blosum_int in blosum_files:
        blosum = blosum_int 

    #O(s), s is the number of the elements in the blosum_files list
    if blosum_int not in blosum_files:
        print("I don't have that file, I am using BLOSUM62 instead\n")
        blosum = "62"



    if alignment == "GLOBAL" :

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("And the parameters are opening score:{}, extension score:{} and blosum_matrix\n".format(opening, extension, blosum))

        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment_nw(seq_list[0], seq_list[1], dna_prot, opening, extension, match, mismatch))
        

    elif alignment == "LOCAL":

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("And the parameters are opening score:{}, extension score:{} and blosum_matrix\n".format(opening, extension, blosum))
        
        (matrix, alignment, score) = alignment_sw(seq_list[0], seq_list[1], dna_prot, opening, extension, match, mismatch)
        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)
        print("The score of the best alignment is", score[0])

    else:
        print("I cannot recognise the alignment\n")
        sys.exit(1)