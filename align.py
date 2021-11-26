#!/usr/bin/env python3
from functions_align import *
import re
import sys


#PSEUDOCODE
    #Call read_files function to get a list with the sequences we are going to align
    #(seq, title) = read_files("dna7.fsa")

    #Choose 2 sequences and call alignment_simple function to create a list of list with the scores
    #To call this function you have to add the 2 sequences and the values for match, mismatch and Indel
    #possible values for match, mismatch and indel:
    #1, -1, -1
    #0, 1, 1
    #0, 1, 10
    # matrix = alignment_simple(seq[0], seq[1], 1, -1, -1)

    #If you want to print this matrix on screen, call print_matrix function
    #print_matrix(matrix)

    #To calculate the best alignment using the scores from the matrix, call traceback matrix
    #As inputs for the function you need the matrix, and the original strings
    #print(traceback(matrix, seq[0], seq[1]))


#The program can get several arguments: file, match, mismatch, indel, extension.
#The first argument is the name of the file
#The second, third, fourth and fifth value are integers

#List of links that contain different blosum matrix
#https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#http://www.cbcb.umd.edu/confcour/CMSC423-materials/BLOSUM80.txt
#https://github.com/noporpoise/seq-align/blob/master/scoring/BLOSUM50.txt

accepted_filetypes = ['.fsa', '.fasta', '.fna']

#Read the file
if len(sys.argv) == 2:
    try:
        infile = sys.argv[1]

    except IOError:
        print("The file you introduce is not accessible")
        sys.exit(1)

else:
    infile = input("Give the name of the infile: ")



# Check that the specified file is a fasta file
if re.search(r'\.\w+$', infile).group(0) in accepted_filetypes:

    #Try to open the file and generates an error message if it fails
    try:
        #O(1), open a file
        infile = open(infile, 'r')

    except IOError as err:
        print("can't open file, reason:", str(err))
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
amino_acids_checker = ["D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "V", "W", "Y"]
is_protein = False
pos = 0


#Saves the titles and the sequences in lists
#O(m), m is the number of lines in the file
for line in infile:
    if line.startswith('>'):
        title.append(line[:-1])

        #If in the next line after ">", sequences is not empty, it is added to seq_list
        if sequences != "":
            seq_list.append(sequences)
        sequences = ""

    else:
        sequences += line[:-1] 

if sequences != "":
    seq_list.append(sequences)


#Create an only string with the sequence
seq = "".join(seq_list)


# Check if sequence from infile is either a DNA or a protein sequence
#O(n), n is the number of sequences (theoretically only 2)
while pos < len(seq):
    
    if seq[pos] in amino_acids_checker:
        is_protein = True

    #O(p), p is the number of elements in amino_acids list
    elif seq[pos] not in amino_acids:
        print("The sequence found in the file", infile.name, "contained an impure sequence.")
        sys.exit(1)
        
    pos = pos + 1


#Once we know if it is DNA or protein, we can decide which alignment methods should the program use

#They are DNA sequences
if is_protein == False:

    print("What you have given me is a DNA sequence \n")
    alignment = input("So do you want to do local or global alignment?\n")
    alignment = alignment.upper()
    dna_prot = "dna"
    blosum = 0

    #O(q), q is the strings we are comparing in the conditional
    if alignment == "GLOBAL":

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        match = int(input("Give me the match value: "))
        mismatch = int(input("Give me the mismatch value: "))
        opening = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        print("\n")

        
        (matrix, alignment) = alignment_nw(seq_list[0], seq_list[1], dna_prot, match, mismatch, opening, extension)
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)

    #O(q), q is the strings we are comparing in the conditional
    elif alignment == "LOCAL":

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        match = int(input("Give me the match value: "))
        mismatch = int(input("Give me the mismatch value: "))
        opening = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        print("\n")



        (matrix, alignment, score) = alignment_sw(seq_list[0], seq_list[1], dna_prot, match, mismatch, opening, extension)

        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)
        print("The score of the best alignment is", score[0])

        


#They are protein sequences
elif is_protein == True:

    print("What you have given me is a protein sequence \n")
    alignment = input("So do you want to do local or global alignment?\n")
    alignment = alignment.upper()
    dna_prot = "protein"
    match = 0 
    mismatch = 0

    #O(q), q is the strings we are comparing in the conditional
    if alignment == "GLOBAL" :

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        opening = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        blosum_int = input("Give me which blosum matrix you want to use (introduce the number only):")
        print("\n")

        if blosum_int == "62":
            blosum = blosum_matrix("BLOSUM62.txt")
        elif blosum_int == "80":
            blosum = blosum_matrix("BLOSUM80.txt")
        elif blosum_int == "50":
            blosum = blosum_matrix("BLOSUM50.txt")
        else:
            print("I don't have that file, I am using BLOSUM62 instead\n")
            blosum = blosum_matrix("BLOSUM62.txt")

        
        (matrix, alignment) = alignment_nw(seq_list[0], seq_list[1], dna_prot, opening, extension, match, mismatch)
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)
        

    #O(q), q is the strings we are comparing in the conditional
    elif alignment == "LOCAL":

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        opening = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        blosum_int = input("Give me which blosum matrix you want to use (introduce the number only): ")
        print("\n")


        if blosum_int == "62":
            blosum = blosum_matrix("BLOSUM62.txt")
        elif blosum_int == "80":
            blosum = blosum_matrix("BLOSUM80.txt")
        elif blosum_int == "50":
            blosum = blosum_matrix("BLOSUM50.txt")
        else:
            print("I don't have that file, I am using BLOSUM62 instead\n")
            blosum = blosum_matrix("BLOSUM62.txt")

    
        (matrix, alignment, score) = alignment_sw(seq_list[0], seq_list[1], dna_prot, opening, extension, match, mismatch)

        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)
        print("The score of the best alignment is", score[0])
