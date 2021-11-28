#!/usr/bin/env python3
from functions_align import *
import re
import sys

"""
PSEUDOCODE
    Initialize the variables.
    Get the name of the file and check if its extension is allowed using a regular expression.
    Open the file
    Read the file and save the header lines and sequences in lists
    Turn the sequences into a unique string

    Check if the string contains nucleotides, amino acids or more characters.
        If the string contain different characters, it is considered an impure sequence.
        Program exists.

    If the sequences contain only nucleotides:
        Ask for the alignment until a valid answer is given.
        Ask for the paramaters.
        If mismatch, indel or opening are positive: 
            it turns them into negative numbers.

        If alignment is global:
            It calculates the global alignment.
            It prints the alignment

        If alignment is local:
            It calculates the local alignment.
            It prints the alignment
    
    If the sequences contain only amino acids:
        Ask for the alignment until a valid answer is given.
        Ask for the paramaters.
        If indel or opening are positive: 
            It turns them into negative numbers.

        If alignment is global:
            It calculates the global alignment.
            It prints the alignment

        If alignment is local:
            It calculates the local alignment.
            It prints the alignment


"""

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
accepted_filetypes = ['.fsa', '.fasta', '.fna']
alignment = ""


#Read the file
if len(sys.argv) == 2:
    infile = sys.argv[1]
else:
    infile = input("Give the name of the infile: ")



# Check that the specified file is a fasta file
if re.search(r'\.\w+$', infile).group(0) in accepted_filetypes:

    #Try to open the file and generates an error message if it fails
    try:
        infile = open(infile, 'r')

    except IOError as err:
        print("Can't open file, reason:", str(err))
        sys.exit(1)

else:
    print('The given file was not a fasta-file.')
    sys.exit(1)



#Saves the header lines and the sequences in lists
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


#Create an only string with the sequences
seq = "".join(seq_list)


#Check if sequence from infile is either a DNA or a protein sequence
while pos < len(seq):
    
    if seq[pos] in amino_acids_checker:
        is_protein = True

    elif seq[pos] not in amino_acids:
        print("The sequence found in the file", infile.name, "contained an impure sequence.")
        sys.exit(1)
        
    pos = pos + 1


#Once we know if it is DNA or protein, we can decide which alignment methods the program should use

# Sequences are DNA
if is_protein == False:
    
    dna_prot = "dna"
    #This variable is not used in DNA
    blosum = 0

    print("The file provided contained DNA sequences.\n")
    
    while alignment not in ['GLOBAL', 'LOCAL']:
        alignment = input("Do you want to perform a local or a global alignment?\n").upper()
    
    print()
    print("Please specify which parameters you want to use:\n")
    match = input_type_limiter("Specify the match value (has to be a number): ", float)
    mismatch = input_type_limiter("Specify the mismatch penalty (has to be a number): ", float)
    opening = input_type_limiter("Specify the indel penalty (has to be a number): ", float)
    extension = input_type_limiter("Specify the extension penalty (has to be a number): ", float)

    if mismatch > 0:
        print("Mismath was positive, I have changed to a negative number")
        mismatch = float("-" + str(mismatch))
    if opening > 0:
        print("Opening penalty was positive, I have changed to a negative number")
        opening = float("-" + str(opening))
    if extension > 0:
        print("Extension penalty was positive, I have changed to a negative number")
        extension = float("-" + str(extension))
    print()



    if alignment == "GLOBAL":

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nNeedleman-Wunsch method (nice guys by the way)\n")
        print("And the parameters are: match:{}, mismatch:{}, opening score:{} and extension score:{}\n".format(match, mismatch, opening, extension))
        
        (matrix, alignment) = alignment_nw(seq_list[0], seq_list[1], blosum, dna_prot, opening, extension, match, mismatch)
        
        print("\n")
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)


    elif alignment == "LOCAL":

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nSmith-Waterman method (nice guys by the way)\n")
        print("And the parameters are: match:{}, mismatch:{}, opening score:{} and extension score:{}\n".format(match, mismatch, opening, extension))

        (matrix, alignment, score) = alignment_sw(seq_list[0], seq_list[1], blosum, dna_prot, opening, extension, match, mismatch)
        
        print()
        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(alignment)
        print("The score of the best alignment is", score)
    


# Sequences are protein
elif is_protein == True:

    dna_prot = "protein"
    #Those variables are not used in proteins
    match = 0 
    mismatch = 0
    
    print("What you have given me is a protein sequence \n")
    while alignment not in ['GLOBAL', 'LOCAL']:
        alignment = input("Do you want to perform a local or a global alignment?\n").upper()

    print("Please tell me which parameters you want to use\n")
    opening = input_type_limiter("Specify the indel penalty (has to be a number): ", float)
    extension = input_type_limiter("Specify the extension penalty (has to be a number): ", float)

    if opening > 0:
        print("Opening penalty was positive, I have changed to a negative number")
        opening = float("-" + str(opening))
    if extension > 0:
        print("Extension penalty was positive, I have changed to a negative number")
        extension = float("-" + str(extension))
    
    print()
    print("The BLOSUM matrix options are: ")
    
    for i in range(len(blosum_files)):
        print("BLOSUM"+ str(blosum_files[i]), end = ", ")
        
    blosum_int = input_type_limiter("Specify the blosum matrix you want to use (introduce the number only): ", int)
    print()

    if blosum_int in blosum_files:
        blosum = blosum_int 

    if blosum_int not in blosum_files:
        print("I don't have that file, I am using BLOSUM62 instead\n")
        blosum = "62"



    if alignment == "GLOBAL" :

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("And the parameters are: opening score:{}, extension score:{} and blosum_matrix\n".format(opening, extension, blosum))
        
        (matrix, alignment) = alignment_nw(seq_list[0], seq_list[1], blosum, dna_prot, opening, extension, match, mismatch)
        
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print("\n")
        
        print(alignment)
        

    elif alignment == "LOCAL":

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("And the parameters are: opening score:{}, extension score:{} and blosum_matrix\n".format(opening, extension, blosum))
        
        (matrix, alignment, score) = alignment_sw(seq_list[0], seq_list[1], blosum, dna_prot, opening, extension, match, mismatch)
        
        print("\n")
        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        
        print(alignment)
        print("The score of the best alignment is", score)
