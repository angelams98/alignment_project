#!/usr/bin/env python3

import sys

#Read the file
if len(sys.argv) == 6:
    try:
        infile = sys.argv[1]
        match = sys.argv[2]
        mismatch = sys.argv[3]
        indel = sys.argv[4]
        extension = sys.argv[5]
    except ValueError:
        print("One or more values weren't correct")
        sys.exit(1)


else:
    infile = input("Give the name of the infile: ")
    match = input("Give me the match value: ")
    mismatch = input("Give me the mismatch value: ")
    indel = input("Give me the indel value: ")
    extension = input("Give me the extension value: ")


#Initialize the variables
seq_list = []
title = []
sequences = ""
nucleotides = ["A", "T", "G", "C"]
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

dna_flag_def = True
protein_flag_def = True

pos = 1

#Try to open the file and generates an error message if it fails
try:
    infile = open(infile, 'r')
except IOError as err:
    print("can't open file, reason:", str(err))
    sys.exit(1)

#Saves the titles and the sequences in lists
for line in infile:
    if line.startswith('>'):
        title.append(line[:-1])
        #If in the next line after ">", sequences is not empty, it is added to seq_list
        if sequences != "":
            seq_list.append(sequences)
        sequences = ""

    else:
        sequences += line[:-1] 

#If after finishing the loop, there is still something in sequences, it is added to seq_list        
if sequences != "":
    seq_list.append(sequences)

#Create an only string with the sequence
seq = "".join(seq_list)

#Assign the initial value for dna_flag and protein_flag according to the first position
if seq[0] in nucleotides:
    dna_flag = True

if seq[0] not in nucleotides:
    dna_flag = False
    dna_flag_def = False

if seq[0] in amino_acids:
    protein_flag = True

if seq[0] not in amino_acids:
    protein_flag = False
    protein_flag_def = False
    

#If none of the flags is True, it skips this part
while pos < len(seq) and (dna_flag == True or protein_flag == True):
    if seq[pos] in nucleotides:
        dna_flag = True

    #When it finds an element that is not a nucleotide, the definitive dna_flag is set to False
    if seq[pos] not in nucleotides:
        dna_flag_def = False

    if seq[pos] in amino_acids:
        protein_flag = True
    
    #When it finds an element that is not an amino acid, the definitive protein_flag is set to False
    if seq[pos] not in amino_acids:
        protein_flag_def = False

    pos += 1


#If the definitive dna_flag hasn't changed in the sequence, it is still True and it is a DNA
#This condition goes first because protein contains the same letters, so protein_flag_def must be True too.
if dna_flag_def == True:
    print("It's DNA")
#If the definitive protein_flag hasn't changed in the sequence, it is still True and it is a protein   
elif protein_flag_def == True:
    print("It's a protein")
#If none of them is True, the file doesn't contain a proper sequence. 
else:
    print("The file doesn't contain a sequence")

