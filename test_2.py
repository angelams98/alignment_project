#!/usr/bin/env python3

#Import libraries
import sys
import re


#Function for creating the BLOSUM matrix (dicts of dicts) from a file  
def blosum_matrix(file):
    """It takes a file with a blosum matrix and saves it as a dictionary of dictionaries"""

    #Initialize variables
    blosum = []
    blosum_dict = dict()
    flag = False


    try:
        infile = open(file, 'r')
    except IOError:
        print("An error ocurred")

    for line in infile:
        line_temp = "".join(line.split())
        first_row = re.match(r'^[A-Z]{23}\**', line_temp)

        if flag :
            blosum.append(line.split())

        if first_row != None:
            line = "-" + line[:-1]
            blosum.append(line.split())
            flag = True

    for row in range(len(blosum)):
        blosum_dict[blosum[row][0]] = {}
        for col in range(len(blosum[0])) :
            if row == 0 or col == 0:
                pass
            else:
                blosum_dict[blosum[row][0]][blosum[0][col]] = int(blosum[row][col])

    return blosum_dict




#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_nw(string1, string2, dna_prot, opening, exten, match, mismatch):
    """Calculates a scoring matrix for two DNA sequences using Needleman-Wunsch's algorithm"""

    #Initialize the variables 
    matrix = []
    matrix_moves = []

    nrow = len(string1) + 1
    ncol = len(string2) + 1


    #Create an empty matrix
    for row in range(nrow):
        matrix.append([])
        matrix_moves.append([])

        for col in range(ncol):
            matrix[row].append(0)
            matrix_moves[row].append(None)


    #Fill out the matrix
    for row in range(nrow):
        for col in range(ncol):
            #In the first row we only calculate the values using the values from the left, so we start in position 1
            if row == 0 and col !=0:
                if matrix_moves[row][col-1] != "diag":
                    matrix[row][col] = matrix[row][col-1] + exten
                else:
                    matrix[row][col] = matrix[row][col-1] + opening
                matrix_moves[row][col] = "gap"

            #In the first column we only calculate the values using the values from the top, so we start in row 1
            elif row != 0 and col == 0:

                if matrix_moves[row-1][col] != "diag":
                    matrix[row][col] = matrix[row-1][col] + exten
                else:
                    matrix[row][col] = matrix[row-1][col] + opening
                matrix_moves[row][col] = "gap"
            
            #When not in the first row and column, the values can be calculated from left, top or diagonal
            elif row != 0 and col != 0:
                if dna_prot == "dna":
                    #We compare the nucleotides in the strings
                    if string1[row-1] == string2[col-1]:
                        value_diag  = matrix[row-1][col-1] + match

                    if string1[row-1] != string2[col-1]:
                        value_diag  = matrix[row-1][col-1] + mismatch

                elif dna_prot == "protein":
                    value_diag  = matrix[row-1][col-1] + int(blosum[string1[row-1]][string2[col-1]])

                if matrix_moves[row][col-1] != "diag":
                    value_left = matrix[row][col-1] + exten

                if matrix_moves[row-1][col] != "diag":
                    value_top = matrix[row-1][col] + exten

                if matrix_moves[row][col-1] == "diag":
                    value_left = matrix[row][col-1] + opening

                if matrix_moves[row-1][col] == "diag":
                    value_top = matrix[row-1][col] + opening
            

                #The correct values is going to be the maximum value from the 3 we have calculated above  
                matrix[row][col] = max(value_diag, value_left, value_top)
                list_values = [value_left, value_top, value_diag]

                
                #Check from which cell we have calculated the score to save the movement
                if list_values.index(max(list_values)) == 0:
                    matrix_moves[row][col] = "diag"
                else:
                    matrix_moves[row][col] = "gap"

                
    matrix_t = [[matrix[col][row] for col in range(len(matrix))] for row in range(len(matrix[0]))]
    return matrix_t



#Function to print the matrix on screen
def print_matrix(matrix):
    """Prints the matrix on screen in a clear way"""
    for row in range(len(matrix)):
        printlist = []
        for column in range(len(matrix[row])):
            printlist.append(str(matrix[row][column]))
        print("\t".join(printlist))
    print("\n")  




#This function calculates the best alignment using the matching scores from the matrix
def traceback_nw (matrix, seq1, seq2):
    """Find the best alignment for two sequences using a scoring matrix"""

    #Initialize variables
    row = len(matrix) - 2
    col = len(matrix[0]) - 2

    score = matrix[row][col]

    align1 = seq1[col]
    align2 = seq2[row]
    total_alignment = ""


    while row > 0 and col > 0:
    
        diag_score = matrix[row-1][col-1]
        left_score = matrix[row][col-1]
        top_score = matrix[row-1][col]
        values = [diag_score, left_score, top_score]


        #Check if the diagonal score is the highest
        if values.index(max(values)) == 0:
            align1 = seq1[col-1] + align1
            align2 = seq2[row-1] + align2

            row = row - 1
            col = col - 1

        #Check if the left score is the highest
        if values.index(max(values)) == 1:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2
            
            col = col - 1

        #Check if the top score is the highest
        if values.index(max(values)) == 2:
            align1 = "-" + align1
            align2 = seq2[row-1] + align2   
            
            row = row - 1


        #If we are in the first row and the second column, it can only go to the left
        if row == 0 and col == 1:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2

            col = col - 1
        
        #If we are in the second row and the first column, it can only go to the top
        if row == 1 and col == 0:
            align1 = "-" + align1
            align2 = seq2[row-1] + align2

            row = row  - 1


    #Save the sequence in lines of 60 characters
    for i in range(0, len(align1), 60):
        total_alignment += align1[i:i+60] + "\n" + align2[i:i+60] + "\n" + "\n"

    return total_alignment


#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_sw(string1, string2, dna_prot, opening, exten, match, mismatch):
    """Calculates a scoring matrix for two DNA sequences using Smith-Waterman's algorithm"""
    
    #Initialize the variables 
    matrix = []
    matrix_moves = []

    nrow = len(string1) + 1
    ncol = len(string2) + 1


    #Create an empty matrix
    for row in range(nrow):
        matrix.append([])
        matrix_moves.append([])

        for col in range(ncol):
            matrix[row].append(0)
            matrix_moves[row].append(None)


    #Fill out the matrix
    for row in range(1, nrow):
        for col in range(1, ncol):

            #We check if the previous value was a gap, so the value score is for gap opening or extension
            if matrix_moves[row][col-1] != "diag":
                value_left = matrix[row][col-1] + exten

            if matrix_moves[row-1][col] != "diag":
                value_top = matrix[row-1][col] + exten

            if matrix_moves[row][col-1] == "diag":
                value_left = matrix[row][col-1] + opening

            if matrix_moves[row-1][col] == "diag":
                value_top = matrix[row-1][col] + opening
        

            if dna_prot == "dna":
                #We compare the nucleotides in the strings
                if string1[row-1] == string2[col-1]:
                    value_diag  = matrix[row-1][col-1] + match

                if string1[row-1] != string2[col-1]:
                    value_diag  = matrix[row-1][col-1] + mismatch

            elif dna_prot == "protein":

                value_diag  = matrix[row-1][col-1] + int(blosum[string1[row-1]][string2[col-1]])



            #The correct value is going to be the maximum from the 3 we have calculated above  
            if max(value_diag, value_left, value_top) >= 0:
                matrix[row][col] = max(value_diag, value_left, value_top)

            if max(value_diag, value_left, value_top) < 0:
                #In Smith-Waterman, the minimum value is always 0    
                matrix[row][col] = 0

        
            list_values = [value_left, value_top, value_diag]

            #We store the cell from which we have calculated the score
            if list_values.index(max(list_values)) == 0:
                matrix_moves[row][col] = "diag"

            else:
                matrix_moves[row][col] = "gap"

                
    matrix_t = [[matrix[col][row] for col in range(len(matrix))] for row in range(len(matrix[0]))]
    return matrix_t





def traceback_sw (matrix, seq1, seq2):

    #Initialize the variables
    matrix_max_value = -1
    matrix_max_position_list = []
    total_align1 = []
    total_align2 = []
    total_alignment_scores = []
    

    for row in range(len(matrix)):
        for col in range(len(matrix[row])):

            #If we find the same value we already have, we save it to do the alignment for that position
            if matrix[row][col] == matrix_max_value:
                matrix_max_position_list.append([row, col])

            #If we find a higher value, we save the value and the position
            if matrix[row][col] > matrix_max_value:
                matrix_max_value = matrix[row][col]
                
                matrix_max_position_list = []
                matrix_max_position_list.append([row, col])
            
    
    #We do the alignments for all the maximum values we have saved
    for i in range(len(matrix_max_position_list)):

        row = matrix_max_position_list[i][0]-1
        col = matrix_max_position_list[i][1]-1

        align1 = seq1[col]
        align2 = seq2[row]
        
        diag_score = -1
        total_alignment = ""

        while diag_score != 0:

            score = matrix[row][col]
            diag_score = matrix[row-1][col-1]
            left_score = matrix[row][col-1]
            top_score = matrix[row-1][col]
            values = [diag_score, left_score, top_score]
                
                    
            #Check if the diagonal score is the highest
            if values.index(max(values)) == 0:
                align1 = seq1[col-1] + align1
                align2 = seq2[row-1] + align2
            
                row = row - 1
                col = col - 1
                score += diag_score
                

            #Check if the left score is the highest
            if values.index(max(values)) == 1:
                align1 = seq1[col-1] + align1
                align2 = "-" + align2
        
                col = col - 1
                score += left_score
            

            #Check if the top score is the highest
            if values.index(max(values)) == 2:
                align1 = "-" + align1
                align2 = seq2[row-1] + align2
        
                row = row - 1
                score += top_score
            

        total_align1.append(align1)
        total_align2.append(align2)
        total_alignment_scores.append(score)

        final_align1 = total_align1[total_alignment_scores.index(max(total_alignment_scores))]
        final_align2 = total_align2[total_alignment_scores.index(max(total_alignment_scores))]


    return final_align1, final_align2, total_alignment_scores



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
while pos < len(seq):
    
    if seq[pos] in amino_acids_checker:
        is_protein = True

    elif seq[pos] not in amino_acids:
        print("The sequence found in the file", infile.name, "contained an impure sequence.")
        sys.exit(1)
        
    pos = pos + 1


#Once we know if it is DNA or protein, we can decide which alignment methods should the program use

#They are DNA sequences
if is_protein == False:

    print("What you have given me is a DNA sequence \n")
    alignment = input("So do you want to do local or global alignment?\n")
    dna_prot = "dna"
    blosum = 0

    if alignment in ("local", "LOCAL", "Local"):

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        match = int(input("Give me the match value: "))
        mismatch = int(input("Give me the mismatch value: "))
        indel = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        print("\n")

        
        matrix = alignment_nw(seq_list[0], seq_list[1], dna_prot, match, mismatch, indel, extension)
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(traceback_nw(matrix, seq_list[0], seq_list[1]))


    elif alignment in ("global", "GLOBAL", "Global"):

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        match = int(input("Give me the match value: "))
        mismatch = int(input("Give me the mismatch value: "))
        indel = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        print("\n")

        matrix = alignment_sw(seq_list[0], seq_list[1], dna_prot, match, mismatch, indel, extension)
        (output1, output2, score) = traceback_sw(matrix, seq_list[0], seq_list[1])

        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        for i in range(0, len(output1), 60):
            print(output1[i:i+60] +"\n" + output2[i:i+60] +"\n" +"\n")
        print("The value for this alignment is", score)
        


#They are protein sequences
elif is_protein == True:

    print("What you have given me is a protein sequence \n")
    alignment = input("So do you want to do local or global alignment?\n")
    dna_prot = "protein"
    match = 0 
    mismatch = 0

    if alignment in ("local", "LOCAL", "Local"):

        print("Perfect, we are going to do a local alignment of your sequence, we will be applying \nNeedleman-Wunshman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        indel = int(input("Give me the indel value: "))
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

        
        matrix = alignment_nw(seq_list[0], seq_list[1], dna_prot, indel, extension, match, mismatch)

        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(traceback_nw(matrix, seq_list[0], seq_list[1]))
        


    elif alignment in ("global", "GLOBAL", "Global"):

        print("Perfect, we are going to do a global alignment of your sequence, we will be applying \nSmith-Waterman method (nice guy by the way)\n")
        print("Please tell me which parameters you want to use\n")

        indel = int(input("Give me the indel value: "))
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

    
        matrix = alignment_sw(seq_list[0], seq_list[1], dna_prot, match, mismatch, indel, extension)

        print("Smith-Waterman alignment for:\n{}\n{}\n".format(title[0], title[1]))
        (output1, output2, score) = traceback_sw(matrix, seq_list[0], seq_list[1])

        for i in range(0, len(output1), 60):
            print(output1[i:i+60] +"\n" + output2[i:i+60] +"\n" +"\n")
        print("The value for this alignment is", score)


