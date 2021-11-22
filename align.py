#!/usr/bin/env python3

#Import libraries
import sys
import re


#Function for creating the BLOSUM matrix (dicts of dicts) from a file  
def blosum_matrix(file):
    """It takes a file with a blosum matrix and saves it as a dictionary of dictionaries"""

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
def alignment_dna_nw(string1, string2, match, mismatch, opening, exten):
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
                if matrix_moves[row][col-1] != "diag":
                    value_left = matrix[row][col-1] + exten
                if matrix_moves[row-1][col] != "diag":
                    value_top = matrix[row-1][col] + exten
                if matrix_moves[row][col-1] == "diag":
                    value_left = matrix[row][col-1] + exten
                if matrix_moves[row-1][col] == "diag":
                    value_top = matrix[row-1][col] + exten
            
                #We compare the nucleotides in the strings
                if string1[row-1] == string2[col-1]:
                    value_diag  = matrix[row-1][col-1] + match

                if string1[row-1] != string2[col-1]:
                    value_diag  = matrix[row-1][col-1] + mismatch


                #The correct values is going to be the maximum value from the 3 we have calculated above  
                matrix[row][col] = max(value_diag, value_left, value_top)
                list_values = [value_left, value_top, value_diag]
                #print(list_values.index(max(list_values)))


                if list_values.index(max(list_values)) == 0:
                    matrix_moves[row][col] = "diag"
                else:
                    matrix_moves[row][col] = "gap"

                
    matrix_t = [[matrix[col][row] for col in range(len(matrix))] for row in range(len(matrix[0]))]
    return matrix_t




#Function for protein sequences
#It calculate the alignment scores using blosum

##ADD EXTENSION AND OPENING
def alignment_protein_nw(string1, string2, dna_prot, opening, exten):
    """Calculates a scoring matrix for two protein sequences using Needleman-Wunsch's algorithm"""

    string1.upper()
    string2.upper()

    #Initialize the variables 
    matrix = []
    matrix_moves = []

    #print(string1, string2)

    nrow = len(string1) + 1
    ncol = len(string2) + 1

    #Create an empty matrix
    for row in range(nrow):
        matrix.append([])
        matrix_moves.append([])
        for col in range(ncol):
            matrix[row].append(0)
            matrix_moves[row].append(None)

    #print(matrix)

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
                if matrix_moves[row][col-1] != "diag":
                    value_left = matrix[row][col-1] + exten

                if matrix_moves[row-1][col] != "diag":
                    value_top = matrix[row-1][col] + exten

                if matrix_moves[row][col-1] == "diag":
                    value_left = matrix[row][col-1] + opening

                if matrix_moves[row-1][col] == "diag":
                    value_top = matrix[row-1][col] + opening

                value_diag  = matrix[row-1][col-1] + int(blosum[string1[row-1]][string2[col-1]])



                #The correct values is going to be the maximum value from the 3 we have calculated above  
                matrix[row][col] = max(value_diag, value_left, value_top)
                list_values = [value_left, value_top, value_diag]
                #print(list_values.index(max(list_values)))


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
    #print(len(matrix))
    #print(matrix[row-1][col-1])
    #print(seq1)
    #print(seq2)

    score = matrix[row][col]

    align1 = seq1[col]
    align2 = seq2[row]
    #print(align1, align2)
    #align_middle = "|"
    total_alignment = ""


    while row > 0 and col > 0:
    
        diag_score = matrix[row-1][col-1]
        left_score = matrix[row][col-1]
        top_score = matrix[row-1][col]
        values = [diag_score, left_score, top_score]


        if values.index(max(values)) == 0:
            align1 = seq1[col-1] + align1
            align2 = seq2[row-1] + align2
            #align_middle = "|" + align_middle
            #print(row, col)
            row = row - 1
            col = col - 1
            score += diag_score


        if values.index(max(values)) == 1:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2
            #align_middle = " " + align_middle
            
            col = col - 1
            score += left_score


        if values.index(max(values)) == 2:
            align1 = "-" + align1
            align2 = seq2[row-1] + align2
            #align_middle = " " + align_middle   
            
            row = row - 1
            score += top_score
        


        else:
            align1 = seq1[col-1] + align1
            align2 = seq2[row-1] + align2
            #align_middle = "|" + align_middle

            row = row - 1
            col = col - 1
            #align1 = "-" + align1
            #align2 = seq2[row-1] + align2
            #align_middle = " " + align_middle   
            
            #row = row - 1
            #score += top_score

        if row == 0 and col == 1:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2
            #align_middle = "|" + align_middle
            col = col - 1
        
        if row == 1 and col == 0:
            align1 = "-" + align1
            align2 = seq2[row-1] + align2
            #align_middle = "|" + align_middle
            row = row  - 1


    for i in range(0, len(align1), 60):
        total_alignment += align1[i:i+60] + "\n" + align2[i:i+60] + "\n" + "\n"

    print(score)
    return total_alignment


#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_dna_sw(string1, string2, match, mismatch, opening, exten):
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

            if matrix_moves[row][col-1] != "diag":
                value_left = matrix[row][col-1] + exten
            if matrix_moves[row-1][col] != "diag":
                value_top = matrix[row-1][col] + exten
            if matrix_moves[row][col-1] == "diag":
                value_left = matrix[row][col-1] + exten
            if matrix_moves[row-1][col] == "diag":
                value_top = matrix[row-1][col] + exten
        
            #We compare the nucleotides in the strings
            if string1[row-1] == string2[col-1]:
                value_diag  = matrix[row-1][col-1] + match

            if string1[row-1] != string2[col-1]:
                value_diag  = matrix[row-1][col-1] + mismatch


            #The correct values is going to be the maximum value from the 3 we have calculated above  
            
            if max(value_diag, value_left, value_top)<0:
                matrix[row][col] = 0
            if max(value_diag, value_left, value_top) >= 0:
                matrix[row][col] = max(value_diag, value_left, value_top)


            list_values = [value_left, value_top, value_diag]
            #print(list_values.index(max(list_values)))


            if list_values.index(max(list_values)) == 0:
                matrix_moves[row][col] = "diag"
            else:
                matrix_moves[row][col] = "gap"

                
    matrix_t = [[matrix[col][row] for col in range(len(matrix))] for row in range(len(matrix[0]))]
    return matrix_t




#Function for protein sequences
#It calculate the alignment scores using blosum

##ADD EXTENSION AND OPENING
def alignment_protein_sw(string1, string2, dna_prot, opening, exten):
    """Calculates a scoring matrix for two protein sequences using Smith-Waterman's algorithm"""

    string1.upper()
    string2.upper()

    #Initialize the variables 
    matrix = []
    matrix_moves = []

    #print(string1, string2)

    nrow = len(string1) + 1
    ncol = len(string2) + 1

    #Create an empty matrix
    for row in range(nrow):
        matrix.append([])
        matrix_moves.append([])
        for col in range(ncol):
            matrix[row].append(0)
            matrix_moves[row].append(None)

    #print(matrix)

    #Fill out the matrix
    for row in range(1, nrow):
        for col in range(1, ncol):

            if matrix_moves[row][col-1] != "diag":
                value_left = matrix[row][col-1] + exten

            if matrix_moves[row-1][col] != "diag":
                value_top = matrix[row-1][col] + exten

            if matrix_moves[row][col-1] == "diag":
                value_left = matrix[row][col-1] + opening

            if matrix_moves[row-1][col] == "diag":
                value_top = matrix[row-1][col] + opening

            value_diag  = matrix[row-1][col-1] + int(blosum[string1[row-1]][string2[col-1]])



            #The correct values is going to be the maximum value from the 3 we have calculated above  
            if max(value_diag, value_left, value_top)<0:
                matrix[row][col] = 0
            if max(value_diag, value_left, value_top) >= 0:
                matrix[row][col] = max(value_diag, value_left, value_top)

            list_values = [value_left, value_top, value_diag]
            #print(list_values.index(max(list_values)))


            if list_values.index(max(list_values)) == 0:
                matrix_moves[row][col] = "diag"
            else:
                matrix_moves[row][col] = "gap"
    
    matrix_t = [[matrix[col][row] for col in range(len(matrix))] for row in range(len(matrix[0]))]
    return matrix_t

def traceback_ws (matrix, seq1, seq2):
    #print(matrix, seq1, seq2)
    #print(len(seq1), len(seq2))

    matrix_max_value = -1
    matrix_max_position_list = []
    total_align1 = []
    total_align2 = []
    total_alignment_scores = []
    
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):

            if matrix[row][col] == matrix_max_value:
                matrix_max_position_list.append([row, col])

            if matrix[row][col] > matrix_max_value:
                matrix_max_value = matrix[row][col]
                
                matrix_max_position_list = []
                matrix_max_position_list.append([row, col])
            
            

    #print(matrix_max_position_list)
    for i in range(len(matrix_max_position_list)):
        row = matrix_max_position_list[i][0]-1
        col = matrix_max_position_list[i][1]-1
        #print(col, len(seq1))
        align1 = seq1[col]
        align2 = seq2[row]
        #print(align1, align2)
        #print(row, col)
        

        diag_score = -1
        total_alignment = ""

        while diag_score != 0:

            score = matrix[row][col]
            diag_score = matrix[row-1][col-1]
            left_score = matrix[row][col-1]
            top_score = matrix[row-1][col]
            values = [diag_score, left_score, top_score]
                
                    
            # Diagonal is highest
            if values.index(max(values)) == 0:
                #print(align1)
                align1 = seq1[col-1] + align1
                align2 = seq2[row-1] + align2
                #align_middle = "|" + align_middle
                #print(row, col)
                #print(align1 +"\n"+ align2)
                row = row - 1
                col = col - 1
                score += diag_score
                
            # left score is highest
            if values.index(max(values)) == 1:
                align1 = seq1[col-1] + align1
                align2 = "-" + align2
                #align_middle = " " + align_middle
                #print(row, col)
                #print(align1 +"\n"+ align2)
                col = col - 1
                score += left_score
            
            # top score is highest
            if values.index(max(values)) == 2:
                align1 = "-" + align1
                align2 = seq2[row-1] + align2
                #align_middle = " " + align_middle   
                #print(row, col)
                #print(align1 +"\n"+ align2)
                row = row - 1
                score += top_score
            
            
        

        total_align1.append(align1)
        total_align2.append(align2)
        total_alignment_scores.append(score)

        final_align1 = total_align1[total_alignment_scores.index(max(total_alignment_scores))]
        final_align2 = total_align2[total_alignment_scores.index(max(total_alignment_scores))]


    return final_align1, final_align2
    #print(align1[j:j+60] + "\n" + align2[j:j+60] + "\n" + "\n")
    #print(total_alignment_list)



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



#Read the file
if len(sys.argv) == 2:
    try:
        infile = sys.argv[1]

    except ValueError:
        print("One or more values weren't correct")
        sys.exit(1)


else:
    infile = input("Give the name of the infile: ")


#Try to open the file and generates an error message if it fails
try:
    infile = open(infile, 'r')
except IOError as err:
    print("can't open file, reason:", str(err))
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
        #print(is_protein)
    elif seq[pos] not in amino_acids:
        #print(is_protein)
        print("The sequence found in the file", infile.name, "contained an impure sequence.")
        sys.exit(1)
        
    pos = pos + 1


#Once we know if it is DNA or protein, we can decide which alignment methods should the program use

#If the definitive dna_flag hasn't changed in the sequence, it is still True and it is a DNA.
#This condition goes first because protein contains the same letters, so protein_flag_def must be True too.
if is_protein == False:
    match = int(input("Give me the match value: "))
    mismatch = int(input("Give me the mismatch value: "))
    indel = int(input("Give me the indel value: "))
    extension = int(input("Give me the extension value: "))

    print("The sequence was DNA")
    matrix1 = alignment_dna_sw("GCATGCG", "GATTACA", match, mismatch, indel, extension)
    print_matrix(matrix1)
    print(traceback_ws(matrix1, "GCATGCG", "GATTACA"))

    matrix2 = alignment_dna_sw(seq_list[0], seq_list[1], match, mismatch, indel, extension)
    #print(matrix2)
    #print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
    print(traceback_ws(matrix2, seq_list[0], seq_list[1]))

#If the definitive protein_flag hasn't changed in the sequence, it is still True and it is a protein  
elif is_protein == True:
    blosum_int = input("Give me which blosum matrix you want to use (introduce the number only):")
    indel = int(input("Give me the indel value: "))
    extension = int(input("Give me the extension value: "))
    print("The sequence was a protein.")
    #matrix1 = alignment_dna("GCATGCG", "GATTACA", match, mismatch, indel, extension)
    #print_matrix(matrix)
    #print(traceback(matrix1, "GCATGCG", "GATTACA"))

    if blosum_int == "62":
        blosum = blosum_matrix("BLOSUM62.txt")
    elif blosum_int == "80":
        blosum = blosum_matrix("BLOSUM80.txt")
    elif blosum_int == "50":
        blosum = blosum_matrix("BLOSUM50.txt")
    else:
        print("I don't have that file, I am using BLOSUM62 instead")
        blosum = blosum_matrix("BLOSUM62.txt")

    
    matrix3 = alignment_protein_sw(seq_list[0], seq_list[1], blosum, indel, extension)
    #print(matrix3)
    
    print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
    
    (output1, output2) = traceback_ws(matrix3, seq_list[0], seq_list[1])

    for i in range(0, len(output1), 60):
        print(output1[i:i+60] +"\n" + output2[i:i+60] +"\n" +"\n")


