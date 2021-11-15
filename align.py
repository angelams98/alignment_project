#!/usr/bin/env python3

#Import libraries
import sys


#This function reads the file and generates a list with the sequences and the names
def read_files(infile):

    #Initialize the variables
    seq = []
    title = []
    sequences = ""

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
            if sequences != "":
                seq.append(sequences)
            sequences = ""
        else:
            sequences += line[:-1] 

    return seq, title




#Function for creating the BLOSUM matrix (dicts of dicts) from a file  
def blosum_matrix(file):

    import re

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


        if first_row != None:
            flag = True
        if flag:
            blosum.append(line[:-1].split())

        

    for row in range(len(blosum)):
        blosum_dict[blosum[row][0]] = {}
        for col in range(len(blosum[0])) :
            blosum_dict[blosum[row][0]][blosum[0][col]] = blosum[row][col] 




#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_simple(string1, string2, match, mismatch, opening, exten):

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
                else:
                    value_left = matrix[row][col-1] + opening
                    value_top = matrix[row-1][col] + opening

            
                #We compare the nucleotides in the strings
                if string1[row-1] == string2[col-1]:
                    value_diag  = matrix[row-1][col-1] + match

                else:
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
    return matrix_t, matrix_moves




#Function for protein sequences
#It calculate the alignment scores using blosum
def alignment_blosum(string1, string2, blosum_file, indel):

    string1.upper()
    string2.upper()

    #Initialize the variables 
    matrix = []

    nrow = len(string1) + 1
    ncol = len(string2) + 1

    #Create an empty matrix
    for row in range(nrow):
        matrix.append([])
        for col in range(ncol):
            matrix[row].append(0)

    #Fill out the matrix
    for row in range(nrow):
        for col in range(ncol):
            #In the first row we only calculate the values using the values from the left, so we start in position 1
            if row == 0 and col !=0:
                #Indel -1
                matrix[row][col] = matrix[row][col-1] + indel

            #In the first column we only calculate the values using the values from the top, so we start in row 1
            elif row != 0 and col == 0:
                #Indel -1
                matrix[row][col] = matrix[row-1][col] + indel
            
            #When not in the first row and column, the values can be calculated from left, top or diagonal
            elif row != 0 and col != 0:
                #Indel -1
                value_left = matrix[row][col-1] + indel
                value_top = matrix[row-1][col] + indel

                try:
                    value_diag  = matrix[row-1][col-1] + blosum62[string1[row-1]][string2[col-1]]

                except KeyError:
                    value_diag  = matrix[row-1][col-1] + blosum62[string2[col-1]][string1[row-1]]
        
                #The correct values is going to be the maximum value from the 3 we have calculated above   
                matrix[row][col] = max(value_left, value_top, value_diag)




#Function to print the matrix on screen
def print_matrix(matrix):
    for row in range(len(matrix)):
        printlist = []
        for column in range(len(matrix[row])):
            printlist.append(str(matrix[row][column]))
        print("\t".join(printlist))
    print("\n")  




#This function calculates the best alignment using the matching scores from the matrix
def traceback (matrix, seq1, seq2):

    #Initialize variables
    row = len(matrix) - 1
    col = len(matrix[0]) - 1

    score = matrix[row][col]

    align1 = seq1[col-1]
    align2 = seq2[row-1]
    #align_middle = "|"
    total_alignment = ""


    while row > 0 and col > 0:
    
        diag_score = matrix[row-1][col-1]
        left_score = matrix[row][col-1]
        top_score = matrix[row-1][col]

        if diag_score > left_score or diag_score > top_score:
            align1 = seq1[col-1] + align1
            align2 = seq2[row-1] + align2
            #align_middle = "|" + align_middle

            row = row - 1
            col = col - 1
            score += diag_score


        if left_score > diag_score or left_score > top_score:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2
            #align_middle = " " + align_middle
            
            col = col - 1
            score += left_score


        if top_score > diag_score or top_score > left_score:
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


    for i in range(0, len(align1), 60):
        total_alignment += align1[i:i+60] + "\n" + align2[i:i+60] + "\n" + "\n"

    
    return total_alignment


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




(seq, title) = read_files(infile)
matrix2 = alignment_simple(seq[0], seq[1], match, mismatch, indel)
print(traceback(matrix2, seq[0], seq[1]))


matrix = alignment_simple("GCATGCG", "GATTACA", 1, -1, -1)
print_matrix(matrix)
print(traceback(matrix, "GCATGCG", "GATTACA"))

blosum_matrix("BLOSUM62.txt")

