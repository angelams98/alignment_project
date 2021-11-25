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




#This function creates a list of lists with the matching scores
def alignment_nw(string1, string2, dna_prot, opening, exten, match, mismatch):
    """
    Calculates a scoring matrix for two DNA/protein sequences using the Needleman-Wunsch algorithm
    PSEUDOCODE:
        Create two matrices (list of lists):
            Dimensions: Rows = length of string1 + 1, columns = length of string2 + 1
            1. Matrix for scores
            2. matrix_moves for keeping track of whether a gap was introduced previously or not
        
        For row in range(length of string1+1):
            For col in range(length of string2+1):
                
                If currently in the first row:
                    check matrix_moves[cell to the left] to determine whether new gap is introduced or one is being continued
                    matrix[row][col] = score from left cell + gap/extend penalty
                    
                If currently in the first column:
                    check matrix_moves[cell abive] to determine whether new gap is introduced or one is being continued
                    matrix[row][col] = score from above cell + gap/extend penalty
                
                If elsewhere:
                    if DNA
                        check if the sequences match at this position
                        if yes:
                            score for diagonal path = diagonal cell score + match
                        if no:
                            score for diagonal path = diagonal cell score + mismatch
                    if protein
                        score for diagonal path = look up score in blosum matrix
                        
                    check matrix_moves[cells above and to the left] whether taking value from top or left would introduce or extend a gap
                        calculate scores for left and top paths
                        
                    value of score matrix[row][col] = max of left/top/diagonal paths
                    value of matrix_moves[row][col]:
                        'gap' if max value came from top or left
                        'diag' if max value came from diagonal

        Transpose matrix

    :string1: The first DNA/protein sequence to be aligned as a string
    :string2: The other DNA/protein sequence to be aligned as a string
    :dna_prot: A string which can be either 'DNA' or 'protein' which defines whether the sequences are protein or DNA
    :match: A score for when the sequences match
    :mismatch: A penalty for mismatching positions of the sequences
    :opening: A penalty for opening a gap in the alignment
    :exten: A penalty for extending a gap in the alignment

    :returns: A matrix of containing the scores of all possible alignments
    """

    #Initialize the variables 
    matrix = []
    matrix_moves = []

    ncol = len(string1)
    nrow = len(string2) 

    string1 = "*" + string1
    string2 = "*" + string2


    #Create an empty matrix
    #O(m), m is the number of rows in the new matrix
    for row in range(nrow + 1):
        matrix.append([])
        matrix_moves.append([])

        #O(n), n is the number of columns in the new matrix
        for col in range(ncol + 1):
            matrix[row].append(0)
            matrix_moves[row].append('diag')


    #Fill out the matrix
    #O(m), m is the number of rows in the matrix
    for row in range(nrow + 1):

        #O(n), n is the number of columns in the matrix
        for col in range(ncol + 1):
            #In the first row we only calculate the values using the values from the left, so we start in position 1
            if row == 0 and col !=0:
                if matrix_moves[row][col-1] != "diag":
                    matrix[row][col] = matrix[row][col-1] + exten
                else:
                    matrix[row][col] = matrix[row][col-1] + opening
                matrix_moves[row][col] = "left"

            #In the first column we only calculate the values using the values from the top, so we start in row 1
            elif row != 0 and col == 0:

                if matrix_moves[row-1][col] != "diag":
                    matrix[row][col] = matrix[row-1][col] + exten
                else:
                    matrix[row][col] = matrix[row-1][col] + opening
                matrix_moves[row][col] = "top"
            
            #When not in the first row and column, the values can be calculated from left, top or diagonal
            elif row != 0 and col != 0:
                if dna_prot == "dna":
                    #We compare the nucleotides in the strings
                    if string1[col] == string2[row]:
                        value_diag  = matrix[row-1][col-1] + match

                    if string1[col] != string2[row]:
                        value_diag  = matrix[row-1][col-1] + mismatch

                elif dna_prot == "protein":
                    value_diag  = matrix[row-1][col-1] + int(blosum[string1[col]][string2[row]])

                if matrix_moves[row][col-1] != "diag":
                    value_left = matrix[row][col-1] + exten

                if matrix_moves[row-1][col] != "diag":
                    value_top = matrix[row-1][col] + exten

                if matrix_moves[row][col-1] == "diag":
                    value_left = matrix[row][col-1] + opening

                if matrix_moves[row-1][col] == "diag":
                    value_top = matrix[row-1][col] + opening
            

                #The correct values is going to be the maximum value from the 3 we have calculated above  
                #O(o), o is the number of elements it has to check
                matrix[row][col] = max(value_left, value_top, value_diag)
                list_values = [value_diag, value_left, value_top]

                
                #Check from which cell we have calculated the score to save the movement
                #O(o), o is the number of elements it has to check
                if list_values.index(max(list_values)) == 0:
                    matrix_moves[row][col] = "diag"
                if list_values.index(max(list_values)) == 1:
                    matrix_moves[row][col] = "left"
                if list_values.index(max(list_values)) == 2:
                    matrix_moves[row][col] = "top"


    print_matrix(matrix_moves)
    print_matrix(matrix)
    return matrix



#Function to print the matrix on screen
def print_matrix(matrix):
    """
        Prints a list of lists matrix on screen in a formatted way
        
        :matrix: A list of lists matrix     
    """
    #O(m), m is the number of rows in the matrix
    for row in range(len(matrix)):
        printlist = []
        #O(n), n is the number of columns in the matrix
        for column in range(len(matrix[row])):
            printlist.append(str(matrix[row][column]))
        print("\t".join(printlist))
    print("\n")  



#Function to print the matrix on screen
def print_matrix(matrix):
    """Prints the matrix on screen in a clear way"""
    for row in range(len(matrix)):
        printlist = []
        for column in range(len(matrix[row])):
            printlist.append(str(matrix[row][column]))
        print("\t".join(printlist))
    print("\n")  



str1 = "QIKDLT"
str2 = "QITDKDVLLV"
blosum = blosum_matrix("BLOSUM62.txt")
result = alignment_nw(str1, str2, "protein", 0, 0, -10, -1)

print_matrix(result)
