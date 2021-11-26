#!/usr/bin/env python3

#Import libraries
import sys
import re

def make_matrix(string1, string2, filler):
    #Initialize the variables 
    matrix = []
    
    
    
    #print(len(string1), len(string2))
    #Create an empty matrix
    #O(m), m is the number of rows in the new matrix
    for row in range(len(string2)+1):
        matrix.append([])
    
        #O(n), n is the number of columns in the new matrix
        for col in range(len(string1)+1):
            matrix[row].append(filler)

    return matrix

#Function for creating the BLOSUM matrix (dicts of dicts) from a file  
def blosum_matrix(file):
    """
        Creates a BLOSUM matrix formatted as a dictionary of dictionaries.
        PSEUDOCODE:
            Open file
        
            # Saving the BLOSUM content of the file as a list of lists (= a "matrix")
            For line in file:
                Join line into a single string
                Regex match first line of BLOSUM matrix
                    
                    If regex match found
                        Save line to "matrix"
                        Set flag to true
                    
                    If flag is true
                        Save line to "matrix"
            
            # Creating dictionary of dictionaries
            For row in range(number of rows in "matrix"):
                For col in range(number of columns in "matrix"):
                    
                    outer key = matrix[row][0] (list of amino acids)
                    value for outer key = inner key 
                    inner key = matrix[0][col] (list of amino acods)
                    value for inner key = matrix[row][col]
        
        :file: A .txt file from NCBI which contains a BLOSUM matrix
        
        :returns: A dictionary of dictionaries
    """

    #Initialize variables
    blosum = []
    blosum_dict = dict()
    flag = False


    try:
        #O(1)
        infile = open(file, 'r')
    except IOError:
        print("An error ocurred")

    #O(k), k is the number of lines in the file
    for line in infile:
        line_temp = "".join(line.split())
        first_row = re.match(r'^[A-Z]{23}\**', line_temp)

        if flag :
            blosum.append(line.split())

        if first_row != None:
            line = "-" + line[:-1]
            blosum.append(line.split())
            flag = True


    #O(n), n is the number of rows
    for row in range(len(blosum)):
        blosum_dict[blosum[row][0]] = {}
        #O(m), m is the number of cols
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
    matrix = make_matrix(string1, string2, 0)
    matrix_moves = make_matrix(string1, string2, 'diag')

    ncol = len(string1)
    nrow = len(string2)
    
    string1 = "*" + string1
    string2 = "*" + string2

    align1 = ""
    align2 = ""
    middle_space = ""
    total_alignment = ""


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

    print_matrix(matrix)

    while nrow > 0 or ncol > 0:

        if matrix_moves[nrow][ncol] == "diag":
            nrow -= 1
            ncol -= 1
        
            align1 = string1[ncol + 1] + align1
            align2 = string2[nrow + 1] + align2

            #Identity
            if string1[ncol + 1] == string2[nrow + 1]:
                middle_space = "|" + middle_space

            #blosum score greater than 0
            elif int(blosum[string1[ncol + 1]][string2[nrow + 1]]) > 0:
                middle_space = ":" + middle_space

            #Blosum score 0 or negative
            elif int(blosum[string1[ncol + 1]][string2[nrow+1]]) <= 0:
                middle_space = "." + middle_space

            #Gap
            else:
                middle_space = " " + middle_space

            
        if matrix_moves[nrow][ncol] == "top":
            nrow -= 1
            
            align1 = "-" + align1
            align2 = string2[nrow + 1] + align2
            middle_space = " " + middle_space


        if matrix_moves[nrow][ncol] == "left":
            ncol -= 1

            align1 = string1[ncol + 1] + align1
            align2 = "-" + align2
            middle_space = " " + middle_space



    for i in range(0, len(align1), 60):
        total_alignment += align1[i:i+60] + "\n" + middle_space[i:i+60] + "\n"+ align2[i:i+60] + "\n" + "\n"

    return matrix, total_alignment



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



#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_sw(string1, string2, dna_prot, opening, exten, match, mismatch):
    print(string1, string2, dna_prot, open, exten, match, mismatch)
    print(opening, exten)
    """
        Calculates a scoring matrix for two DNA/protein sequences using the Smith-Waterman algorithm
        PSEUDOCODE:
            Create two matrices (list of lists):
                Dimensions: Rows = length of string1 + 1, columns = length of string2 + 1
                1. Matrix for scores
                2. matrix_moves for keeping track of whether a gap was introduced previously or not
            
            For row in range(1, length of string1+1):
                For col in range(1, length of string2+1):
                    
                    check matrix_moves[cells above and to the left] whether taking value from top or left would introduce or extend a gap
                    calculate scores for left and top paths
                    
                    if DNA:
                        check if the sequences match at this position
                        if yes:
                            score for diagonal path = diagonal cell score + match
                        if no:
                            score for diagonal path = diagonal cell score + mismatch
                    
                    if protein:
                        score for diagonal path = look up score from blosum matrix
                        
                    find max score of left, above or diagonal. 
                    If this is below 0, matrix[row][col] = 0
                    If this is => 0, matrix[row][col] = max score

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
    matrix = make_matrix(string1, string2, 0)
    matrix_moves = make_matrix(string1, string2, "zero")
    ncol = len(string1)
    nrow = len(string2)
    
    string1 = "*" + string1
    string2 = "*" + string2


    #Fill out the matrix
    #O(m), m is the number of rows of the matrix
    for row in range(1, nrow + 1):
        #O(n), n is the number of rows of the matrix
        for col in range(1, ncol + 1):

            if dna_prot == "dna":
                #We compare the nucleotides in the strings
                if string1[col] == string2[row]:
                    value_diag  = matrix[row-1][col-1] + match

                if string1[col] != string2[row]:
                    value_diag  = matrix[row-1][col-1] + mismatch

            elif dna_prot == "protein":
                #We calculate the score using the blosum matrix
                value_diag  = matrix[row-1][col-1] + int(blosum[string1[col]][string2[row]])


            #We check if the previous value was a gap, so the value score is for gap opening or extension
            if matrix_moves[row][col-1] == "gap":
                value_left = matrix[row][col-1] + exten

            if matrix_moves[row-1][col] == "gap":
                value_top = matrix[row-1][col] + exten

            if matrix_moves[row][col-1] != "gap":
                value_left = matrix[row][col-1] + opening

            if matrix_moves[row-1][col] != "gap":
                value_top = matrix[row-1][col] + opening


            #The correct value is going to be the maximum from the 3 we have calculated above  
            #O(o), o is the number of elements it has to check

            zero = 0

            matrix[row][col] = max(value_diag, value_left, value_top, zero)
 
        
            list_values = [value_diag, value_left, value_top, zero]

            #We store the cell from which we have calculated the score
            #O(o), o is the number of elements it has to check
            if list_values.index(max(list_values)) == 0:
                matrix_moves[row][col] = "diag"

            elif list_values.index(max(list_values)) == 3:
                matrix_moves[row][col] = "zero"

            else:
                matrix_moves[row][col] = "gap"


    #Initialize the variables
    matrix_max_value = -1
    matrix_max_position_list = []
    total_align1 = []
    total_align2 = []
    total_alignment_scores = []
    final_alignment = ""
    

    #O(m), m is the number of rows in the matrix
    for row in range(len(matrix)):
        #O(n), n is the number of columns in the matrix
        for col in range(len(matrix[row])):

            #If we find the same value we already have, we save it to do the alignment for that position
            if matrix[row][col] == matrix_max_value:

                #O(1), it is appending elements into a list
                matrix_max_position_list.append([row, col])

            #If we find a higher value, we save the value and the position
            if matrix[row][col] > matrix_max_value:
                matrix_max_value = matrix[row][col]
                matrix_max_position_list = []

                #O(1), it is appending elements to a list
                matrix_max_position_list.append([row, col])
            
    
    #We do the alignments for all the maximum values we have saved
    #O(p), p is the number of elements we have saved in matrix_max_position_list
    for i in range(len(matrix_max_position_list)):

        row = matrix_max_position_list[i][0] 
        col = matrix_max_position_list[i][1] 

        align1 = ""
        align2 = ""
        
        diag_score = -1
        total_alignment = ""
        middle_space = ""


                #O(p*q), p is the number of rows we are checking, q is the number of columns we are checking
        while diag_score != 0:

            score = matrix[row][col]
            diag_score = matrix[row-1][col-1]
            left_score = matrix[row][col-1]
            top_score = matrix[row-1][col]
            values = [diag_score, left_score, top_score]

            #print(align1)
            #print(align2)
                
                    
            #Check if the diagonal score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if values.index(max(values)) == 0:
                row = row - 1
                col = col - 1
                score += diag_score

                align1 = string1[col + 1] + align1
                align2 = string2[row + 1] + align2

                

                if string1[col + 1] == string2[row + 1]:
                    middle_space = "|" + middle_space

                #blosum score greater than 0
                elif int(blosum[string1[col+1]][string2[row+1]]) > 0:
                    middle_space = ":" + middle_space

                #Blosum score 0 or negative
                elif int(blosum[string1[col+1]][string2[row+1]]) <= 0:
                    middle_space = "." + middle_space

                #Gap
                else:
                    middle_space = " " + middle_space
            
                

            #Check if the left score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if values.index(max(values)) == 1:
                col = col - 1
                score += left_score

                align1 = string1[col + 1] + align1
                align2 = "-" + align2
                middle_space = " " + middle_space
        

            #Check if the top score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if values.index(max(values)) == 2:
                row = row - 1
                score += top_score

                align1 = "-" + align1
                align2 = string2[row] + align2
                middle_space = " " + middle_space


        total_align1.append(align1)
        total_align2.append(align2)
        

        #O(1), it is appending elements to a list
        total_alignment_scores.append(score)

        #O(p), p is the number of elements in total_alignment_scores, the same as in matrix_max_position_list
        final_align1 = total_align1[total_alignment_scores.index(max(total_alignment_scores))]
        final_align2 = total_align2[total_alignment_scores.index(max(total_alignment_scores))]
    #print(final_align1, final_align2)

    for i in range(0, len(final_align1), 60):
        final_alignment += final_align1[i:i+60] + "\n" + middle_space[i:i+60] + "\n" + final_align2[i:i+60] + "\n"  

    return matrix, final_alignment, total_alignment_scores

