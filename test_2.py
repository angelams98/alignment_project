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


        """
        #O(p*q), p is the number of rows we are checking, q is the number of columns we are checking
        while diag_score != 0:

            score = matrix[row][col]
            diag_score = matrix[row-1][col-1]
                                
            #Check if the diagonal score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if matrix_moves[row][col] == "diag":

                row -= 1
                col -= 1

                align1 = string1[col + 1] + align1
                align2 = string2[row + 1] + align2

                #Identity
                if string1[col + 1] == string2[row + 1]:
                    middle_space = "|" + middle_space
                
                #Blosum score greater than 0
                elif int(blosum[string1[col+1]][string2[row+1]]) > 0:
                    middle_space = ":" + middle_space

                #Blosum score 0 or less than 0
                elif int(blosum[string1[col + 1]][string2[row + 1]]) <= 0:
                    middle_space = "." + middle_space
                
                #Gaps
                else:
                    middle_space = " " + middle_space
                

            #Check if the left score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if matrix_moves[row][col] == "top":
                row -= 1
                align1 = "-" + align1
                align2 = string2[row + 1] + align2
                middle_space = " " + middle_space

            if matrix_moves[row][col] == "left":
                col -= 1
                align1 = string1[col + 1] + align1
                align2 = "-" + align2
                middle_space = " " + middle_space

            
        """

def alignment_sw1(string1, string2, blosum_number, dna_prot, opening, exten, match, mismatch):
    """
        Calculates a scoring matrix for two DNA/protein sequences using the Smith-Waterman algorithm

        PSEUDOCODE:
            Initialize two list of lists matrices, one for keeping track of scores (MATRIX) 
            and one for keeping track of which cell scores came from (MATRIX_MOVES).
            Create a BLOSUM matrix if sequences are protein
            
            # Creating score and moves matrices
            Double-loop over rows and columns of matrices
                Look at the scores of the cells to let, top and top-left of the current cell. 
                Calculate scores each of these would give, then select highest one to record to MATRIX.
                If highest score is a negative number, set score to 0.
                
                Record to MATRIX_MOVES whether the score was positive and came from top-left 
                cell (introduce match in alignment), positive and came from the cell above or to 
                the left (introduce gap in alignment), or if the score was negative/zero (end 
                alignment if encountered).
            
            # Finding highest scores in MATRIX
            Double-loop over rows and columns of MATRIX
                If current cell has same score as previous highest score
                    Append position to list
                If current cell has higher score than previous highest score
                    Set as new highest score
                    Clear list
                    Append position to list
            
            # Traceback
            For i highest scoring positions, the traceback has to run i times:
                Start in first of the highest scoring positions
                While a zero is not in the cell top-left of the current cell:
                    
                    Get scores from cells to the left, top and top-left of current cell
                    
                    If top-left score highest:
                        Position is a match, assign amino acids / nucleotides for this position to alignment
                        Move up and to the left
                    
                    If left score highest:
                        A gap is opened in sequence 1 of the alignment
                        Move up
                    
                    If above score highest:
                        A gap is opened in sequence 2 of the alignment
                        Move left
            
            Format and return the output


        :string1: The first DNA/protein sequence to be aligned as a string
        :string2: The other DNA/protein sequence to be aligned as a string
        :blosum_number: User-specified BLOSUM matrix to use for scoring amino acid substitutions
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

    # Call blosum_matrix function to create the user-requested BLOSUM substitution matrix
    if dna_prot == "protein":
        blosum = blosum_matrix(str("BLOSUM" + str(blosum_number) + ".txt"))
    
    #Fill out the matrix
    for row in range(1, nrow + 1):
        for col in range(1, ncol + 1):

            if dna_prot == "dna":
                #We compare the nucleotides in the strings
                if string1[col] == string2[row]:
                    value_diag  = matrix[row - 1][col - 1] + match

                if string1[col] != string2[row]:
                    value_diag  = matrix[row - 1][col - 1] + mismatch

            elif dna_prot == "protein":
                #We calculate the score using the blosum matrix
                value_diag  = matrix[row - 1][col - 1] + int(blosum[string1[col]][string2[row]])

            #We check if the previous value was a gap, so the value score is for gap opening or extension
            if matrix_moves[row][col - 1] == "gap":
                value_left = matrix[row][col - 1] + exten

            if matrix_moves[row - 1][col] == "gap":
                value_top = matrix[row - 1][col] + exten

            if matrix_moves[row][col - 1] != "gap":
                value_left = matrix[row][col - 1] + opening

            if matrix_moves[row - 1][col] != "gap":
                value_top = matrix[row - 1][col] + opening

            #The correct value is going to be the maximum from the 3 we have calculated above or 0
            list_values = [value_diag, value_left, value_top, 0]
            matrix[row][col] = max(list_values)
 

            #We store the cell from which we have calculated the score
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
    total_align_middle = []
    total_alignment_scores = []
    final_alignment = ""


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

        row = matrix_max_position_list[i][0] 
        col = matrix_max_position_list[i][1] 

        align1 = ""
        align2 = ""
        
        diag_score = -1
        middle_space = ""


        while diag_score != 0:

            score = matrix[row][col]
            diag_score = matrix[row-1][col-1]
            left_score = matrix[row][col-1]
            top_score = matrix[row-1][col]
            values = [diag_score, left_score, top_score]
                    
            #Check if the diagonal score is the highest
            if values.index(max(values)) == 0:
                row = row - 1
                col = col - 1
                score += diag_score

                align1 = string1[col + 1] + align1
                align2 = string2[row + 1] + align2
                
                if dna_prot == "protein":

                    if string1[col + 1] == string2[row + 1]:
                        middle_space = "|" + middle_space
    
                    #blosum score greater than 0
                    elif int(blosum[string1[col + 1]][string2[row + 1]]) > 0:
                        middle_space = ":" + middle_space
    
                    #Blosum score 0 or negative
                    elif int(blosum[string1[col + 1]][string2[row + 1]]) <= 0:
                        middle_space = "." + middle_space
    
                    #Gap
                    else:
                        middle_space = " " + middle_space
                        
                if dna_prot == "dna":
                    if string1[col + 1] == string2[row + 1]:
                        middle_space = "|" + middle_space
                    
                    else:
                        middle_space = " " + middle_space
                    
            #Check if the left score is the highest
            if values.index(max(values)) == 1:
                col = col - 1
                score += left_score

                align1 = string1[col + 1] + align1
                align2 = "-" + align2
                middle_space = " " + middle_space

            #Check if the top score is the highest
            if values.index(max(values)) == 2:
                row = row - 1
                score += top_score

                align1 = "-" + align1
                align2 = string2[row + 1] + align2
                middle_space = " " + middle_space


        total_align1.append(align1)
        total_align2.append(align2)
        total_align_middle.append(middle_space)
        total_alignment_scores.append(score)


        final_align1 = total_align1[total_alignment_scores.index(max(total_alignment_scores))]
        final_align2 = total_align2[total_alignment_scores.index(max(total_alignment_scores))]
        final_align_middle = total_align_middle[total_alignment_scores.index(max(total_alignment_scores))]

    
    for i in range(0, len(final_align1), 60):
        final_alignment += final_align1[i:i+60] + "\n" + final_align_middle[i:i+60] + "\n" + final_align2[i:i+60] + "\n"  + "\n"

    return matrix, final_alignment, max(total_alignment_scores)


def alignment_nw1(string1, string2, blosum_number, dna_prot, opening, exten, match, mismatch):
    """
    Calculates a scoring matrix for two DNA/protein sequences using the Needleman-Wunsch algorithm,
    followed by traceback to return the optimal global alignment.
    
    PSEUDOCODE:
        Initialize two list of lists matrices, one for keeping track of scores (MATRIX) 
        and one for keeping track of where the score was calculated from (MATRIX_MOVES)
        Create a BLOSUM matrix if sequences are protein
        
        # Creating score and moves matrices
        Double-loop over rows and columns of matrices
            Look at the scores of the cells to let, top and top-left of the current cell. 
            Calculate scores each of these would give, then select highest one to record to MATRIX.
            Record which direction the current score came from to MATRIX_MOVES
            
        # Traceback to get final alignment
        Start in the last cell of MATRIX_MOVES
        While not in the first row or column of MATRIX_moves:
            If MATRIX_MOVES = 'diag':
                Position is a match, assign amino acids / nucleotides for this position to alignment
                Move up and to the left
            If MATRIX_MOVES = 'top':
                A gap is opened in sequence 1 of the alignment
                Move up
            If MATRIX_MOVES = 'left':
                A gap is opened in sequence 2 of the alignment
                Move left
                
        Format and return the output
        

    :string1: A string containing the first DNA/protein sequence to be aligned
    :string2: A string containing the other DNA/protein sequence to be aligned
    :blosum_number: User-specified BLOSUM matrix to use for scoring amino acid substitutions
    :dna_prot: A string which can be either 'DNA' or 'protein' to define whether the sequences are protein or DNA
    :match: A score for when the sequences match
    :mismatch: A penalty for mismatching positions of the sequences
    :opening: A penalty for opening a gap in the alignment
    :exten: A penalty for extending a gap in the alignment

    :returns: The optimal global alignment of two sequences
    """

    #Initialize variables 
    ncol = len(string1)
    nrow = len(string2)
    string1 = "*" + string1
    string2 = "*" + string2
    align1 = ""
    align2 = ""
    middle_space = ""
    global_alignment = ""
    
    # Initialize score and movement matrices
    matrix = make_matrix(string1, string2, 0)
    matrix_moves = make_matrix(string1, string2, 'diag')
    
    # Call blosum_matrix function to create the user-requested BLOSUM substitution matrix
    if dna_prot == "protein":
        blosum = blosum_matrix(str("BLOSUM" + str(blosum_number) + ".txt"))

    # Create matrix of scores and matrix of moves
    for row in range(nrow + 1):

        for col in range(ncol + 1):
            
            #In the first row we only calculate the values using the values from the left, so we start in position 1
            if row == 0 and col !=0:
                if matrix_moves[row][col - 1] != "diag":
                    matrix[row][col] = matrix[row][col - 1] + exten
                else:
                    matrix[row][col] = matrix[row][col - 1] + opening
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
                
                # For DNA look at whether the position is a match or a mismatch
                if dna_prot == "dna":
                    
                    if string1[col] == string2[row]:
                        value_diag  = matrix[row-1][col-1] + match

                    if string1[col] != string2[row]:
                        value_diag  = matrix[row-1][col-1] + mismatch
                
                # For proteins look at the blosum matrix to see if a gap is opened in the alignment or not,
                # then check substituion matrix (BLOSUM) for score in case of not a gap
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
                list_values = [value_diag, value_left, value_top]
                matrix[row][col] = max(list_values)
                

                #Check from which cell we have calculated the score to save the movement
                if list_values.index(max(list_values)) == 0:
                    matrix_moves[row][col] = "diag"
                if list_values.index(max(list_values)) == 1:
                    matrix_moves[row][col] = "left"
                if list_values.index(max(list_values)) == 2:
                    matrix_moves[row][col] = "top"


    # If we reach last position from the end for each string (first position contains "*"), program stops
    while nrow > 1 or ncol > 1:
        
        # Match between sequences
        #nrow and ncol are the length of the strings and matrix is bigger
        if matrix_moves[nrow + 1][ncol + 1] == "diag":
            nrow -= 1
            ncol -= 1
        
            align1 = string1[ncol + 1] + align1
            align2 = string2[nrow + 1] + align2
        
        # Lines below creates a middle line that shows the similarity of positions in the output
            if dna_prot == "protein":
                # Check whether sequences are the same in this position
                if string1[ncol + 1] == string2[nrow + 1]:
                    middle_space = "|" + middle_space
    
                # If sequences are not the same in this position, check the score to determine how great the similarity is
                # If score above 0 = high similarity
                elif int(blosum[string1[ncol + 1]][string2[nrow + 1]]) > 0:
                    middle_space = ":" + middle_space
    
                # If score is below 0 = low similarity
                elif int(blosum[string1[ncol + 1]][string2[nrow + 1]]) <= 0:
                    middle_space = "." + middle_space
    
                # If there is a gap
                else:
                    middle_space = " " + middle_space
                    
            elif dna_prot == "dna":
                if string1[ncol + 1] == string2[nrow + 1]:
                    middle_space = "|" + middle_space
                else:
                    middle_space = " " + middle_space
                    
        # Gap in sequence
        if matrix_moves[nrow + 1][ncol + 1] == "top":
            nrow -= 1
            
            align1 = "-" + align1
            align2 = string2[nrow + 1] + align2
            middle_space = " " + middle_space

        # Gap in sequence
        if matrix_moves[nrow + 1][ncol + 1] == "left":
            ncol -= 1

            align1 = string1[ncol + 1] + align1
            align2 = "-" + align2
            middle_space = " " + middle_space

    # Format output into a single string
    for i in range(0, len(align1), 60):
        global_alignment += align1[i : i + 60] + "\n" + middle_space[i : i + 60] + "\n"+ align2[i : i + 60] + "\n" + "\n"
    
    
    return matrix, global_alignment