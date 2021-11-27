#!/usr/bin/env python3

#Import libraries
import re



def input_type_limiter(string, input_type):
    """
        A function which limits what type of input a user is allowed to input
        
        PSEUDOCODE:
            While loop until user inputs something which can be returned
                try:
                    Return what the user input as the required data type
                except:
                    If the user input cannot be converted to the required data type, print error and keep looping
                    
        :string: Message to be provided with user-input request
        :input_type: The datatype required
        
        :return: User-input as the required data type
    """
    
    type_dict = {float: "number", int: "number", str: "text string"}
    
    while True:
        try:
            return input_type(input(string))
        except:
            print("The given input could not be returned as a", type_dict[input_type] + ".", " Please try again.")
            pass


def make_matrix(string1, string2, filler):
    """
        Creates a list of lists matrix with dimensions n*m, where n the length of
        string2+1 and m is the length of string1+1.
        The matrix has the same value in all entrie, and the value is defined by the user.
        
        PSEUDOCODE:
            Make empty list
            Double for-loop over length of input strings.
                For each row, append a new list to make list of lists
                For each column, fill with filler values.
        
            Return output
        
    
        :string1: Text string which determines number of columns
        :string2: Text string which determines the number of rows
        :filler: User-defined value with which to fill out the matrix
        
        :return: A list of lists matrix
    """
    
    #Initialize the list
    matrix = []

    #Create an empty matrix
    #O(n), n is the number of rows in the new matrix
    for row in range(len(string2)+1):
        matrix.append([])
    
        #O(m), m is the number of columns in the new matrix
        for col in range(len(string1)+1):
            matrix[row].append(filler)

    return matrix



 
def blosum_matrix(blosum_file):
    """
        Creates a BLOSUM matrix as a dictionary of dictionaries with values for amino acid substitutions
        
        PSEUDOCODE:
            Open file
        
            # Saving the BLOSUM content of the file as a list of lists matrix
            For line in file:
                Join lines into a single string
                Regex match first line of BLOSUM matrix
                    If regex match found, save line and all following lines to list of lists matrix (BLOSUM_LIST)
            
            # Creating dictionary of dictionaries
            Double-loop over rows and columns of BLOSUM_LIST
                Outer keys = first entire row of BLOSUM_LIST
                Values for outer key = inner keys
                Inner keys = first entire column of BLOSUM_LIST
                Values for inner keys = score from BLOSUM_LIST[amino_acid1][amino_acid2]
        
            Return the output
        
        
        :blosum_file: A .txt file (downloaded from NCBI) which contains a BLOSUM matrix
        
        :return: A dictionary of dictionaries
    """

    #Initialize variables
    blosum_list = []
    blosum_dict = dict()
    flag = False

    # Open file
    try:
        infile = open(blosum_file, 'r')
    except IOError:
        print("An error ocurred")

    # Create list of lists from lines in infile containing BLOSUM matrix
    #O(k), k is the number of lines in the file
    for line in infile:
        line_temp = "".join(line.split())
        first_row = re.match(r'^[A-Z]{23}\**', line_temp)

        if flag:
            blosum_list.append(line.split())

        if first_row != None:
            line = "-" + line[:-1]
            blosum_list.append(line.split())
            flag = True

    # Create a dictionary of dictionaries
    #O(n), n is the number of rows
    for row in range(len(blosum_list)):
        blosum_dict[blosum_list[row][0]] = {}
        #O(m), m is the number of cols
        for col in range(len(blosum_list[0])) :
            if row == 0 or col == 0:
                pass
            else:
                blosum_dict[blosum_list[row][0]][blosum_list[0][col]] = int(blosum_list[row][col])

    return blosum_dict




def alignment_nw(string1, string2, blosum_number, dna_prot, opening, exten, match, mismatch):
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
        blosum = blosum_matrix(str("BLOSUM"+str(blosum_number)+".txt"))

    # Create matrix of scores and matrix of moves
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
                elif int(blosum[string1[ncol + 1]][string2[nrow+1]]) <= 0:
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
        global_alignment += align1[i:i+60] + "\n" + middle_space[i:i+60] + "\n"+ align2[i:i+60] + "\n" + "\n"
    
    
    return matrix, global_alignment




def print_matrix(matrix):
    """
        Prints a list of lists matrix on screen in a formatted way
        
        PSEUDOCODE:
            For each row of list of lists matrix:
                Make/reset printlist
                    Loop over columns of the given row    
                        Fill printlist with the contents of each row
                Print the printlist as a tab-separated string
            
            
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




def alignment_sw(string1, string2, blosum_number, dna_prot, opening, exten, match, mismatch):
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
        blosum = blosum_matrix(str("BLOSUM"+str(blosum_number)+".txt"))
    
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

            #The correct value is going to be the maximum from the 3 we have calculated above or 0
            #O(o), o is the number of elements it has to check
            matrix[row][col] = max(value_diag, value_left, value_top, 0)
 
            list_values = [value_diag, value_left, value_top, 0]

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
    total_align_middle = []
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
        middle_space = ""

        #O(p*q), p is the number of rows we are checking, q is the number of columns we are checking
        while diag_score != 0:

            score = matrix[row][col]
            diag_score = matrix[row-1][col-1]
            left_score = matrix[row][col-1]
            top_score = matrix[row-1][col]
            values = [diag_score, left_score, top_score]
                    
            #Check if the diagonal score is the highest
            #O(r), r is the number of elements we are comparing, only 3
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
                        #print(align1)
                        #print(middle_space)
                        #print(align2)
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
                #MODIFIED
                align2 = string2[row + 1] + align2
                middle_space = " " + middle_space

        total_align1.append(align1)
        total_align2.append(align2)
        total_align_middle.append(middle_space)

        #O(1), it is appending elements to a list
        total_alignment_scores.append(score)

        #O(p), p is the number of elements in total_alignment_scores, the same as in matrix_max_position_list
        final_align1 = total_align1[total_alignment_scores.index(max(total_alignment_scores))]
        final_align2 = total_align2[total_alignment_scores.index(max(total_alignment_scores))]
        final_align_middle = total_align_middle[total_alignment_scores.index(max(total_alignment_scores))]

    for i in range(0, len(final_align1), 60):
        final_alignment += final_align1[i:i+60] + "\n" + final_align_middle[i:i+60] + "\n" + final_align2[i:i+60] + "\n"  + "\n"

    return matrix, final_alignment, max(total_alignment_scores)
