#!/usr/bin/env python3

#Import libraries
import sys
import re


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
                list_values = [value_left, value_top, value_diag]

                
                #Check from which cell we have calculated the score to save the movement
                #O(o), o is the number of elements it has to check
                if list_values.index(max(list_values)) == 2:
                    matrix_moves[row][col] = "diag"
                else:
                    matrix_moves[row][col] = "gap"


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




#This function calculates the best alignment using the matching scores from the matrix
def traceback_nw (matrix, seq1, seq2):
    """
        Perform a traceback through a scoring matrix using the Needleman-Wunsch algorithm to find the best global alignment
        PSEUDOCODE:
            
            Initialize row and col to the last cell of the matrix to start looping bottom to top
            
            While row > 0 or col > 0:
                Look at scores in cells to the left, above and diagonal
                
                if diagonal score is highest:
                    append letter from seq1 to align1
                    append letter from seq2 to align2
                    row = row-1
                    col = col-1
                
                if left score is highest:
                    append letter from seq1 to align1
                    introduce gap "-" to align2
                    col = col-1
                    
                if above score is highest:
                    introduce gap "-" to align1
                    append letter from seq2 to align2
                    row = row-1
                
            for i in range(length of alignment):
                save alignments as a string, 60 characters at a time separated by \n
                
        
        :matrix: A list of lists containing a scoring matrix
        :seq1: The first of the sequences to be aligned
        :seq2: The second of the sequences to be aligned
    
        :return: \n-formatted string with total alignment
    """

    #Initialize variables
    row = len(seq2) 
    col = len(seq1) 
    
    align1 = ""
    align2 = ""
    total_alignment = ""

    while row > 0 or col > 0:
    
        diag_score = matrix[row-1][col-1]
        left_score = matrix[row][col-1]
        top_score = matrix[row-1][col]
        values = [diag_score, left_score, top_score]

        if values.index(max(values)) == 0:
            align1 = seq1[col-1] + align1
            align2 = seq2[row-1] + align2

            row = row - 1
            col = col - 1

        if values.index(max(values)) == 1:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2

            col = col - 1

        if values.index(max(values)) == 2:
            align1 = "-" + align1
            align2 = seq2[row-1] + align2 

            row = row - 1
        
    for i in range(0, len(align1), 60):
        total_alignment += align1[i:i+60] + "\n" + align2[i:i+60] + "\n" + "\n"

    return total_alignment


#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_sw(string1, string2, dna_prot, opening, exten, match, mismatch):
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
    matrix = []
    matrix_moves = []

    nrow = len(string1) 
    ncol = len(string2) 

    string1 = "*" + string1
    string2 = "*" + string2

    #Create an empty matrix
    #O(m), m is the number of rows of the new matrix
    for row in range(nrow + 1):
        matrix.append([])
        matrix_moves.append([])

        #O(n), n isthe number of columns of the new matrix
        for col in range(ncol + 1):
            matrix[row].append(0)
            matrix_moves[row].append('diag')


    #Fill out the matrix
    #O(m), m is the number of rows of the matrix
    for row in range(nrow + 1):
        #O(n), n is the number of rows of the matrix
        for col in range(ncol + 1):

            if row != 0 and col != 0:

                if dna_prot == "dna":
                    #We compare the nucleotides in the strings
                    if string1[row-1] == string2[col-1]:
                        value_diag  = matrix[row-1][col-1] + match

                    if string1[row-1] != string2[col-1]:
                        value_diag  = matrix[row-1][col-1] + mismatch

                elif dna_prot == "protein":
                    #We calculate the score using the blosum matrix
                    value_diag  = matrix[row-1][col-1] + int(blosum[string1[row-1]][string2[col-1]])


                #We check if the previous value was a gap, so the value score is for gap opening or extension
                if matrix_moves[row][col-1] != "diag":
                    value_left = matrix[row][col-1] + exten

                if matrix_moves[row-1][col] != "diag":
                    value_top = matrix[row-1][col] + exten

                if matrix_moves[row][col-1] == "diag":
                    value_left = matrix[row][col-1] + opening

                if matrix_moves[row-1][col] == "diag":
                    value_top = matrix[row-1][col] + opening


                #The correct value is going to be the maximum from the 3 we have calculated above  
                #O(o), o is the number of elements it has to check
                if max(value_diag, value_left, value_top) >= 0:
                    matrix[row][col] = max(value_diag, value_left, value_top)

                #O(o), o is the number of elements it has to check
                if max(value_diag, value_left, value_top) < 0:
                    #In Smith-Waterman, the minimum value is always 0    
                    matrix[row][col] = 0

            
                list_values = [value_left, value_top, value_diag]

                #We store the cell from which we have calculated the score
                #O(o), o is the number of elements it has to check
                if list_values.index(max(list_values)) == 2:
                    matrix_moves[row][col] = "diag"

                else:
                    matrix_moves[row][col] = "gap"


    return matrix


def traceback_sw (matrix, seq1, seq2):
    """
        Perform a traceback through a scoring matrix using the Smith-Waterman algorithm to find the best local alignment
        PSEUDOCODE:
            for row in matrix:
                for col in matrix:
                    loop over all positions, save coordinates of highest score(s) encountered to a list
            
            for i in range(number of places where the highest score of the matrix can be found):
                row, col = coordinates saved in position i of high score list
                
                while a 0 is not encountered in diagonal:
                    total score = score for starting point i
                    
                    find max of left, above and diagonal scores
                    if diagonal max:
                        no gap added
                        row = row-1
                        col = col-1
                        total score = score + diagonal score
                    
                    if left max:
                        gap in align2
                        col = col-1
                        total score = score + left score
                        
                    if above max:
                        gap in align1
                        row = row-1
                        total score = score + above score
                
                Append align1 and align2 to list
            
            Final alignment is highest scoring (when there is more than one starting point)
        
        :matrix: A list of lists containing a scoring matrix
        :seq1: The first of the sequences to be aligned
        :seq2: The second of the sequences to be aligned
        
        :return: a tuple containing the aligned sequences
    """
    
    #Initialize the variables
    matrix_max_value = -1
    matrix_max_position_list = []
    total_align1 = []
    total_align2 = []
    total_alignment_scores = []
    

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

        row = matrix_max_position_list[i][0]-1
        col = matrix_max_position_list[i][1]-1

        align1 = seq1[col]
        align2 = seq2[row]
        
        diag_score = -1
        total_alignment = ""

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
                align1 = seq1[col-1] + align1
                align2 = seq2[row-1] + align2
            
                row = row - 1
                col = col - 1
                score += diag_score
                

            #Check if the left score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if values.index(max(values)) == 1:
                align1 = seq1[col-1] + align1
                align2 = "-" + align2
        
                col = col - 1
                score += left_score
            

            #Check if the top score is the highest
            #O(r), r is the number of elements we are comparing, only 3
            if values.index(max(values)) == 2:
                align1 = "-" + align1
                align2 = seq2[row-1] + align2
        
                row = row - 1
                score += top_score
            

        total_align1.append(align1)
        total_align2.append(align2)
        #O(1), it is appending elements to a list
        total_alignment_scores.append(score)

        #O(p), p is the number of elements in total_alignment_scores, the same as in matrix_max_position_list
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
        indel = int(input("Give me the indel value: "))
        extension = int(input("Give me the extension value: "))
        print("\n")

        
        matrix = alignment_nw(seq_list[0], seq_list[1], dna_prot, match, mismatch, indel, extension)
        print("Needleman-Wunsch alignment for:\n{}\n{}\n".format(title[0], title[1]))
        print(traceback_nw(matrix, seq_list[0], seq_list[1]))

    #O(q), q is the strings we are comparing in the conditional
    elif alignment == "LOCAL":

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
        
        #O(r), r is the length of the output1 sequence
        for i in range(0, len(output1), 60):
            print(output1[i:i+60] +"\n" + output2[i:i+60] +"\n" +"\n")
        print("The value for this alignment is", score)
        


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
        

    #O(q), q is the strings we are comparing in the conditional
    elif alignment == "LOCAL":

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

        #O(r), r is the length of the output1 sequence
        for i in range(0, len(output1), 60):
            print(output1[i:i+60] +"\n" + output2[i:i+60] +"\n" +"\n")
        print("The value for this alignment is", score)

