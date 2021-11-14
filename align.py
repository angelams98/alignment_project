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



#Function for DNA samples
#This function creates a list of lists with the matching scores
def alignment_simple(string1, string2, match, mismatch, indel):

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
            
                #We compare the nucleotides in the strings
                if string1[row-1] == string2[col-1]:
                    #Match +1
                    value_diag  = matrix[row-1][col-1] + match

                else:
                    #Mismatch -1
                    value_diag  = matrix[row-1][col-1] + mismatch


                #The correct values is going to be the maximum value from the 3 we have calculated above  
                matrix[row][col] = max(value_left, value_top, value_diag) 
                
    matrix_t = [[matrix[col][row] for col in range(len(matrix))] for row in range(len(matrix[0]))]
    return matrix_t



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
    #print(row)
    col = len(matrix[0]) - 1
    #print(col)

    score = matrix[row][col]

    align1 = seq1[col-1]
    align2 = seq2[row-1]
    align_middle = "|"
    total_alignment = ""


    while row > 0 and col > 0:
    
        diag_score = matrix[row-1][col-1]
        left_score = matrix[row][col-1]
        top_score = matrix[row-1][col]

        if diag_score > left_score or diag_score > top_score:
            align1 = seq1[col-1] + align1
            align2 = seq2[row-1] + align2
            align_middle = "|" + align_middle

            row = row - 1
            col = col - 1
            score += diag_score


        elif left_score > diag_score or left_score > top_score:
            align1 = seq1[col-1] + align1
            align2 = "-" + align2
            align_middle = " " + align_middle
            
            col = col - 1
            score += left_score


        elif top_score > diag_score or top_score > left_score:
            align1 = "-" + align1
            align2 = seq2[row-1] + align2
            align_middle = " " + align_middle
            
            row = row - 1
            score += top_score

    for i in range(0, len(align_middle), 60):
        total_alignment += align1[i:i+60] + "\n" + align_middle[i:i+60] + "\n" + align2[i:i+60] + "\n"
    
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

#Read the file
if len(sys.argv) != 2:
    infile = input("Give the name of the infile: ")
else:
    infile = sys.argv[1]



(seq, title) = read_files(infile)
matrix2 = alignment_simple(seq[0], seq[1], 1, -1, -1)
#print(matrix2)
print(traceback(matrix2, seq[0], seq[1]))



#matrix = alignment_simple("GCATGCG", "GATTACA", 1, -1, -1)
#print_matrix(matrix)
#print(traceback(matrix, "GCATGCG", "GATTACA"))
#print_matrix(matrix2)
#print(traceback(matrix2, seq[0], seq[1]))




"""

blosum62 = {"C": {"C": 9, "S": -1, "T": -1,"P": -3, "A": 0, "G": -3,"N": -3, "D": -3, 
                "E": -4, "Q":  -3, "H": -3, "R":-3, "K": -3, "M": -1, "I": -1, "L": -1, 
                "V": -1, "F": -2, "Y": -2, "W": -2},
            "s":{"S": 4, "T": 1, "P": -1, "A": 1, "G": 0, "N": 1, "D": 0, "E": 0, "Q": 0, 
                "H": -1, "R": -1, "K": 0, "M": -1, "I": -2, "L": -2, "V": -2, "F": -2,
                "Y": -2, "W": -3},
            "T":{"T": 5, "P": -1, "A": 0, "G": -2, "N": 0, "D": -1, "E": -1, "Q": -1, 
                "H": -2, "R": -1, "K": -1, "M": -1, "I": -1, "L": -1, "V": 0, "F": -2,
                "Y": -2, "W": -2},
            "P":{"P": 7, "A": -1, "G": -2, "N": -2, "D": -1, "E": -1, "Q": -1, "H": -2, 
                "R": -2, "K": -1, "M": -2, "I": -3, "L": -3, "V": -2, "F": -4, "Y": -3, 
                "W": -4},
            "A":{"A": 4, "G": 0, "N": -2, "D": -2, "E": -1, "Q": -1, "H": -2, "R": -1, 
                "K": -1, "M": -1, "I": -1, "L": -1, "V": 0, "F": -2, "Y": -2, "W": -3},
            "G":{"G": 6, "N": 0, "D": -1, "E": -2, "Q": -2, "H": -2, "R": -2, "K": -2, 
                "M": -3, "I": -4, "L": -4, "V": -3, "F":-3, "Y": -3, "W": -2},
            "N":{"N": 6, "D": 1, "E": 0, "Q": 0, "H": 1, "R": 0, "K": 0, "M": -2, 
                "I":-3, "L": -3, "V": -3, "F": -3, "Y": -2, "W": -4},
            "D":{"D": 6, "E": 2, "Q": 0, "H": -1, "R": -2, "K": -1, "M": -3, "I": -3,
                "L": -4, "V": -3, "F": -3, "Y": -3, "W": -4},
            "E":{"E": 5, "Q": 2, "H": 0, "R": 0, "K": 1, "M": -2, "I": -3, "L": -3, 
                "V": -2, "F": -3, "Y": -2, "W": -3},
            "Q":{"Q": 5, "H": 0, "R": 1, "K": 1, "M": 0, "I": -3, "L": -2, "V": -2, 
                "F": -3, "Y": -1, "W": -2},
            "H":{"H": 8, "R": 0, "K": -1, "M": -2, "I": -3, "L": -3,  "V": -3, "F": -1,
                "Y": 2, "W": -2},
            "R":{"R": 5, "K": 2, "M": -1, "I": -3, "L": -2, "V": -3, "F": -3,"Y": -2, 
                "W": -3},
            "K":{"K": 5, "M": -1, "I": -3, "L": -2, "V": -2, "F": -3, "Y": -2, "W": -3},
            "M":{"M": 5, "I": 1, "L": 2, "V": 1, "F": 0, "Y": -1, "W": -1},
            "I":{"I": 4, "L": 2, "V": 3, "F": 0, "Y": -1, "W": -3},
            "L":{"L": 4, "V": 1, "F": 0,"Y": -1, "W": -2},
            "V":{"V": 4, "F": -1, "Y": -1, "W": -3},
            "F":{"F": 6, "Y": 3, "W": 1},
            "Y":{"Y": 7, "W": 2},
            "Y":{"W": 11}}


"""



            











#matrix = alignment_blosum("TNDEQHA", "HNPEVMA", -8)
#print_matrix(matrix)

#Function for creating a Blosum matrix out of a file
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

        
    #print(blosum)
    dictionary = {blosum[row][0]:{blosum[0][col]:blosum[row][col]} for col in range(len(blosum[0])) for row in range(len(blosum)) }
    #print(dictionary)
                

#blosum_matrix("BLOSUM62.txt")