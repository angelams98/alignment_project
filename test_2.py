#!/usr/bin/env python3

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


#Function to print the matrix on screen
def print_matrix(matrix):
    for row in range(len(matrix)):
        printlist = []
        for column in range(len(matrix[row])):
            printlist.append(str(matrix[row][column]))
        print("\t".join(printlist))
    print("\n")  

(matrix, matrix_moves) = alignment_simple("GCATGCG", "GATTACA", 1, -1, -5, -1)
print_matrix(matrix_moves)
print_matrix(matrix)