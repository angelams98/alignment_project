#!/usr/bin/env python3

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

        

    for row in range(len(blosum)):
        blosum_dict[blosum[row][0]] = {}
        for col in range(len(blosum[0])) :
            blosum_dict[blosum[row][0]][blosum[0][col]] = blosum[row][col] 


blosum_matrix("BLOSUM62.txt")