#!/usr/bin/env python3

#Function for creating a Blosum matrix out of a file
def blosum_matrix(file):
    blosum = []
    blosum_dict = dict()

    try:
        infile = open(file, 'r')
    except IOError:
        print("An error ocurred")

    for line in infile:
        line = "".join(line.split())
        print(line)
        if line.startswith("  "):
            blosum.append(line[:-1].split("  "))
        elif line[0:1] != "#":
            blosum.append(line[:-1].split())

    
    dictionary = {blosum[row][0]:{blosum[0][col]:blosum[row][col]} for col in range(len(blosum[0])) for row in range(len(blosum)) }

                

blosum_matrix("BLOSUM62.txt")