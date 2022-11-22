import numpy as np
import os,sys

def read_matrix(file):
    f = open(file, 'r',encoding='utf-8')
    lines = f.readlines()
    matrix = []
    m = len(lines)
    l = lines[0].strip('\n').split(' ')
    n = len(l)
    i= 0    
    matrix = np.zeros((m, n), dtype=float)
    for line in lines:
        l = line.strip('\n').split(' ')
        matrix[i: ] = l[0:n]
        i += 1
    f.close()
    return matrix

def read_vector(file):
    f = open(file, 'r',encoding='utf-8')
    lines = f.readlines()
    line = lines[0].strip('\n').split(' ')
    n = len(line)
    vector = np.zeros(n,dtype=float)
    vector[0:] = line[0:n]
    f.close()
    return(vector)