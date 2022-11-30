


import numpy as np
import os,sys
import configparser
import argparse

from matrix import *
from utils import *


parser = argparse.ArgumentParser(description='test')
parser.add_argument('--m', type=str, choices=['LU','GS','H','G','URV'],default='LU')
parser.add_argument('--input', type=str, default='../data/LU.txt')
parser.add_argument('--solve',type=int, default=0)
parser.add_argument('--vector',type=str, default='../data/V_LU.txt')

args = parser.parse_args()



mat = read_matrix(args.input)
if(mat.size <= 0):
    print("Input is not a matrix\n")
    sys.exit(1)
print("Input Matrix:")
print(mat)
#print(np.linalg.det(mat))
#print(np.linalg.qr(mat))
if (args.m == 'LU'):
    print("Your choice is : LU factorization")
    m,n = mat.shape
    if (m != n):
        print("Your input matrix should be a square matrix\n")
        sys.exit(1)
    P,L,U,det = LU_Fact(mat)
    print("P:")
    print(P)
    print("L:")
    print(L)
    print("U:")
    print(U)
    print("det：")
    print(det)
    if(args.solve == 1):
        b = read_vector(args.vector)
        b = P.dot(b)
        y = Solver_L(L,b)
        x = Solver_U(U,y)
        print("The reslut of Ax = b is: ")
        print(x)

if (args.m == 'GS'):
    print("Your choice is : Gram-Schmidt QR factorization")
    rank = np.linalg.matrix_rank(mat)
    m,n = mat.shape
    if(rank !=m and rank !=n):#Gram-Schmidt得不到正定矩阵
        print("This matrix cannot be QR factorized by Gram-Schmidt")
    Q,R,det = GramSchmidt(mat)
    print("Q:")
    print(Q)
    print("R:")
    print(R)
    if(det == "Error"):
        print("Your matrix must be square to calculate determinant\n")
    else:
        print("det：")
        print(det)
    if(args.solve == 1):
        b = read_vector(args.vector)
        b = Q.T.dot(b)
        x = Solver_U(R,b)
        print("The reslut of Ax = b is: ")
        print(x)

elif (args.m == 'H'):
    print("Your choice is : Householder QR factorization")
    m,n = mat.shape
    Q,R,det = Householder(mat)
    print("Q:")
    print(Q)
    print("R:")
    print(R)
    if(det == "Error"):
        print("Your matrix must be square to calculate determinant\n")
    else:
        print("det：")
        print(det)
    if(args.solve == 1):
        b = read_vector(args.vector)
        b = Q.T.dot(b)
        x = Solver_U(R,b)
        print("The reslut of Ax = b is: ")
        print(x)

elif (args.m == 'G'):
    print("Your choice is : Givens QR factorization")
    m,n = mat.shape
    Q,R,det = Givens(mat)
    print("Q:")
    print(Q)
    print("R:")
    print(R)
    if(det == "Error"):
        print("Your matrix must be square to calculate determinant\n")
    else:
        print("det：")
        print(det)
    if(args.solve == 1):
        b = read_vector(args.vector)
        b = Q.T.dot(b)
        x = Solver_U(R,b)
        print("The reslut of Ax = b is: ")
        print(x)

elif (args.m == 'URV'):
    print("Your choice is : URV factorization")
    m,n = mat.shape
    U,R,V,det = URV(mat)
    print("U:")
    print(U)
    print("R:")
    print(R)
    print("V: ")
    print(V)
    if(det == "Error"):
        print("Your matrix must be square to calculate determinant\n")
    else:
        print("det：")
        print(det)
    if(args.solve == 1):
        b = read_vector(args.vector)
        b = U.T.dot(b)
        R = R.dot(V.T)
        x = Solver_U(R,b)
        print("The reslut of Ax = b is: ")
        print(x)

