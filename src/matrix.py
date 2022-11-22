import numpy as np
import os,sys

#solve det(Q) in QR factorization
def detQ(Q):
    n = Q.size()
    I = np.eye(n)
    r = np.linalg.matrix_rank(-I-Q)#to judge det(Q)
    n = np.linalg.matrix_rank(Q)
    if ((n-r) % 2 == 1):
        det = -1
    else:
        det = 1 
    return det

#solve det(R) in QR factorization
def detR(R):
    n = R.size()
    det_R = 1
    for i in range(n):
        det_R = det_R * R[i,i]
    return det_R

#solve Ax=y when A is a upper triangular matrix
def Solver_U(A,y):
    m,n = A.shape
    if m!=n:#consistency
        print("This system does not have a certain solution.\n")
        x = []
        return x
    x = np.zeros(n)
    for i in range(n-1,-1,-1):
        sum = 0
        for j in range(i+1,n):
            sum += A[i][j] * x[j]
        x[i] = (y[i] - sum) / A[i][i] 
    return x

#solve Ax=y when A is a lower triangular matrix
def Solver_L(A,y):
    m,n = A.shape
    x = np.zeros(n)
    for i in range(0,n):
        sum = 0
        for j in range(0,i):
            sum += A[i][j] * x[j]
        x[i] = (y[i] - sum) / A[i][i]
    
    return x

#solve PA=LU factorization
def LU_Fact(mat):
    n = len(mat)
    interchange = 1
    P = np.eye(len(mat))

    for k in range(n):
        m_abs = np.abs(mat[k:,k])
        pivot = np.max(m_abs)
        if (pivot != mat[k,k]):
            p_index = np.where(m_abs == pivot)[0] + k
            P[p_index],P[k] = P[k],P[p_index]
            mat[p_index],mat[k] = mat[k],mat[p_index]
            interchange = -interchange
        for i in range(k+1,n):
            factor = mat[i][k] / mat[k][k]
            mat[i,k] = factor 
            for j in range(k+1,n):
                mat[i,j] = mat[i,j] - factor * mat[k,j]

    L = np.tril(mat)
    U = np.triu(mat)

    det_U = 1
    for i in range(n):
        L[i,i] = 1;
        det_U = det_U * U[i,i]
    
    det = det_U * interchange
    return P,L,U,det

#solve A=QR factorization by GramSchmidt
def GramSchmidt(mat):
    m,n = mat.shape
    Q = np.zeros((m,min(m,n)))
    R = np.zeros((min(m,n),n))
    for k in range(n):
        Q[:,k] = mat[:,k]
        for i in range(0,k):
            R[i,k] = Q[:,i].T.dot(mat[:,k])
            #print("R:",i,k,R[i,k],"\n")
            Q[:,k] = Q[:,k] - R[i,k] * Q[:,i]
        R[k,k] = np.linalg.norm(Q[:,k])
        #print("R:",k,k,R[k,k],"\n")
        Q[:,k] = Q[:,k] / R[k,k]
        #print("Q:",k,Q[:,k],"\n")
    
    if m == n:
        det = detR(R) * detQ(Q)
    else:
        det = "Error"
    return Q,R,det

#solve A=QR factorization by Householder
def Householder(mat):
    m,n = mat.shape
    Q = np.eye(m)
    R = mat.copy()
    if (m <= n):
        size = m - 1
    else:
        size = n
    for i in range(0,size):
        e = np.zeros_like(R[i:,i])
        I = np.eye(m-i)
        a= np.eye(m)
        e[0] = 1
        mod = np.linalg.norm(R[i:,i])
        u = R[i:,i] - mod * e;
        #print(u)
       
        v1 = u.reshape(len(u),1)
        v2 = u.reshape(1,len(u))
        frac = v2.dot(v1)  
        a[i:,i:] = I - (2/frac) * v1.dot(v2)
        #print(a)
        R = a.dot(R)
        Q = Q.dot(a)
        #print("Q: ", Q,"\n")

    if m == n:
        det = detR(R) * detQ(Q)
    else:
        det = "Error"
    return Q,R,det

#solve A=QR factorization by Givens
def Givens(mat):
    m,n = mat.shape
    Q = np.eye(m)
    R = mat.copy()
    for i in range(0,n):
        for j in range(i+1,m):
            P = np.eye(m)
            x = np.hypot(R[j,i],R[i,i])
            print(R[i,i],R[j,i])
            print(x)
            cos = R[i,i]/x
            sin = R[j,i]/x
            P[j,i] = -sin
            P[i,j] = sin
            P[i,i] = cos
            P[j,j] = cos
            #print(a)
            R = P.dot(R)
            Q = P.dot(Q)
            #print("Q: ", Q,"\n")

    Q = Q.T
    if m == n:
        det = detR(R) * detQ(Q)
    else:
        det = "Error"
    return Q,R,det

#solve A=URV^t factorization
def URV(mat):
    m,n = mat.shape
    r = np.linalg.matrix_rank(mat)
    Q,B_1,_ = Householder(mat)
    B = B_1[:r,:]
    print(B.T)
    q,B_2,_ = Householder(B.T)
    T = B_2[:r,:r]
    V = q
    U = Q
    R = np.zeros_like(mat)
    R[:r,:r] = T.T
    if m == n:
        det = detR(R) * detQ(U) * detQ(V.T)
    else:
        det = "Error"
    return U,R,V,det