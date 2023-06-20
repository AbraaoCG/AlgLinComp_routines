
import numpy as np
def LU_Solver(A,B):
    L,U = LU_decomposition(A)
    Y = forward_substitution(L, B)
    X = back_substitution(U, Y)

    return X

def LU_decomposition(A):
    n = len(A)
    L = [[0.0] * n for i in range(n)]
    U = [[0.0] * n for i in range(n)]

    for j in range(n):
        L[j][j] = 1.0

        for i in range(j + 1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = A[i][j] - s1

        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (A[i][j] - s2) / U[j][j]

    return L, U


def forward_substitution(L, B):
    n = len(L)  # Dimensão da matriz L
    m = len(B[0])  # Número de colunas de B
    # Inicializa matriz solução Y com zeros
    Y = []  
    for i in range(n):
        Y.append([])
        for j in range(m):
            Y[i].append(0)
    for j in range(m):  # Loop pelas colunas de B
        for i in range(n):  # Loop pelas linhas de L e Y
            # Somatorio dos termos com coeficiente já calculados ( à esquerda do pivo atual i)
            
            s = sum(L[i][k] * Y[k][j] for k in range(i))
            if (L[i][i] != 0): # Verifica Divisão por 0 ( matrix singular ).
                Y[i][j] = (B[i][j] - s) / L[i][i] # Calculo do coeficiente na linha e coluna atuais.
            else:
                raise ZeroDivisionError
    return Y  # Retorna a matriz solução Y


def back_substitution(U, Y):
    m = len(U)  # Dimensão da matriz U
    n = len(U[0])  # Número de colunas de Y
    yl = len(Y)
    yc = len(Y[0])
    # Inicializa matriz solução X com zeros
    X= []
    for i in range(yl):
        X.append([])
        for j in range(yc):
            X[i].append(0)
    for k in range(yc):  # Loop pelas colunas de Y
        for i in range(m-1, -1, -1):  # Loop pelas linhas de U e X (de trás para frente)
            s = 0
            for j in range(i+1,m):
                s += U[i][j] * X[j][k]
            X[i][k] = (Y[i][k] - s) / U[i][i]

    return X  # Retorna a matriz solução X

def transpose(X):
    Xt = []
    for i in range(X.shape[1]):
        Xt.append([])
        for j in range(X.shape[0]):
            Xt[i].append(X[j][i])
    return Xt
            
def matrix_multiplication(A, B):
    m = len(A)
    n = len(B)
    p = len(B[0])
    if len(A[0]) != n:
        raise ValueError("Matrizes não compatíveis para multiplicação")
    C = [[0 for j in range(p)] for i in range(m)]
    for i in range(m):
        for j in range(p):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

A = [[64.09606402,64.09606402,4.006004],[2.001     ,  2.001      , 2.001  ],[ 3.003001 ,  -1.       ,   1.  ]]
B = [[1,0,0],[0,1,0],[0,0,1]]
print(LU_Solver(A,B))