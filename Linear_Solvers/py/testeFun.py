from functions1 import *

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime

Y = []
X = np.arange(12)

for x in X:
    #value = 2 + 5 * x
    value = 15 + 22 * x + 30 * (x**2) + 123 * (x**20) # POLI
    # value = np.log(x) # * ( 1 + np.random.random() * 0.15 * (-1 ** x)) # LN ou EXP
    Y.append( value)

janela = 5


# medMovExp = media_movel_exponencial(Y,janela)
# medMovArt = media_movel_aritmetica(Y,janela)

# RegExp = exponential_regression(X,Y)
# RegLN = logarithmic_regression(X,Y)
print(X)
RegPoli_pred = polyfit(X,Y,20)

plt.plot(X,Y, label="Input")

# #Plot de Medias Móveis

# plt.plot(X,medMovExp, label="Exponencial")
# plt.plot(X[janela - 1:len(X)],medMovArt, label="Aritm")

# Plot de Regressões

# plt.plot(X, RegExp,label="Reg Exp")
# plt.plot(X, RegLN,label="Reg LN")

plt.plot(X, RegPoli_pred,label="Reg Polinomial")

# plt.ylim(5,15)
# plt.xlim(115,120)
plt.legend()
plt.show()



# A = np.random.random((3,3))
# B = np.identity(3)


# for x in A:
#     print(x)

# print('\n ----------------------------')

# for x in B:
#     print(x)

# print('\n ----------------------------')

# time0 = datetime.now()

# A_inv = LU_Solver(A,B)
# print('\n ----------------------------')

# time1 = datetime.now()

# A_inv_np = np.linalg.inv(A)

# time2 = datetime.now()

# test = np.dot(A,A_inv)
# test2 = np.dot(A,A_inv_np)
# print('Solucao propria')
# for x in test:
#     print(x)

# # print('\n ----------------------------')

# print('Solucao Numpy')

# for x in test2:
#     print(x)

# print(f'\n\n{time1 - time0}       {time2 - time1}')