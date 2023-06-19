import numpy as np

# # def func(x):
# #     g = 9.806
# #     k = 0.00341
# #     f = np.log10(np.cosh(x*np.sqrt(g*k)))- 50
# #     f2 = ( ( 4 * np.cos(x) ) - ( np.e**(2*x) ) )
# #     return f2

# # print(func(0))

# def getFunctions(path):
#     with open(path, 'r') as arquivo:
#         linhas = arquivo.readlines()

#     for linha in linhas:
#         nome, expressao = linha.strip().split('=')

#     funcoes = []
#     for linha in linhas:
#         nome, expressao = linha.strip().split('=')
#         args = nome.split('(')[1].split(')')[0].split(',')
#         expressao = expressao.replace('log','np.log10')
#         expressao = expressao.replace('sen','np.sin')
#         expressao = expressao.replace('cos','np.cos')
#         expressao = expressao.replace('e','np.e')
        
#         print('lambda ' + ','.join(args) + ':' + expressao)
#         funcao = eval('lambda ' + ','.join(args) + ':' + expressao)
#         funcoes.append(funcao)

#     return [np.array(funcoes).T, args]

# # import numpy as np

# funcoes,args = getFunctions('funcoes.txt')

# # # print(funcoes[0](*[1,2]))
# # from sympy import symbols, diff

# # f = funcoes[0]
# # # derv = diff(funcoes[0], symbols('x',real = True))
# # x, y, z = symbols('x y z')
# # print(diff(x**2, x))

# # # print(funcoes[0](*[3,2]))

# def partialDer(f,h,xArray,n): # f
#     xArrDif = xArray.copy()
#     xArrDif[n-1] += h # n é a n-ésima variável.
#     fa = f(*xArray)
#     fi = f(*xArrDif) 
#     return ((fi - fa) / h)

# f = funcoes[0]
# print(f)
# #print(f(1,2,3))
# print(partialDer(f,0.00001,[1,2,3],2))

# numArgsMatrix = 3
# xArray = np.ones(shape=(numArgsMatrix,))
# print(xArray)

# functionsVector,numArgsMatrix = getFunctions('funcoes.txt')
# xArray = np.ones(shape=(numArgsMatrix,))
# jacobianMatrix = getJacobianM(functionsVector,numArgsMatrix,xArray)

# print(jacobianMatrix)
