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

import numpy as np

def integrate_polynomial(f, a, b, n):
    # Função para realizar a integração polinomial
    x = np.linspace(a, b, n+1)  # Pontos de integração igualmente espaçados
    print(x)
    w = np.ones_like(x) / n  # Pesos iguais para a integração polinomial
    print(w)
    integral = np.dot(f(x), w) * (b - a)  # Cálculo da integral
    return integral

def integrate_gauss(f, a, b, n, weight_func):
    # Função para realizar a integração por quadratura de Gauss
    x, w = weight_func(n)  # Obtenção dos pontos de integração e pesos usando a função fornecida
    x_scaled = 0.5 * (b - a) * x + 0.5 * (b + a)  # Escalonamento dos pontos de integração para o intervalo [a, b]
    w_scaled = 0.5 * (b - a) * w  # Escalonamento dos pesos
    integral = np.dot(f(x_scaled), w_scaled)  # Cálculo da integral
    return integral

# Exemplo de função a ser integrada
def f(x):
    return np.exp(- (x**2) / 2) / np.sqrt(2 * np.pi)# x ** 3 + 2 * x ** 2 + 3 * x + 4

# Exemplo de função para cálculo dos pesos na quadratura de Gauss
def gauss_weights(n):
    x, w = np.polynomial.legendre.leggauss(n)  # Função do NumPy para obter os pontos de integração e pesos da quadratura de Gauss
    return x, w

def construct_vandermonde_matrix(points):
    n = len(points) - 1
    V = np.vander(points, increasing=True)
    return V

def construct_k_vector(a, b, n):
    k = np.arange(n+1)
    K = (b**k - a**k) / k
    K[0] = 1.0  # Tratamento especial para k = 0 para evitar divisão por zero
    return K

def solve_integration_system(points, a, b):
    n = len(points) - 1
    V = construct_vandermonde_matrix(points)
    K = construct_k_vector(a, b, n)
    W = np.linalg.solve(V, K)
    return W



# Exemplo de uso da função de integração polinomial
integral_polynomial = integrate_polynomial(f, 0, 1, 3)
print("Integração Polinomial:", integral_polynomial)

# Exemplo de uso da função de integração por quadratura de Gauss
integral_gauss = integrate_gauss(f, 0, 1, 3, gauss_weights)
print("Integração Gauss:", integral_gauss)
