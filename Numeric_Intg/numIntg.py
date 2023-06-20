
import numpy as np

def integrate_polynomial(f, a, b, n):
    # Função para realizar a integração polinomial
    x = np.linspace(a, b, n)  # Pontos de integração igualmente espaçados
    # w = np.ones_like(x) / n  # Pesos iguais para a integração polinomial
    W = getW(x,a,b)
    integral = np.dot(f(x), W) # Cálculo da integral
    return integral[0]

def integrate_gauss(f, a, b, n, weight_func):
    # Função para realizar a integração por quadratura de Gauss
    x, w = weight_func(n)  # Obtenção dos pontos de integração e pesos usando a função fornecida
    x_scaled = 0.5 * (b - a) * x + 0.5 * (b + a) 
    w_scaled = 0.5 * (b - a) * w 
    integral = np.dot(f(x_scaled), w_scaled)  # Cálculo da integral
    return integral

# Exemplo de função a ser integrada
def f(x):
    # return (x**3 + 2)
    return np.exp(- (x**2) / 2) / np.sqrt(2 * np.pi)# x ** 3 + 2 * x ** 2 + 3 * x + 4

# Função para cálculo dos pesos na quadratura de Gauss
def gauss_weights(n):
    x, w = np.polynomial.legendre.leggauss(n)  # Obtem os pontos de integração e pesos da quadratura de Gauss tabelados
    return x, w

def construct_vandermonde_matrix(xArr):
    # V = np.vander(xArr, increasing=True)
    n = len(xArr)
    V = np.ndarray(shape=(n,n))
    for i in range(n):
        for j in range(n):
            V[i][j] = xArr[j]**i
    return V

def construct_k_vector(a, b, n):
    K = np.ndarray(shape=(n,1))
    for i in range(1,n+1):
        K[i-1][0] = (b**i - a**i) / i
    return K

def getW(xArr, a, b):
    n = len(xArr)
    VandM_Matrix = construct_vandermonde_matrix(xArr)
    KMatrix = construct_k_vector(a, b, n)

    W = np.linalg.solve(VandM_Matrix, KMatrix)
    return W


def getFunctions(path):
    with open(path, 'r') as arquivo:
        linhas = arquivo.readlines()

    for linha in linhas:
        nome, expressao = linha.strip().split('=')

    funcoes = []
    numArgsF = []
    for linha in linhas:
        nome, expressao = linha.strip().split('=')
        args = nome.split('(')[1].split(')')[0].split(',')
        expressao = expressao.replace('log','np.log10')
        expressao = expressao.replace('ln','np.log')
        expressao = expressao.replace('sen','np.sin')
        expressao = expressao.replace('cos','np.cos')
        expressao = expressao.replace('e','np.e')
        # print('lambda ' + ','.join(args) + ':' + expressao)
        funcao = eval('lambda ' + ','.join(args) + ':' + expressao)
        funcoes.append(funcao)
        numArgsF.append(len(args))

    return [np.array(funcoes).T,max(numArgsF)]


def Poly_Integ_adaptative(f, a, b, tol):
    it = 2 ; itMax = 10000
    Ik = integrate_polynomial(f=f, a=a, b=b, n=1 ) # Calcula integral com 1 ponto.
    while(it < itMax):
        Ik_proximo = integrate_polynomial(f=f, a=a, b=b, n=it) # Calcula integral 
        tolk = abs(Ik_proximo - Ik) / abs(Ik_proximo)
        print(f'TolK com {it} pontos = {tolk}')
        Ik = Ik_proximo
        if (tolk < tol):
            print(f'Tolerância desejada alcançada com n = {it}')
            return Ik_proximo
        it+= 1
    print(f'Não houve convergência para a tolerância desejada em {itMax} iterações.')

def QuadGauss_Integ_adaptative(f, a, b, tol):
    it = 2 ; itMax = 10000
    Ik = integrate_gauss(f=f, a=a, b=b, n=1, weight_func=gauss_weights)  # Calcula integral com 1 ponto.
    while(it < itMax):
        Ik_proximo = integrate_gauss(f=f, a=a, b=b, n=it, weight_func=gauss_weights)  
        tolk = abs(Ik_proximo - Ik) / abs(Ik_proximo)
        print(f'TolK com {it} pontos = {tolk}')
        Ik = Ik_proximo
        if (tolk < tol):
            print(f'Tolerância desejada alcançada com n = {it}')
            return Ik_proximo
        it+= 1
    print(f'Não houve convergência para a tolerância desejada em {itMax} iterações.')

ICOD = int(input("Programa de Integração numérica. Defina a função no arquivo 'funcaoInt.txt' e selecione o tipo de integração desejado:\n 1 : Integração Polinomial\n 2 : Integração por Quadratura de Gauss\n Seleção = "))

funcs,numF = getFunctions('funcaoInt.txt')

if(numF > 0):
    f = funcs[0]
    a = float(input("Seja a e b os limites de integração: a --> b\nInsira a: "))
    b = float(input("Insira b: "))
    tol = float(input("Insira a tolerância aceita: "))

    if (ICOD == 1):
        # Uso da função de integração polinomial
        integral_polynomial = Poly_Integ_adaptative(f, a, b, tol)
        print("Integração Polinomial:", integral_polynomial)
    elif (ICOD == 2):
        # Uso da função de integração por quadratura de Gauss
        integral_gauss = QuadGauss_Integ_adaptative(f, a, b, tol)
        print("Integração Gauss:", integral_gauss)
    
else:
    print("Certifique-se de inserir a função f(x) no arquivo 'funcaoInt.txt', no formato f(x) = ...função...")