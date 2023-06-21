
import numpy as np

def integrate_polynomial(f, a, b, n):
    # Função para realizar a integração polinomial
    x = np.linspace(a, b, n)  # Pontos de integração igualmente espaçados
    print(x)
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

def simpson_Rule(f,a,b,L):
    return (L/6) * (f(a) + 4*f((a+b)/2) + f(b))

def simpson_Rule_adaptative(f, a, b, tol):
    itMax = 50
    L = b - a
    Ik = simpson_Rule(f,a,b,L)  # Calcula integral com 3 pontos e 1 sub-intervalo.
    it = 1 ; n = 3
    while(it < itMax):
        # Parâmetros da iteração.
        numS = 2**it # Número de subintervalos.
        tSub = L*(1/2)**it # Tamanho de um subintevalo.
        n = int(3 + (2 * ( ((2**(it)) - 1 ) / 1 ))) # PA de segunda ordem que indica número de pontos.
        # É preciso iterar no número de subintervalos e somar as integras associadas a cada um.
        Ik_proximo = 0
        b_new = a # 'b_new' é setado como a ( inicio do intervalo de integração )
        for j in range(1,numS+1):
            a_new = b_new # Novo 'a' para regra de simpson = recebe inicio do intervalo de int ou fim do subintervalo anterior.
            b_new = a_new + tSub # Novo 'b' para regra de simpson. = é o inicio do subintervalo + tamanho do subintervalo.
            Ik_proximo += simpson_Rule(f,a_new,b_new,tSub) # L novo é o tamanho do subintervalo.
        # Supoe-se que f seja uma função de uma variável e realiza a regra de simpson.
        tolk = abs(Ik_proximo - Ik) / abs(Ik_proximo)
        print(f'{numS} Sub-intervalos e {n} pontos: \ntolk = {tolk}  ;  I{it+1} = {Ik_proximo}')
        Ik = Ik_proximo
        
        if (tolk < tol):
            print(f'Tolerância desejada alcançada com n = {n}')
            return Ik_proximo
        it+= 1
    print(f'Não houve convergência para a tolerância desejada em {itMax} iterações.')


ICOD = print("Programa de Integração numérica adaptativa utilizando regra de simpson em um número crescente de subintervalos.")

funcs,numF = getFunctions('funcaoInt.txt')

if(numF > 0):
    f = funcs[0]
    a = float(input("Seja a e b os limites de integração: a --> b\nInsira a: "))
    b = float(input("Insira b: "))
    tol = float(input("Insira a tolerância aceita: "))

    If = simpson_Rule_adaptative(f, a, b, tol)
    print(f'\nIntegral = {If}')
else:
    print("Certifique-se de inserir a função f(x) no arquivo 'funcaoInt.txt', no formato f(x) = ...função...")