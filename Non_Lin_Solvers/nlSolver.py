# from rootSolver import *
import numpy as np
# Método para obter funcoes de um arquivo de texto no formato do arquivo 'funcoes.txt'
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
        expressao = expressao.replace('pi','np.pi')
        # print('lambda ' + ','.join(args) + ':' + expressao)
        funcao = eval('lambda ' + ','.join(args) + ':' + expressao)
        funcoes.append(funcao)
        numArgsF.append(len(args))

    return [np.array(funcoes).T,max(numArgsF)]


# Função para obtenção de derivadas parciais com diferenças finitas.

def partialDer(f,h,xArray,n):
    xArrDif = xArray.copy()
    xArrDif[n-1] += h # n é a n-ésima variável.
    fa = f(*xArray)
    fi = f(*xArrDif) 
    return ((fi - fa) / h)

# Montar Jacobiana
def getJacobianM(functionsVector,numArgsMax,xArray,hList):
    m = len(functionsVector)
    n = numArgsMax
    jacobianM = np.ndarray(shape=(m,n))
    for i in range(m):
        for j in range(n): # Itero sobre cada variável
            f = functionsVector[i]
            jacobianM[i][j] = partialDer(f,hList[j],xArray,j+1)
    return jacobianM

# Método de Newton para resolver sistema de equações não lineares.
def newtonMethod(funcPath,tol):
    # Argumentos de execução
    itMax = 10000
    # Obter funcoes do usuário a partir de um TXT em 'funcPath'
    functionsVector,numArgsMax = getFunctions(funcPath)
    numFunc = len(functionsVector)
    # Prepara Valores de partida
    dxi = 0.001
    hList = [dxi] * numArgsMax
    xArray = np.ones(numArgsMax) # X0 = Vetor coluna preenchido com '1's
    jacobianMatrix = getJacobianM(functionsVector,numArgsMax,xArray,hList) # Jacobiana Inicial
    fX_Vector = np.ndarray(shape=(numFunc,1))
    it = 0
    # Algorítimo principal
    while(it < itMax):
        # Computar inversa da Jacobiana
        invJacM = np.linalg.inv(jacobianMatrix)
        # Computar f(Xk)
        for i in range(numFunc):
            fX_Vector[i][0] = functionsVector[i](*xArray)
        # Computar DX
        dX_Matrix = np.dot(-invJacM,fX_Vector)
        dX_Vector = np.reshape(dX_Matrix,newshape=(numArgsMax,))
        # Computar Xk+1
        xArray += dX_Vector
        # Computar as normas de Xk+1 e dXk
        normX = 0 ; normDx = 0
        for i in range(numArgsMax):
            normX += xArray[i]**2
            normDx += dX_Vector[i]**2
        normX = np.sqrt(normX) ; normDx = np.sqrt(normDx)
        tolK = normDx / normX
        print(f'Solução na iteração {it} : {xArray}')
        print(f'Tolerância na iteração {it} : {tolK}\n')
        # Verificar condição de parada.
        if (tolK < tol):
            print(f'Convergência alcançada na iteração {it}!')
            return xArray
        # Recalcular Jacobiana
        jacobianMatrix = getJacobianM(functionsVector,numArgsMax,xArray,dX_Vector)
        it += 1

    print('\nConvergência não alcançada!')
    return np.nan

# Método de Broyden para resolver sistema de equações não lineares.
def broydenMethod(funcPath,tol):
    # Argumentos de execução
    itMax = 10000
    # Obter funcoes do usuário a partir de um TXT em 'funcPath'
    functionsVector,numArgsMax = getFunctions(funcPath)
    numFunc = len(functionsVector)
    # Prepara Valores de partida
    dxi = 0.001
    hList = [dxi ] * numArgsMax
    xArray = np.ones(numArgsMax) # X0 = Vetor coluna preenchido com '1's
    jacobianMatrix = getJacobianM(functionsVector,numArgsMax,xArray,hList) # Jacobiana Inicial
    fX_Vector = np.ndarray(shape=(numFunc,1))
    it = 0
    # Algorítimo principal
    while(it < itMax):
        # Computar inversa da Jacobiana
        invJacM = np.linalg.inv(jacobianMatrix)
        # Computar f(Xk)
        for i in range(numFunc):
            fX_Vector[i][0] = functionsVector[i](*xArray)
        # Computar DX
        dX_Matrix = np.dot(-invJacM,fX_Vector)
        dX_Vector = np.reshape(dX_Matrix,newshape=(numArgsMax,))
        # Computar Xk+1
        xArray += dX_Vector
        # Computar as normas de Xk+1 e dXk
        normX = 0 ; normDx = 0
        for i in range(numArgsMax):
            normX += xArray[i]**2
            normDx += dX_Vector[i]**2
        normX = np.sqrt(normX) ; normDx = np.sqrt(normDx)
        tolK = normDx / normX
        print(f'Solução na iteração {it} : {xArray}')
        print(f'Tolerância na iteração {it} : {tolK}\n')
        # Verificar condição de parada.
        if (tolK < tol):
            print(f'Convergência alcançada na iteração {it}!')
            return xArray
        it = it + 1
        # Recalcular Jacobiana com método de Broyden
        Yk = np.ndarray(shape=(numFunc,1))
        for i in range(numFunc): # Computar Y a partir de f(Xk+1) - f(Xk)
            Yk[i][0] = functionsVector[i](*xArray) - fX_Vector[i][0] 
        dX_MatrixT = dX_Matrix.T
        jacobianMatrix += ( np.dot(Yk - ( np.dot(jacobianMatrix,dX_Matrix) ) , dX_MatrixT)  / np.dot(dX_MatrixT,dX_Matrix))
    print('\nConvergência não alcançada!')
    return np.nan



print("Insira no arquivo funcoes.txt as funções fi(x) tal que fi(x) = 0\n Faça isso no formato f(x1,...,xn) = ... , onde n é o número máximo de variáveis no sistema de equações.\n")
funcPath = 'funcoes.txt'
tol = float(input("Insira a tolerância aceita: "))
ICOD = int(input("Insira o método de solução: \n Método de Newton -> 1\n Método de Broyden -> 2\nICOD : "))

xArr = np.nan
if (ICOD == 1):
    xArr = newtonMethod(funcPath=funcPath,tol=tol)

elif(ICOD == 2):
    xArr = broydenMethod(funcPath=funcPath,tol=tol)


if(type(xArr) != type(np.nan)):
    print(f'\nO Vetor solução é{xArr}')
else:
    print('Não foi encontrado resultado válido.\n Caso tenha convergido, revise a entrada das funções e os argumentos e tente novamente')