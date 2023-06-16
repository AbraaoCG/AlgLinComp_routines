# from rootSolver import *
import numpy as np
# Método para obter funcoes de um arquivo de texto no formato do arquivo 'funcoes.txt'
def getFunctions(path):
    import math

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
        expressao = expressao.replace('sen','np.sin')
        expressao = expressao.replace('cos','np.cos')
        expressao = expressao.replace('e','np.e')
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

def newtonMethod(funcPath):
    # Argumentos de execução
    tol = 1e-05
    itMax = 10000
    # Obter funcoes do usuário a partir de um TXT em 'funcPath'
    functionsVector,numArgsMax = getFunctions(funcPath)
    numFunc = len(functionsVector)
    # Prepara Valores de partida
    hList = [0.001 ] * numArgsMax
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
        # Verificar condição de parada.
        if (tolK < tol):
            print('Convergência alcançada!')
            return xArray
        # Recalcular Jacobiana
        jacobianMatrix = getJacobianM(functionsVector,numArgsMax,xArray,dX_Vector)
        
    print('Convergência não alcançada!')

# Método de Broyden para resolver sistema de equações não lineares.

def broydenMethod(funcPath):
    # Argumentos de execução
    tol = 1e-05
    itMax = 10000
    # Obter funcoes do usuário a partir de um TXT em 'funcPath'
    functionsVector,numArgsMax = getFunctions(funcPath)
    numFunc = len(functionsVector)
    # Prepara Valores de partida
    Y = [ 1 ] * numFunc
    hList = [0.001 ] * numArgsMax
    xArray = np.ones(numArgsMax) # X0 = Vetor coluna preenchido com '1's
    jacobianMatrix = getJacobianM(functionsVector,numArgsMax,xArray,hList) # Jacobiana Inicial
    
    fX_Vector = np.ndarray(shape=(numFunc,1))
    it = 0
    # Algorítimo principal
    while(it < itMax):
        # Computar inversa da Jacobiana
        # print(jacobianMatrix)
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
        # Verificar condição de parada.
        if (tolK < tol):
            print('Convergência alcançada!')
            return xArray
        # Recalcular Jacobiana com método de Broyden
        Yk = np.ndarray(shape=(numFunc,1))
        for i in range(numFunc): # Computar Y a partir de f(Xk+1) - f(Xk)
            Yk[i][0] = functionsVector[i](*xArray) - fX_Vector[i][0] 
        dX_MatrixT = dX_Matrix.T
        jacobianMatrix += ( np.dot(Yk - ( np.dot(jacobianMatrix,dX_Matrix) ) , dX_MatrixT)  / np.dot(dX_MatrixT,dX_Matrix))
        
    print('Convergência não alcançada!')


funcPath = 'funcoes.txt'
xArr = broydenMethod(funcPath)
print(xArr)

functionsVector,numArgsMax = getFunctions(funcPath)
for f in functionsVector:
    print(f(*xArr))


# xArray = [2,1,1]
# functionsVector,numArgsMatrix = getFunctions('funcoes.txt')
# jacobianMatrix = getJacobianM(functionsVector,numArgsMatrix,xArray)

# print(jacobianMatrix)