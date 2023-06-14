# Esse programa pergunta ao usuário qual o método desejado para encontrar
# a raiz de uma função e realiza o cálculo da raiz, imprimindo resultado na tela.

import numpy as np

# Definir funções complementares.

# Cálculo de X
def func(x):
    g = 9.806
    k = 0.00341
    f1 = np.log10(np.cosh(x*np.sqrt(g*k)))- 50 # Método bissecante: a,b = -600,-650 ; Método de Newton x0 = -600
    #f2 = ( ( 4 * np.cos(x) ) - ( np.e**(2*x) ) ) # Método bissecante: a,b = 1,0 ; Método de Newton x0 = 0-10
    return f1

# Cálculo da derivada primeira de x de forma literal ( com expressão dada ).
def d1_func(x):
    g = 9.806
    k = 0.00341
    # Derivada de f1 e f2 para método de Newton Original
    df1 = (np.tanh(x*np.sqrt(g*k)) * np.sqrt(g*k)) / np.log(10)
    #df2 = ( -4 * np.sin(x) )- ((np.e**(2*x)) * 2)

    return df1

# Cálculo da derivada primeira de x de forma numérica ( com diferenças finitas ).

def d1_func_numericMethod(xi,xa):
    h = xi - xa
    return ((func(xi) - func(xa)) / h)

# ------------------------------------------------------------------------------------------------------
# Definir algorítimo método da bisseção
def bissecMethod():
    # Definir Tolerância e número máximo de iterações
    tol = 1e-5
    itMax = 10000
    print(f'Tolerância utilizada: {tol}\n Para mudar, realize alteração no código.')
    a = b = ''
    while( (type(a) != type(1.0)) and ( type(a) != type(1.0))):# Enquanto a e b não forem numéricos.
        a = float(input('Insira a: Primeiro valor de x, tal que f(a) < 0\n a = '))
        b = float(input('Insira b: Segundo valor de x, tal que f(b) > 0\n b = '))
        if((type(a) == type('')) or ( type(a) == type(''))): # Se inserido valor não-numérico
            print('Valor inserido invalido. Deve ser um valor numérico! \n')
    # Algorítimo do método da biseção
    xi = np.nan ; fi = np.nan
    it = 0
    while(abs(a - b) > tol):
        if (it > itMax):
            print('Convergência não alcançada')
            return np.nan    
        xi = (a+b) / 2
        fi = func(xi)
        it += 1
        if (fi < 0):
            a = xi
        elif(fi > 0):
            b = xi
        if(fi == 0.0):
            break
    print(f'Convergência alcançada em {it} iterações.')
    return xi    

# ------------------------------------------------------------------------------------------------------
# Definir Algorítimo do método de Newton

def newtonMethod(method):
    # Valores de dx, número máximo de iterações e tolerância definidos em código.
    dx = 0.001
    itMax = 10000
    tol = 1e-5
    print(f'Método de Newton utilizando: \ndx = {dx}\nNúmero máximo de iterações = {itMax}\n')
    print("Tenha certeza de definir f(x) em código antes de executar o programa.\n Caso Método Original seja selecionado, também definda f'(x) em código.")
    
    # xo é definido pelo usuário, e é equivalente à xa na primeria iteração.
    xa = float(input('Insira x0( o valor inicial de x ) para execução do método: \n xo = '))

    # Algorítimo de aproximação da raiz de f(x)
    xi = xa + dx # x na iteração atual.
    it = 0
    while ( it < itMax ):
        # Cálculo de f'(x), em função do método ( Original ou Secante ).
        if (method == 0): # Original
            df_xi = d1_func(xi)
        elif(method == 1): # Secante
            df_xi = d1_func_numericMethod(xi = xi, xa = xa)
        
        # Cálculo do x próximo
        x_next = xi - (func(xi) / df_xi)

        if(abs(x_next - xi) <= tol):
            print(f'O método de Newton convergiu em {it} iterações!')
            return x_next # x_next é a raiz da função.
        else:
            xa = xi
            xi = x_next
        it += 1

    # Se não retornar é porque não convergiu
    print('O método de Newton não convergiu!')

# ------------------------------------------------------------------------------------------------------
# Definir algorítimo do método da interpolação Inversa.

def invInterpolMethod():
    tol = 1e-5
    itMax = 10000
    print(f'Método da Interpolação Inversa utilizando: \nNúmero máximo de iterações = {itMax}\n')
    # Obter Valores de incialização
    x1 = float(input('Insira x1: '))
    x2 = float(input('Insira x2: '))
    x3 = float(input('Insira x3: '))    
    # Algorítimo principal
    xi = 1e+36 # Inicializo 'xi' como 'teto'.
    it = 0
    while(it < itMax):

        y1 = func(x1) ; y2 = func(x2) ; y3 = func(x3)
        xi_next = (y2 * y3 * x1) / ((y1 - y2) * (y1 - y3)) + (y1 * y3 * x2) / ((y2 - y1) * (y2 - y3) )+ (y1 * y2 * x3) / ((y3 - y1) * (y3 - y2))

        if (abs(xi_next - xi) < tol):
            print(f'Convergência alcançada em {it} iterações!')
            return xi_next
        # Loop para identificar maior elemento entre os f(xk), k = 1,2,3
        cont = 1 ; max = 0 ; imax = -1
        for x in [y1,y2,y3]:
            if (abs(x) > max):
                max = abs(x)
                imax = cont
            cont += 1
        # Trocar xk por xi_next calculado.
        match imax:
            case 1:
                x1 = xi_next
            case 2:
                x2 = xi_next
            case 3:
                x3 = xi_next
        xi = xi_next
        it += 1
    print('Convergência não alcançada')
# ------------------------------------------------------------------------------------------------------
# Apresentação do programa
print("Esse programa é destinado ao cálculo de raizes de funções utilizando três diferentes métodos")
print("Tenha certeza de programar sua função corretamente dentro do código.\n")

# Obter ICOD --> método de resolução
ICOD = -1
ICOD = int(input('Qual método desejado para encontrar a raiz?\n   ICOD 1: Método da Bisseção\n   ICOD 2: Método de Newton\n   ICOD 3: Método da Interpolação inversa.\n'))

match ICOD:
    case 1: # Método da Bisseção.
        root = bissecMethod()
        print(f'A raiz da função é {root}')
    case 2: # Método de Newton.
        method = int(input('Insira 0 para método de newton original e 1 para método de newton Secante:\n'))
        root = newtonMethod(method)
        print(f'A raiz da função é {root}')
    case 3:
        root = invInterpolMethod()
        print(f'A raiz da função é {root}')


