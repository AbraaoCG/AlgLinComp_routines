# Esse programa pergunta ao usuário qual o método desejado para encontrar
# a raiz de uma função.

import numpy as np

# Definir f(x)

def func(x):
    g = 9.806
    k = 0.00341
    f = np.log(np.cosh(x*np.sqrt(g*k)))- 50
    f2 = x
    return f

# Definir algorítimo método da bisseção
def bissec_method():
    # Definir Tolerância
    tol = 1e-5
    print(f'Tolerância utilizada: {tol}\n Para mudar, realize alteração no código.')
    # # Definir a e b de inicio. Supoe-se que f(a) < 0 e f(b) > 0 , e que f(x) tem ao menos raiz.
    # a = -1
    # b = 1
    # Perguntar Parâmetros ao usuário
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

# Apresentação do programa
print("Esse programa é destinado ao cálculo de raizes de funções utilizando três diferentes métodos")
print("Tenha certeza de programar sua função corretamente dentro do código.\n")

# Obter ICOD --> método de resolução
ICOD = -1
ICOD = int(input('Qual método desejado para encontrar a raiz?\n   ICOD 1: Método da Bisseção\n   ICOD 2: Método de Newton\n   ICOD 3: Método da Interpolação inversa.\n'))

match ICOD:
    case 1:
        root = bissec_method()
        print(f'A raiz da função é {root}')
    case 2:
        pass
    case 3:
        pass



