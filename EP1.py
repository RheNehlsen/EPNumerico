import numpy as np


def calcula_coefs(alpha, beta):
  '''
    Calcula os coeficentes da matriz de Givens pelo método estável
    Inputs :
      Alpha(float) : Elemento da diagonal principal
      Beta (float) : Elmento da subdiagonal
    Output:
     (ck(float), sk(float)): Uma tupla com o cosseno e seno associados
  '''
  if abs(alpha)>abs(beta):
    tau = -beta/alpha
    ck = 1/np.sqrt(1 + tau**2)
    sk = ck*tau
  else:
    tau = -alpha/beta
    sk = 1/np.sqrt(1 + tau**2)
    ck = sk*tau
  return ck, sk

def calcula_R(A, B, n):
  '''
    Calcula a digonal e a sobrediagonal da matriz R e todos os coencientes da matrizes Q
    Input:
      A(list) : diagonal da matriz A
      B(list) : subdiagonal da matriz A
      n(int)  : Tamanho da matriz A
    Outputs:
      A(list)              : diagonal da matriz R
      sobrediagonal (list) : sobrediagonal da matriz R
      Givens(list)         : Vetor de tuplas dos coeficentes de Q utilizados

  '''
  Givens = []
  sobrediagonal = B.copy()
  for i in range(n):
    ck, sk = calcula_coefs(A[i], B[i])
    Givens.append((ck, sk))
    ai = ck*A[i] -sk*B[i] #A[i] novo
    bi = ck*sobrediagonal[i] - sk*A[i+1] #B[i] novo
    aii = sk*sobrediagonal[i] + ck*A[i+1] #A[i+1] novo
    bii = ck*sobrediagonal[i+1] #B[i+1] novo
    A[i], sobrediagonal[i], A[i+1], sobrediagonal[i+1] = ai, bi, aii, bii
  return A,sobrediagonal, Givens


def calcula_A(A,sobrediagonal, Givens,n):
  '''
    Calcula a diagonal e subdiagonal da nova matriz A
    Inputs:
      A(list)             : diagonal da matriz R
      sobrediagonal(list) : sobrediagonal da matriz R
      Givens(list)        : Vetor de tuplas dos coeficentes de Q utilizados
      n(int)              : Tamanho da matriz A
    Outputs:
      A(list)             : diagonal da matriz A
      subdiagonal(list)   : subdiagonal da matriz A
  '''

  subdiagonal = sobrediagonal.copy()
  for i in range(n):
    ck, sk = Givens[i]
    ai = A[i]*ck - sk*subdiagonal[i]
    bi = -A[i+1]*sk
    aii = A[i+1]*ck
    A[i], subdiagonal[i], A[i+1] = ai, bi, aii
  return A,subdiagonal


def calcula_autovetores(V,givens):
  """
    Faz a multiplicação da matriz de autovetores por Q^t
    Inputs:
      V(np.array) : matriz de autovetores
      Q(list)     : vetor de coeficentes de Givens
    Outputs:
      V(np.array) : matriz de autovetores calculada
  """
  for i in range(len(givens)):  
      for j in range(len(V[0])):
          Vji = (V[j][i]*givens[i][0]) - (V[j][i+1]*givens[i][1])
          Vji1 = (V[j][i]*givens[i][1]) + (V[j][i+1]*givens[i][0])
          V[j][i] = Vji
          V[j][i+1] = Vji1
  return(V)

def wilkinson(a_anterior, a_atual, beta):
  """
    Calcula o coefiente mu a partir da heurística de wilkinson
  """
  d = (a_anterior - a_atual)/2
  return a_atual + d -(np.sign(d)*np.sqrt((d**2) + (beta**2)))

def QR(A,B,V):
  """
    Aplica o algoritmo QR para encontrar autovalores e autovetores
    Inputs:
      A(list)      : diagonal da matriz A
      B(list)      : subdiagonal da matriz A
      V(np.array)  : Matriz H^t
    Outpus:
      A(list)      : autovalores
      V(np.array)  : matriz de autovetores
      k(int)       : numero de iteracoes
  """

  k = 0

  for i in range(len(A)-1, 0 , -1):  #loop principal
    convergiu = False
    while not convergiu:
      if k>0:
       mu = wilkinson(A[i-1], A[i], B[i-1]) #aplica a heurística de wilkinson para calcular mu
      else:
        mu = 0
      A = [(elem - mu) for elem in A] # subtrai mu da diagonal principal
      A,R,Q = calcula_R(A,B, i) # calcula as diagonais de R e os coeficentes de Q
      A,B = calcula_A(A,R,Q, i) # calcula a nova A a partir de R
      A = [(elem + mu) for elem in A] # adiciona mu a diagonal principal
      V = calcula_autovetores(V, Q) #calcula os autovetores a partir de Q
      k+=1
      if abs(B[i-1]) < 10**-6:
        convergiu= True
  return A, V, k