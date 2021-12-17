'''
Métodos Numéricos e Aplicações - MAP3121
Exercício Programa 2
Rhenan Silva Nehlsen - Turma 3 - 11374871
Vinícius Maalouli Vinha - Turma 7 -11257421

'''

import operator
from EP1 import QR
import numpy as np
import re

def produto_interno(x,y):
  #Retorna o produto interno entre os vetores x e y 
  return sum([x[i]*y[i] for i in range(len(x))])

def norma(x):
  #Retorna a norma do vetor x
  return np.sqrt(produto_interno(x,x))

def monta_omega(a):
  #Monta o vetor w das transformações de Householder a partir
  # de um vetor a 
  a0 = a[0] + np.sign(a[0])*norma(a)

  return  [a0] + a[1:] #Apenas o primeiro elemento é modificado

def househoulder_vetor(x,w,prod_int):
  #Aplica a transfromação de Householder a um vetor x
  escalar = 2*produto_interno(x,w)/prod_int
  vetor =[x[i] - escalar*w[i] for i in range(len(x))]
  return vetor

def tridiagonalizacao_householder(A):
  #Tridiagonaliza uma matriz simétrica A
  # por meio das Trasformações de Householder

  H_transposto = np.eye(len(A)) #Inicia H^t como a identidade
  for i in range(0,len(A)-2):
    w = monta_omega(A[i+1:,i].tolist()) #Monta o vetor omega a partir das colunas
    norma2w = produto_interno(w,w) # Calcula a norma^2 de w

    # Aplica a trasformação na coluna i da matriz, aproveitando que os outros elementos são nulos
    primeiro_elemento = A[i+1,i] - (2*(produto_interno(A[i+1:, i],w))/norma2w)*w[0] 
    A[i+1:, i] = [primeiro_elemento] + [0]*(len(A)-(2+i))


    #Faz a multiplicação por Hw a esquerda
    for j in range(i+1, len(A)):
      A[i+1:,j] = househoulder_vetor(A[i+1:,j],w, norma2w)
    
    #Copia a coluna i para a linha i, aproveitando a simetria
    A[i, i+1:] = A[i+1:, i]

    #Faz a multiplicação por Hw a direita
    for j in range(i+1, len(A)):
      A[j,i+1:] = househoulder_vetor(A[j,i+1:],w, norma2w)
    
    #Faz a multiplicação de H^t por Hw
    for j in range(len(A)):
      H_transposto[j, i+1:] = househoulder_vetor( H_transposto[j, i+1:], w, norma2w)
  
  return A, H_transposto #Retorna a matriz tridiagonal e H^t

def separa_diagonais(A):
  #Separa a diagonal e a subdiagonal da matriz tridiagonal para
  # aplicação do algoritmo QR

  diagonal = []
  subdiagonal = []
  for i in range(len(A)):
    if i==len(A)-1:
      subdiagonal.append(0)
    else:
      subdiagonal.append(A[i+1][i])
    diagonal.append(A[i][i])
  
  return diagonal,subdiagonal


def calcula_valores(A):

  #Calcula os autovalores e autovetores de uma matriz simétrica

  tridiagonal, H_t = tridiagonalizacao_householder(A)


  diagonal, subdiagonal = separa_diagonais(tridiagonal)
  
  autovalores, autovetores, _ = QR(diagonal, subdiagonal, H_t)
  return autovalores, autovetores

def verifica_autovalor(A,autovalores,autovetores):
  for i in range(len(autovalores)):
    autovetor = autovetores[:,i]
    AV = np.matmul(A, autovetor)
    lv = autovalores[i]*autovetor
    print('AxV ', AV)
    print('lxV ', lv )
    print('----------------')
 


def ler_arquivoT1(filename):
  # Lê o arquivo de input do primeiro teste e retorna a matriz A
  A = []
  with open(filename, 'r') as arquivo:
    next(arquivo)
    for linha in arquivo:
      A.append(list(map(float, re.findall(r'\d+', linha))))
  
  return np.array(A)


def teste1(filename):

  A = ler_arquivoT1(filename)
  A_original = A.copy()
  autovalores, autovetores = calcula_valores(A)
  if filename =='input-a':
    autovalores_esperados =   [7, -2, 2, -1]
  else:
    autovalores_esperados = [0.5*((1 - np.cos((2*i -1)*np.pi/(2*len(autovalores) + 1)))**-1) for i in range(1,len(autovalores)+1)]
    autovalores_esperados[-2], autovalores_esperados[-1] = autovalores_esperados[-1], autovalores_esperados[-2]
  
  for i in range(len(autovalores)):
    print('λ ', i+1 )
    print('λ esperado: ', autovalores_esperados[i])
    print('λ obtido: ', autovalores[i])
    print('Erro:' , abs(autovalores[i]- autovalores_esperados[i]))
    print('Autovetor associado a λ: ', autovetores[:, i])
    print('AV - λV = 0? :', (A_original @ autovetores[:, i]) - autovalores[i]*autovetores[:,i])
    print('-----------------------------------------------------------------------------')
  
  return autovalores, autovalores_esperados, A_original, autovetores 


# ------------------------------------------TRELIÇAS--------------------------------------------------------#

def contribuicaoK(K,i,j, angulo, tamanho, AE):
  # Calcula as contribuições de cada barra na matriz de ridigez global

  C = np.cos(angulo*np.pi/180)
  S = np.sin(angulo*np.pi/180)
  C2 = C*C
  S2 = S*S
  CS = C*S

  valores = [C2*AE/tamanho,CS*AE/tamanho, -C2*AE/tamanho, -CS*AE/tamanho,S2*AE/tamanho,  -S2*AE/tamanho]
  posicoes = {
    0  : [(2*(i-1), 2*(i-1)), (2*(j-1), 2*(j-1))],
    1  : [(2*(i-1), 2*i-1), (2*i -1, 2*(i-1)), (2*j -1, 2*(j-1)), (2*(j-1), 2*j -1)],
    2 : [(2*(j-1), 2*(i-1)), (2*(i-1), 2*(j-1))],
    3 : [(2*j -1, 2*(i-1)), (2*(j-1), 2*i -1), (2*i -1, 2*(j-1)), (2*(i-1), 2*j -1)],
    4 : [(2*i -1, 2*i -1),(2*j -1, 2*j -1)],
    5 : [(2*j -1, 2*i -1), (2*i -1, 2*j-1)]
  }
  
  for valor in posicoes:
      for indice in posicoes[valor]:
        try:
          K[indice[0], indice[1]] += valores[valor]
        except IndexError:
          #As colunas e linhas referentes aos nos 13 e 14 não entram na matriz K
          pass
  
  return K



def ler_arquivoT2(filename):
  # Lê os arquivos de input do segundo teste
  #Retorna a matriz K e a diagonal da matriz de Massas

  with open(filename, 'r') as arquivo:
      n_nos, n_nos_livres, n_barras = map(int, arquivo.readline().strip().split(' '))
      densidade, area, elasticidade = map(float, arquivo.readline().strip().split(' '))
      AE = area*elasticidade*(10**9)
      K = np.zeros((2*n_nos_livres, 2*n_nos_livres))
      M = [0]*(n_nos_livres)
      for linha in arquivo:
        i, j , angulo, tamanho = map(float, linha.strip().split(' '))
        i, j = int(i), int(j)
        K = contribuicaoK(K, i, j, angulo, tamanho, AE)
        massa = tamanho*area*densidade/2
        if i<= n_nos_livres:
          M[i-1] += massa
        if j<= n_nos_livres:
          M[j-1] += massa
  return K, M


def raiz_M(M):
  #Calcula M^(-1/2)

  M_raiz = []
  for elem in M:
    M_raiz += [1/np.sqrt(elem)]*2 
  return M_raiz

def monta_K_til(K,M_raiz):
  #Calcula a matriz K~

  #Multiplica a direita por M^(-1/2)
  for i in range(len(K)):
    for j in range(len(K)):
      K[i,j] = K[i,j]*M_raiz[j]
  
  #Multiplica a esquerda por M^(-1/2)
  for i in range(len(K)):
    for j in range(len(K)):
      K[i,j] = M_raiz[i]*K[i,j]
  
  return K


def trelicas(filename, printa=True):
  #Calcula as 5 menores frequências de vibração da treliça e seus respectivos modos de vibração

  K, M = ler_arquivoT2(filename)
  K_original = K.copy()
  M_raiz = raiz_M(M)
  K = monta_K_til(K, M_raiz)
  autovalores,  modos_de_vibracao = calcula_valores(K)

  # w = raiz(autovalor)
  frequencias = [np.sqrt(autovalor) for autovalor in autovalores]

  # Z = M^(-1/2)*y
  for i in range(len( modos_de_vibracao)):
    for j in range(len( modos_de_vibracao)):
       modos_de_vibracao[i,j] = M_raiz[i]* modos_de_vibracao[i,j]
  
  #Ordena os menores w mantendo seus indices
  indices = {i : frequencias[i] for i in range(len(frequencias))}
  indices = sorted(indices.items(), key = operator.itemgetter(1))

  
  menores_frequencias = []
  menores_modos = np.zeros((24,5))
  for i in range(5):
    menores_frequencias.append(indices[i][1])
    menores_modos[:,i] = modos_de_vibracao[:, indices[i][0]]

    if printa:
      print(f"w{i+1} = {indices[i][1]} (rad/s)")
      print(f'Modos de vibração: {modos_de_vibracao[:, indices[i][0]]}')
      print('---------------------------------------------------------')
    
  return menores_frequencias, menores_modos




def main():
  
  continua = True 
  while  continua:
    print('Selecione um teste')
    print('1) Autovalores e autovetores')
    print('2) Aplicação da treliça')
    teste = int(input('Digite 1 ou 2 para o teste desejado '))

    if teste==1:
      
      print('Selecione um arquivo de input')
      print('1) input-a')
      print('2) input-b')
      print('3) Outro')
      escolha = int(input('Digite 1 ou 2 para arquivo desejado ou 3 para digitar o nome de outro arquivo: '))
      if escolha==3:
        arquivo = input('Digite o nome do arquivo: ')
      else:
        arquivo = 'input-a' if escolha==1 else 'input-b'
      teste1(arquivo)

    elif teste ==2:

      print('Selecione um arquivo de input')
      print('1) input-c')
      print('2) Outro')
      escolha = int(input('Digite 1 para arquivo padrão ou 2 para digitar o nome de outro arquivo: '))

      if escolha == 2:
        arquivo = input('Digite o nome do arquivo: ')
      else:
        arquivo = 'input-c'

      trelicas(arquivo, True)
    else:
      print('Opção inválida')
      continue
    
    continua = input('Deseja continuar(S/N): ').lower() == 's'




if __name__=="__main__":  
    main()








