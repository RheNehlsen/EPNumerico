from EP2 import trelicas
import numpy as np

import pygame
from pygame.locals import *

from OpenGL.GL import *
from OpenGL.GLU import *

import os
import imageio


def read_file():
    #Le o arquivo de coordenadas 
    vertices = []
    barras =[]
    with open('coordenadas', 'r') as file:
        for _ in range(14):
            coord_ponto = tuple(map(int,file.readline().strip().split(' ')))
            vertices.append([1.5*coord_ponto[0]/15, 1.5*(coord_ponto[1]-20)/40])
        for linha in file:
            barra = tuple(map(lambda x: int(x)-1, linha.strip().split(' ')))
            barras.append(barra)
    return vertices, barras


def calcula_descolacamentos(tempo, n_frequencia):
    # Calcula os deslocamentos dos nós num intervalo de tempo dada uma frequencia
    frequencias, modos = trelicas('input-c', False)
    modo = modos[:, n_frequencia]
    frequencia = frequencias[n_frequencia]
    escala = [1.5, 2, 2,1, 0.8]
    x = []
    for valor in modo:
        deslocamentos = [escala[n_frequencia]*valor*np.cos(t*frequencia) for t in tempo]
        x.append(deslocamentos)
    return x
  
def monta_XY(tempo, n_frequencia):
    #Separa os deslocamentos horizontais e verticais de cada no
    deslocamentos = calcula_descolacamentos(tempo, n_frequencia)
    X =[]
    Y =[]
    for i in range(len(deslocamentos)):
        if i%2==0: 
            X.append(deslocamentos[i])
        else:
            Y.append(deslocamentos[i])
    return X, Y


def desenha_trelica(vertices, barras):
    #Desenha a treliça na tela
    glBegin(GL_LINES)
    for b in barras:
        for no in b:
            glColor3i(0,0,0)
            glVertex2dv(vertices[no])
    glEnd()

def atualiza_vertices(vertices, X,Y,n):
    #Atualiza as coordenadas dos nós
    for i in range(len(X)):
        vertices[i][0] += X[i][n]
        vertices[i][1] += Y[i][n]
    return vertices


def main():
    t0 = float(input('Digite o valor de t0(tempo inicial): '))
    tf = float(input('Digite o valor de tf(tempo final): '))
    passoT = float(input('Digite o valor do passo temporal: '))

    tempo =[]
    t = t0
    while t<= tf:
        tempo.append(t)
        t+= passoT

    n_frequencia = input('De 1 a 5, qual frequência deseja utilizar? (Sendo 1 a menor e 5 a quinta menor): ')
    snaps = input('Deseja salvar as imagens(S/N): ')
    
    if not os.path.exists(r'./snaps'):
        #Cria a pasta auxiliar 
        os.makedirs(r'./snaps')

    if not os.path.exists(r'./gifs'):
        #Cria a pasta auxiliar 
        os.makedirs(r'./gifs')
    
    i =0
    vertices, barras = read_file()
    X, Y = monta_XY(tempo, int(n_frequencia) -1)

    #Configurações pygame
    pygame.init()
    display = (600,600)
    screen = pygame.display.set_mode(display, OPENGL)

    #Configurações OpenGL
    gluPerspective(40, display[0]/display[1], 1, 10)
    glTranslatef(0.0,0.0,-5)
    glRotatef(180, 1, 0, 0 )

    fim = False
    filenames = []
    while not fim:
        
        glClearColor(255, 255,255,1) #fundo branco
        glClear( GL_COLOR_BUFFER_BIT) #Limpa o desenho
        desenha_trelica(vertices, barras) #Desenha novo frame

        
        buffer = glReadPixels(0, 0, *display, GL_RGBA, GL_UNSIGNED_BYTE) #Le o frame desenhado
        pygame.display.flip()
        
        screen_surf = pygame.image.fromstring(buffer , display, "RGBA") #Le o frame desenhado
        
        #Arquivo de cada frame
        filename = f"./snaps/{i}.png"
        filenames.append(filename)

        pygame.image.save(screen_surf, filename)# Salva o arquivo


        vertices = atualiza_vertices(vertices, X, Y, i)# Desloca os nós
        
        i+=1
        if i== len(X[0]):
            pygame.quit()
            fim = True

    #Monta o GIF
    images =[]
    for file in filenames:
        images.append(imageio.imread(file))
    imageio.mimsave(f'./gifs/trelica_frequencia{n_frequencia}.gif', images)
    
    #Apaga os arquivos temporários
    if not snaps.lower()=='s':
        for file in filenames:
            os.remove(file)
    
    print("GIF pronto na pasta /gifs")
main()

