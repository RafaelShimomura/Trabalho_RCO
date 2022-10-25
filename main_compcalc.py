import numpy as np
from plot_tsai_hill import plot_tsai_hill
from plot_tsai_wu import plot_tsai_wu
from regmist_input import regmist_input
from matriz import matriz
from deformacao import deformacao
import matplotlib.pyplot as plt
from tensao import tensao
from Falha import *
'''
OBJETIVO: Plotar a superfície de falha pelos 3 critérios e mostrar o ponto de cada lâmina. Printar se cada lâmina falhou ou não e o tipo da falha.
    -> Entrada pelas propriedades diretas ou regra das misturas;
    
CONSIDERAÇÕES:
    -> Todas as lâminas do mesmo material e mesma espessura;
'''
#---------------- 1- INÍCIO DOS INPUTS ----------------
# Opções de entrada: 0-> Propriedades diretas do material; 1->Regra das Misturas
ent_opt = 0

if ent_opt == 0:
    E1 = 76e3 # MPa
    E2 = 5.5e3 # MPa
    v12 = 0.25 
    G12 = 2e3 # MPa
    props = np.array([E1, E2, v12, G12])
elif ent_opt == 1:
    Ef = 1
    Em = 1
    Vf = 0.5
    props = regmist_input(Ef, Em, Vf)
else:
    print('Opção Inválida de entrada')
# Carregamentos
Nx = 1000 # N/mm
Ny = 0 # N/mm
Nz = 0 # N/mm
Mx = 0 # N/mm
My = 0 # N/mm
Mz = 0 # N/mm

# Dados de material
Xt = 1380      # Resistência à tração X
Xc = -280     # Resistência à compressão X
Yt = 28        # Resistência à tração Y
Yc = -140      # Resistência à compressão Y
S12 = 55       # Resistência ao cisalhamento no plano 1-2  
pos_lam = [0, 45, -45, -45, 45, 0]
h = 0.5 # mm (espessura de cada lâmina)
#---------------------- FIM DOS INPUTS -----------------------
F = [Nx, Ny, Nz, Mx, My, Mz]
F = np.array(F)
material = np.array([Xt, Xc, Yt, Yc, S12])
pos_lam_rad = np.multiply(np.pi/180, pos_lam)
n_lam = np.size(pos_lam) # número de camadas
h_lam = np.zeros(n_lam+1)
E1 = props[0]
E2 = props[1]
v12 = props[2]
G12 = props[3]

if n_lam % 2 == 0: # se o numero de lâminas for PAR entra aqui
    for i in range(0 ,n_lam+1, 1):
        h_lam.itemset((i), -((n_lam/2)-i)*h)

else: # se o numero de lâminas for IMPAR entra aqui
    for i in range(0,n_lam+1, 1):
        h_lam.itemset((i), -((n_lam/2)-i)*h) 

h_lam_dist = np.zeros(3*n_lam) # h para calcular 3 pontos por lamina para o caso de flexão
j = 0
for i in range(0,np.size(h_lam_dist),3):
    
    h_lam_dist.itemset(i,h_lam[j])
    h_lam_dist.itemset(i+1,(h_lam[j+1]+h_lam[j])/2)
    h_lam_dist.itemset(i+2,h_lam[j+1])
    j = j+1

#---------------------- 2 - CALCULO DAS MATRIZES DE RIGIDEZ -----------------------
Q_lam, ABBD, ABBD_inv = matriz(E1, E2, G12, v12, pos_lam_rad, h_lam, n_lam)

#---------------------- 3 - CALCULO DAS DEFORMAÇÕES E TENSÕES -----------------------
def_global, k_global, def_lamina, def_local= deformacao(n_lam, h_lam_dist, pos_lam_rad, ABBD_inv, F)
tensao_global, tensao_local= tensao(n_lam, h_lam_dist, pos_lam_rad, k_global, def_global, Q_lam)

#---------------------- 4 - PRINT DOS RESULTADOS DAS ETAPAS 2,3-----------------------
print('ABBD :\n', ABBD)
print('ABBD invertida:\n', ABBD_inv)
for fiona in range(0,n_lam,1):
    print('Matriz de rigidez local da lâmina %.0f: \n' %(fiona), Q_lam[fiona])
print('Deformações_global :\n', def_global)
print('Curvatura :\n', k_global)
print('Deformações_lâmina :\n', def_lamina)
print('Deformações_local :\n', def_local)
for i in range(n_lam):
    print('Tensões Globais na lâmina %.0f:\n' %(i+1), tensao_global[i])
    print('Tensões Local %.0f:\n' %(i+1), tensao_local[i])

MS_maxT_v = []
MS_TsaiH_v = []
MS_TsaiW_v = []


for i in range(0,n_lam):
    MS_maxT = MaxTensao(tensao_global[i], material)
    MS_TsaiH = TsaiHill(tensao_global[i], material)
    MS_TsaiW = TsaiWu(tensao_global[i], material)
    MS_maxT_v.append(MS_maxT)
    MS_TsaiH_v.append(MS_TsaiH)
    MS_TsaiW_v.append(MS_TsaiW)

MS_maxT_v = np.array(MS_maxT_v)
MS_TsaiH_v = np.array(MS_TsaiH_v)
MS_TsaiW_v = np.array(MS_TsaiW_v)
minor_maxT = 1
minor_TsaiH = 1
minor_TsaiW = 1
for i in range(0,n_lam-1):
    if MS_maxT_v[i+1]<=MS_maxT_v[i]:
        minor_maxT = i+1
        MS_minor_maxT = MS_maxT_v[i+1]
for i in range(0,n_lam-1):
    if MS_TsaiH_v[i+1]<=MS_TsaiH_v[i]:
        minor_TsaiH = i+1
        MS_minor_TsaiH = MS_TsaiH_v[i+1]
for i in range(0,n_lam-1):
    if MS_TsaiW_v[i+1]<=MS_TsaiW_v[i]:
        minor_TsaiW = i+1
        MS_minor_TsaiW = MS_TsaiW_v[i+1]

print(i)
MS_maxT_v = np.array(MS_maxT_v)
MS_TsaiH_v = np.array(MS_TsaiH_v)
MS_TsaiW_v = np.array(MS_TsaiW_v)

if np.amin(MS_maxT_v) < 0:
        print('Falha pelo critério da máxima tensão na lâmina %.0f, com ângulo de %.1f, tal que MS=%.2f' %(i+1,pos_lam[i],MS_minor_maxT))
if np.amin(MS_TsaiH_v) < 0:
        print('Falha pelo critério de Tsai-Hill na lâmina %.0f, com ângulo de %.1f, tal que MS=%.2f' %(i+1,pos_lam[i],MS_minor_TsaiH))
if np.amin(MS_TsaiW_v) < 0:
        print('Falha pelo critério de Tsai-Wu na lâmina %.0f, com ângulo de %.1f, tal que MS=%.2f' %(i+1,pos_lam[i],MS_minor_TsaiW))

plt.figure(1)
#plt.axis('equal')

#plotando Tensão Máxima
plt.plot([Xt,Xc,Xc,Xt, Xt], [Yt,Yt,Yc,Yc,Yt], label="Máxima Tensão", color='yellow')

#plotando Tsai-Hill
tsai_hill, sigc, sigt = plot_tsai_hill(Xt,Yt,S12,Xc,Yc,tensao_global[minor_TsaiH][2][0])
plt.plot(sigt, tsai_hill[0:499], label = "Tsai-Hill", color='blue')
plt.plot(sigc, tsai_hill[500:999], label = False, color='blue')
plt.plot(sigc, tsai_hill[1000:1499], label = False, color='blue')
plt.plot(sigt, tsai_hill[1500:1999], label = False, color='blue')

#plotando Tsai-Wu
sig1, sig2, sig3 = plot_tsai_wu(Xt, Yt, S12, Xc, Yc, tensao_global[minor_TsaiW][2][0])
plt.plot(sig1, sig2, label = "Tsai-Wu", color='green')
plt.plot(sig1, sig3, label = False, color='green')

# plotando os pontos

plt.plot(tensao_global[minor_maxT][0][0],tensao_global[minor_maxT][1][0], 'bo', color = 'yellow')
plt.plot(tensao_global[minor_TsaiH][0][0],tensao_global[minor_TsaiH][1][0], 'bo', color = 'blue')
plt.plot(tensao_global[minor_TsaiW][0][0],tensao_global[minor_TsaiW][1][0], 'bo', color = 'green')
plt.show()