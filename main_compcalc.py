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
    E1 = 54.87e3 # MPa
    E2 = 18.32e3 # MPa
    v12 = 0.25 
    G12 = 8.9e3 # MPa
    props = np.array([E1, E2, v12, G12])
elif ent_opt == 1:
    Ef = 1
    Em = 1
    Vf = 0.5
    props = regmist_input(Ef, Em, Vf)
else:
    print('Opção Inválida de entrada')
# Carregamentos
Nx = 500 # N/mm
Ny = 100 # N/mm
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
pos_lam = [45, 0, 45]
h = 0.5 # mm (espessura de cada lâmina)
#---------------------- FIM DOS INPUTS -----------------------
F = [Nx, Ny, Nz, Mx, My, Mz]
F = np.array(F)
pos_lam_rad = np.multiply(np.pi/180, pos_lam)
n_lam = np.size(pos_lam) # número de camadas
h_lam = np.zeros(n_lam+1)
E1 = props[0]
E2 = props[1]
v12 = props[2]
G12 = props[3]

if n_lam % 2 == 0: # se o numero de lâminas for PAR entra aqui
    for i in range(0 ,n_lam+1, 1):
        h_lam.itemset((i), +((n_lam/2)-i)*h)

else: # se o numero de lâminas for IMPAR entra aqui
    for i in range(0,n_lam+1, 1):
        h_lam.itemset((i), +((n_lam/2)-i)*h) 
h_lam = np.array([-6, -3, 3, 6])
print(h_lam)
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
for i in range(3*n_lam):
    print('Tensões Globais na lâmina %.0f:\n' %(i), tensao_global[i])
    #print('Tensões Local %.0f:\n' %(i), tensao_local[i])

sigma_12 = 0 #isso deve ser retirado, é so teste


plt.figure(1)
#plt.axis('equal')

#plotando Tensão Máxima
plt.plot([Xt,Xc,Xc,Xt, Xt], [Yt,Yt,Yc,Yc,Yt], label="Máxima Tensão", color='yellow')

#plotando Tsai-Hill
tsai_hill, sigc, sigt = plot_tsai_hill(Xt,Yt,S12,Xc,Yc,sigma_12)
plt.plot(sigt, tsai_hill[0:499], label = "Tsai-Hill", color='blue')
plt.plot(sigc, tsai_hill[500:999], label = False, color='blue')
plt.plot(sigc, tsai_hill[1000:1499], label = False, color='blue')
plt.plot(sigt, tsai_hill[1500:1999], label = False, color='blue')

#plotando Tsai-Wu
sig1, sig2, sig3 = plot_tsai_wu(Xt, Yt, S12, Xc, Yc, sigma_12)
plt.plot(sig1, sig2, label = "Tsai-Wu", color='purple')
plt.plot(sig1, sig3, label = False, color='purple')

# plotando a ditribuição de tensões e deformações
dist_sigma1 = np.zeros(3*n_lam)
dist_sigma2 = np.zeros(3*n_lam)
dist_sigma12 = np.zeros(3*n_lam)
dist_def1 = np.zeros(3*n_lam)
dist_def2 = np.zeros(3*n_lam)
dist_def12 = np.zeros(3*n_lam)

for i in range(0,np.size(h_lam_dist),1):
    dist_sigma1.itemset(i,tensao_global[i][0])
    dist_sigma2.itemset(i,tensao_global[i][1])
    dist_sigma12.itemset(i,tensao_global[i][2])
    dist_def1.itemset(i,def_lamina[i][0])
    dist_def2.itemset(i,def_lamina[i][0])
    dist_def12.itemset(i,def_lamina[i][0])

fig, axs = plt.subplots(2, 3)
axs[0,0].plot(dist_sigma1,h_lam_dist)
axs[0,1].plot(dist_sigma2,h_lam_dist)
axs[0,2].plot(dist_sigma12,h_lam_dist)
axs[1,0].plot(dist_def1,h_lam_dist)
axs[1,1].plot(dist_def2,h_lam_dist)
axs[1,2].plot(dist_def12,h_lam_dist)
plt.show()