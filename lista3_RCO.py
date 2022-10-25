#--------------------Lista 2 - Relações Constitutivas-------------------
#----------------------Rafael Yukio Kono Shimomura----------------------

import numpy as np
import math
from matriz import matriz
from deformacao import deformacao
from tensao import tensao
#------------------------------Exercicio 1-------------------------------
def ex1():
    Nx = 500000
    Ny = 100000
    esf = np.array([[Nx, Ny, 0, 0, 0, 0]])

    num = 3
    E1 = 54.87e9
    E2 = 18.32e9
    v12 = 0.25
    G12 = 8.9e9
    teta1 = 45
    teta2 = 0
    teta3 = 45

    teta1_rad = teta1*math.pi/180
    teta2_rad = teta2*math.pi/180
    teta3_rad = teta3*math.pi/180

    lam_ori = np.array([teta1, teta2, teta3])
    lam_ori_rad = np.array([teta1_rad, teta2_rad, teta3_rad])

    h = np.array([-6e-3,-3e-3,3e-3,6e-3])
    [Q_lam, Qg, Qg_inv] = matriz(E1, E2, G12, v12, lam_ori_rad, h, num)
    [def_global, k_global, def_lamina, def_local] = deformacao(num, h, lam_ori_rad, Qg_inv, esf)
    [tensao_global, tensao_local] = tensao(num, h, lam_ori_rad, k_global, def_global, Q_lam)
    print(tensao_local)

    return(deform_g_1, tensao_kg_1, tensao_kl_1, deform_kl_1, lam_ori)

[deform_g_1, tensao_kg_1, tensao_kl_1, deform_kl_1, lam_ori] = ex1()

'''
File = open("Respostas3.txt",'w')
print('Exercício 1:', file=File)
print('\nSistema Global:', file=File)
print('\nTensão da lâmina %i°: \n' % lam_ori[0],  tensao_kg_1[0], file=File)
print('\nTensão da lâmina %i°: \n' % lam_ori[1], tensao_kg_1[1], file=File)
print('\nTensão da lâmina %i°: \n' % lam_ori[2], tensao_kg_1[2], file=File)
print('\nDeformação do laminado: \n', deform_g_1, file=File)

print('\nSistema Local:', file=File)
print('\nTensão da lâmina %i°: \n' % lam_ori[0],  tensao_kl_1[0], file=File)
print('\nDeformação da lâmina %i°: \n' % lam_ori[0], deform_kl_1[0], file=File)
print('\nTensão da lâmina %i: \n' % lam_ori[1],  tensao_kl_1[1], file=File)
print('\nDeformação da lâmina %i°: \n' % lam_ori[1], deform_kl_1[1], file=File)
print('\nTensão da lâmina %i°: \n' % lam_ori[2],  tensao_kl_1[2], file=File)
print('\nDeformação da lâmina %i°: \n' % lam_ori[2], deform_kl_1[2], file=File)

File.close()
'''