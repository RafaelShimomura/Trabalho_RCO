import numpy as np
import math

def deformacao(num, h, lam_ori_rad, Qg_inv, esf):

    def_global = np.zeros(num)
    def_global = def_global.ravel().tolist()

    def_local = np.zeros(num)
    def_local = def_local.ravel().tolist()

    #-----------------------Deformação Global------------------------
    dk_medio = np.dot(Qg_inv, esf.T)
    #Deformações
    def_medio = np.array([[dk_medio[0]], [dk_medio[1]], [dk_medio[2]]])
    def_medio = np.reshape(def_medio, (3,1))
    #Curvatura
    k_medio = np.array([[dk_medio[3]], [dk_medio[4]], [dk_medio[5]]])
    k_medio = np.reshape(k_medio, (3,1))

    #--------------------Deformação da Lamina---------------------   
    for lamina in range(num):
        def_global[lamina] = def_medio + h[lamina]*k_medio   

    #---------------------------Matriz R----------------------------
    R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    R_inv = np.linalg.inv(R)

    for lamina in range(num):
        m = math.cos(lam_ori_rad[lamina])
        n = math.sin(lam_ori_rad[lamina])
        T = np.array([[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, (m**2 - n**2)]])
        def_local[lamina] = np.linalg.multi_dot([R, T, R_inv, def_global[lamina]])

    return(def_medio, k_medio, def_global, def_local)

    #def plotdeformacao 