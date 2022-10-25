import numpy as np
import math

def tensao(num, h, lam_ori_rad, k_global, def_global, Q_lam):
    
    tensao_global = np.zeros(3*num)
    tensao_global = tensao_global.ravel().tolist()

    tensao_local = np.zeros(3*num)
    tensao_local = tensao_local.ravel().tolist()

#----------------------------Tensão Global-----------------------------
    def_des = np.zeros(3*num)
    def_des = def_des.ravel().tolist()
    des = np.zeros(3*num)
    des = des.ravel().tolist()

    for lamina in range(3*num):
        for j in range(num):
            des[lamina] = h[lamina]*k_global
            def_des[lamina] = def_global - des[lamina]
            tensao_global[lamina] = np.dot(Q_lam[j], def_des[lamina])

#----------------------------Tensão Local-----------------------------
    for lamina in range(3*num):
        for j in range(num):
            m = math.cos(lam_ori_rad[j])
            n = math.sin(lam_ori_rad[j])
            T = np.array([[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, (m**2 - n**2)]])
            T_inv= np.linalg.inv(T)
            tensao_local[lamina] = np.dot(T, tensao_global[lamina])
        
    return(tensao_global, tensao_local)