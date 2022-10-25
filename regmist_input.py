import numpy as np

def regmist_input(Ef, Em, Vf, vf, vm, Gf, Gm): # (Mod. elasticidade; Mod Elasticidade, Fração volumética de fibra; Poisson fibra; Poisson Matriz)
    Vm = 1 - Vf
    E1 = Ef*Vf + Em*Vm
    E2 = 1/((Vf/Ef)+(Vm/Em))
    v12 = vf*Vf + vm*Vm
    v21 = v12*E2/E1
    G12 = 1/((Vf/Gf)+(Vm/Gm))
    inputs = np.array([E1, E2, v12, G12])
    return inputs