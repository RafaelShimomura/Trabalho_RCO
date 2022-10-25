import numpy as np

def TsaiHill(sigma, material):
    sigma_1 = sigma[0]
    sigma_2 = sigma[1]
    sigma_12 = sigma[2]
    print(sigma_1)
    Xt = material[4]
    Yt = material[5]
    Xc = material[6]
    Yc = material[7]
    S12 = material[8]

    if sigma_1 > 0 and sigma_2 > 0:
        f = (sigma_1/Xt)**2 + (sigma_2/Yt)**2 - (sigma_1*sigma_2/Xt**2) + (sigma_12/S12)**2
        FS = np.sqrt(f)
    elif sigma_1 > 0 and sigma_2 < 0:
        f = (sigma_1/Xt)**2 + (sigma_2/Yc)**2 + (sigma_1*sigma_2/Xt**2) + (sigma_12/S12)**2
        FS = np.sqrt(f)
    elif sigma_1 < 0 and sigma_2 < 0:
        f = (sigma_1/Xc)**2 + (sigma_2/Yc)**2 - (sigma_1*sigma_2/Xc**2) + (sigma_12/S12)**2
        FS = np.sqrt(f)
    else:
        f = (sigma_1/Xc)**2 + (sigma_2/Yt)**2 + (sigma_1*sigma_2/Xc**2) + (sigma_12/S12)**2
        FS = np.sqrt(f)

    MS_TH = (1/FS) - 1
    return(MS_TH)

def TsaiWu(sigma, material):
    sigma_1 = sigma[0]
    sigma_2 = sigma[1]
    sigma_12 = sigma[2]

    Xt = material[4]
    Yt = material[5]
    Xc = material[6]
    Yc = material[7]
    S12 = material[8]

    F1 = 1/Xt + 1/Xc
    F11 = - 1/(Xt*Xc)
    F2 = 1/Yt + 1/Yc
    F22 = - 1/(Yt*Yc)
    F66 = (1/S12)**2
    F12 = - np.sqrt(F11*F22)/2
    
    A = F11*sigma_1**2 + F22*sigma_2**2 + F66*sigma_12**2 + 2*F12*sigma_1*sigma_2
    B = F1*sigma_1 + F2*sigma_2

    Sf_mais = (-B + np.sqrt(B**2 + 4*A))/(2*A)
    Sf_menos = abs((-B - np.sqrt(B**2 + 4*A))/(2*A))
    MS_mais = Sf_mais - 1
    MS_menos = Sf_menos - 1
    
    return(MS_mais, MS_menos)

def MaxTensao(sigma, material):
    sigma_1 = sigma[0]
    sigma_2 = sigma[1]
    sigma_12 = sigma[2]

    Xt = material[4]
    Yt = material[5]
    Xc = material[6]
    Yc = material[7]
    S12 = material[8]
    
    if sigma_1 > 0:
        MSx_max = Xt/sigma_1
    else:
        MSx_max = Xc/sigma_1
    if sigma_2 > 0:
        MSy_max = Yt/sigma_2
    else:
        MSy_max = Yc/sigma_2
    MSc_max = S12/sigma_12

    MS_max = max([MSx_max, MSy_max, MSc_max])
    
    return(MS_max)