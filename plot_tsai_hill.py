import matplotlib.pyplot as plt
import numpy as np

def plot_tsai_hill(Xt,Yt,S12,Xc,Yc,sigma_12):

    plot1 = np.zeros(2000)
    sigc = np.linspace(Xc,0,499)
    sigt = np.linspace(0,Xt,499)

    for k in range(0,499,1):

        plot1[k] = (Yt*np.sqrt(-4*Xt**4*sigma_12**2+(S12**2*Yt**2-4*S12**2*Xt**2)*sigt[k]**2+4*S12**2*Xt**4)-S12*Yt**2*sigt[k])/(2*S12*Xt**2)
        plot1[k+500] = (Yt*np.sqrt(-4*Xc**4*sigma_12**2+(S12**2*Yt**2-4*S12**2*Xc**2)*sigc[k]**2+4*S12**2*Xc**4)+S12*Yt**2*sigc[k])/(2*S12*Xc**2)
        plot1[k+1000] = (Yc*np.sqrt(-4*Xc**4*sigma_12**2+(S12**2*Yc**2-4*S12**2*Xc**2)*sigc[k]**2+4*S12**2*Xc**4)-S12*Yc**2*sigc[k])/(2*S12*Xc**2)
        plot1[k+1500] = (Yc*np.sqrt(-4*Xt**4*sigma_12**2+(S12**2*Yc**2-4*S12**2*Xt**2)*sigt[k]**2+4*S12**2*Xt**4)+S12*Yc**2*sigt[k])/(2*S12*Xt**2)
    

    return plot1, sigc, sigt