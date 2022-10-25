import matplotlib.pyplot as plt
import numpy as np

def plot_tsai_wu(Xt, Yt, S12, Xc, Yc, sigma_12):

    x_inter1 = -(2*(Yc*Yt*(S12**2*Xc*Xt*Yc**2 + S12**2*Xc*Xt*Yt**2 + S12**2*Xc**2*Yc*Yt + S12**2*Xt**2*Yc*Yt + S12**2*Xc*Xt*Yc*Yt + 3*Xc*Xt*Yc*Yt*sigma_12**2
                - S12**2*Xc*Xt**2*Yc*Yt**2*(1/(Xc*Xt*Yc*Yt))**(1/2) - S12**2*Xc*Xt**2*Yc**2*Yt*(1/(Xc*Xt*Yc*Yt))**(1/2) - S12**2*Xc**2*Xt*Yc*Yt**2*(1/(Xc*Xt*Yc*Yt))**(1/2) 
                - S12**2*Xc**2*Xt*Yc**2*Yt*(1/(Xc*Xt*Yc*Yt))**(1/2)))**(1/2) - 2*S12*Xc*Yc*Yt - 2*S12*Xt*Yc*Yt + S12*Xc*Xt*Yc*Yt**2*(1/(Xc*Xt*Yc*Yt))**(1/2) + 
                    S12*Xc*Xt*Yc**2*Yt*(1/(Xc*Xt*Yc*Yt))**(1/2))/(3*S12*Yc*Yt)

    x_inter2 = (2*(Yc*Yt*(S12**2*Xc*Xt*Yc**2 + S12**2*Xc*Xt*Yt**2 + S12**2*Xc**2*Yc*Yt + S12**2*Xt**2*Yc*Yt + S12**2*Xc*Xt*Yc*Yt + 3*Xc*Xt*Yc*Yt*sigma_12**2 
                - S12**2*Xc*Xt**2*Yc*Yt**2*(1/(Xc*Xt*Yc*Yt))**(1/2) - S12**2*Xc*Xt**2*Yc**2*Yt*(1/(Xc*Xt*Yc*Yt))**(1/2) - S12**2*Xc**2*Xt*Yc*Yt**2*(1/(Xc*Xt*Yc*Yt))**(1/2) 
                - S12**2*Xc**2*Xt*Yc**2*Yt*(1/(Xc*Xt*Yc*Yt))**(1/2)))**(1/2) + 2*S12*Xc*Yc*Yt + 2*S12*Xt*Yc*Yt - S12*Xc*Xt*Yc*Yt**2*(1/(Xc*Xt*Yc*Yt))**(1/2) 
                - S12*Xc*Xt*Yc**2*Yt*(1/(Xc*Xt*Yc*Yt))**(1/2))/(3*S12*Yc*Yt)

    sig1 = np.linspace(x_inter1,x_inter2, 1999)
    sig2 = np.zeros(1999)
    sig3 = np.zeros(1999)

    for k in range(0,1999,1):

        sig2[k] = (np.sqrt(abs(4*Xc**2*Xt**2*Yc*Yt*sigma_12**2-3*S12**2*Xc*Xt*Yc*Yt*sig1[k]**2+(((4*S12**2*Xc*Xt**2+4*S12**2*Xc**2*Xt)*Yc-2*S12**2*Xc**2*Xt**2*Yc**2*np.sqrt(abs(1/(Xc*Xt*Yc*Yt))))*Yt
                    -2*S12**2*Xc**2*Xt**2*Yc*np.sqrt(1/(Xc*Xt*Yc*Yt))*Yt**2)*sig1[k]+S12**2*Xc**2*Xt**2*Yt**2-2*S12**2*Xc**2*Xt**2*Yc*Yt+S12**2*Xc**2*Xt**2*Yc**2))
                    -S12*Xc*Xt*Yc*np.sqrt(abs(1/(Xc*Xt*Yc*Yt)))*Yt*sig1[k]+S12*Xc*Xt*Yt+S12*Xc*Xt*Yc)/(2*S12*Xc*Xt)

        sig3[k] = -(np.sqrt(abs(4*Xc**2*Xt**2*Yc*Yt*sigma_12**2-3*S12**2*Xc*Xt*Yc*Yt*sig1[k]**2+(((4*S12**2*Xc*Xt**2+4*S12**2*Xc**2*Xt)*Yc-2*S12**2*Xc**2*Xt**2*Yc**2*np.sqrt(abs(1/(Xc*Xt*Yc*Yt))))*Yt
                  -2*S12**2*Xc**2*Xt**2*Yc*np.sqrt(1/(Xc*Xt*Yc*Yt))*Yt**2)*sig1[k]+S12**2*Xc**2*Xt**2*Yt**2-2*S12**2*Xc**2*Xt**2*Yc*Yt+S12**2*Xc**2*Xt**2*Yc**2))
                  +S12*Xc*Xt*Yc*np.sqrt(abs(1/(Xc*Xt*Yc*Yt)))*Yt*sig1[k]-S12*Xc*Xt*Yt-S12*Xc*Xt*Yc)/(2*S12*Xc*Xt)

    return sig1, sig2, sig3