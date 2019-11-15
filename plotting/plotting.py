import matplotlib.pyplot as plt
import numpy as np

def plotMagPhase(field, Lmax=1):
    fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(15,15))
    ax1.imshow(np.square(np.square(np.abs(field))), extent=[-Lmax/2,Lmax/2,-Lmax/2,Lmax/2])
    ax1.set_xlabel('x(mm)')
    ax1.set_ylabel('y(mm)')

    ax2.imshow(np.angle(field),extent=[-Lmax/2,Lmax/2,-Lmax/2,Lmax/2])
    ax2.set_xlabel('x(mm)')
    ax2.set_ylabel('y(mm)')

