import numpy as np
import matplotlib.pyplot as plt

def load_airfoil(filename): #filename is airfoil number as STRING
    loc = 'airfoils/'+filename+'.dat'
    data = np.genfromtxt(fname=loc,skip_header=1)
    x = data[:,0]
    y = data[:,1]
    Np = len(x)
    return x, y, Np









if __name__ is '__main__':
    x,y, Np = load_airfoil('0008')
    plt.plot(x,y)
    plt.show()