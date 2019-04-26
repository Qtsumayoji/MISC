import numpy as np
import pylab as plt

if __name__ == '__main__':
    dim = 100
    datax = np.array([i for i in range(dim)])
    datay = np.array([i*i-i for i in range(dim)])
    plt.plot(datax, datay)