
import numpy as np
import matplotlib.pyplot as plt

SIG=np.loadtxt("signal.dat",delimiter=",")
X = np.array(SIG[0])
Y = np.array(SIG[1])
plt.xlabel("Tiempo [s]")
plt.ylabel("Amplitud")
plt.plot(X,Y)
