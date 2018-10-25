
import numpy as np
import matplotlib.pyplot as plt

# Almacene los datos de signal.dat y de incompletos.dat

print("************signal.dat*****************")
SIG=np.loadtxt("signal.dat",delimiter=",").T
X = np.array(SIG[0])
Y = np.array(SIG[1])


#Haga una grafica de los datos de signal.dat y guarde dicha grafica sin mostrarla en ApellidoNombre_signal.pdf
plt.figure(figsize=[14,6])
plt.xlabel("Tiempo [s]")
plt.ylabel("Amplitud")
plt.plot(X,Y)
plt.savefig("JulianPedroza_signal.pdf")
plt.close()


#=========IMPLEMENTANDO TRANSFORMACIONES===========
def frec(x):
	tam=len(x)
	sp=x[1]-x[0]
	if(tam%2==0):
        	a=np.arange(0,int(tam/2))
        	b=np.arange(-int(tam/2),0)
       		f=np.concatenate((a,b))
        	return f/(sp*tam)

