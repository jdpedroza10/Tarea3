#!/bin/bash
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
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
	else:
        	a=np.arange(0,int((tam-1)/2))
       		b=np.arange(-int((tam-1)/2),0)
        	f=np.concatenate((a,b))
        	return f/(sp*tam)
def DFT(y):
	tam=len(y)
	Four=np.zeros(tam,dtype=np.complex_)
	for k in range(tam):
		s=[]
		for n in range(tam):
			s.append(y[n]*np.exp(-2*np.pi*complex(0,1)*k*n/tam))
		Four[k]=np.sum(s)
	return Four
def IDFT(y):
	tam=len(y)
	Four=np.zeros(tam,dtype=np.complex_) #T fourier
	for n in range(tam):
        	s=[]
		for k in range(tam):
            		s.append(y[k]*np.exp(2*np.pi*complex(0,1)*k*n/tam))
		Four[n]=np.sum(s)/tam
	return Four
def IDFT(y):
	tam=len(y)
    	Four=np.zeros(tam,dtype=np.complex_) #T fourier
	for n in range(tam):
        	s=[]
		for k in range(tam):
            		s.append(y[k]*np.exp(2*np.pi*complex(0,1)*k*n/tam))
		Four[n]=np.sum(s)/tam
	return Four
##=======Comprobar errores de DFT a IDFT========
def prueba():
	YT=DFT(Y)
    	Y_b=IDFT(YT)
    	Y_r=np.sqrt(np.mean((Y-np.real(Y_b))**2))
    	print("errores medios = ",Y_r)

















