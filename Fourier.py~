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
plt.savefig("PedrozaJulian_signal.pdf")
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
#===============================================
XX=frec(X)
YY=DFT(Y)
print("NO se usa el paquete fftfreq")
print("La transformada discreta de Fourier se implementa manualmente.")
print("No se usan librerias de python para este proposito.")
#gr=int(len(XX)/2)
XX,YY=XX,(np.real(YY))
#Se elimina la refleccion en x que el metodo tiene
plt.figure(figsize=[14,6])
plt.title("DFT signal.dat")
plt.xlabel("Frecuencia [Hz]")
plt.ylabel("Amplitud")

plt.plot(XX,YY)
plt.savefig("PedrozaJulian_TF.pdf")
plt.close()

def encontrarFmax(xx,yy,n_max=5):
   	 s=sorted(zip(yy,xx))
   	 y,x=map(list,zip(*s))
    	 for i in range(1,n_max+1):
        	print( " frecuencia principal",i," = "+str(round(x[-i],2))+" Hz")
def FiltrarBajo(xx,yy,Fc,enc=False):
    	XT=frec(xx)
    	YT=DFT(yy)
    	if(enc):
		encontrarFmax(XT[XT>0],YT[XT>0])
	cut=np.abs(XT)>Fc
    	YT[cut]=0
    	X_n,Y_n=xx,IDFT(YT)
    	return X_n,Y_n
new_x,new_y=FiltrarBajo(X,Y,1000,enc=True)
plt.figure(figsize=[14,6])
plt.plot(X,Y,alpha=0.8,label="Signal original")
plt.plot(new_x,np.real(new_y),label="Signal filtrada ")
plt.legend()
plt.xlabel("Tiempo [s]")
plt.ylabel("Amplitud")
plt.title("Signal.dat Filtro Pasabajos a 1000 Hz")
plt.savefig("PedrozaJulian_filtrada.pdf")
plt.close()

#********************************************
#=============INCOMPLETOS==================
INC=np.loadtxt("incompletos.dat",delimiter=",").T
X_i,Y_i=np.array(INC[0]),np.array(INC[1])
print("************incompletos.dat*****************")
print("No se puede hacer la transfomrada de fourier en los datos incompletos.dat")
print("puesto que sus datos no se encuentran a intervalos regulares de tiempo.")
print("  ")

x_inter=np.linspace(np.amin(X_i),np.amax(X_i),512)

f1=interp1d(X_i,Y_i,kind="quadratic")
f2=interp1d(X_i,Y_i,kind="cubic")
y_q=f1(x_inter)
y_c=f2(x_inter)

x_frec=frec(x_inter)
Ty_q=DFT(y_q)
Ty_c=DFT(y_c)
#===GRAFICAR TRANSFORMADAS=========
plt.figure(figsize=[12,16])
plt.subplot(3,1,1)
plt.plot(x_frec,np.real(Ty_q),c="g")
plt.xlabel("Frecuencia [Hz]")
plt.ylabel("Amplitud")
plt.title("Interpolacion cuadratica")
plt.subplot(3,1,2)
plt.plot(x_frec,np.real(Ty_c),c="r")
plt.xlabel("Frecuencia [Hz]")
plt.ylabel("Amplitud")
plt.title("Interpolacion cubica")
plt.subplot(3,1,3)
plt.plot(frec(X),np.real(DFT(Y)))
plt.xlabel("Frecuencia [Hz]")
plt.ylabel("Amplitud")
plt.title("Signal")
plt.savefig("PedrozaJulian_TF_interpola.pdf")
plt.close()
#===================================

print("La transformada de fourier para la senial original es mas ruidosa hacia las frecuencias altas")
print("mientras que las interpolaciones disminuyen este ruido")
x_inter,y_q1000=FiltrarBajo(x_inter,y_q,1000)
x_inter,y_q500=FiltrarBajo(x_inter,y_q,500)

x_inter,y_c1000=FiltrarBajo(x_inter,y_c,1000)
x_inter,y_c500=FiltrarBajo(x_inter,y_c,500)


nx,ny500=FiltrarBajo(X,Y,500)
nx,ny1000=FiltrarBajo(X,Y,1000)
#==============GRAFICAR 2 FILTROS========================
plt.figure(figsize=[14,10])
plt.subplot(2,1,1)
plt.title("Frecuencia de corte = 500 Hz")
plt.plot(x_inter,np.real(y_q500),c="g",label="interpolacion cuadratica")
plt.plot(x_inter,np.real(y_c500),c="r",label="interpolacion cubica")
plt.plot(nx,np.real(ny500),c="b",label="signal")
plt.legend()

plt.subplot(2,1,2)
plt.title("Frecuencia de corte = 1000 Hz")
plt.plot(x_inter,np.real(y_q1000),c="g",label="interpolacion cuadratica")
plt.plot(x_inter,np.real(y_c1000),c="r",label="interpolacion cubica")
plt.plot(nx,np.real(ny1000),c="b",label="signal")
plt.legend()
plt.savefig("PedrozaJulian_2Filtros.pdf")
plt.close()














