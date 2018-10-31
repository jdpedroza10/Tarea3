
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg

datos=pd.read_csv("WDBC.dat",header=None)
_=np.ones(len(datos))
_[datos[1]=="B"]=0
datos[1]=_
datos=np.array(datos).T
datos=datos[1:,:]

#====CALCULANDO LA MATRIZ DE COVARIANZA==========
def cov(x1,x2):
	tam=len(x1)
    	m_x1,m_x2=np.mean(x1),np.mean(x2)
    	s=0
    	for i in range(tam):
        	s=s+((x1[i]-m_x1)*(x2[i]-m_x2))
    	return s/tam
def COV(X):
    	tam=np.shape(X)[0]
    	s=np.zeros([tam,tam])
    	for i in range (tam):
		for j in range (tam):
			s[i,j]=cov(X[i],X[j])
	return s





