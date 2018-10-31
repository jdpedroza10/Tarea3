
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
COV_M=COV(datos)
print("Matriz de covarianza asociada a los datos (sin incluir ID number) : ")
print(COV_M)
#******************************************+
#headers=['ID', 'Diagnosis', 'Radius', 'Texture', 'Perimeter', 'Area','Smoothness', 'Compactness', 'Concavity', 'Concave.points','Symmetry', 'Fractal.dimension']

#Autovalores y auto vectores
vals,vecs=linalg.eig(COV_M)
print("Aoutovalores y Autovectores de la matriz de covarianza")
for i in range(len(vals)):
	print("******************************************")
    	print ("El autovalor ",vals[i])
    	print("Esta asociado con el vector ")
    	print(vecs[i])

print("\n Para determinar los parametros mas importantes en un diagnostico 'B' o 'M' es necesario hallar los valores propios mas altos")
print("Y proyectar los datos sobre los autovectores asociados")

#proy=np.dot(vecs,datos)

ii=list((-vals).argsort()[:2]) #Extrayendo los indices de las componentes principales

COMP1=vecs[ii][0]
COMP2=vecs[ii][1]

print("\n Las componentes principales son : ")
print("PC1 = ")
print(COMP1)
print("\nPC2 = ")
print(COMP2)


xb=np.dot(COMP1,datos)[datos[0]==0]
yb=np.dot(COMP2,datos)[datos[0]==0]
xm=np.dot(COMP1,datos)[datos[0]==1]
ym=np.dot(COMP2,datos)[datos[0]==1]

#GRAFICANDO#################################
print("En la figura se ve que claramente hay una linea de separacion entre los dos diagnosticos.")
print("Estos e debe a que la proyeccion sobre los ejes de las componentes principales maximiza")
print("la clasificacion de los datos segun determinadas caracterisitcas. El PCA es una herramienta")
print("Muy util a la hora de establecer relevancia para los factores involucrados")
print(" y hasta para predecir diagnosticos, en este caso.")
plt.figure(figsize=[10,10])
plt.title("Principal Component Analysis",fontsize=18)
plt.scatter(xm,ym,c="r",label="Maligno")
plt.scatter(xb,yb,c="g",label="Benigno")
plt.legend()
plt.xlabel("Componente principal 1",fontsize=14)
plt.ylabel("Componente principal 2",fontsize=14)
plt.savefig("PedrozaJulian_PCA.pdf")
plt.close()





#,p1y,p2x,p2y=np.dot






