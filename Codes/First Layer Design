#First layer Design

import numpy as np
import math as m
import cmath
import matplotlib.pyplot as plt

#Defind parameter values
air = 1.0
n1 = 1.46 #SiO2 Average of 4 "Thin Film" values 
n2 = 2.32 #TiO2 Average of 5 "Thin Film" values 
tar = 656.279 #target wavelength (hydrogen alpha) up to 3 decimal paces
d1 = tar/4/n1 
d2 = tar/4/n2 

#Define boundary/propagation calcultion functions
def Boundary(n1, n2):
    return np.array([[n2+n1, n2-n1], [n2-n1, n2+n1]])/(2*n2)

def Propagation(n, d, w):
    k = 2*m.pi/w
    p = n*d*k
    return np.array([[cmath.exp(-1j*p), 0], [0, cmath.exp(1j*p)]])


#Function that calculates transmittance/reflection for inputted wavelength w
def Calculate(w):
    
    # Definition
    B_AtoL = Boundary(air, n1) # air to first material (SiO2)
    P_L = Propagation(n1, d1, w) # Propagation through first material
    B_LtoH = Boundary(n1, n2) # first to second material (TiO2)
    P_H = Propagation(n2, d2, w) # Propagation through second material
    B_HtoL = Boundary(n2, n1) # second to first material
    B_HtoA = Boundary(n2, air) # second material to air


    #Details of filter design
    # (LH)^4 / H / (LH)^14 / H / (LH)^5
    TMM = np.matmul(B_HtoA, np.linalg.pinv(B_HtoL)) #Transfer Matrix Method (TMM) calculation 
    
    sub = np.matmul(B_HtoL, P_H)
    sub = np.matmul(sub, B_LtoH)
    sub = np.matmul(sub, P_L)

    sub2 = np.linalg.matrix_power(sub, 4)
    sub1 = np.linalg.matrix_power(sub, 14)
    sub = np.linalg.matrix_power(sub, 5)
    
    experim = np.matmul(sub, sub2)
    experim = np.matmul(experim, P_H)
    experim = np.matmul(experim, sub1)
    experim = np.matmul(experim, P_H)
    
    sub = np.matmul(experim, sub)
    
    #TMM = TMM*sub*B_AtoL
    TMM = np.matmul(TMM, sub)
    TMM = np.matmul(TMM, B_AtoL)
    
    l =[]
    
    # Transmittance Calculation
    t = 1/TMM[1,1]
    T = abs(t)*abs(t)
    l.append(T)
    
    # Reflection Calculation
    r = t*TMM[0, 1]
    R = abs(r)*abs(r)
    l.append(R)
    
    return l


#Creating list with wavelength values - from 400 to 800 with 0.01 increments 
x = []
y = []
y2 = []
precision = []
i = 400
for _ in range (40001):
    precision.append(i)
    i +=0.01

#Calculating transmittance/reflectance for each wavelengths
for i in precision:  #range (400, 801): 
    l = []
    l = Calculate(i)
    
    x.append(i)
    y.append(l[0])
    y2.append(l[1])

X = np.arange(400, 801)

#Plot transmittance (+reflection) graph for this design
plt.plot(x, y) #Transmittance
#plt.plot(x, y2) #Reflection
plt.xlabel("Wavelength  (nm)")
plt.ylabel("Spectral Power ")
plt.title("Filter Transmittance for First Layer")
plt.gca().legend(('Transmittance', 'Reflecttance'))
plt.ylim(-0.1, 1.1)
plt.show()
