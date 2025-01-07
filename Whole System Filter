#Filter Whole System
import numpy as np
import math as m
import cmath
import matplotlib.pyplot as plt

#Define Parameters
air = 1.0
n1 = 1.46 #SiO2 Average of 4 / "Thin Film"
n2 = 2.32 #TiO2 Average of 5 / "Thin Film"
n3 = 1.5234
d1 = 656.279/4/n1
d2 = 656.279/4/n2 #target wavelength (hydrogen alpha) in 3 decimal
d3 = 2000 #Soda lime substrate

d21 = 458.3/4/n1
d22 = 458.3/4/n2

d31 = 522.8/4/n1
d32 = 522.8/4/n2

#Define boundary/propagation calcultion functions
def Boundary(n1, n2):
    return np.array([[n2+n1, n2-n1], [n2-n1, n2+n1]])/(2*n2)

def Propagation(n, d, w):
    k = 2*m.pi/w
    p = n*d*k
    return np.array([[cmath.exp(-1j*p), 0], [0, cmath.exp(1j*p)]])

#Function that calculates transmittance/reflection for inputted wavelength w

def Calculate(w):
    
    # Definition Part
    P_Air = Propagation(air, 1000, w)
    P_Air2 = Propagation(air, 1000, w)
    B_AtoL = Boundary(air, n1)
    P_L = Propagation(n1, d1, w)
    B_LtoH = Boundary(n1, n2)
    P_H = Propagation(n2, d2, w)
    B_HtoL = Boundary(n2, n1)
    B_HtoA = Boundary(n2, air)
    #P_Substrate = Propagation(n3, d3, w)
    #B_HtoSub = Boundary(n2, n3)
    #B_SubtoL = Boundary(n3, n1)
    P_L2 = Propagation(n1, d21, w)
    P_H2 = Propagation(n2, d22, w)
    P_L3 = Propagation(n1, d31, w)
    P_H3 = Propagation(n2, d32, w)


    #First layer
    TMM = np.matmul(B_HtoA, np.linalg.pinv(B_HtoL))
    
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
    
    TMM = np.matmul(TMM, sub)
    TMM = np.matmul(TMM, B_AtoL)
    
    #between first/second layer
    TMM = np.matmul(TMM, P_Air)
    
    #second layer
    TMM2 = np.matmul(B_HtoA, np.linalg.pinv(B_HtoL))
    
    temp = np.matmul(B_HtoL, P_H2)
    temp = np.matmul(temp, B_LtoH)
    temp = np.matmul(temp, P_L2)
    temp = np.linalg.matrix_power(temp, 15)
    
    TMM2 = np.matmul(TMM2, temp)
    TMM2 = np.matmul(TMM2, B_AtoL)
    
    TMM = np.matmul(TMM, TMM2)
    
    #between second/third layer
    TMM = np.matmul(TMM, P_Air2)
    
    #third layer
    TMM3 = np.matmul(B_HtoA, np.linalg.pinv(B_HtoL))
    
    temp2 = np.matmul(B_HtoL, P_H3)
    temp2 = np.matmul(temp2, B_LtoH)
    temp2 = np.matmul(temp2, P_L3)
    temp2 = np.linalg.matrix_power(temp2, 15)
    
    TMM3 = np.matmul(TMM3, temp2)
    TMM3 = np.matmul(TMM3, B_AtoL)
    
    TMM = np.matmul(TMM, TMM3)
    
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

#For zoomed in graph at peak at H-alpha
#i = 570
#for _ in range (21000):
    #precision.append(i)
    #i +=0.01

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
plt.title("Filter Transmittance for Whole System")
plt.gca().legend(('Transmittance', 'Reflecttance'))
plt.ylim(-0.1, 1.1)
plt.show()
