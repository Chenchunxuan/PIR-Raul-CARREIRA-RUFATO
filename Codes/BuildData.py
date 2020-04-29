# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:19:41 2020

@author: RaulCarreira
"""


"""
 This routine calculates the aerodynamic characteristics of airfoils using 
 the airfoil software and stores them in a text file
"""


import StructCharacteristics as Struct
from xfoil import XFoil

from sklearn import linear_model
import numpy as np

def transposeMatrix(M):
    """Calculates the transpose of a matrix."""
    aux=[]
    for j in range(len(M[0])):
        lin=[]
        for i in range(len(M)):
            lin.append(M[i][j])
        aux.append(lin)
    return aux


archive1 = open("Airfoils.txt","r")
text     = []  
matrix   = [] 
text     = archive1.readlines()    # breaks the lines of the file into vectors 

for i in range(len(text)):          
    matrix.append(text[i].split())

b        = []                      # local variable for airfoil's coordinates y
a        = []                      # local variable for airfoil's coordinates x
v        = []
airfoils = []                      # stores all generated airfoils


for i in range(len(text)):
    v = matrix[i]
    airfoils.append(v)
    
airfoils = transposeMatrix(airfoils)
j = 0

g = np.shape(airfoils)

Xstruct = np.zeros((g[1],2))

archive2 = open('CgCla.txt', 'w')
archive2.write(" Clalpha"+"            "+"Xcg"+"                 "+"Ycg"+'\n')
archive2.write("-----------------------------------------------------"+'\n')

for i in range(int(len(airfoils)/2)):  
    xf = XFoil()
    a = airfoils[2*j]
    b = airfoils[2*j+1]
    archive3 = open('load.txt', 'w')
   
    for line in range(len(a)):       
        aux = [a[line],b[line]] 
        string = str(aux).replace("["," ").replace(","," ").replace("]"," ").replace("'"," ")
        archive3.writelines(string+ '\n')
        aux.clear
   
    archive3.close()
    xf.load('load.txt')
    an,cl,cd,cm,co = xf.aseq(0, 6, 0.5)
    
    lm = linear_model.LinearRegression() # Linear regression to take Clalpha
    model = lm.fit(an.reshape(12,1),cl)
    clalpha = lm.coef_
    
    Xstruct[:,0] = a
    Xstruct[:,1] = b
    Surf,Iben,Itor,Xcg,Ycg = Struct.ArfoilProperties(1.0,Xstruct)
    archive2.write(str(clalpha).replace("["," ").replace("]"," ")+" "+str(Xcg)+" "+str(Ycg)+'\n')
    j =j+1 
archive2.close()

    
    
    
    