# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:19:41 2020

@author: RaulCarreira

 This routine calculates the aerodynamic characteristics of airfoils using 
 the Xfoil software and stores them in a txt file.
 
1) To change the .txt paths:  lines 28,53,54
2) To change the alpha analisys: lines 79,80

"""


import StructCharacteristics as Struct
from xfoil import XFoil
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

archive1 = open("Database/Airfoils995dif.txt","r")
text     = []  
matrix   = [] 
text     = archive1.readlines()    # breaks the lines of the file into vectors 

for i in range(len(text)):          
    matrix.append(text[i].split())

b        = []                      # local variable for airfoil's coordinates y
aa       = []                      # local variable for airfoil's coordinates x
v        = []
airfoils = []                      # stores all generated airfoils


for i in range(len(text)):
    v = matrix[i]
    airfoils.append(v)
    
airfoils = transposeMatrix(airfoils)
j = 0

g = np.shape(airfoils)
Xstruct = np.zeros((g[1],2))

archive2 = open('Training.txt', 'w')
archive4 = open("Database/Parameters995dif.txt","r")

text2     = archive4.readlines() 

try:
    for i in range(int(len(airfoils)/2)):  
        xf = XFoil()
        aa = airfoils[2*j]
        b = airfoils[2*j+1]
        archive3 = open('load.txt', 'w')
       
        for line in range(len(aa)):       
            aux = [aa[line],b[line]] 
            string = str(aux).replace("["," ").replace(","," ").replace("]"," ").replace("'"," ")
            archive3.writelines(string+ '\n')
            aux.clear
       
        # Calculates the aerodynamic coeff using Xfoil
        archive3.close()
        xf.load('load.txt')
        xf.Re = 1e3
        xf.max_iter = 100
        ann, cl, cd, cm, co = xf.aseq(1, 11, 1)
        an = np.linspace(1,10,10)
    
        # Calculates the geometric position of the center of gravity for an airfoil
        Xstruct[:,0] = aa
        Xstruct[:,1] = b
        Surf,Iben,Itor,Xcg,Ycg = Struct.ArfoilProperties(1.0,Xstruct)
        
        # Stores the data in an external file
        if i==0:
            for z in range(len(an)):
                archive2.write( "Cl(A="+str(an[z]).replace("["," ").replace("]"," ")
                               +");")
            for z in range(len(an)):
                archive2.write( "Cd(A="+str(an[z]).replace("["," ").replace("]"," ")
                               +");")
            for z in range(len(an)):
                archive2.write( "Cm(A="+str(an[z]).replace("["," ").replace("]"," ")
                               +");")
            archive2.write("Xcg;xt;yt;xc;yc;Bte;Yle;Ate;dzte;zte;rle;b8;b0;b2;b15;b17"
                           +'\n')
    
        for z in range(len(cl)):
            archive2.write(str(cl[z]).replace("["," ").replace("]"," ")+";")
        for z in range(len(cl)):
            archive2.write(str(cd[z]).replace("["," ").replace("]"," ")+";")
        for z in range(len(cl)):
            archive2.write(str(cm[z]).replace("["," ").replace("]"," ")+";")
        archive2.write(str(Xcg)+";"+str(text2[i])+'\n')
    
        j = j+1 
except:
    print('error')
archive2.close()