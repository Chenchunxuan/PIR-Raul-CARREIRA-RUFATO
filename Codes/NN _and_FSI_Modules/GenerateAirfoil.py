# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 17:08:13 2020

@author: RaulCarreira

    This routine generates the coordinates of an airfoil using the 
    BP3434 parametrization.
    
"""
import numpy as np
import matplotlib.pyplot as plt
import math as m
import mpmath as mp

def subs(m1, m2):
    #Calculates the term-to-term difference m1-m2 of matrices with the same size
    matrix_subs = []
    # Assuming the two arrays are the same size
    number_columns = len(m1) # Counts how many elements are in a row

    #Creates a new row in the sum_array
    for j in range(number_columns):
        # Somando os elementos que possuem o mesmo índice
        summ = m1[j]- m2[j]
        matrix_subs.append(summ)
    return matrix_subs

def CurvePoints(xt,yt,xc,yc,Bte,Yle,Ate,dzte,zte,rle,b8,b0,b2,b15,b17):
    # Calculation of bezier control points 3434 parametrization
    # Leading edge thickness curve [LET]
    # LET : Matrix 4 columns [points] by 2 rows [x, y] of bezier points
    LET = np.zeros((2,4))
    LET[0,2] = -3*b8**2/(2*rle)
    LET[0,3] = xt
    LET[1,1] = b8
    LET[1,2] = yt
    LET[1,3] = yt
    
    # Leading edge camber curve [LEC]
    # LEC: Matrix 4 columns [points] by 2 rows [x, y] of bezier points
    LEC = np.zeros((2,4))
    LEC[0,1] = b0
    LEC[0,2] = b2
    LEC[0,3] = xc
    LEC[1,1] = b0*m.tan(Yle)
    LEC[1,2] = yc
    LEC[1,3] = yc
    
    # Trailing edge thickness curve [TET]
    # TET: Matrix 5 columns [points] by 2 rows [x, y] of bezier points
    TET = np.zeros((2,5))
    TET[0,0] = xt
    TET[0,1] = (7*xt+9*b8**2/(2*rle))/4
    TET[0,2] = 3*xt-5*LET[0,2]/2
    TET[0,3] = b15
    TET[0,4] = 1
    TET[1,0] = yt
    TET[1,1] = yt
    TET[1,2] = (yt+b8)/2
    TET[1,3] = dzte+(1-b15)*m.tan(Bte)
    TET[1,4] = dzte

    # Trailing edge camber curve [TEC]
    # TEC: Matrix 5 columns [points] by 2 rows [x, y] of bezier points
    TEC = np.zeros((2,5))
    TEC[0,0] = xc
    TEC[0,1] = (3*xc-yc*mp.cot(Yle))/2
    TEC[0,2] = (-8*yc*mp.cot(Yle)+13*xc)/6
    TEC[0,3] = b17
    TEC[0,4] = 1
    TEC[1,0] = yc
    TEC[1,1] = yc
    TEC[1,2] = 5*yc/6
    TEC[1,3] = zte+(1-b17)*m.tan(Ate)
    TEC[1,4] = zte
    
    return LET,LEC,TET,TEC

def BuildAirfoil(Vector):
    # Using BP3434 control points, generates airfoil curve
    xt = Vector[0]
    yt = Vector[1]
    xc = Vector[2]
    yc = Vector[3]
    Bte = Vector[4] 
    Yle = Vector[5] 
    Ate = Vector[6] 
    dzte =Vector[7] 
    zte = Vector[8]
    rle = Vector[9]
    b8 = Vector[10]
    b0 = Vector[11]
    b2 = Vector[12]
    b15 = Vector[13]
    b17 = Vector[14]
    
    
    LET,LEC,TET,TEC = CurvePoints(xt,yt,xc,yc,Bte,Yle,Ate,dzte,zte,rle,b8,b0,
                                  b2,b15,b17)
    
    
    """ Generating 4 curves to describe the airfoil: 2 third degree curves to 
        describe the leading edge and 2 fourth degree curves to describe the
        trailling edge
    """
    
    Plet=[-1*LET[:,0]+3*LET[:,1]-3*LET[:,2]+1*LET[:,3],3*LET[:,0]+-6*LET[:,1]
          +3*LET[:,2],-3*LET[:,0]+3*LET[:,1],LET[:,0]]
   
    Plec=[-1*LEC[:,0]+3*LEC[:,1]-3*LEC[:,2]+1*LEC[:,3],3*LEC[:,0]+-6*LEC[:,1]
          +3*LEC[:,2],-3*LEC[:,0]+3*LEC[:,1],LEC[:,0]]
   
    Ptet=[1*TET[:,0]-4*TET[:,1]+6*TET[:,2]-4*TET[:,3]+1*TET[:,4],-4*TET[:,0]
          +12*TET[:,1]-12*TET[:,2]+4*TET[:,3],6*TET[:,0]-12*TET[:,1]+6*TET[:,2],
          -4*TET[:,0]+4*TET[:,1],1*TET[:,0]]
    
    Ptec=[1*TEC[:,0]-4*TEC[:,1]+6*TEC[:,2]-4*TEC[:,3]+1*TEC[:,4],-4*TEC[:,0]
          +12*TEC[:,1]-12*TEC[:,2]+4*TEC[:,3],6*TEC[:,0]-12*TEC[:,1]+6*TEC[:,2],
          -4*TEC[:,0]+4*TEC[:,1],1*TEC[:,0]]
    
    PXlet,PYlet = [Plet[0][0],Plet[1][0],Plet[2][0],Plet[3][0]],[Plet[0][1],Plet[1][1],Plet[2][1],Plet[3][1]]
   
    PXlec,PYlec = [Plec[0][0],Plec[1][0],Plec[2][0],Plec[3][0]],[Plec[0][1],Plec[1][1],Plec[2][1],Plec[3][1]]
   
    PXtet,PYtet = [Ptet[0][0],Ptet[1][0],Ptet[2][0],Ptet[3][0],Ptet[4][0]],[Ptet[0][1],Ptet[1][1],Ptet[2][1],Ptet[3][1],Ptet[4][1]]
   
    PXtec,PYtec = [Ptec[0][0],Ptec[1][0],Ptec[2][0],Ptec[3][0],Ptec[4][0]],[Ptec[0][1],Ptec[1][1],Ptec[2][1],Ptec[3][1],Ptec[4][1]]
    
    n = 300     # Desired number of discretization points 
    X  = np.zeros(n)
    X[0:int(n/2)] = np.linspace(1,0,int(n/2))**2
    X[int(n/2)::] = np.linspace(0,1,int(n/2))**2
    a = 0
    n = len(X)
    T = np.zeros(n)
    C = np.zeros(n)
    Yn = np.zeros(n) 
    
    for i in range (n):
        # thickness
        if X[i]<=xt: # leading edge
            A = np.roots(subs(PXlet, [0,0,0,X[i]]))
            for j in range (3):
                if A[j].imag==0 and A[j].real<=1 and A[j].real>=0:
                    a=A[j]
            T[i]= abs(np.polyval(PYlet,a))
            if i==0:
                T[i]=0
        else: # trailing edge
            A = np.roots(subs(PXtet, [0,0,0,0,X[i]]))
            for j in range (4):
                if A[j].imag==0 and A[j].real<=1 and A[j].real>=0:
                    a = A[j].real        
            T[i] = abs(np.polyval(PYtet,a))
            if i==0:
                T[i]=0
        # camber
        if X[i]<=xc: # leading edge
            A = np.roots(subs(PXlec, [0,0,0,X[i]]))
            for j in range (3):
                if A[j].imag==0.0 and A[j].real<=1.0 and A[j].real>=0.0:
                    a=A[j]
            C[i] = abs(np.polyval(PYlec,a))
            if i==0:
                C[i]=0
        else: # trailing edge
            A = np.roots(subs(PXtec, [0,0,0,0,X[i]]))
            for j in range (4):
                if A[j].imag==0.0 and A[j].real<=1.0 and A[j].real>=0.0:
                    a=A[j]
            C[i] = abs(np.polyval(PYtec,a))
            if i==0:
                C[i]=0
        if i<=n/2:
            Yn[i]=C[i]+T[i]/2
        else:
            Yn[i]=C[i]-T[i]/2
    return X,Yn


"""Small exemple to compare 2 airfoils with their respectives BP3434 parameters"""
# Vector1 = [ 3.04399524e-01,  1.53153153e-01,  4.88588589e-01,
#          5.79729730e-02,  2.00000000e-16,  3.11855800e-01,
#          1.41371669e+00,  3.21621622e-03,  8.04054054e-04,
#         -2.56546453e-02,  7.37611157e-03,  0.00000000e+00,
#          2.20185991e-01,  9.27584445e-01,  1.00000000e+00]


# Vector2 = [ 0.30403334,  0.14194477,  0.47115624,  0.05529083,  0.09542751,
#         0.26981494,  0.5853721 ,  0.00344465, -0.00338946, -0.02265443,
#         0.01030791,  0.02749763,  0.24010088,  0.8920514 ,  0.9746181 ]
# X1,Yn1 = BuildAirfoil(Vector1)
# X2,Yn2 = BuildAirfoil(Vector2)
# plt.plot(X1, Yn1,label='Real') #Criando o gráfico
# plt.plot(X2,Yn2,label='Prediction') #Criando o gráfico
# plt.title('Airfoil') #adicionando o título
# plt.xlabel('x/c')
# plt.legend()
# plt.ylabel('y')
# plt.axis('equal')
# plt.grid('on')
# plt.show()



