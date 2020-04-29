# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:38:11 2020

@author: NG78D39
"""

import numpy as np
import matplotlib.pyplot as plt
import AeroCharacteristics as Aero
import StructCharacteristics as Struct
import time

def CalculateRBFmatrix(Xs,Cs):
    M=np.zeros((Xs.shape[0],Cs.shape[0]))
    sig2=1.5*Xs.shape[1]
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            r2=np.dot(Xs[i,:]-Cs[j,:],Xs[i,:]-Cs[j,:])
            M[i,j]=np.exp(-r2/(2*sig2)) 
    return M

def CalculateHmatrix(Xm,Xs):
    #Xm: volume grid points
    #Xs: surface grid points
    Cs=Xs*1
    M=CalculateRBFmatrix(Xs,Cs)
    invM=np.linalg.inv(M)
    A=CalculateRBFmatrix(Xm,Cs)
    H=np.dot(A,invM)
    return H
    
#1. Define Aero and Structural grid
Nel_struct = 6 #Number of elements per surface in struct grid
Nel_aero = 30 #Number of elements per surface in aero grid
  
Xwing_aero=Aero.NACA4(1.0,0.12,0.03,0.4,Nel_aero) #create aero grid nodes
Xwing_struct=Aero.NACA4(1.0,0.12,0.03,0.4,Nel_struct) #create aero grid nodes

#2.Initialize Aero and Structural Models

#2.a. Aero Model
Elements,Xc,Xn,Xt,DL,A=Aero.InitializeAeroModel(Xwing_aero)
AoA,Vnorm=0.0,6.0
Vinf=np.array([np.cos(AoA/57.3),np.sin(AoA/57.3)])*Vnorm #wind speed
#Vinf=Vinf*Vnorm
Vc=Xc*0.0

#2.b. Struct Model
ispring=7 #structural spring location at node 14
X=np.zeros(6) #struct model state vector
# =============================================================================
X[0:2]=0.5*Xwing_struct[2,:]+0.5*Xwing_struct[8,:] #CG location
# =============================================================================
dt=0.1 #time step
Vels=Xwing_struct*0 #velocities of structural nodes

#3. Create Interpolation Matrix H
H=CalculateHmatrix(Xc,Xwing_struct)
H2=CalculateHmatrix(Xwing_aero,Xwing_struct)
   

#Plot Geometry
plt.figure(figsize=(8,6))
plt.plot(Xwing_aero[:,0],Xwing_aero[:,1]+0.2,'kx',lw=2)
plt.plot(Xwing_struct[:,0],Xwing_struct[:,1]+0.2,'ks',markersize=8,markerfacecolor='w',markeredgewidth=2)
plt.plot(Xwing_aero[:,0],Xwing_aero[:,1]+0.2,'-k',lw=2)
plt.axis('equal')
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.legend(['Aero Grid','Struct Grid'],loc='best',fontsize=10)
plt.grid()


"""
#Test H matrix Regression
Xc2=np.dot(H,Xwing_struct)
plt.figure()
plt.plot(Xc[:,0],Xc[:,1],'-r.')
plt.plot(Xc2[:,0],Xc2[:,1],'-b.')
F2=np.random.rand(Xc.shape[0],2)
Fs=np.dot(H.T,F2)
print(np.sum(F2,axis=0))
print(np.sum(Fs,axis=0))

Xwing2=np.dot(H2,Xwing_struct)
plt.figure()
plt.plot(Xwing_aero[:,0],Xwing_aero[:,1],'-r.')
plt.plot(Xwing2[:,0],Xwing2[:,1],'-b.')
"""



fig,ax=plt.subplots()
ax.set_xlim([-0.2,1.2])
Nsteps=100
Xhist=np.zeros((Nsteps,5))
for i in range(Nsteps):
    #print(i,Vinf)
    #4. Calculate Aero Forces
    #4.a Convert Aero Node Locations and Speeds Using H matrix
    if i>0:
        Xwing_aero=np.dot(H2,Xwing_struct)
        Vc=np.dot(H,Vels)
    
    #4.b. Perform Aero Calculation
    Elements,Xc,Xn,Xt,DL,A=Aero.InitializeAeroModel(Xwing_aero)
    V_el,Cp,CL,CD,CF,CFt=Aero.GetAeroSolution(Elements,Xc,Xn,Xt,DL,A,Vinf,Vc)
    Forces_aero=0.5*1.225*1.0*np.dot(Vinf,Vinf)*CF
    
    #5. Convert Forces to structural Grid
    Forces_struct=np.dot(H.T,Forces_aero) #forces at structural nodes
    Moments=Xwing_struct[:,0]*0 #concentrated moments are zero
    
    Xhist[0:i,0:2]+=Vinf*dt
    Xhist[0:i,2:4]+=Vinf*dt
    Xhist[i,:]=np.array([X[0],X[1],Xwing_struct[0,0],Xwing_struct[0,1],X[4]])
    #Plot Geometry
    ax.clear()
    ax.plot(Xwing_aero[:,0],Xwing_aero[:,1],'-k')
    ax.plot([Xwing_struct[ispring,0]]*2,[Xwing_struct[ispring,1],0.0],'-bo',lw=2)
    ax.plot(X[0],X[1],'kx')
    ax.plot(Xwing_struct[:,0],Xwing_struct[:,1],'k.')
    i1=np.max([i-20,0])
    ax.plot(np.clip(Xhist[0:i+1,0],-0.2,1.6),Xhist[0:i+1,1],'-m',lw=3,alpha=0.5)
    ax.plot(np.clip(Xhist[0:i+1,2],-0.2,1.6),Xhist[0:i+1,3],'-g',lw=3,alpha=0.5)
    
    for ik in range(Xc.shape[0]):
        ax.arrow(Xc[ik,0],Xc[ik,1],Forces_aero[ik,0]*1,Forces_aero[ik,1]*1,fc='r',ec='r')
    
    ax.plot([-0.1,1.1],[-0.6,0.6],'w.')
    
    plt.axis('equal')
    plt.grid()
    plt.draw()
    #plt.savefig(str(i)+'.png',bbox='tight')
    plt.pause(1e-17)
    time.sleep(0.001)
    
    #6. Run Struct Model to Get Displacements and Velocities
    #6.a Add spring force
    Fs,Ms=Struct.GetSpringForcesMoments(Xwing_struct[ispring,:],Vels[ispring,:],X[4],X[5])
    Forces_struct[ispring,:]+=Fs
    Moments[ispring]+=Ms
    #6.b. Get Forces, Moments, Displacements and Velocities
    SF,SM=Struct.GetForceAndMomentatReferencePoint(Xwing_struct,Forces_struct,Moments,X[0:2])
    X,Xwing_struct,Vels=Struct.GetDisplacementsVelocities(Xwing_struct,SF,SM,X,dt)
    
#    plt.plot(Xwing_struct[:,0],Xwing_struct[:,1],'-g.')
#    plt.plot(Xc[:,0],Xc[:,1],'-rx')
#    plt.plot(i*dt,X[1],'r.')
#    plt.plot(i*dt,X[4],'b.')

plt.axis('equal')   
plt.grid()
fig,ax = plt.subplots(nrows=2,ncols=1,sharex=True)



# =============================================================================
ax[0].plot(np.linspace(0,Nsteps*dt,Nsteps),Xhist[:,1],'-k')
ax[0].set_ylabel('CG Y position [m]')
ax[0].grid('on')

ax[1].plot(np.linspace(0,Nsteps*dt,Nsteps),Xhist[:,4]*57.3,'-k')
ax[1].set_ylabel('Theta [deg]')
ax[1].set_xlabel('Time [s]')
ax[1].grid('on')
# =============================================================================
plt.show()





