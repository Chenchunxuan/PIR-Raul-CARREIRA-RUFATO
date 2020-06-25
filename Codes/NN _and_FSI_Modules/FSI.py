# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:38:11 2020

@author: NG78D39
"""
import NeuralNetwork as NN
import numpy as np
import matplotlib.pyplot as plt
import AeroCharacteristics as Aero
import StructCharacteristics as Struct
from PIL import Image
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from xfoil import XFoil

def func(x, a, b):
     return a*3.14*np.exp(-b*x)

def Verify_instability(x,y):
    p, pcov = curve_fit(func, x, y)
    return p



def Determine_Airfoil(desired,n):
    X,Y = NN.Find_Airfoil(desired,n)
    return X,Y

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

def fig2data ( fig ):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA
    channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )
 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.tostring_argb(), dtype=np.uint8 )
    buf.shape = ( w, h,4 )
 
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel 
    # to have it in RGBA mode
    buf = np.roll ( buf, 3, axis = 2 )
    return buf

def fig2img (fig):
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data (fig)
    w, h, d = buf.shape
    return Image.frombytes( "RGBA", ( w ,h ), buf.tostring( ) )
    
#1. Define Aero and Structural grid
Nel_struct = 6 # Number of elements per surface in struct grid
Nel_aero   = 30 # Number of elements per surface in aero grid
  
Xwing_aero=Aero.NACA4(1.0,0.12,0.03,0.4,Nel_aero) #create aero grid nodes
Xwing_struct=Aero.NACA4(1.0,0.12,0.03,0.4,Nel_struct) #create aero grid nodes


Xwing_aero1=Aero.NACA4(1.0,0.13,0.03,0.3,Nel_aero) #create aero grid nodes
Xwing_struct1=Aero.NACA4(1.0,0.13,0.03,0.3,Nel_struct) #create aero grid nodes
xf = XFoil()
aa = Xwing_aero1[:,0]
b = Xwing_aero1[:,1]
archive3 = open('lt.txt', 'w')
   
for line in range(len(aa)):       
    aux = [aa[line],b[line]] 
    string = str(aux).replace("["," ").replace(","," ").replace("]"," ").replace("'"," ")
    archive3.writelines(string+ '\n')
    aux.clear
   
# Calculates the aerodynamic coeff using Xfoil
archive3.close()
xf.load('lt.txt')
xf.Re = 1e3
xf.max_iter = 100
ann, cl, cd, cm, co = xf.aseq(1, 11, 1)
an = np.linspace(1,10,10)
desired = np.zeros(31)
desired[0:10] = cl
desired[10:20] = cd
desired[20:30] = cm
Surf,Iben,Itor,Xcg,Ycg = Struct.ArfoilProperties(1.0,Xwing_struct1)
desired[30] = Xcg
#2.Initialize Aero and Structural Models

#2.a. Aero Model
Elements,Xc,Xn,Xt,DL,A=Aero.InitializeAeroModel(Xwing_aero)
AoA,Vnorm=0.0,6.0
Vinf=np.array([np.cos(AoA/57.3),np.sin(AoA/57.3)])*Vnorm #wind speed
#Vinf=Vinf*Vnorm
Vc=Xc*0.0

# =============================================================================
#2.b. Struct Model
Surf,Iben,Itor,Xcg,Ycg = Struct.ArfoilProperties(1.0,Xwing_struct)
X = np.zeros(8) #struct model state vector
X[0] = 0.5
X[1] = 0.0297
Xcg = 0.3
Xec =  X[0]
x_theta = (Xcg-Xec)
mass = 7.0
I_theta = 7.0 + mass*x_theta**2    ## I_t = I_0 + mx_t**2
# =============================================================================
dt = 0.1 #time step
Vels=Xwing_struct*0 #velocities of structural nodes

#3. Create Interpolation Matrix H
H=CalculateHmatrix(Xc,Xwing_struct)
H2=CalculateHmatrix(Xwing_aero,Xwing_struct)
   

#Plot Geometry
plt.figure(figsize=(8,6))
plt.plot(Xwing_aero[:,0],Xwing_aero[:,1]+0.2,'kx',lw=2)
plt.plot(Xwing_struct[:,0],Xwing_struct[:,1]+0.2,'ks',markersize=8,
         markerfacecolor='w',markeredgewidth=2)
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
images = []
ax.set_xlim([-0.2,1.2])
Nsteps = 400
Xhist=np.zeros((Nsteps,5))
CL_curve = np.zeros((Nsteps+1,1))
realTime = 0.0
cont = 0;
re = True

for i in range(Nsteps):
    print('Iteration: ',i,'/',Nsteps)
    realTime = dt*i
    print('Real Time: ',realTime,'s')
    print('_________________________')
    
    if i==20:
        Xcg = 0.6
        x_theta = (Xcg-Xec)
        I_theta = 7.0 + mass*x_theta**2    ## I_t = I_0 + mx_t**2
        
    
    if i>0:
        Xwing_aero=np.dot(H2,Xwing_struct)
        Vc=np.dot(H,Vels)
    
    #4.b. Perform Aero Calculation
    Elements,Xc,Xn,Xt,DL,A=Aero.InitializeAeroModel(Xwing_aero)
    V_el,Cp,CL,CD,CF,CFt=Aero.GetAeroSolution(Elements,Xc,Xn,Xt,DL,A,Vinf,Vc)
    Forces_aero=0.5*1.225*1.0*np.dot(Vinf,Vinf)*CF
    CL_curve[i+1] = CL

    
    #5. Convert Forces to structural Grid
    Forces_struct=np.dot(H.T,Forces_aero) # forces at structural nodes
    Moments=Xwing_struct[:,0]*0           # concentrated moments are zero
    
    Xhist[0:i,0:2]+=Vinf*dt
    Xhist[0:i,2:4]+=Vinf*dt
    Xhist[i,:]=np.array([X[0],X[1],Xwing_struct[0,0],Xwing_struct[0,1],X[4]])
   
    #Plot Geometry
    ax.clear()
    ax.plot(Xwing_aero[:,0],Xwing_aero[:,1],'-k')
    ax.plot([Xec]*2,[X[1],0.0],'-bo',lw=2)
    ax.plot(Xcg,X[1],'kx')
    # ax.plot(Xwing_struct[:,0],Xwing_struct[:,1],'k.')
    i1=np.max([i-20,0])
    ax.plot(np.clip(Xhist[0:i+1,0],-0.2,1.6),Xhist[0:i+1,1],'-m',lw=3,alpha=0.5)
    ax.plot(np.clip(Xhist[0:i+1,2],-0.2,1.6),Xhist[0:i+1,3],'-g',lw=3,alpha=0.5)
    
    for ik in range(Xc.shape[0]):
        ax.arrow(Xc[ik,0],Xc[ik,1],Forces_aero[ik,0]*1,Forces_aero[ik,1]*1,
                 fc='r',ec='r')
    
    ax.plot([-0.1,1.1],[-0.6,0.6],'w.')
    plt.axis('equal')
    plt.grid()
    im = fig2img (fig)
    images.append(im)

    
    #6. Run Struct Model to Get Displacements and Velocities
    
    #6.a Add spring force
    Fs,Ms=Struct.GetSpringForcesMoments(X[0:2],X[2:4],X[6],X[4],X[5],X[7],
                                        x_theta,mass,I_theta)

    #6.b. Get Forces, Moments, Displacements and Velocities
    SF,SM=Struct.GetForceAndMomentatReferencePoint(Xwing_struct,Forces_struct,
                                                   Moments,X[0:2],Fs,Ms)
    X,Xwing_struct,Vels=Struct.GetDisplacementsVelocities(Xwing_struct,SF,SM,X,
                                                          dt,x_theta,mass,I_theta)
    if i*dt%7<0.0001 and i>0:
        print('Here')
        try:
            # Find Peaks Method
            peaks, _ = find_peaks(Xhist[cont:i,4]*57.3)
            auxx = np.zeros(len(peaks))
            cont+=i
            x = peaks*dt
            
            for j in range(len(peaks)):
                auxx[j] = Xhist[peaks[j],4]*57.3
            p = Verify_instability(x,auxx)
            if p[1]<0 and cont==i:
                X_new,Y_new = Determine_Airfoil(desired.reshape(1,31),30)
                Xwing_aero = np.zeros((len(X_new),2))
                X_new1,Y_new1 = Determine_Airfoil(desired.reshape(1,31),6)
                Xwing_struct = np.zeros((len(X_new1),2))
                Xwing_aero[:,0] = X_new
                Xwing_aero[:,1] = Y_new
                Xwing_struct[:,0] = X_new1
                Xwing_struct[:,1] = Y_new1
                Surf,Iben,Itor,Xcg,Ycg = Struct.ArfoilProperties(1.0,Xwing_aero)
                Xcg = Xcg
                #Updating vectors with current airfoil's position in y
                #2.Initialize Aero and Structural Models
                x_theta = (Xcg-Xec)
                print(x_theta)
                I_theta = 7.0 + mass*x_theta**2    ## I_t = I_0 + mx_t**2 
                #2.a. Aero Model
                Elements,Xc,Xn,Xt,DL,A=Aero.InitializeAeroModel(Xwing_aero)
                AoA = 0.0
                Vinf=np.array([np.cos(AoA/57.3),np.sin(AoA/57.3)])*Vnorm #wind speed
                #Vinf=Vinf*Vnorm
                Vc=Xc*0.0
                Vels = Xwing_struct*0.0 #velocities of structural nodes
                
                
                X = np.zeros(8) #struct model state vector
                X[0] = 0.5
                X[1] = 0.0297
                # #3. Create Interpolation Matrix H
                H  = CalculateHmatrix(Xc,Xwing_struct)
                H2 = CalculateHmatrix(Xwing_aero,Xwing_struct)

                
        except:
            print('there are no peaks to calculate')
            
            
images[0].save('out.gif', save_all=True, append_images=images[1:],duration=10,
               loop=0)


plt.show()
plt.grid()
fig,ax = plt.subplots(nrows=2,ncols=1,sharex=True)


ax[0].plot(np.linspace(0,Nsteps*dt,Nsteps),Xhist[:,1],'-k')
ax[0].set_ylabel('CG Y position [m]')
ax[0].grid('on')

ax[1].plot(np.linspace(0,Nsteps*dt,Nsteps),Xhist[:,4]*57.3,'-k')
ax[1].set_ylabel('Theta [deg]')
ax[1].set_xlabel('Time [s]')
ax[1].grid('on')
plt.show()
