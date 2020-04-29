# -*- coding: utf-8 -*-

import numpy as np

def ArfoilProperties(chord,Xwing_struct):
    x = np.zeros(len(Xwing_struct))
    y = np.zeros(len(Xwing_struct))
            
    for i in range(len(Xwing_struct)-1):
        x[i] = chord*Xwing_struct[i,0]
        y[i] = chord*Xwing_struct[i,1]
            
    Ixx    = 0.0
    Iyy    = 0.0
    Surf   = 0.0
    sumXcg = 0.0
    sumYcg = 0.0
    
    for i in range(len(Xwing_struct)-1):
                
        Surf   = Surf + (y[i]+y[i+1])*(x[i]-x[i+1])/2

        bqx    = (x[i]+x[i+1])/2
        bqy    = y[i]/2
        aq     = (x[i]-x[i+1])*y[i]
        btx    = x[i+1]+(x[i]-x[i+1])/3
        bty    = y[i]+(y[i+1]-y[i])/3
        at     = (x[i]-x[i+1])*(y[i+1]-y[i])/2
        sumXcg = sumXcg+bqx*aq+btx*at
        sumYcg = sumYcg+bqy*aq+bty*at

        Ixq = (x[i]-x[i+1])*y[i]**3/12
        Iyq = (x[i]-x[i+1])**3 *y[i]/12

        Ixt = ((x[i]-x[i+1])*(y[i+1]-y[i])**3)/36
        Iyt = ((x[i]-x[i+1])**3 *(y[i+1]-y[i]))/36

        Ixx = Ixx + (Ixq + aq*(bqy**2))+(Ixt+at*(bqy**2))
        Iyy = Iyy + (Iyq + aq*(bqx**2)) + (Iyt+at*(btx**2))                
                

    Xcg = sumXcg/Surf
    Ycg = sumYcg/Surf

    IxxCG = Ixx - Surf*Ycg**2
    IyyCG = Iyy - Surf*Xcg**2

    Iben = IxxCG
    Itor  = IxxCG + IyyCG
              
 
    #arfldata[1] = Surf
    #arfldata[2] = Iben
    #arfldata[3] = Itor
    #arfldata[4] = Xcg
    
    return Surf,Iben,Itor,Xcg,Ycg


def GetForceAndMomentatReferencePoint(Coords,Forces,Moments,RefPoint):
    #Coords: structural grid coordinates
    #Forces: Forces at structural grid nodes
    #Moments: Moments at struct grid nodes
    #RefPoint: Coordinates of reference point

    SF,SM=0.0,0.0
    for i in range(Coords.shape[0]):
        SF+=Forces[i]
        XR=Coords[i,:]-RefPoint
        SM+=np.cross(np.array([XR[0],XR[1],0.0]),np.array([Forces[i,0],Forces[i,1],0.0]))[2]
        SM+=Moments[i]
    return SF,SM

def GetSpringForcesMoments(Point,V,Theta,Omega):
    #Point: attachment point coordinates
    #V: attachment point local speed
    #Theta: pitch angle
    #Omega: rotational rate
    Ks,F0=30.0,0.0
    Cs=6.0
    Fs=np.array([0.0,-Ks*Point[1]+F0-Cs*V[1]])
    Km,M0=50.0,0.0
    Cm=10.0
    Ms=-Km*Theta+M0-Cm*Omega
    
    return Fs,Ms

def Rotate(Xs,Ys,theta):
    Xs2=Xs*np.cos(theta)-Ys*np.sin(theta)
    Ys2=Xs*np.sin(theta)+Ys*np.cos(theta)
    return Xs2,Ys2
    
def GetDisplacementsVelocities(Coords,SF,SM,X,dt):
    
    mass = 10.0
    Inertia = 10.0
    
    # unpack X vector
    RefPoint=X[0:2]
    Vref=X[2:4]
    Theta=X[4]
    Om=X[5]
    
    #prepare grid
    Coords2=Coords*1
    Coords2[:,0]-=RefPoint[0]
    Coords2[:,1]-=RefPoint[1]    
    
    #calculate derivatives
    Vd=np.array([0.0,SF[1]/mass])
    Omd=SM/Inertia
    Xd=Vref*1
    Thetad=Om*1
    
    #update states
    Vref+=Vd*dt
    Om+=Omd*dt
    RefPoint+=Xd*dt
    Theta+=Thetad*dt
    
    #update Xvector
    X2=np.array([RefPoint[0],RefPoint[1],Vref[0],Vref[1],Theta,Om])
    #Calculate grid Velocities
    Vels=Coords*0.0
    Vels[:,0]=Vref[0]
    Vels[:,1]=Vref[1]
    for i in range(Coords.shape[0]):
        Xr=Coords[i,:]-RefPoint
        Vels[i,:]+=np.cross(np.array([0.0,0.0,Om]),np.array([Xr[0],Xr[1],0.0]))[0:2]
    
    #calculate grid displacements
    Coords2[:,0],Coords2[:,1]=Rotate(Coords2[:,0],Coords2[:,1],Thetad*dt)
    Coords2[:,0]+=RefPoint[0]
    Coords2[:,1]+=RefPoint[1]

    return X2,Coords2,Vels
