# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 09:26:11 2020

@author: NG78D39
"""

import numpy as np
import matplotlib as plt
import time

def NACA4(c,t,m,p,NeL):
    #NACA4: Outputs NACA4 Coordinates
    #yc: y due to camber
    #yt: y due to thickness
    #c: chord
    #t: maximum thickness
    #m: maximum camber
    #p: location of maximmum camber
    import numpy as np
    x=np.linspace(0,1,NeL)**2
    #x=0.5*(1-np.cos(x))
    yt,yc=x*0.0,x*0.0
    for i in range(len(x)):
        yt[i]=5*t*(0.2969*np.sqrt(x[i])-0.126*x[i]-0.3516*x[i]**2+0.2843*x[i]**3\
                    -0.1015*x[i]**4)*c
        if x[i]<p:
            yc[i]=m/p**2*(2*p*x[i]-x[i]**2)*c
        else:
            yc[i]=m/(1-p)**2*((1-2*p)+2*p*x[i]-x[i]**2)*c
    x=x*c
    #create points surrounding shape
    x2=np.zeros((2*len(x)-1,2))
    x2[0:len(x),0]=np.array(list(reversed(x)))
    y2=yc+yt
    #y2[-1]=0.0
    x2[0:len(x),1]=np.array(list(reversed(y2)))
    x2[len(x)::,0]=x[1:len(x)]
    y2=yc-yt
    #y2[-1]=0.0
    x2[len(x)::,1]=y2[1:len(x)]
    return x2



def Rotate(Xs,Ys,theta):
    Xs2=Xs*np.cos(theta)-Ys*np.sin(theta)
    Ys2=Xs*np.sin(theta)+Ys*np.cos(theta)
    return Xs2,Ys2

def InducedVelocity(Gamma,Xwing,Xt,Xn,iel,Elements,Xp):
    #Gamma: vorticity
    #Xp: point coordinates
    Pt1=Xwing[int(Elements[iel,0]),:]
    Pt2=Xwing[int(Elements[iel,1]),:]
    
    #print('Pt1=',Pt1)
    #print('Pt2=',Pt2)
    #print('Xp=',Xp)
    Xr=Xp-Pt1 #calculate relative position vector
    #convert to panel coordinates
    Xr_pc=np.array([np.dot(Xr,Xt[iel,:]),np.dot(Xr,Xn[iel,:])])
    Pt1_pc=np.array([0.0,0.0])
    Pt2_pc=np.array([np.dot(Pt2-Pt1,Xt[iel,:]),0.0])
    
    
    #print('Xp_pc=',Xr_pc)
    #print('X2_pc=',Pt2_pc)
            
    R1=np.linalg.norm(Xr_pc-Pt1_pc)
    R2=np.linalg.norm(Xr_pc-Pt2_pc)#np.sqrt((X[0]-X2[0])**2+X[1]**2)
    
    
    TH1=np.arctan2(Xr_pc[1],Xr_pc[0])
    TH2=np.arctan2(Xr_pc[1],Xr_pc[0]-Pt2_pc[0])
    
    
    UL=0.15916*(TH2-TH1)
    WL=0.15916*np.log(R2/R1)
    
    #print(UL,WL)
    Vind=UL*Xt[iel,:]+WL*Xn[iel,:]

    return Vind


def CreateAirfElements(Xwing):
    #create grid
    Elements=np.zeros((Xwing.shape[0]-1,2))
    Xc=Elements*0
    Xn=Elements*0
    Xt=Elements*0
    DL=np.zeros(Elements.shape[0])
    for i in range(len(Elements)):
        Elements[i,:]=np.array([i+1,i])
        
        #calculate normal and tangential vectors
        Xr=Xwing[i+1,:]-Xwing[i,:]
        DL[i]=np.linalg.norm(Xr)
        Xr=Xr/DL[i]
        Xt[i,:]=-Xr
        Xn[i,:]=np.cross(np.array([Xr[0],Xr[1],0.0]),np.array([0,0,1]))[0:2]
        #calculate centroids
        Xc[i,:]=0.5*(Xwing[i,:]+Xwing[i+1,:])-0.0*DL[i]*Xn[i,:] #centroids slightly inside panel
    return Elements,Xc,Xn,Xt,DL

def CreateInfluenceMatrix(Xwing,Elements,Xc,Xn,Xt):
    N_el=Elements.shape[0]
    A=np.zeros((N_el,N_el)) #initialize influence matrix
    for i in range(N_el):
        for j in range(N_el):
            if i==j:
                A[i,j]=-0.5
            else:
                Vind=InducedVelocity(1.0,Xwing,Xt,Xn,j,Elements,Xc[i,:])
                A[i,j]=np.dot(Vind,Xt[i,:])
    #add equation for Kutta condition
    m=int(3*N_el/4)
    A[m,:]=A[m,:]*0
    A[m,0],A[m,-1]=1.0,1.0
    return A

def CreateRHSMatrix(Vinf,Xt,Xc,Vc):
    N_el=Xt.shape[0]
    B=np.zeros(N_el)
    for i in range(N_el):
        Vinfi=Vinf-Vc[i,:]
        B[i]=-np.dot(Vinfi,Xt[i,:])
    #equation for Kutta condition
    m=int(3*N_el/4)
    B[m]=0.0
    return B

def GetFlowField(Elements,Xc,Xt,Xn,Vinf,Gamma,DL,A,Vc):
    N_el=Xc.shape[0]
    V_el=Elements*0
    Cp=np.zeros(N_el)
    CF=np.zeros((N_el,2))
    normVinf=np.linalg.norm(Vinf)
    for i in range(N_el):
        Vinfi=Vinf-Vc[i,:]
        A2=A*1
        A2[i,i]=0.5
        V_el[i,:]=np.dot(A2[i,:],Gamma[:])*Xt[i,:]+Xt[i,:]*np.dot(Vinfi,Xt[i,:])
        Cp[i]=1-(np.linalg.norm(V_el[i,:])/normVinf)**2
    for i in range(N_el):
        if i==int(3*N_el/4):
            Cp[i]=np.interp(Xc[i,0],[Xc[i-1,0],Xc[i+1,0]],[Cp[i-1],Cp[i+1]])
        CF[i,:]=-Cp[i]*DL[i]*Xn[i,:]
    CFt=np.sum(CF,axis=0)
    CD=np.dot(CFt,Vinf/normVinf)
    CL=np.linalg.norm(CFt-CD*Vinf/normVinf)*np.sign(CFt[1])
    return V_el,Cp,CL,CD,CF,CFt

def InitializeAeroModel(Xwing):
    #create elements, centroids, normal, tangent and length vectors
    Elements,Xc,Xn,Xt,DL=CreateAirfElements(Xwing)
    #create influence matrix
    A=CreateInfluenceMatrix(Xwing,Elements,Xc,Xn,Xt)
    return Elements,Xc,Xn,Xt,DL,A
    
def GetAeroSolution(Elements,Xc,Xn,Xt,DL,A,Vinf,Vc):
    B=CreateRHSMatrix(Vinf,Xt,Xc,Vc) #create Right-Hand side matrix
    Gam=np.linalg.solve(A,B) #solve linear system
    V_el,Cp,CL,CD,CF,CFt=GetFlowField(Elements,Xc,Xt,Xn,Vinf,Gam,DL,A,Vc) #get flowfield
    return V_el,Cp,CL,CD,CF,CFt
