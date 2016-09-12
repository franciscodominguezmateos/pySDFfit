'''
Created on Sept 11, 2016

@author: Francisco Dominguez
A GaussNewtong version in 2D
Theta is alwais at 0 then change points on each iteration
'''
from visual import *
import math
import numpy as np
import numpy.linalg as la

r=10

def SDF(p):
    offset=np.array([0,10])
    offset1=np.array([-5,-5])
    d0=la.norm(p)-r
    d1=la.norm(p-offset)-r/1.25
    d2=la.norm(p+offset)-r/1.35
    d3=la.norm(p+offset1)-r/1.15
    return min(min(min(d0,d1),d2),d3)
    #return d0

def JSDF(p):
    delta=0.0001
    ix=np.array([delta,0])
    iy=np.array([0,delta])
    dx=(SDF(p+ix)-SDF(p))/delta
    dy=(SDF(p+iy)-SDF(p))/delta
    return np.array([dx,dy])

def getT(theta):
    alpha=theta[0,0]
    tx=theta[1,0]
    ty=theta[2,0]
    cs=math.cos(alpha)
    sn=math.sin(alpha)
    T=np.matrix([[cs,-sn,tx],
                 [sn, cs,ty],
                 [0 ,  0, 1]])
    return T


def getParam(theta):
    return theta[0,0],theta[1,0],theta[2,0]

def JF(i,Q,P,theta):
    alpha,tx,ty=getParam(theta)
    cs=math.cos(alpha)
    sn=math.sin(alpha)
    pxi=P[0,i]
    pyi=P[1,i]
    qxi=Q[0,i]
    qyi=Q[1,i]
    q=np.array([qxi,qyi])
    Dq=SDF(q)
    JDqx,JDqy=JSDF(q)
    #Jacobian is a row vector not a columne vector
    JF=np.matrix([Dq*(JDqx*(-pxi*sn-pyi*cs)+JDqy*(pxi*cs-pyi*sn)),
                 Dq*JDqx,
                 Dq*JDqy])
    return JF

def JEls(P,theta):
    t=0
    T=getT(theta)
    Q=T*P
    n=shape(P)[1]
    for i in range(n):
        t+=JF(i,Q,P,theta)
    return t/n

def Els(P,theta):
    T=getT(theta)
    Q=T*P
    e=0
    n=Q.shape[1]
    for i in range(n):
        qxi=Q[0,i]
        qyi=Q[1,i]
        q=np.array([qxi,qyi])
        e+=SDF(q)**2
    return e/n


def Hg(P,theta):
    k=theta.shape[0]
    H=np.matrix(np.eye(k,k))
    g=np.matrix(np.zeros((k,1)))
    T=getT(theta)
    Q=T*P
    n=shape(P)[1]
    for i in range(n):
        J=JF(i,Q,P,theta)
        JT=J.T
        D=SDF(P[:,i])
        H+=JT*J
        g+=JT*D;
    return H,g
#Draw SDF
for x in linspace(-20, 20, 60):
    for y in linspace(-20,20,60):
        p=np.array([x,y])
        d=SDF(p)/20
        if d<0:
            c=(-d,-d,-d)
        else:
            c=(0,d,d)
        box(pos=p,color=c)
#transformation to fit
alpha=-math.pi/4
cs=math.cos(alpha)
sn=math.sin(alpha)
tx=5
ty=-2;
T=np.matrix([[cs,-sn,tx],
             [sn, cs,ty],
             [0 ,  0, 1]])
#Create point to transform
lpts=[]
for a in np.arange(0,math.pi/4,math.pi/20):
    print a
    x=r*math.cos(a)+np.random.randn()*0.0125
    y=r*math.sin(a)+np.random.randn()*0.0125
    p=np.array([x,y,1])
    lpts.append(p)
pts=np.array(lpts)
vpts1=points(pos=pts,color=(1,0,0))
Q=np.matrix(pts).T
P=T*Q
vpts2=points(pos=P.T,color=(0,1,0))

print P
theta=np.matrix([[0.0],[0.0],[0.0]])
print "Els=",Els(P,theta)
Pp=P
thetaZero=theta
for i in range(20000):
    #Just need Jacobian at 0
    H,g=Hg(Pp,thetaZero)
    Hi=np.linalg.inv(H)
    dTheta=Hi*g
    thetap=-0.02*dTheta
    theta=theta+thetap
    Tp=getT(thetap)
    Pp=Tp*Pp
    
    e=Els(P,theta)
    if e<0.05:
        break
    if i % 10==0:
        #print "JEls=",JE
        #print "theta =" , theta
        T=getT(theta)
        R=T*P
        pR=points(pos=R.T, color=color.yellow)
        print "Els=",e,"i=",i
#print theta
T=getT(theta)
R=T*P
pR=points(pos=R.T, color=color.blue)
print np.linalg.inv(T)
print "Els=",e,"i=",i



