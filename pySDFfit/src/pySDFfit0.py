'''
Created on Mar 17, 2016

@author: Francisco Dominguez

'''
from visual import *
import math
import numpy as np
import numpy.linalg as la

r=10

def SDF(p):
    return la.norm(p)-r

def JSDF(p):
    n=la.norm(p)
    return p/n

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
    JF=np.matrix([[Dq*(JDqx*(-pxi*sn-pyi*cs)+JDqy*(pxi*cs-pyi*sn))],
                 [Dq*JDqx],
                 [Dq*JDqy]])
    return JF

def Els(P,theta):
    T=getT(theta)
    Q=T*P
    e=0
    n=Q.shape[1]
    for i in range(n):
        qxi=Q[0,i]
        qyi=Q[1,i]
        q=np.array([qxi,qyi])
        e+=SDF(q)
    return e

def JEls(P,theta):
    t=0
    T=getT(theta)
    Q=T*P
    for i in range(shape(P)[1]):
        t+=JF(i,Q,P,theta)
    return t


#transformation to fit
theta=math.pi/4
cs=math.cos(theta)
sn=math.sin(theta)
tx=5
ty=5;
T=np.matrix([[cs,-sn,tx],
             [sn, cs,ty],
             [0 ,  0, 1]])
#Create point to transform
lpts=[]
for a in np.arange(0,math.pi/2,math.pi/20):
    print a
    x=r*math.cos(a)
    y=r*math.sin(a)
    p=np.array([x,y,1])
    lpts.append(p)
pts=np.array(lpts)
vpts1=points(pos=pts,color=(1,0,0))
varr=arrow(pos=(0,0,0),axis=(20,0,0))
Q=np.matrix(pts).T
P=T*Q
vpts2=points(pos=P.T,color=(0,1,0))

print P
theta=np.matrix([[0.0],[0.0],[0.0]])
print "Els=",Els(P,theta)
for i in range(700):
    JE=JEls(P,theta)
    theta=theta-0.001*JE
    e=Els(P,theta)
    if e<0.1:
        break
    if i % 10==0:
        #print "JEls=",JE
        #print "theta =" , theta
        T=getT(theta)
        R=T*P
        pR=points(pos=R.T, color=color.yellow)
        print "Els=",e,theta
print theta
T=getT(theta)
R=T*P
pR=points(pos=R.T, color=color.yellow)
print np.linalg.inv(T)

print SDF(np.array([[1,0],[2,0],[3,0]]))



