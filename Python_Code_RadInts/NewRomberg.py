import numpy as np
from scipy.special import *
import scipy.integrate as integ
from NewHydrogenic import *

#Romberg Algorithm radial integrals 
#Takes an open continuum kinetic energy in atomic units and an effective quantum number (subject to chenge is not integer since there is possibility of divergence).
#d wave regular function
def dcont_RadIntRombReg(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1.5e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,2,xx)*xx*W(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#g wave regular function
def gcont_RadIntRombReg(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1.5e-7,5*nu**2,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,4,xx)*xx*W(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#d wave irregular function
def dcont_RadIntRombIreg(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1.5e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,2,xx)*xx*W(nu,3,xx)*SB(3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#g wave irregular function
def gcont_RadIntRombIreg(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1.5e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,4,xx)*xx*W(nu,3,xx)*SB(4,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

#Adaptative quadrature integrals.
def dRadIntReg(EcontAU,nu):
    return integ.quadrature(lambda x: x*Ucl(EcontAU,2,x)*W(nu,3,x),0,5*(nu)**2,tol=1.49e-10, rtol=1.49e-10, maxiter=1000)[0]
def gRadIntReg(EcontAU,nu):
    return integ.quadrature(lambda x: x*Ucl(EcontAU,4,x)*W(nu,3,x),0,5*(nu)**2,tol=1.49e-10, rtol=1.49e-10, maxiter=1000)[0]

def dRadIntIreg(EcontAU,nu):
    return integ.quadrature(lambda x: x*Ccl(EcontAU,2,x)*W(nu,3,x)*SB(3,x),0,5*(nu)**2,tol=1.49e-10, rtol=1.49e-10, maxiter=1000)[0]
def gRadIntIreg(EcontAU,nu):
    return integ.quadrature(lambda x: x*Ccl(EcontAU,4,x)*W(nu,3,x)*SB(4,x),0,5*(nu)**2,tol=1.49e-10, rtol=1.49e-10, maxiter=1000)[0]

#d wave regular function with the hydrogenic function
def dcont_RadIntRombRegHyd(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,2,xx)*xx*Unl(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#g wave regular function
def gcont_RadIntRombRegHyd(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,4,xx)*xx*Unl(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#d wave irregular function
def dcont_RadIntRombIregHyd(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,2,xx)*xx*Unl(nu,3,xx)*SB(3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#g wave irregular function
def gcont_RadIntRombIregHyd(EcontAU,nu):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-6,5*nu**2,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,4,xx)*xx*Unl(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

#d whittaker function
def dWhit_RadIntRombIregHyd(EnegAU,nu):
    dBW=2*np.pi/np.sqrt(-2*EnegAU)
    nuw=1/np.sqrt(-2*EnegAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*nu**2,Nromb,endpoint=True)
    y =  [W(nuw,2,xx)*xx*Unl(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
def dHyd_RadIntRombIregHyd(EnegAU,nu):
    dBW=2*np.pi/np.sqrt(-2*EnegAU)
    nuw=int(1/np.sqrt(-2*EnegAU))
    print(nuw)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*nu**2,Nromb,endpoint=True)
    y =  [Unl(nuw,2,xx)*xx*Unl(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint
#g wave irregular function
def gWhit_RadIntRombIregHyd(EnegAU,nu):
    dBW = 2*np.pi/np.sqrt(-2*EnegAU)
    nuw = 1/np.sqrt(-2*EnegAU)
    NdBW = max(int((5.*nu**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-6,k(EnegAU)*5*nu**2,Nromb,endpoint=True)
    y =  [W(nuw,4,xx)*xx*Unl(nu,3,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

#Arbitrary orbital angular momentum
def lcont_RadIntRombRegHyd(EcontAU,l,n,ln):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*n**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*n**2,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,l,xx)*xx*Unl(n,ln,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

def lcont_RadIntRombIregHyd(EcontAU,l,n,ln):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    NdBW = max(int((5.*n**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    mesh = np.linspace(1e-5,5*n**2,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,l,xx)*xx*Unl(n,ln,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    return radint

#Arbitrary orbital angular momentum

# coulomb f and whittaker function
def lcont_RadIntRombRegWhit(EcontAU,l,EnegAU,ln, plt, name=None):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    nuw = 1/np.sqrt(-2*EnegAU)
    NdBW = max(int((5.*nuw**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    nlim = max(nuw,10)
    lm = max(l,ln)
    mesh = np.linspace(1e-5,5*nlim**2,Nromb,endpoint=True)
    y =  [Ucl(EcontAU,l,xx)*xx*Wstich(nuw,ln,xx)*SBb(lm,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    if(plt):
        if(name==None):
            name="%i-%i_sanity_check%.3e.dat"%(l,ln,EcontAU)
        with open(name, mode='w') as f:
            for xx, yy in zip(mesh, y):
                f.write("%.8e   %.8e \n"%(xx,yy))
    return radint

#Coulomb g and whittaker function
def lcont_RadIntRombIregWhit(EcontAU,l,EnegAU,ln, plt, name=None):
    dBW = 2*np.pi/np.sqrt(2*EcontAU)
    nuw = 1/np.sqrt(-2*EnegAU)
    NdBW = max(int((5.*nuw**2.)//dBW),10)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    lm = max(l,ln)
    nlim = max(nuw,10)
    mesh = np.linspace(1e-5,5*nlim**2,Nromb,endpoint=True)
    y =  [Ccl(EcontAU,l,xx)*xx*Wstich(nuw,ln,xx)*SBb(lm,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    if(plt):
        if(name==None):
            name="%i-%i_Iregsanity_check%.3e.dat"%(l,ln,EcontAU)
        with open(name, mode='w') as f:
            for xx, yy in zip(mesh, y):
                f.write("%.8e   %.8e \n"%(xx,yy))
    return radint

#Two Whittakers
def lcont_RadIntRombWhitWhit(Eneg1,l,Eneg2,ln, plt, name=None):
    nu1 = 1/np.sqrt(-2*Eneg1)
    nu2 = 1/np.sqrt(-2*Eneg2)
    nuw = min(nu1, nu2)
    NdBW = max(int((5.* nuw**2.)),30)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    lm = max(l,ln)
    nlim = max(nuw,10)
    mesh = np.linspace(1e-5,5*nlim**2,Nromb,endpoint=True)
    y =  [Wstich(nu1,l,xx)*xx*Wstich(nu2,ln,xx)*SBb(lm,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.romb(y,dx=dmesh)
    if(plt):
        if(name==None):
            name="%i-%i_Iregsanity_check%.3e.dat"%(l,ln,nu1)
        with open(name, mode='w') as f:
            for xx, yy in zip(mesh, y):
                f.write("%.8e   %.8e \n"%(xx,yy))
    return radint


#Two Whittakers, overlap
def OVERLAP_lcont_RadIntRombWhitWhit(Eneg1,l,Eneg2,ln, plt, name=None):
    nu1 = 1/np.sqrt(-2*Eneg1)
    nu2 = 1/np.sqrt(-2*Eneg2)
    nuw = min(nu1, nu2)
    NdBW = max(int((5.* nuw**2.)),30)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    lm = max(l,ln)
    nlim = max(nuw,10)
    mesh = np.linspace(1e-5,5*nlim**2,Nromb,endpoint=True)
    y =  [Wstich(nu1,l,xx)*Wstich(nu2,ln,xx)*SB(lm,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.trapezoid(y,dx=dmesh)
    if(plt):
        if(name==None):
            name="%i-%i_Iregsanity_check%.3e.dat"%(l,ln,nu1)
        with open(name, mode='w') as f:
            for xx, yy in zip(mesh, y):
                f.write("%.8e   %.8e \n"%(xx,yy))
    return radint

def OVERLAP_lcont_RadIntRombWhitWhitst(Eneg1,l,Eneg2,ln, plt, name=None):
    nu1 = 1/np.sqrt(-2*Eneg1)
    nu2 = 1/np.sqrt(-2*Eneg2)
    nuw = min(nu1, nu2)
    NdBW = max(int((5.* nuw**2.)),30)
    Nromb = 2**(int(np.log2(200*NdBW)))+1
    lm = max(l,ln)
    nlim = max(nuw,10)
    mesh = np.linspace(1e-5,5*nlim**2,Nromb,endpoint=True)
    y =  [Wstich(nu1,l,xx)*Wstich(nu2,ln,xx)*SBb(lm,xx) for xx in mesh]
    dmesh = mesh[1]-mesh[0]
    radint=integ.trapezoid(y,dx=dmesh)
    if(plt):
        if(name==None):
            name="%i-%i_Iregsanity_check%.3e.dat"%(l,ln,nu1)
        with open(name, mode='w') as f:
            for xx, yy in zip(mesh, y):
                f.write("%.8e   %.8e \n"%(xx,yy))
    return radint