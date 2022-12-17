import coul92 as cl
import couln as cln
import gensub as gs
import numpy as np
import scipy.special as spc
import math 
import time

def k(e):
    return np.sqrt(2*e)
def kappa(e):
	return np.sqrt(-2*e)
def etaf(e,l):
	return 1/k(e) * np.log(2*k(e)) + np.angle(spc.gamma(l+1-1j/k(e)))-1/2. * l *np.pi
def nuf(e):
	return np.sqrt(1/(-2*e))
def D(e,l):
	if(spc.gamma(nuf(e)-l)>0):
		return np.sqrt(np.pi) * (2/nuf(e))**(nuf(e)) * 1/np.sqrt(spc.gamma(nuf(e)-l)*spc.gamma(l+1+nuf(e)))
	else:
		return np.sqrt(np.pi) * (2/nuf(e))**(nuf(e)) * 1/np.sqrt((-1*spc.gamma(nuf(e)-l))*spc.gamma(l+1+nuf(e)))

def F(l,eta,rho):
    [f,*args] = cl.coul90(rho,eta,0,l,0,0)
    return f[-1]

def Fp(l,eta,rho):
    [f,g,fp,*args] = cl.coul90(rho,eta,0,l,0,0)
    return fp[-1]

def G(l,eta,rho):
    [f,g,*args] = cl.coul90(rho,eta,0,l,0,0)
    return -1*g[-1]

def Gp(l,eta,rho):
    [*args,gp] = cl.coul90(rho,eta,0,l,0,0)
    return -1*gp[-1]



#Regular Coulomb function for an attractive potential and positive energy
def Ucl(e,l,r):
	if(k(e)*r > 1):
		return np.sqrt(2/(np.pi*k(e)))*F(l,-1/k(e),r*k(e))
	else:
		return gs.seaton(l,e*2.0,r,1.0)[0]

#Long Range behavior
def LRUcl(e,l,r):
	return np.sqrt(2/(np.pi*k(e))) * np.sin(k(e)*r+1/k(e) * np.log(r) + etaf(e,l)) 

#And its derivative
def Uclp(e,l,r):
	if(k(e)*r>1):
		return np.sqrt(2/(np.pi*k(e)))*Fp(l,-1/k(e),r*k(e))*k(e)
	else:
		return gs.seaton(l,e*2.0,r,1.0)[1]
 

#Regular Coulomb function for an attractive potential and negative energy
def Ucln(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[0]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1)[0]

#Analytic coulomb function for negative energies and nu<l
def Uclno(e,l,r):
	return gs.coulfg(l,e*2.0,r,1e-14)[0]

def Uclnop(e,l,r):
	return gs.coulfg(l,e*2.0,r,1e-14)[1]

def LRUclno(e,l,r):
	nu = nuf(e)
	prefac = np.exp(r/nu)*r**(-nu)
	
	r8 = (2**(-15 - nu)*nu**(9 + l + nu)*(-l + nu)*(1 - l + nu)*(2 - l + nu)*(3 - l + nu)*(4 - l + nu)*(5 - l + nu)*(6 - l + nu)*(7 - l + nu)*(1 + l + nu)*(2 + l + nu)*(3 + l + nu)*(4 + l + nu)*(5 + l + nu)*(6 + l + nu)*(7 + l + nu)*(8 + l + nu))/(315.*r**8*spc.gamma(1 + l - nu))
	r7 = (2**(-11 - nu)*nu**(8 + l + nu)*(-l + nu)*(1 - l + nu)*(2 - l + nu)*(3 - l + nu)*(4 - l + nu)*(5 - l + nu)*(6 - l + nu)*(1 + l + nu)*(2 + l + nu)*(3 + l + nu)*(4 + l + nu)*(5 + l + nu)*(6 + l + nu)*(7 + l + nu))/(315.*r**7*spc.gamma(1 + l - nu))
	r6 = (2**(-10 - nu)*nu**(7 + l + nu)*(-l + nu)*(1 - l + nu)*(2 - l + nu)*(3 - l + nu)*(4 - l + nu)*(5 - l + nu)*(1 + l + nu)*(2 + l + nu)*(3 + l + nu)*(4 + l + nu)*(5 + l + nu)*(6 + l + nu))/(45.*r**6*spc.gamma(1 + l - nu))
	r5 = (2**(-8 - nu)*nu**(6 + l + nu)*(-l + nu)*(1 - l + nu)*(2 - l + nu)*(3 - l + nu)*(4 - l + nu)*(1 + l + nu)*(2 + l + nu)*(3 + l + nu)*(4 + l + nu)*(5 + l + nu))/(15.*r**5*spc.gamma(1 + l - nu))
	r4 = (2**(-7 - nu)*(-1 + l/nu)*(-1 - 3/nu + l/nu)*(-1 - 2/nu + l/nu)*(-1 - 1/nu + l/nu)*(1 + 1/nu + l/nu)*(1 + 2/nu + l/nu)*(1 + 3/nu + l/nu)* (1 + 4/nu + l/nu)*nu**(12 - 2*(-0.5 - l/2. - nu/2.))*spc.gamma(1 + 2*(0.5 + l)))/(3.*r**4*spc.gamma(2 + 2*l)*spc.gamma(1 + l - nu))
	r3 = -1/3.*(2**(-4 - nu)*(-1 + l/nu)*(-1 - 2/nu + l/nu)*(-1 - 1/nu + l/nu)*(1 + 1/nu + l/nu)*(1 + 2/nu + l/nu)*(1 + 3/nu + l/nu)*nu**(10 - 2*(-0.5*l - nu/2.))*spc.gamma(1 + 2*(0.5 + l)))/(r**3*spc.gamma(2 + 2*l)*spc.gamma(1 + l - nu))
	r2 = (2**(-3 - nu)*(-1 + l/nu)*(-1 - 1/nu + l/nu)*(1 + 1/nu + l/nu)*(1 + 2/nu + l/nu)*nu**(6 - 2*(-0.5 - l/2. - nu/2.))*spc.gamma(1 + 2*(0.5 + l)))/(r**2*spc.gamma(2 + 2*l)*spc.gamma(1 + l - nu))
	r1 = -((2**(-1 - nu)*(-1 + l/nu)*(1 + 1/nu + l/nu)*nu**(4 - 2*(-0.5*l - nu/2.))*spc.gamma(1 + 2*(0.5 + l)))/(r*spc.gamma(2 + 2*l)*spc.gamma(1 + l - nu)))
	r0 = spc.gamma(1 + 2*(0.5 + l))/(2**nu*nu**(2*(-0.5 - l/2. - nu/2.))*spc.gamma(2 + 2*l)*spc.gamma(1 + l - nu))

	return prefac*(r0+r1+r2+r3+r4+r5+r6+r7+r8)

def Uclno_wLR(e,l,r):
	cfg = gs.coulfg(l,e*2.0,r,1e-14)
	if(cfg[-2]==2):
		return LRUclno(e,l,r)
	else: 
		return cfg[0]

#Long Range behavior
def LRUcln(e,l,r):
	nu = nuf(e)
	A = spc.gamma(nu+1+l)/(spc.gamma(nu-l) * nu**(2*l+1))
	return np.sqrt(A) * LRUclno(e,l,r)

# Long range behavior for the regular function. Normalized in energy. 
def Ucln_wLR(e,l,r):
	nu = nuf(e)
	A = spc.gamma(nu+1+l)/(spc.gamma(nu-l) * nu**(2*l+1))
	return np.sqrt(A) * Uclno_wLR(e,l,r)

#And its derivative
def Uclpn(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[1]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1.0)[1]

#Irregular Coulomb function for an attractive potential and positive energy
def Ccl(e,l,r):
	if(k(e)*r>1):
		return np.sqrt(2/(np.pi*k(e)))*G(l,-1/k(e),r*k(e))
	else:
		return gs.seaton(l,e*2.0,r,1.0)[2]

#And its derivative.
def Cclp(e,l,r):
	if(k(e)*r>1): 
		return np.sqrt(2/(np.pi*k(e)))*Gp(l,-1/k(e),r*k(e))*k(e)
	else:
		return gs.seaton(l,e*2.0,r,1.0)[3]


#Long Range behavior
def LRCcl(e,l,r):
	return -np.sqrt(2/(np.pi*k(e))) * np.cos(k(e)*r+1/k(e) * np.log(r) + etaf(e,l))

#Regular Coulomb function for an attractive potential and negative energy
def Ccln(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[2]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1.0)[2]

#Long range behavior
def LRCcln(e,l,r):
	nu = nuf(e)
	prefac =  (1.*(-1)**l*np.exp(-0.6931471805599453*nu + (1.*r)/nu)*nu**(-0.5 + nu)*np.cos(np.pi*nu)*spc.gamma(1. + l - 1.*nu)*np.sqrt(spc.gamma(-1.*l + nu))*np.sqrt(spc.gamma(1. + l + nu)))/(r**nu*spc.gamma(1. + l - nu))

	r8 = (-3.083824387166251e-8*nu**9.*(-1.*l + nu)*(1. - 1.*l + nu)*(2. - 1.*l + nu)*(3. - 1.*l + nu)*(4. - 1.*l + nu)*(5. - 1.*l + nu)*(6. - 1.*l + nu)*(7. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu)*(3. + l + nu)*(4. + l + nu)*(5. + l + nu)*(6. + l + nu)*(7. + l + nu)*(8. + l + nu))/r**8
	r7 = (-4.934119019466002e-7*nu**8.*(-1.*l + nu)*(1. - 1.*l + nu)*(2. - 1.*l + nu)*(3. - 1.*l + nu)*(4. - 1.*l + nu)*(5. - 1.*l + nu)*(6. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu)*(3. + l + nu)*(4. + l + nu)*(5. + l + nu)*(6. + l + nu)*(7. + l + nu))/r**7
	r6 =         (-6.907766627252403e-6*nu**7.*(-1.*l + nu)*(1. - 1.*l + nu)*(2. - 1.*l + nu)*(3. - 1.*l + nu)*(4. - 1.*l + nu)*(5. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu)*(3. + l + nu)*(4. + l + nu)*(5. + l + nu)*(6. + l + nu))/r**6

	r5 =         (-0.00008289319952702883*nu**6.*(-1.*l + nu)*(1. - 1.*l + nu)*(2. - 1.*l + nu)*(3. - 1.*l + nu)*(4. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu)*(3. + l + nu)*(4. + l + nu)*(5. + l + nu))/r**5

	r4 = (-0.0008289319952702883*nu**5.*(-1.*l + nu)*(1. - 1.*l + nu)*(2. - 1.*l + nu)*(3. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu)*(3. + l + nu)*(4. + l + nu))/r**4

	r3 = (-0.006631455962162306*nu**4.*(-1.*l + nu)*(1. - 1.*l + nu)*(2. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu)*(3. + l + nu))/r**3

	r2 = (-0.039788735772973836*nu**3.*(-1.*l + nu)*(1. - 1.*l + nu)*(1. + l + nu)*(2. + l + nu))/r**2

	r1 = (-0.15915494309189535*nu**2.*(-1.*l + nu)*(1. + l + nu))/r

	r0 = -0.3183098861837907*nu**1.
	return prefac*(r0+r1+r2+r3+r4+r5+r6+r7+r8)

#Long range form for the irregular function. Energy normalized. 
def Ccln_wLR(e,l,r):
	cfg = gs.coulfg(l,e*2.0,r,1e-14)
	[a, g] = gs.ganda(l,2.0*e,1,99)
	if(cfg[-2]==2):
		return LRCcln(e,l,r)
	else: 
		return cfg[0]*g/np.sqrt(a)+cfg[2]*1/np.sqrt(a)

#And its derivative
def Cclpn(e,l,r):
	nu = nuf(e)
	if(nu>l):
		return gs.seaton(l,e*2.0,r,1.0)[3]
	else:
		return -1*gs.seaton1(l,e*2.0,r,1.0)[3]
#Hydrogenic orbitals
def Unl(n,l,r):
    return r*np.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))*spc.assoc_laguerre(2*r/n,n-l-1,2*l+1)*np.exp(-r/n)*(r**l)*((2/n)**(l+3./2.))
#Closed channel function. 

#Scaled Whittaker function, not energy normalized.
def W_scal(nu,l,r):
	wa1 = np.zeros((20,10))
	wa2 = np.zeros((20,10))
	wa3 = np.zeros((20,10))
	wa4 = np.zeros((20,10))
	W = cln.couln(l,1,-1/nu**2,r,1e-12,wa1,wa2,wa3,wa4,0,70)
	return W[0]
#Long range behavior.
def LRW_scal(nu,l,r):
	prefac = np.exp(-r/nu) * r**nu
	r0 = 2**nu *(1/nu)**nu
	r4  = (2**(-7 + nu)*(1 + l - nu)*(2 + l - nu)*(3 + l - nu)*(4 + l - nu)*nu**(4 - nu)*(-3 + l + nu)*(-2 + l + nu)*(-1 + l + nu)*(l + nu))/(3.*r**4)
	r3 = (2**(-4 + nu)*(1 + l - nu)*(2 + l - nu)*(3 + l - nu)*nu**(3 - nu)*(-2 + l + nu)*(-1 + l + nu)*(l + nu))/(3.*r**3)
	r2 = (2**(-3 + nu)*(1 + l - nu)*(2 + l - nu)*nu**(2 - nu)*(-1 + l + nu)*(l + nu))/r**2
	r1 = (2**(-1 + nu)*(1 + l - nu)*nu**(1 - nu)*(l + nu))/r
	return prefac*(r0+r1+r2+r3+r4)
#Derivative of the scaled function.
def W_scalp(nu,l,r):
	wa1 = np.zeros((20,10))
	wa2 = np.zeros((20,10))
	wa3 = np.zeros((20,10))
	wa4 = np.zeros((20,10))
	W = cln.couln(l,1,-1/nu**2,r,1e-12,wa1,wa2,wa3,wa4,0,70)
	return W[1]

#Whittaker from Coulumb
def Wseat(e,l,r):
	nu = np.sqrt(-0.5/e)
	beta = np.pi*(nu-l)
	f, fp, g, gp = gs.coulfg(l,e*2.0,r,1e-14)[0:4]
	return -np.sin(beta)*g-np.cos(beta)*f

#Closed channel function normalized to unit norm at zero quantum defect
def Wstich(nu,l,r):
	if((2*r/nu)>=100 ):
		gam = spc.gamma(nu-l)
		prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * abs(gam))
		return prefac*LRW_scal(nu,l,r)
	elif((2*r/nu)>=5 ):
		gam = spc.gamma(nu-l)
		prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * abs(gam))
		return prefac*W_scal(nu,l,r)
	else:
		e=-0.5/nu**2
		beta = np.pi*(nu-l)
		f = Ucln(e,l,r)
		g = Ccln(e,l,r)
		return (-np.sin(beta)*SB(l,r)*g-np.cos(beta)*f)
#Regular Whittaker
def W(nu,l,r):
	gam = spc.gamma(nu-l)
	prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * gam)
	return prefac*W_scal(nu,l,r)
#Long range behavior of the Whittaker function.
def LRW(nu,l,r):
	gam = spc.gamma(nu-l)
	prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * gam)
	return prefac*LRW_scal(nu,l,r)
#Derivative of the Whittaker function. 
def Wp(nu,l,r):
	gam = spc.gamma(nu-l)
	prefac = nu**(3/2.) * 1/np.sqrt(nu**2 * spc.gamma(nu+l+1) * gam)
	return prefac*W_scalp(nu,l,r)



#Seaton-Burgess cutoff
def SB(l,r):
    return (1-np.exp(-10*r/(l*(l+1))))**(2*l+1) if l!=0 else (1-np.exp(-10*r/1.5))

def SBb(l,r):
    return (1-np.exp(-5*r/(l*(l+1))))**(2*l+1) if l!=0 else (1-np.exp(-5*r/1.5))

'''
#to test the new module compare with hydogenic orbitals for 1s, 4f and 6p
#r = np.linspace(1e-6,5*4**2,900)
#tot = np.zeros((900,4))
#for i in range(900):
#	tot[i] = [r[i],Unl(4,3,r[i]),W(4,3,r[i]),Wp(4,3,r[i])]
#np.save("1scom.npy",tot)
#print(nu)
#r = np.zeros(401)
#wtest = np.zeros((401,2))
#for i in range(401):	
#	wtest[i,0]=(Ucln(eau,l,i*5*8**2/400)*Wp(nu,l,i*5*8**2/400)-Uclpn(eau,l,i*5*8**2/400)*W(nu,l,i*5*8**2/400))-(-np.sin(np.pi*(nu-l))*2/np.pi)
#	wtest[i,1]=Ucln(eau,l,i*5*8**2/400)*Cclpn(eau,l,i*5*8**2/400)-Uclpn(eau,l,i*5*8**2/400)*Ccln(eau,l,i*5*8**2/400)-(2/np.pi)
#np.save("wtest.npy",np.array(wtest))
def Wtest():
	l = 5
	zion = 1.0
	eau = -0.5/(4**2)
	print(nuf(eau))
	nu = nuf(eau)
	r = np.linspace(1e-5,150,601)
	fogt = np.zeros((601,4))
	for i in range(601):
		[f,fp,g,gp] = [Uclno(eau,l,r[i]),Uclnop(eau,l,r[i]),W_scal(nuf(eau),l,r[i]),W_scalp(nuf(eau),l,r[i])]
		#[f,fp,g,gp,*args]
		fogt[i] = [r[i],f*gp-fp*g,-nu**l * 2/spc.factorial(l-nu),np.sqrt(nuf(eau))]
	np.savetxt("TestOutput/ftest.dat",fogt)

#Wtest()
#l=1
#eau = -0.5*(-1/1**2)
#np.sqrtnu = np.sqrt(-1/(2*eau))
#r = np.linspace(0,50,401)
#fogt = np.zeros((401,6))
#fgt = np.zeros((401,4))
#for i in range(401):
#	fogt[i] = gs.fogo(eau*2.0,l,r[i],1e-12)
#	fgt[i] = gs.seaton(l,eau*2.0,r[i],1)
#np.save("fogotbnu.npy",fogt)
#np.save("fgtbnu.npy",fgt)
# 

#Lets test a little bit the precision of our functions. We can do the test for the regular function to compare it with the analytic one. 
# Generate in a mesh for large r, define it at energy e = -0.5/n**2-0.0380 and we can test it for d waves and g waves.
n = 8
en = -0.5/n**2-0.0380
ldl = 2


#And for smaller r where the fortran code is known to be convergent
rtest = np.linspace(5,800,1200)
file = open("TestOutput/comparingfandg_lr.dat",mode='w',buffering=1)

for r in rtest:
	file.write("%10.8e   %10.8e   %10.8e \n"%(r ,Ucln_wLR(en,ldl,r), Ccln_wLR(en,ldl,r)))

file.close() 
# 
# 
# 
# '''