import numpy as np
import matplotlib.pyplot as plt




## Variables globales et constantes

n = 1000
l = 1
p0 = np.pi



kmax = 100
eps = 1e-3


h = 2*l/n

# Les vecteurs sont de taille n+1






## Fonctions de base


def long(f):
	l = 0
	for i in range(n):
		l += np.sqrt((f[i+1]-f[i])**2+h**2)
	return l




def m(f):
	return np.sum(np.sqrt((f[1:]+f[:-1])**2 + h**2))

def aire(f):
	return h * (np.sum(f) + (f[0]+f[-1])/2)



## Fonctions de l'algorithme d'Uzawa

def f_uza(f):
	return -aire(f)

def h_uza(f):
	return np.array([long(f)-p0, f[0], f[-1]])

def g_uza(f):
	return -f



## Lagrangien et gradient

def L(f,lam,mu):
	return f_uza(f) + np.sum(lam*h_uza(f)) + np.sum(mu*g_uza(f))

def gradL(f,lam,mu):
	u = np.ones_like(f)
	u[0] = .5
	u[-1] = .5
	e0 = np.zeros_like(f)
	e0[0] = 1
	en = np.zeros_like(f)
	en[-1] = 1

	gradH = np.zeros_like(f)
	gradH[0] = -(f[1]-f[0])/np.sqrt((f[1]-f[0])**2+h**2)
	for i in range(1,n):
		gradH[i] = (f[i]-f[i-1])/np.sqrt((f[i]-f[i-1])**2+h**2) - (f[i+1]-f[i])/np.sqrt((f[i+1]-f[i])**2+h**2)
	gradH[-1] = (f[-1]-f[-2])/np.sqrt((f[-1]-f[-2])**2+h**2)

	return -h*u + lam[0]*gradH + lam[1]*e0 + lam[2]*en - mu




## Algorithme d'Uzawa


def grad_pas_fixe(f0,lam,mu, pas=0.01):
	f = f0.copy()
	
	crt = 1+eps
	k = 0

	while (crt >= eps and k<kmax):
		f_tmp = f - pas * gradL(f,lam,mu)
		crt = np.linalg.norm(f-f_tmp)
		f[:] = f_tmp

		k+=1


	# print(k)

	return f


def uzawa(f0,lam0,mu0, rho=0.01):

	f = f0.copy()
	lam = lam0.copy()
	mu = mu0.copy()
	
	crt = 1+eps
	k = 0

	while (crt >= eps and k<kmax):

		# rhos = np.linspace(-1,1,200)
		# vals = np.zeros_like(rhos)
		# for i in range(len(rhos)):
		# 	vals[i] = L(f - rhos[i] * gradL(f,lam,mu), lam, mu)

		# plt.plot(rhos, vals)


		crt = eps/2

		f_tmp = grad_pas_fixe(f,lam,mu)
		lam_tmp = lam + rho*h_uza(f)
		mu_tmp = np.maximum(mu + rho*g_uza(f), np.zeros_like(f))

		crt = np.linalg.norm(f-f_tmp)
		k+=1

		f[:] = f_tmp
		lam[:] = lam_tmp
		mu[:] = mu_tmp

		print('Aire:', aire(f), '\tPer-p0:', long(f)-p0, end='\r')

	# plt.show()
	return f


x = np.linspace(-l,l,n+1)

# f0 = np.ones(n+1)
f0 = np.sqrt(l**2-x**2)
p0 = long(f0)
lam0 = np.ones(3)
mu0 = np.ones_like(f0)

sol = uzawa(f0,lam0,mu0)

plt.plot(x,sol)
plt.plot(x,np.sqrt(l**2-x**2),'r.')

plt.show()




def fonc():
	lgs = []
	for i in range(10,10000,100):
		x = np.linspace(-l,l,i)
		y = np.sqrt(l**2-x**2)

		lg = 0
		for i in range(i-1):
			lg += np.sqrt((y[i+1]-y[i])**2+h**2)

		lgs.append(lg)
	plt.plot(lgs)
	plt.show()


