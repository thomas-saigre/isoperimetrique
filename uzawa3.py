import numpy as np
import matplotlib.pyplot as plt






## Variables globales et constantes

n  = 1000
l  = 1
p0 = np.pi
a0 = np.pi / 2



kmax = 100
eps = 1e-3


h = 2*l/n

# Les vecteurs sont de taille n+1



## Classe fonction pour gérer les vecteurs de R^{n+1} avec des valeurs nulles au bord


class fonction:
	"""docstring for fonction"""
	def __init__(self, n):
		self.n = n
		self.vals = np.empty(n-1)	# les valeurs au bord sont 0

	def init_random(self):
		self.vals = np.random.random(self.n-1)

	def init_cst(self, c=1):
		self.vals = c * np.ones_like(self.vals);

	def init_vals(self, tab):
		self.vals = np.copy(tab)

	def aire(self):
		return h * np.sum(self.vals)

	def long(self):
		l = np.sqrt(self.vals[0]**2 + h**2)
		for i in range(n-2):
			l += np.sqrt((self.vals[i+1]-self.vals[i])**2+h**2)
		l += np.sqrt(self.vals[-1]**2 + h**2)
		return l

	def proj(self):
		self.vals = np.maximum(self.vals, np.zeros_like(self.vals))

	def plot(self):
		p = np.concatenate(([0],self.vals,[0]))
		x = np.linspace(-l,l,n+1)
		plt.plot(x,p)







## Fonctions de l'algorithme d'Uzawa

def f_uza(f):
	return f.long()
#	return -aire(f)

def h_uza(f):
	return f.aire() - a0




## Lagrangien du problème, et gradient


def L(f, lam):
	return f.long() + lam*(f.aire() - a0)
def L2(f, lam):
	fc = fonction(n)
	fc.init_vals(f)
	return fc.long() + lam*(fc.aire() - a0)

def gradL(f, lam):
	x = f.vals
	u = np.ones_like(f.vals)
	
	gradlg = np.zeros_like(f.vals)
	gradlg[0] = x[0] / np.sqrt(x[0]**2 + h**2) - (x[1]-x[0]) / np.sqrt((x[1]-x[0])**2 + h**2)
	for i in range(1,n-2):
		gradlg[i] = (x[i]-x[i-1]) / np.sqrt((x[i]-x[i-1])**2 + h**2) - (x[i+1]-x[i]) / np.sqrt((x[i+1]-x[i])**2 + h**2)
	gradlg[-1] = (x[-1]-x[-2]) / np.sqrt((x[-1]+x[-2])**2 + h**2) + x[-1] / np.sqrt(x[-1]**2 + h**2)

	return gradlg + lam * u




## Algorithme d'Uzawa

def grad_pas_fixe(f,lam,pas=0.01):

	x = fonction(n)
	x.init_vals(f.vals)
	
	crt = 1+eps
	k = 0

	while (crt >= eps and k<kmax):
		# x_tmp n'est juste qu'on vecteur, alors que x est une fonction
		x_tmp = x.vals - pas * gradL(x,lam)
		crt = np.linalg.norm(x.vals-x_tmp)
		x.init_vals(x_tmp)

		k+=1

	# print(k)

	return x


def uzawa(f0,lam,rho=0.01):

	x = fonction(n)
	x.init_vals(f0.vals)
	
	crt = 1+eps
	k = 0

	while (crt >= eps and k<5):

		# rhos = np.linspace(-.1,.1,200)
		# vals = np.zeros_like(rhos)
		# for i in range(len(rhos)):
			# vals[i] = L2(x.vals - rhos[i] * gradL(x,lam), lam)

		# plt.plot(rhos, vals, label=str(k))
		# crt = eps / 2


		x_tmp = grad_pas_fixe(x,lam)
		print(x_tmp.vals)
		x_tmp.proj()
		lam_tmp = lam + rho*h_uza(x)

		crt = np.linalg.norm(x.vals-x_tmp.vals)
		k+=1


		x.init_vals(x_tmp.vals)
		lam = lam_tmp

		print(k, "Per :", x.long(), "Aire - a0 :", np.abs(x.aire() - a0))

	# plt.legend()
	# plt.show()
	return x








## Exécution

x_cercle = np.linspace(-l,l,n+1)
y_cercle = np.sqrt(l**2 - x_cercle**2)

f0 = fonction(n)
f0.init_vals(y_cercle[1:-1])
# f0.init_cst()
lam0 = 1

sol = uzawa(f0,lam0)


# plt.plot(x_cercle, y_cercle, 'r.')
# sol.plot()
plt.show()
