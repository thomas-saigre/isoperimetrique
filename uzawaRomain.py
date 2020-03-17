import numpy as np
import matplotlib.pyplot as plt


class grid:
    
    def __init__(self, n, l=1, a0=np.pi):
        self.n = n
        self.l = l
        self.a0 = a0
        
        self.x = np.linspace(-self.l, self.l, n+2)[1:-1]
        self.h = 2*self.l / (n+2)
    
    def display(self, f):
        plt.plot(self.x, f)
        plt.title("Résultat pour n={0}".format(self.n))
        plt.show()

#------------------------------------------------------------------------------
# Gradient à pas variable
#------------------------------------------------------------------------------

def lg(grid, f):
    somme = np.sqrt(f[0]**2+grid.h**2)
    for i in range(1, len(f)-1):
        somme += np.sqrt((f[i+1]-f[i])**2 + grid.h**2)
    somme += np.sqrt(f[-1]**2+grid.h**2)
    
    return somme

def A(grid, f):
    return grid.h * sum(f)

def L(grid, f, lamb):
    return lg(grid, f) + lamb*(A(grid, f)-grid.a0)

"""
n = 60
grille = grid(n)
f = np.ones(n)
lamb = 1
print(grille.h)
print(lg(grille, f))
print(A(grille, f))
print(L(grille, f, lamb))
"""

def gradL(grid, f, lamb):
    n = grid.n
    h = grid.h
    
    grad_lg = np.zeros(n)
    
    grad_lg[0] = f[0] / np.sqrt(f[0]**2+h**2) - (f[1]-f[0])/np.sqrt((f[1]-f[0])**2+h**2)
    
    for i in range(1, n-1):
        grad_lg[i] = (f[i]-f[i-1])/np.sqrt((f[i]-f[i-1])**2+h**2) - (f[i+1]-f[i])/np.sqrt((f[i+1]-f[i])**2+h**2)
    
    grad_lg[-1] = (f[-1]-f[-2])/np.sqrt((f[-1]-f[-2])**2+h**2) + f[-1] / np.sqrt(f[-1]**2+h**2)

    return grad_lg + lamb * h*np.ones(n)
"""
n = 60
grille = grid(n)
#f = np.ones(n)
f = np.sqrt(1 - np.linspace(-1, 1, n+2)[1:-1]**2)
lamb = 1
gradient = gradL(grille, f, lamb)
print(gradient)
grille.display(gradient)
"""

def proj(f):
    return np.maximum(f, 0)

def rech_rho1(grid, f, lamb, kmax=100):
    rho = 1
    Deltalag = L(grid, f-rho*gradL(grid, f, lamb), lamb) - L(grid, f, lamb)
    
    k = 0
    while (Deltalag > 0 and k < kmax):
        rho /= 1.3
        Deltalag = L(grid, f-rho*gradL(grid, f, lamb), lamb) - L(grid, f, lamb)
        k += 1
    #print(k)
    return rho
"""
n = 60
grille = grid(n)
f = 0.25*np.ones(n)
lamb = 1
rho1 = rech_rho1(grille, f, lamb)
print(rho1)
"""
def rech_rho2(grid, f, lamb, kmax=100):
    rho = 1
    
    k = 0
    while (lamb+rho*(A(grid, f)-grid.a0) < lamb and k < kmax):
        rho /= 1.3
        k += 1
    
    return rho

temps = [0, 10, 50, 100, 500, 1000, 2000, 5000, 7500, 9990]

def grad(grid, f0, lamb, eps=0.01, kmax=1000):
    
    Rho = []
    compt = 0
    
    crt = 1+eps
    k = 0
    while (crt >= eps and k<kmax):
        rho1 = rech_rho1(grid, f, lamb)
        Rho.append(rho1)
        f1 = f0 - rho1 * gradL(grid, f0, lamb)
        
        crt = np.linalg.norm(f1-f0)
        f0 = f1
        """
        if compt<len(temps) and k == temps[compt]:
            grid.display(f0)
            compt += 1
        """
        k += 1
    """
    plt.plot(np.arange(k), Rho)
    plt.title("rho")
    plt.show()
    #print(Rho)
    #print("k :", k)
    """
    #print("gradient k :", k)
    return f0
"""
n = 60
grille = grid(n)
f = np.ones(n)
lamb = 1
liste = grad(grille, f, lamb)
#print("deuxieme fois \n \n")
#liste = grad(grille, liste, lamb)
"""
#------------------------------------------------------------------------------
# Uzawa
#------------------------------------------------------------------------------

temps = [0, 5, 10, 20, 50, 100]

def uzawa(grid, f, eps=0.01, kmax=1000):
    #Rho = []
    f0 = f.copy()
    lamb0 = 0
    
    crt = 1+eps
    k = 0
    compt = 0
    while (crt >= eps and k<kmax):
        f1 = grad(grid, f0, lamb0)
        #f1 = proj(f1)
        rho2 = rech_rho2(grid, f, lamb0)
        lamb1 = lamb0 + rho2*(A(grid, f)-grid.a0)
        
        crt = np.linalg.norm(f1-f0)
        f0 = f1
        lamb0 = lamb1
        """
        if k == temps[compt]:
            grid.display(f0)
            compt += 1
        """
        k += 1
    #print("uzawa k :", k)
    return f0

n = 60
grille = grid(n)
f = np.ones(n)
#f = np.sqrt(1-np.linspace(-1,1,n+2)[1:-1]**2)
lamb = 1
liste = uzawa(grille, f)
grille.display(liste)
