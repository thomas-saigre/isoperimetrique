import numpy as np
import matplotlib.pyplot as plt








## Classe fonction pour gérer les vecteurs de R^{n+1} avec des valeurs nulles au bord


class fonction:

    def __init__(self, n):
        self.n = n
        self.vals = np.empty(n-1)	# les valeurs au bord sont 0

    def init_random(self):
        self.vals = np.random.random(self.n-1)

    def init_cst(self, c=1):
        self.vals = c * np.ones_like(self.vals);

    def init_vals(self, tab):
        self.vals = np.copy(tab)
    
    def get_vals(self, f):
        self.vals[:] = f.vals
        
    def __sub__(self,vct):
        assert(self.n == 1+len(vct))
        dif = fonction(self.n)
        dif.init_vals(self.vals - vct)
        return dif

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

    def plot(self,fen):
        p = np.concatenate(([0],self.vals,[0]))
        x = np.linspace(-l,l,n+1)
        fen.plot(x,p)







## Fonctions de l'algorithme d'Uzawa

def f_uza(f):
    return f.long()

def h_uza(f):
    return f.aire() - a0




## Lagrangien du problème, et gradient


def L(f, lam):
    return f.long() + lam*h_uza(f)


def gradL(f, lam):
    x = f.vals
    u = np.ones_like(f.vals)

    gradlg = np.zeros_like(f.vals)
    gradlg[0] = x[0] / np.sqrt(x[0]**2 + h**2) - (x[1]-x[0]) / np.sqrt((x[1]-x[0])**2 + h**2)
    for i in range(1,n-2):
        gradlg[i] = (x[i]-x[i-1]) / np.sqrt((x[i]-x[i-1])**2 + h**2) - (x[i+1]-x[i]) / np.sqrt((x[i+1]-x[i])**2 + h**2)
    gradlg[-1] = (x[-1]-x[-2]) / np.sqrt((x[-1]+x[-2])**2 + h**2) + x[-1] / np.sqrt(x[-1]**2 + h**2)

    return gradlg + lam * u


## Lagrangien augmenté

def L_aug(f, lam, b, c):
    return f.long() + lam * (h_uza(f) + c) + b/2 * (h_uza(f) + c)**2

def gradL_aug(f, lam, b, c):
    x = f.vals
    u = np.ones_like(f.vals)

    gradlg = np.zeros_like(f.vals)
    gradlg[0] = x[0] / np.sqrt(x[0]**2 + h**2) - (x[1]-x[0]) / np.sqrt((x[1]-x[0])**2 + h**2)
    for i in range(1,n-2):
        gradlg[i] = (x[i]-x[i-1]) / np.sqrt((x[i]-x[i-1])**2 + h**2) - (x[i+1]-x[i]) / np.sqrt((x[i+1]-x[i])**2 + h**2)
    gradlg[-1] = (x[-1]-x[-2]) / np.sqrt((x[-1]+x[-2])**2 + h**2) + x[-1] / np.sqrt(x[-1]**2 + h**2)


    return gradlg + (lam + b*(h_uza(f)+c)) * u




## Algorithme d'Uzawa


def uzawa(f0,lam,pas_grad=0.2,pas_uzawa=0.1):

    f = fonction(n)
    f.init_vals(f0.vals)

    f_m1 = fonction(n)
    f_m1.init_vals(2 * f0.vals)
    

    crt = 1+eps
    k = 0
    
    aires = []
    per = []
    lag = []

    
    while (crt >= eps and k<kmax):
        
        f_tmp = fonction(n)
        
        # Recherche du pas pour la variable primale

        pas = pas_grad
        cpt_rho = 0
        Delta_L = 1

        while (cpt_rho < 1e2) and (Delta_L >= 0):
            pas /= 1.3
            f_tmp.get_vals(f - pas * gradL(f, lam))
            Delta_L = L(f_tmp, lam) - L(f, lam)
            # f_tmp.get_vals(f - pas * gradL_aug(f, lam, b, c))
            # Delta_L = L_aug(f_tmp, lam, b, c) - L_aug(f, lam, b, c)
            
            cpt_rho += 1
        if (cpt_rho == 1e2): print("Nombre maximal d'itération atteint")

        f_tmp.proj()
        

        # Recherche du pas pour la variable duale

        cpt_pas = 0
        pasM = pas_uzawa
        Delta_L = -1

        while (cpt_pas < 1.e2) and (Delta_L <= 0):
            pasM /= 1.2
            lam_tmp = lam + pasM * h_uza(f_tmp)
            Delta_L = L(f_tmp,lam_tmp) - L(f_tmp,lam)
            cpt_pas += 1

        if (cpt_pas == 1e2): print("Nombre maximal atteint")



        
        f.get_vals(f_tmp)
        lam_tmp = lam + pasM*h_uza(f)

        crt = np.linalg.norm(f.vals-f_m1.vals)

        
        k += 1

        

        f_m1.get_vals(f)
        # f.get_vals(f_tmp)
        lam = lam_tmp
        
        aires.append(f.aire())
        per.append(f.long())
        lag.append(L(f,lam))
        # lag.append(L_aug(f,lam,b,c))

        print(k, "Per :", f.long(), "Aire - a0 :", np.abs(f.aire() - a0), (cpt_rho, cpt_pas), crt)

    return f,aires,per,lag





## Variables globales et constantes

n  = 1000
l  = 1
p0 = np.pi
# a0 = np.pi / 2
a0 = 0.5



kmax = 1000
eps = 1e-3


h = 2*l/n

b = 1
c = 1

# Les vecteurs sont de taille n+1




## Exécution

x_cercle = np.linspace(-l,l,n+1)
y_cercle = np.sqrt(l**2 - x_cercle**2)

f0 = fonction(n)
f0.init_vals(y_cercle[1:-1])
# f0.init_cst()
lam0 = 1

sol,aires,long,lag = uzawa(f0,lam0)


fig, ax = plt.subplots(2,2)

ax[0,0].plot(x_cercle, y_cercle, 'r--')
ax[0,0].set_title("Solution")
sol.plot(ax[0,0])

ax[0,1].plot(aires)
ax[0,1].set_title("Aires")
ax[0,1].plot([0,len(aires)],[a0,a0])

ax[1,0].plot(long)
ax[1,0].set_title("Long")

ax[1,1].plot(lag)
ax[1,1].set_title("Lagrangien")

fig.tight_layout()

plt.show()
