# Résolution du problème de Didon
# Thomas Saigre & Romain Vallet
# M1 CSMI Année 2019-2020


import numpy as np
import matplotlib.pyplot as plt
from time import time
import sys
import getopt



## Variables globales et constantes

n  = 60
l  = 1



kmax = 2e5
eps = 1e-3


h = 2*l/(n+1)



# Les vecteurs sont de taille n+1

output = 'aff'
a0 = np.pi/2
nom = 'uzawa'
verbose = False





# Gestion des arguments

try:
    opts,args = getopt.getopt(sys.argv[1:],"hva:",["a0=","help","output="])
except getopt.GetoptError as err:
    print(err)
    print("usage : {} [options]".format(sys.argv[0]))
    sys.exit(2)


for o,a in opts:

    if o == "-a" or o == "--a0":
        a0 = float(a)

    elif o == "-h" or o == "--help":
        print("Résolution du problème iso-aire avec l'algorithme d'Uzawa")
        print("usage : {} [options]".format(sys.argv[0]))
        print("\t--a0=Val | -aVal\tmodifier la valeur de a0 (défaut pi/2)")
        print("\t-h | --help\t\tafficher l'aide")
        print("\t-v\t\t\tafficher les itérations")
        print("\t--output=nom\t\tPréfixe des fichiers pdf générés (par défaut, affichage dans la console)")
        sys.exit(0)

    elif o == "--output":
        output = "exp"
        nom = a

    elif o == "-v":
        verbose = True

    else:
        assert False, "Option non valide"












## Classe fonction pour gérer les vecteurs de R^{n+1} avec des valeurs nulles au bord


class fonction:

    def __init__(self, n):
        self.n = n
        self.vals = np.empty(n)	# les valeurs au bord sont 0

    def init_cst(self, c=1):
        self.vals = c * np.ones_like(self.vals)

    def init_vals(self, tab):
        self.vals = np.copy(tab)
    
    def get_vals(self, f):
        self.vals[:] = f.vals
        
    def __sub__(self,vct):
        assert(self.n == len(vct))
        dif = fonction(self.n)
        dif.init_vals(self.vals - vct)
        return dif

    def aire(self):
        return h * np.sum(self.vals)

    def long(self):
        l = np.sqrt(self.vals[0]**2 + h**2)
        for i in range(self.n-1):
            l += np.sqrt((self.vals[i+1]-self.vals[i])**2+h**2)
        l += np.sqrt(self.vals[-1]**2 + h**2)
        return l

    def proj(self):
        self.vals = np.maximum(self.vals, np.zeros_like(self.vals))
    


    def plot(self,fen):
        p = np.concatenate(([0],self.vals,[0]))
        x = np.linspace(-l,l,self.n+2)
        fen.plot(x,p,label="Solution calculée")

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
    for i in range(1,f.n-1):
        gradlg[i] = (x[i]-x[i-1]) / np.sqrt((x[i]-x[i-1])**2 + h**2) - (x[i+1]-x[i]) / np.sqrt((x[i+1]-x[i])**2 + h**2)
    gradlg[-1] = (x[-1]-x[-2]) / np.sqrt((x[-1]+x[-2])**2 + h**2) + x[-1] / np.sqrt(x[-1]**2 + h**2)

    return gradlg + lam * u






## Algorithme d'Uzawa


def uzawa(f0,lam,pas_grad=0.1,pas_uzawa=0.1):

    N = f0.n

    f = fonction(N)
    f.init_vals(f0.vals)

    f_m1 = fonction(N)
    f_m1.init_vals(2 * f0.vals)
    

    crt = 1+eps
    k = 0
    
    aires = []
    per = []
    lag = []

    
    while (crt >= eps and k<kmax):
        
        f_tmp = fonction(N)
        
        # Recherche du pas pour la variable primale

        pas = pas_grad
        cpt_rho = 0
        Delta_L = 1

        while (cpt_rho < 1e2) and (Delta_L >= 0):
            pas /= 1.3

            f_tmp.get_vals(f - pas * gradL(f, lam))
            Delta_L = L(f_tmp, lam) - L(f, lam)
            
            cpt_rho += 1

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



        f.get_vals(f_tmp)
        lam_tmp = lam + pasM*h_uza(f)

        crt = np.linalg.norm(f.vals-f_m1.vals)+np.linalg.norm(h_uza(f))
        
        k += 1

        

        f_m1.get_vals(f)
        lam = lam_tmp
        
        aires.append(f.aire())
        per.append(f.long())
        lag.append(L(f,lam))

        if verbose:
            print(k, "Per :", f.long(), "\tAire - a0 :", np.abs(f.aire() - a0), (cpt_rho, cpt_pas), "\t", crt)

    return f,aires,per,lag






## Exécution



x_cercle = np.linspace(-l,l,n+2)
y_cercle = np.sqrt(l**2 - x_cercle**2)


f0 = fonction(n)
# f0.init_vals(y_cercle[1:-1])
f0.init_cst(a0/(2*l))
lam0 = 1

t0 = time()
sol,aires,lg,lag = uzawa(f0,lam0)
t1 = time()

# Affichage des résultats :
chn = "(pi/2)" if a0==np.pi/2 else ''
print("Problème isoaire, avec a0 = {}".format(a0),chn)
print("Temps d'exécution : {}".format(t1-t0), "en {} itérations".format(len(aires)))
print("Aire finale : {}".format(aires[-1]))
print("Périmètre final : {}".format(lg[-1]))


if output == 'aff':

    fig, ax = plt.subplots(2,2)

    ax[0,0].plot(x_cercle, y_cercle, 'r--',label=r"Solution avc $a_0=\pi/2$")
    ax[0,0].set_title("Solution")
    sol.plot(ax[0,0])
    ax[0,0].legend()

    ax[0,1].plot(aires)
    ax[0,1].set_title("Aires")
    ax[0,1].plot([0,len(aires)],[a0,a0],label=r"$a_0$")
    ax[0,1].legend()

    ax[1,0].plot(lg)
    ax[1,0].set_title("Longueurs")

    ax[1,1].plot(lag)
    ax[1,1].set_title("Lagrangien")

    fig.tight_layout()
    plt.show()

elif output == "exp":

    figSol = plt.figure()
    plt.plot(x_cercle, y_cercle, 'r--',label=r"Solution avc $a_0=\pi/2$")
    plt.title("Solution")
    sol.plot(plt)
    plt.legend()
    figSol.savefig(nom+"_sol.pdf", bbox_inches='tight')

    figAires = plt.figure()
    plt.plot(aires)
    plt.title("Aires")
    plt.plot([0,len(aires)],[a0,a0],label=r"$a_0$")
    plt.legend()
    figAires.savefig(nom+"_aires.pdf", bbox_inches='tight')

    figLon = plt.figure()
    plt.plot(lg)
    plt.title("Longueurs")
    figLon.savefig(nom+"_longueurs.pdf", bbox_inches='tight')

    figLag = plt.figure()
    plt.plot(lag)
    plt.title("Lagrangien")
    figLag.savefig(nom+"_lagrangien.pdf", bbox_inches='tight')
