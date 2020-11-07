# Résolution du problème isopérimétrique

Travail effectué par [Thomas Saigre](https://github.com/thomas-saigre) et [Romain Vallet](https://github.com/ValletRomain) dans le cadre de l'UE optimisation en M1 CSMI (année 2019 - 2020)


## Rapport et support de présentation

* Résumé du travail effectué : [sur ce lien](LaTeX/resume.pdf)
* Support de présentation : [ici](LaTeX/diapo.pdf)


## Codes

* Le programme Python [didon.py](didon.py) utilise l'algorithme d'Uzawa projeté avec le Lagrangien vu en cours
* Le programme Python [aug.py](aug.py) utilise aussi cet algorithme, mais utilise une version augmentée du Lagrangien.

La commande à entrer dans le terminal pour faire tourner ces programmes est la suivante :

```
python3 {didon|aug}.py [options]
```

avec ces options :
```
--a0=Val | -aVal    modifier la valeur de a0 (défaut pi/2)
--b=Val | -bVal     modifier la valeur de b (défaut 1
-h | --help		    afficher l'aide
-v			        afficher les itérations
--output=nom		Préfixe des fichiers pdf générés (par défaut, affichage dans la console)
```
