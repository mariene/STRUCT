# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:03:31 2017

@author: 3202002
"""
import math
import numpy as np
import pylab
import matplotlib.pyplot as plt
import ForceField


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#A. Lecture d’un fichier PDB
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


def lire_pdb(nom_fichier="3pdz.pdb"):
    """ Permet de lire un fichier pdb et de stocker les donnees 
    
    Parameters
    ----------
    nom_fichier : str, nom du fichier ou chemin

    Returns
    -------
    liste : list, liste compose de liste contenant chaque mot de chaque ligne 
    """
    fichier = open(nom_fichier, "r")
    liste = []
    for i in fichier:
        #l = ' '.join(i.split(' ')).split()
        liste.append(i)
    fichier.close()
    return liste 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def details(liste):
    """
    
    Parameters
    ----------
    liste : list, liste obtenue avec la fonction lire_pdb()

    Returns
    -------
    dico :  dict, 
    
    """
    dico = dict()
    
    for l in liste:
        
        if l[0:6] == "NUMMDL":
            dico["nbmodele"] = int(l[10:14])
        if l[0:6] == "EXPDTA":
            #dico["methode"] = ' '.join(l[1:len(l)])
            dico["methode"] = l[10:80]
        if l[0:6] == "REMARK" : 
            if "RESOLUTION." in l :
                if 'NOT APPLICABLE' in l:
                    dico["resolution"] = None
                else :
                    pass
                    dico["resolution"] = l[26:40]
    return dico

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def lecture_dbref(liste):
    dico = dict()
    for l in liste:
        if l[0:5] == "DBREF":
            i = ' '.join(l.split(' ')).split()
            dico[l[7:11]] = i[2:len(i)]
    return dico
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def coordonnees(liste):
    modele = 0
    res = dict()
    dico = dict()
    dico2={}
    string=''
    for j in range(len(liste)):
        l=liste[j]
        
        if l[0:6]=="ENDMDL" :
            res["MODEL_"+str(modele)] = dico
            modele +=1
            dico = dict()
            dico2={}
        
        if l[0:4] == "ATOM" :
                i = ' '.join(l.split(' ')).split()
                
                if l[17:20] + "-"+l[22:26] == string :
                    dico2[i[2]]= map(float,i[6:9])
                    
                    
                if l[17:20] + "-"+l[22:26] != string :
                    
                    dico[string]=dico2
                    
                    string = l[17:20] + "-"+l[22:26]
                    dico2={}
                    dico2[i[2]]= map(float,i[6:9])
   
        res["MODEL_"+str(modele)] = dico

    return res


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def atome(liste):
    modele = 0
    res = dict()
    dico = dict()
    liste_res=[]

    for l in liste :

        if l[0:6]=="ENDMDL" :
            res["MODEL_"+str(modele)] = liste_res
            modele +=1
            dico = dict()
            liste_res=[]
            
            
        if l[0:4] == "ATOM" :
            dico={}
            dico["nomAtome"] =  l[12:16].replace(" ","") 
            dico["nomResidus"] = l[17:20]
            dico["numeroResidus"] = int(l[22:26])
            liste_res.append(dico)
  
        res["MODEL_"+str(modele)] = liste_res
       
    return res
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def main(nom_fichier = "1cll.pdb") :
    fichier = lire_pdb(nom_fichier)
    #print fichier
    #print details(fichier)  
    #print lecture_pdb_bis(fichier)    
    print coordonnees(fichier)
    #print atome(fichier)

#main()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def trie(atome, liste,sel):
    """Permet de recuperer les residus qui nous interessent
    
    Parameters
    ----------
    atome : str, nom de l'atome
    liste : list, liste de dictionnaire 
    sel : list, liste des positions 
    
    Returns
    -------
    res : list, liste des residus 
    
    """
    res = list()
    for i in liste : 
        if (i["nomAtome"] == atome) and (i["numeroResidus"] in sel): # and (i not in res) : 
            res.append(i)
    return res

      
def coord(dico,liste_coord):
    liste_residus = liste_coord.keys()
    for i in liste_residus : 
        tmp = i.split("-")
        #print tmp
        if tmp[0] == dico['nomResidus'] and int(tmp[1]) == dico['numeroResidus']:
            #print dico['nomResidus']
            #print i, liste_coord[i]
            n=dico['nomAtome']
            return liste_coord[i][n]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#B. Comparaison de deux structures protéiques
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


def RMSD(coordonnees_1, liste_atomes_1,coordonnees_2, liste_atomes_2,sel_p1,sel_p2, atome = "CA"):
    """Permet de calculer le RMSD
    
    Parameters
    ----------
    coordonnees_1 : dict, dictionnaire contenant les coordonnees des atomes de 
                    chaque residus pour une premiere proteine donnee. Ce dico 
                    est obtenu avec la fonction coordonnees() qui retourne le 
                    dico suivant {aa : {atome : coord}}
    liste_atomes_1 : list, liste de dictionnaire où chaque dictionnaire correspond 
                    à : {nomAtome : str,  nomResidus : str, numeroResidus : int}
                    pour la premieme proteine
    coordonnees_2 : dict, dictionnaire contenant les coordonnees des atomes de 
                    chaque residus pour une deuxieme proteine donnee. Ce dico 
                    est obtenu avec la fonction coordonnees() qui retourne le 
                    dico suivant {aa : {atome : coord}}
    liste_atomes_2 : list, liste de dictionnaire où chaque dictionnaire correspond 
                    à : {nomAtome : str,  nomResidus : str, numeroResidus : int}
                    pour la deuxieme proteine
    sel_p1 : list, liste des positions qui nous interesse dans la sequence 1
    sel_p2 : list, liste des positions qui nous interesse dans la sequence 2
    atome : str, nom de l'atome (par defaut ca sera CA)
    
    
    Return
    ------
        float, le RMSD 
    
    
    Test
    ----
    >>> sel3PDZ = range(21,25) + [26] + range(28,52) + range(53,69)
    >>> sel1FCF = range(159, 164) + range(165, 179) + range(184, 210)  
    
    >>> fichier1 = lire_pdb("3pdz.pdb")
    >>> fichier2 = lire_pdb("1fcf_aliSeq.pdb")
    
    >>> coord1 = coordonnees(fichier1)
    >>> coord1 = coord1["MODEL_0"]
    >>> at1 = atome(fichier1)
    >>> at1=at1["MODEL_0"]

    >>> coord2 = coordonnees(fichier2)
    >>> coord2 = coord2["MODEL_0"]
    >>> at2 = atome(fichier2)
    >>> at2=at2["MODEL_0"]
    >>> round(RMSD (coord1, at1,coord2,at2,sel3PDZ,sel1FCF),2)
    9.17
    
    """
    def calcul(P1,P2):
        """Fonction intermediaire de calcul
        
        Parameters
        ----------
        P1 : list, coordonnees d'un point
        P2 : list, coordonnees d'un point
        
        Returns
        -------
            float
        """
        
        return pow((P1[0]-P2[0]),2) + pow((P1[1]-P2[1]),2) + pow((P1[2]-P2[2]),2)
    
    def trie(atome, liste,sel):
        """Recupere les elements qui nous interessent
        
        """
        res = list()
        for i in liste : 
            if (i["nomAtome"] == atome) and (i["numeroResidus"] in sel): # and (i not in res) : 
                res.append(i)
        return res
                 
    somme = 0.0
    trier1 = trie(atome,liste_atomes_1,sel_p1)
    trier2 = trie(atome,liste_atomes_2,sel_p2)

    
    if len(trier1) == len(trier2):
        
        for i in range (len(trier1)):
            
            coord1 = coord(trier1[i],coordonnees_1)
            coord2 = coord(trier2[i],coordonnees_2)
            #print coord1,coord2
            #somme = somme + math.sqrt(calcul(coord1,coord2))
            somme = somme + calcul(coord1,coord2)
            #print i, somme
            #print calcul(coord1,coord2)
            
    #print math.sqrt(somme/len(trier1))
    #print somme  
    return math.sqrt(somme/len(trier1))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#RMSD (coord1, at1,coord2,at2,sel3PDZ,sel1FCF)

# Trouver la plus petite distance


#for t in range(30):
#    coord1 = coordonnees(fichier1)
#    coord1 = coord1["MODEL_"+str(t)]
#    at1 = atome(fichier1)
#    at1=at1["MODEL_"+str(t)]
#    Score = RMSD (coord1, at1,coord2,at2,sel3PDZ,sel1FCF)
#
#    if Score < s:
#        s=Score
#        i=t


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def calcul_distance(P1,P2):    
    """Permet de connaitre la distance entre deux points
    
    Parameters
    ----------
    P1 : list, coordonnees d'un point
    P2 : list, coordonnees d'un point
    
    Returns
    -------    
    float, distance entre 2 points
    
    """
    return math.sqrt(pow((P1[0]-P2[0]),2) + pow((P1[1]-P2[1]),2) + pow((P1[2]-P2[2]),2))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  C. Cartes de contacts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
def carte_contact(coordonnees):
    """
    retourne une matrice symétrique, 
    des distances entre les 'CA' de tous les aa d'une meme proteine.
    
    Parameters
    ----------
    coordonnees : dictionnary, {aa : { atome : [ x, y, z ]  }} 
    
    Returns
    -------
    matrix, distance entre aa
    
    
    """
    dim=(len(coordonnees.keys())+10)
    matrice=np.zeros((dim,dim))
    
    liste_aa_ord= [0] *(len(coordonnees.keys())+10) 


    for aa1 in coordonnees.keys() :
        a = aa1.split()
        if len(a)==2 :
            #print(type(int(a[1])))

            liste_aa_ord[int(a[1])]= aa1
    
    for i in range(1, len(liste_aa_ord)):
        for j in range(1, len(liste_aa_ord)):
            aa1=liste_aa_ord[i]
            aa2=liste_aa_ord[j]
            if aa1 != 0 :
                if aa2 != 0:
                    p1 = coordonnees[aa1]['CA']
                    p2 = coordonnees[aa2]['CA']
                    
                    matrice[i,j]=calcul_distance(p1,p2)
                                  
                                  
                                  
    # on enlève les colonnes et les lignes vides
    liste_vide=[]
                  
    for i in range(0,dim)  :
        if all(matrice[:,i]== [0] * dim ):
            liste_vide.append(i)
        
        
    m= np.delete(matrice, liste_vide, axis=0)
    m= np.delete(m, liste_vide, axis=1)
    
    
    return m


def dissimilarite():
    """
    Difference entre 2 cartes de contacts
    """
    pass



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  D. Variance circulaire
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def norme(P):
    """Permet de calculer la norme
    
    Parameters
    ----------
    P : list, coordonnees de vecteur 
    
    Returns
    -------
    float, la norme du vecteur
    """
    return math.sqrt(pow(P[0],2) + pow(P[1],2) + pow(P[2],2))
    
def vect(P1,P2):
    """ Permet de calculer les coord d'un vecteur 
    
    Parameters
    ----------
    P1 : list, coordonnees d'un point
    P2 : list, coordonnees d'un point
    """
    return [P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]]
    

def CV_i(rc,coord_i,coordonnees):
    """Permet de calculer la variance ciculaire d'un atome
    
    Parameters
    ----------
    coord_i : list, coordonnees d'un atome, par ex : [11.09, 1.768, 0.092]
    coordonnees : dict, {aa : {atome : coord}}
    rc : float, rayon rc
    
    Returns
    -------
    float, la variance circulaire 
    """
    #Recup tous les atomes ayant une distance inf a rc
    res_at = []
    
    for liste_at in coordonnees.keys() : 
        for atome in coordonnees[liste_at].keys():
            c = coordonnees[liste_at][atome]
            #print calcul_distance(c,coord_i)
            if calcul_distance(c,coord_i) < rc and (coord_i[0]!=c[0] and coord_i[1]!=c[1] and coord_i[2]!=c[2]) : 
                #print c
                res_at.append(c) 
            
            
    #Avec la liste d'atome recupere je calcule la variance circulaire
    somme_x = 0.
    somme_y = 0.
    somme_z = 0.     

    for i in res_at : 
        #print i
        tmp = vect(coord_i,i) #vecteur ij
        somme_x += (tmp[0]/norme(tmp))
        somme_y += (tmp[1]/norme(tmp))
        somme_z += (tmp[2]/norme(tmp))
    return 1 - (norme([somme_x,somme_y,somme_z])/len(res_at))
        

def CV_residus(rc,coordonnees):
    """Permet de calculer la variance ciculaire de chaque residus
    
    Cela correspond a la moyenne des variances circulaires de chaque atome de 
    chaque residus
    
    Parameters
    ----------
    coordonnees : dict, {aa : {atome : coord}}
    rc : float, rayon rc
    
    Returns
    -------
    cv_res : dict, {residus : Cv} 
    """
    res = dict()
    for i in coordonnees.keys() :
        tmp = dict()
        for atome in coordonnees[i].keys():
            tmp[atome] = CV_i(rc,coordonnees[i][atome],coordonnees)
        res[i]=tmp
    #print res
    
    cv_res = dict()

    for i in res.keys() :
        somme = 0.
        if i != '' : 
            for atome in res[i].keys():
                #print somme
                somme+=res[i][atome]
            
            cv_res[i] = somme/len(res[i].keys())
            #print "sortie : ", somme/len(res[i].keys())
    #print cv_res
    return cv_res
        
def pourcentage(residus):
    """Calcule le pourcentage de residus les plus enfouis et les plus protuberant
    
    Parameters
    ----------
    residus : dict
    
    Returns
    -------
    tuple, couple contenant le pourcentage de residus les plus enfouis et 
            les plus protuberant    
    
    Comments
    --------
    J'ai suppose le 0.5, en vrai je ne sais pas donc -> A REVOIR ! 
    
    """
    enfoui = 0.0
    protuberant = 0.0
    for i in residus.keys() :
        c = residus[i]
        if c > 0.5 : 
            protuberant +=1
        else :
            enfoui +=1
    print enfoui/len(residus.keys())
    print protuberant /len(residus.keys()) 
    
    return (enfoui/len(residus.keys()),protuberant /len(residus.keys()))
        
    
def ecriture_pdb(valeur,liste):
    """Permet d'ecrire un fichier pdb avec les valeurs donnees en entree
    
    """
    
    
    
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  E. Champ de force
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# dr : rayon
# dE : epsilon
# dcharge : charges




    
    
    
    
    
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  
sel3PDZ = range(21,25) + [26] + range(28,52) + range(53,69)
sel1FCF = range(159, 164) + range(165, 179) + range(184, 210)  

fichier1 = lire_pdb("3pdz.pdb")
fichier2=lire_pdb("1cll.pdb")
fichier3 = lire_pdb("1fcf_aliSeq.pdb")

coord1 = coordonnees(fichier1)
coord1 = coord1["MODEL_0"]
at1 = atome(fichier1)
at1=at1["MODEL_0"]

coord2 = coordonnees(fichier2)
coord2 = coord2["MODEL_0"]
at2 = atome(fichier2)
at2=at2["MODEL_0"]


    
d1=carte_contact(coord1)
plt.figure()
plt.axis(  [0,len(d1) ,0,len(d1)]    )
pylab.pcolor(d1)
plt.colorbar()
#plt.show()

d2=carte_contact(coord2)
plt.figure()
plt.axis(  [0,len(d2) ,0,len(d2)]    )
pylab.pcolor(d2)
plt.colorbar()


#distance (coord1, at1,coord2,at2,sel3PDZ,sel1FCF)
print CV_i(20.0,[11.09, 1.768, 0.092],coord1)

cv = CV_residus(20,coord1)
pourcentage(cv)
#
#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()








