# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:03:31 2017

@author: 3202002
"""
#Wesh

def lire_pdb(nom_fichier="3pdz.pdb"):
    """
    @param nom_fichier : str, nom du fichier ou chemin
    @return liste : list, liste compose de liste contenant chaque mot de chaque 
    ligne 
    """
    fichier = open(nom_fichier, "r")
    liste = []
    for i in fichier:
        l = ' '.join(i.split(' ')).split()
        liste.append(l)
    fichier.close()
    return liste 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def details(liste):
    """
    @param liste
    """
    dico = dict()
    
    for l in liste:
        
        if l[0] == "NUMMDL":
            dico["nbmodele"] = int(l[1])
        if l[0] == "EXPDTA":
            dico["methode"] = ' '.join(l[1:len(l)])
        if l[0] == "REMARK" : 
            if "RESOLUTION." in l :
                if ' '.join(l[l.index("RESOLUTION.")+1:len(l)]) == 'NOT APPLICABLE.':
                    dico["resolution"] = None
                else :
                    dico["resolution"] = float(l[l.index("RESOLUTION.")+1])
    return dico

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def lecture_pdb_bis(liste):
    dico = dict()
    for l in liste:
        if l[0] == "DBREF":
            dico[l[1]] = l[2:len(l)]
    return dico
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def coordonnees(liste):
    modele = 0
    res = dict()
    j =0
    while j != (len(liste)):
        dico = dict()
        
        if liste[j][0] == "ATOM" :
            while liste[j][0] == "ATOM" :   
                #dico [' '.join(liste[j][0:2])] = map(float,liste[j][6:9])
                string = liste[j][3] + "-"+liste[j][5]
                dico [string] = map(float,liste[j][6:9])
                j+=1
            res["modele"+str(modele)] = dico
            modele +=1
        j+=1
    
    return res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def atome(liste):
    res = []
    for l in liste :
        dico = dict()
        if l[0] == "ATOM" :
            dico["nomAtome"] = l[2]
            dico["nomResidus"] = l[3]
            dico["numeroResidus"] = l[5]
            res.append(dico)
    return res
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def main(nom_fichier = "1cll.pdb") :
    fichier = lire_pdb(nom_fichier)
    #print fichier
    #print details(fichier)  
    #print lecture_pdb_bis(fichier)    
    print coordonnees(fichier)
    #print atome(fichier)

main()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RMSD(coordonnees, liste_atomes):
    """
    J en sais rien !
    """
    def calcul(P1,P2):
        return pow((P1[0]-P2[0]),2)+pow((P1[1]-P2[1]),2)+pow((P1[2]-P2[2]),2)
    
    def trie(atome, liste):
        res = list()
        for i in liste : 
            if i["nomAtome"] == atome : 
                res.append(i)
        return res
    
    trier = trie("CA",liste_atomes)
    
        
    
fichier = lire_pdb("1cll.pdb")
coord = coordonnees(fichier)
at = atome(fichier)
RMSD (coord, at)
       
