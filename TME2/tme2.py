# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:03:31 2017

@author: 3202002
"""
import math

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

#main()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RMSD(coordonnees_1, liste_atomes_1,coordonnees_2, liste_atomes_2,sel_p1,sel_p2):
    
    def calcul(P1,P2):
        return pow((P1[0]-P2[0]),2)+pow((P1[1]-P2[1]),2)+pow((P1[2]-P2[2]),2)
    
    def trie(atome, liste,sel):
        res = list()
        for i in liste : 
            if (i["nomAtome"] == atome) and (eval(i["numeroResidus"]) in sel) and (i not in res) : 
                res.append(i)
        return res
        
    def coord(dico,liste_coord):
        liste_residus = liste_coord.keys()
        for i in liste_residus : 
            tmp = i.split("-")
            if tmp[0] == dico['nomResidus'] and tmp[1] == dico['numeroResidus']:
                return liste_coord[i]
                
    somme = 0
    trier1 = trie("CA",liste_atomes_1,sel_p1)
    trier2 = trie("CA",liste_atomes_2,sel_p2)
    
    if len(trier1) == len(trier2):
        for i in range (len(trier1)):
            coord1 = coord(trier1[i],coordonnees_1)
            coord2 = coord(trier2[i],coordonnees_2)
            somme = somme + calcul(coord1,coord2)
            #print calcul(coord1,coord2)
    print math.sqrt(somme/len(trier1))
            
            
        
sel3PDZ = range(21,25) + [26] + range(28,52) + range(53,69)
sel1FCF = range(159, 164) + range(165, 179) + range(184, 210)  

fichier1 = lire_pdb("3pdz.pdb")
fichier2 = lire_pdb("1fcf_aliSeq.pdb")

coord1 = coordonnees(fichier1)
coord1 = coord1["modele0"]
at1 = atome(fichier1)

coord2 = coordonnees(fichier2)
coord2 = coord2["modele0"]
at2 = atome(fichier2)
RMSD (coord1, at1,coord2,at2,sel3PDZ,sel1FCF)
       
