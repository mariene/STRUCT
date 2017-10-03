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
        #l = ' '.join(i.split(' ')).split()
        liste.append(i)
    fichier.close()
    return liste 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def details(liste):
    """
    @param liste
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
                #if ' '.join(l[l.index("RESOLUTION.")+1:len(l)]) == 'NOT APPLICABLE.':
                if 'NOT APPLICABLE' in l:
                    dico["resolution"] = None
                else :
                    dico["resolution"] = float(l[l.index("RESOLUTION.")+1])
    return dico

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def lecture_pdb_bis(liste):
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
    for j in range(len(liste)):
        
        #dico = dict()
        l=liste[j]
        if l[0:6]=="ENDMDL" :
            res["MODEL_"+str(modele)] = dico
            modele +=1
            dico = dict()
            
            
        if l[0:4] == "ATOM" :
            
                i = ' '.join(l.split(' ')).split()
                string = l[17:20] + "-"+l[22:26]
                dico[string] = map(float,i[6:9])
                

            

                
        res["MODEL_"+str(modele)] = dico
        #modele +=1
        #dico = dict()

    
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


def coord(dico,liste_coord):
    liste_residus = liste_coord.keys()
    for i in liste_residus : 
        tmp = i.split("-")
        #print tmp
        if tmp[0] == dico['nomResidus'] and int(tmp[1]) == dico['numeroResidus']:
            print i, liste_coord[i]
            return liste_coord[i]




def RMSD(coordonnees_1, liste_atomes_1,coordonnees_2, liste_atomes_2,sel_p1,sel_p2):
    
    def calcul(P1,P2):
        return pow((P1[0]-P2[0]),2) + pow((P1[1]-P2[1]),2) + pow((P1[2]-P2[2]),2)
    
    def trie(atome, liste,sel):
        res = list()
        for i in liste : 
            if (i["nomAtome"] == atome) and (i["numeroResidus"] in sel): # and (i not in res) : 
                res.append(i)
        return res
    
#    def dico_correp(numero,liste_dico):
#        for i in liste_dico : 
#            if eval(i['numeroResidus']) == numero : 
#                return i
                
    def coord(dico,liste_coord):
        liste_residus = liste_coord.keys()
        for i in liste_residus : 
            tmp = i.split("-")
            #print tmp
            if tmp[0] == dico['nomResidus'] and int(tmp[1]) == dico['numeroResidus']:
                print i, liste_coord[i]
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
    
    return trier1, trier2
            
            
        
sel3PDZ = range(21,25) + [26] + range(28,52) + range(53,69)
sel1FCF = range(159, 164) + range(165, 179) + range(184, 210)  

fichier1 = lire_pdb("3pdz.pdb")
fichier2 = lire_pdb("1fcf_aliSeq.pdb")

coord1 = coordonnees(fichier1)
coord1 = coord1["MODEL_0"]
at1 = atome(fichier1)
at1=at1["MODEL_0"]



coord2 = coordonnees(fichier2)
coord2 = coord2["MODEL_0"]
at2 = atome(fichier2)
at2=at2["MODEL_0"]

t1,t2 = RMSD (coord1, at1,coord2,at2,sel3PDZ,sel1FCF)
#       
