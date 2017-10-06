# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 11:14:04 2017

@author: 3301269
"""
import re
import numpy as np
import pylab
import matplotlib.pyplot as plt


def methodexp(f):
    """
    Recupere la methode experimentale, le nombre de modeles (pour RMN) et la resolution (quand pas RMN ni solution) d'un fichier pdb.
    
    f -  flux ouvert de fichier pdb
    """
    line = f.readline()
    exp = re.compile("^EXPDTA[ ]*([A-Z\-]+ [A-Z]+)")
    model = re.compile("^NUMMDL[ ]*([0-9]+)")
    resol = re.compile("^REMARK[0-9 ]*RESOLUTION\.[ ]*([0-9]+\.[0-9]+) ANGSTROMS\.")
    experience = ""
    resolution = "Not applicable"
    modele = "Not applicable"
    while not(experience) and line != "" :
    	hits = exp.search(line)
    	if hits:
    		experience = hits.group(1)	
    		if "NMR" in hits.group(1):
    			line = f.readline()			
    			hi = model.search(line)
    			if hi:
    				modele = int(hi.group(1))
    	line = f.readline()			
    if not("SOLUTION" in experience) and not("NMR" in experience):
    	while resolution == "Not applicable" and line != "":
    		hit2 = resol.search(line)
    		if hit2:
    			resolution = float(hit2.group(1))
    		line = f.readline()
    return (experience, modele, resolution) #/!\ Nombre de modeles est une chaine de caracteres




def DBREFparser(f) :
	"""
	Recupere les informations de la ligne DBREF des fichiers pdb (informations de la sequence de reference de chaque chaine) et rend un dictionnaire
	avec les infos de DBREF {codePDB_chaine : (donnees DBREF)}.
    
	f - flux ouvert du fichier pdb
	"""
	dico = {}
	line = f.readline()
	dbref = re.compile("^DBREF")
	seqres = re.compile("^SEQRES") #SEQRES et SEQADV permettent d'arrêter le programme avant la fin (marqueur qui dit qu'on n'est plus sur DBREF)
	seqadv = re.compile("^SEQADV")
	espace = re.compile("[ ]+")
	ReadEnd = None # Arret de lecture tant que l'on ne trouve pas SEQRES et SEQADV
	while not(ReadEnd) and line != "":
		hits = dbref.search(line)
		if hits:
			chaine = line[12]
			if chaine == " ":
				chaine = "-"
			cle = line[7:11] + "_" + chaine #ID code de l'entree
			pdbini = espace.sub("", line[14:19]) #nombre de la seq initial pour ce segment de sequence PDB
			pdbend = espace.sub("", line[20:25]) #nombre final de la seq pour ce segment de sequence PDB
			db = espace.sub("", line[26:32]) #database
			dbcode = espace.sub("", line[33:41]) #dbAccession : code d'acces de la seq ds la base de donnees
			dbID = espace.sub("", line[42:54]) #ID de la seq ds la base de donnes
			seqini = espace.sub("", line[55:61]) #nombre initial du segment de la seq pour la base de donnees
			seqend = espace.sub("", line[62:68]) #nombre final du segment de la seq pour la base de donnees
			dico[cle] = ((int(pdbini), int(pdbend)), db, dbcode, dbID, (int(seqini), int(seqend)))
		elif seqres.search(line) or seqadv.search(line):
			ReadEnd = 1
		line = f.readline()
	return dico
  
  


def ATOMparser(modeles, chaines, f):
	"""
	Recupere donnees ATOM d'un fichier pdb. Rend un dictionnaire {codePDB_chaine : {numeroModele : [[infos d'atome1], [atome2], ...]}}.
 
	modeles - nombre de modeles dans le fichier
	chaines - liste avec le nom des chaines de la forme codePDB_chaine (cles )
  	f - flux ouvert du fichier pdb
	"""
      # Initialisation du dictionnaire
	datoms = {}
	for a in chaines:
		if modeles == "Not applicable":
			datoms[a] = {1:[]}
		else: #Liste vide pour chaque modele
			datoms[a] = {}
			for m in xrange(int(modeles)):
				datoms[a][m+1] = []
    
	line = f.readline()
	pdb1 = datoms.keys()[0].split("_")[0] #Code PDB
	model = re.compile("^MODEL")
	espace = re.compile("[ ]+")
	atom = re.compile("^ATOM")
	while line != "":
		hits = model.search(line)
		hit = atom.search(line)
		if hits:
			mod = int(espace.sub("",line[10:14])) #Numero du modele
		if hit:
			chaine = line[21]
			if chaine == " ":
				chaine = "-"
			atome = espace.sub("", line[6:11]) #numero atome
			nom = espace.sub("", line[13:16]) #nom atome CA...
			loca = espace.sub("-", line[16]) #localisation alternative
			resName = espace.sub("", line[17:20]) #nom residu a 3 lettres
			resNum = int(espace.sub("", line[22:26])) #numero de sequence du résidu
			insert = espace.sub("-", line[26]) #a quel numero de résidu d'aa cet atome appartient
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[46:54])
			occup = float(line[54:60]) #occupancy
			temp = float(line[60:66]) #facteur de temperature
			ele = espace.sub("", line[76:78]) #nom de l'element
			charge = espace.sub("", line[78:80]) #charge
			champAtoms = [nom, loca, resName, resNum, insert, x, y, z, occup, temp, ele, charge]
			cle = pdb1 + "_" + chaine
			if len(datoms[cle]) == 1:
				mod = 1
			datoms[cle][mod].append(champAtoms)
		line = f.readline()
	return datoms


##B
#def Selection_CA_Residus (P_1chaine, sel_P) :
#    """
#    Selection des x, y, z des CA (carbone alpha) des residus de la liste donnee en entree. Rend un liste de tuples des coord. des residus 
#    par ordre croissant d'index dans sel_P (qui n'est pas forcement l'ordre de ces indices dans sel_P)
#    
#    P_1chaine - pour 1 seul modele, 1 seule chaine, liste de listes d'infos de chaque atome pour P [[infos atome1], [infos atome2], ...]
#    sel_P - liste des indices des residus dans P selectionnes
#    """
#    CoordRes = []
#    ResIndex = 0
#    atom = 0
#    while ResIndex <= max(sel_P) and atom <= len(P_1chaine) :
#        print "RES", ResIndex
#        while P_1chaine[atom][0] != 'CA': #Cherche le carbone alpha du residu
#            atom += 1
#        print "ATOM3",P_1chaine[atom]
#        if ResIndex in sel_P : #Si on est dans la liste de residus selectionnes
#            CoordRes.append((P_1chaine[atom][5], P_1chaine[atom][6], P_1chaine[atom][7])) #Recupere les coordonnees du residu
#        atom += 1
#        ResIndex += 1
#    return CoordRes
            

def Selection_CA_Residus (P_1chaine, sel_P) :
    """
    Selection des x, y, z des CA (carbone alpha) des residus de la liste donnee en entree. 
    Rend un dictionnaire {numéro_résidu : tuples des coord. du CA du residu}
    
    P_1chaine - pour 1 seul modele, 1 seule chaine, liste de listes d'infos de chaque atome pour P [[infos atome1], [infos atome2], ...]
    sel_P - liste des indices des residus dans P selectionnes
    """
    for i in sel_P:
        assert P_1chaine[-1][3] >= i and i >= P_1chaine[0][3], "Numero de residu " + str(i) + " incorrect"
    CoordRes = {}
    ResIndex = 0
    atom = 0
    while len(CoordRes) < len(sel_P) and atom <= len(P_1chaine) :           
        while atom < len(P_1chaine) and P_1chaine[atom][0] != 'CA': #Cherche le carbone alpha du residu
            atom += 1
        ResIndex = P_1chaine[atom][3]
        if int(ResIndex) in sel_P : #Si on est dans la liste de residus selectionnes
            CoordRes[ResIndex] = (P_1chaine[atom][5], P_1chaine[atom][6], P_1chaine[atom][7]) #Recupere les coordonnees du residu
        atom += 1
    return CoordRes
            
  


def RMSD (P1, P2, sel_P1, sel_P2, modeleNum, modeleNum2, nomChaine, nomChaine2) :
    """
    Calcul du RMSD entre 2 structures de proteines ALIGNEES a partir d'une selection de residus.
    Le calcul ne se fait que sur les carbones alpha. 1 chaine 1 modele VS 1 chaine 1 modele
    
    P1 - structure de proteine P1 (dictionnaire de ATOMparser())
    P2 - structure de proteine P2 (dictionnaire de ATOMparser())
    sel_P1 - liste des indices des residus dans P1 selectionnes (resID)
    sel_P2 - liste des indices des residus dans P2 selectionnes (resID)
    modeleNum - numero du modele pour P1 pour faire le calcul de RMSD
    modeleNum2 - numero du modele pour P2 pour faire le calcul de RMSD
    nomChaine - nom de la chaine a comparer codePDB_chaineNum pour P1
    nomChaine2 - nom de la chaine a comparer codePDB_chaineNum pour P2
    """
    assert len(sel_P1) == len(sel_P2), "Nombre de residus incorrect"
    
    # Selection des bonnes coordonnees
    P1_1chaine = P1[nomChaine][modeleNum]
    print len(P1_1chaine)
    P2_1chaine = P2[nomChaine2][modeleNum2]
    CoordResP1 = Selection_CA_Residus (P1_1chaine, sel_P1)
    CoordResP2 = Selection_CA_Residus (P2_1chaine, sel_P2)
    print CoordResP1
    #Calcul du RMSD
    rmsd = 0.
    for i in xrange(len(sel_P1)):
        x1 = CoordResP1[sel_P1[i]][0]
        y1 = CoordResP1[sel_P1[i]][1]
        z1 = CoordResP1[sel_P1[i]][2]
        x2 = CoordResP2[sel_P2[i]][0]
        y2 = CoordResP2[sel_P2[i]][1]
        z2 = CoordResP2[sel_P2[i]][2]
        rmsd += ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)/len(sel_P1) #Prend que CA donc nombre d'atomes == nombre de residus
    return rmsd**0.5




#C
def Distance2Atomes (x1, y1, z1, x2, y2, z2):
    return ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5

def CarteContacts (P1, P2, sel_P1, sel_P2, modeleNum, modeleNum2, nomChaine, nomChaine2):
    
    assert len(sel_P1) == len(sel_P2), "Nombre de residus incorrect"
    
    CarteContact = np.zeros((len(sel_P1), len(sel_P1)))
    
    # Selection des bonnes coordonnees
    P1_1chaine = P1[nomChaine][modeleNum]
    P2_1chaine = P2[nomChaine2][modeleNum2]
    CoordResP1 = Selection_CA_Residus (P1_1chaine, sel_P1)
    CoordResP2 = Selection_CA_Residus (P2_1chaine, sel_P2)
    
    #Calcul des distances
    for i, res1 in enumerate(sel_P1):
        x1 = CoordResP1[res1][0]
        y1 = CoordResP1[res1][1]
        z1 = CoordResP1[res1][2]
        for j, res2 in enumerate(sel_P2):
            x2 = CoordResP2[res2][0]
            y2 = CoordResP2[res2][1]
            z2 = CoordResP2[res2][2]
            CarteContact[i, j] = Distance2Atomes(x1, y1, z1, x2, y2, z2)
    return CarteContact

def Dissimilarite (Contact1, Contact2):
    return


#D
def Rij(xi, yi, zi, xj, yj, zj):
	"""
	Calcul des coordonnees d'un vecteur Rij a partir des coordonnees tridimensionnelles de 2 points.
    
	xi, yi, zi - coord d'un point
	xj, yj, zj - coord d'un autre point
	"""
	xij = xj - xi
	yij = yj - yi
	zij = zj - zi
	return xij, yij, zij
 
def Norme(x, y, z):
    """
    Calcul de la norme d'un vecteur a partir de ses coordonnees.
    
    x, y, z - coordonnees du vecteur
    """
    return (x*x + y*y + z*z)**0.5
    
#D
def ListesCoord (datoms, nomChaine, modeleNum):
	"""
	Recupere du dictionnaire datoms les listes des X, Y, Z des atomes de tous les residus, ainsi que le numéro des residus dans la protéine et le nom de ces résidus.
    
	datoms - dictionnaire avec l'info des lignes ATOM des fichiers PDB (sortie de la fonction ATOMparser())
	nomChaine - nom de la chaine pour laquelle on applique la fonction (clé du dico datoms)
	modeleNum - numero du modele pour lequel on applique la fonction (clé du dico datoms[nomChaine])
	"""
	lesX = [] # x
	lesY = [] # y
	lesZ = [] # z
	AAnum = [] #numero des acides amines
	AA = [] # nom des acides amines
	P1_1chaine = datoms[nomChaine][modeleNum] 
	for j in xrange(len(P1_1chaine)) :
         x = P1_1chaine[j][5] 
         y = P1_1chaine[j][6]
         z = P1_1chaine[j][7]
         aanum = P1_1chaine[j][3]
         aa = P1_1chaine[j][2]
         lesX.append(x)
         lesY.append(y)
         lesZ.append(z)
         AAnum.append(aanum)
         AA.append(aa)
	return lesX, lesY, lesZ, AAnum, AA
 
 
def Env_i(rc, x, y, z, datoms, nomChaine, modeleNum):
    """
    Fonction qui donne les coordonnees des atomes qui sont a moins de rc de l'atome (x,y,z), c-a-d dans son environnement, sous forme de 3 listes.
    
    rc - rayon de l'enveloppe limite de l'atome, en Angstrom
    x, y, z - coordonnees d'un atome donné
    datoms - dictionnaire avec l'info des lignes ATOM des fichiers PDB (sortie de la fonction ATOMparser())
    nomChaine - nom de la chaine pour laquelle on applique la fonction (clé du dico datoms)
    modeleNum - numero du modele pour lequel on applique la fonction (clé du dico datoms[nomChaine])
    """
    Xenv = []
    Yenv = []
    Zenv = []
    X, Y, Z = ListesCoord(datoms, nomChaine, modeleNum)[0:3]
    for j in range(len(X)):
        if x != X[j] or y != Y[j] or z != Z[j]: # Atome donné pas pris en compte
            if Distance2Atomes(x, y, z, X[j], Y[j], Z[j]) < rc:
                Xenv.append(X[j])
                Yenv.append(Y[j])
                Zenv.append(Z[j])
    return Xenv, Yenv, Zenv


def CVi (rc, x, y, z, datoms, nomChaine, modeleNum) :
	"""
	Calcul de la variance circulaire d'un atome a partir de ses coordonnees 
	et de la liste des coordonnees de tous les atomes de la proteine.
    
	rc - rayon de l'enveloppe limite de l'atome, en Angstrom
	x, y, z - coordonnees d'un atome donné (i)
	datoms - dictionnaire avec l'info des lignes ATOM des fichiers PDB (sortie de la fonction ATOMparser())
	nomChaine - nom de la chaine pour laquelle on applique la fonction (clé du dico datoms)
	modeleNum - numero du modele pour lequel on applique la fonction (clé du dico datoms[nomChaine])
    	"""
	envX, envY, envZ = Env_i(rc, x, y, z, datoms, nomChaine, modeleNum)
	xVS = 0
	yVS = 0
	zVS = 0 #Coordonnees du vecteur de population d'atomes dans l'enveloppe de l'atome i
	for j in range(len(envX)):
		#Calcul du vecteur Rij allant de l'atome i à l'atome j
		xRij, yRij, zRij = Rij(x, y, z, envX[j], envY[j], envZ[j]) #Coordonnes du vecteur Rij
		normeRij = Norme(xRij, yRij, zRij)

		#Actualisation du vecteur de population avec le nouveau vecteur  
		xVS = xVS + xRij/normeRij
		yVS = yVS + yRij/normeRij
		zVS = zVS + zRij/normeRij
  
	CVi = 1 - Norme(xVS, yVS, zVS)/len(envX)
	return CVi

 
def CV (rc, datoms, nomChaine, modeleNum):
    """
    Calcul de CV pour chaque atome d'une protéine, étant donné rc, rend une liste de CV.
    
    rc - rayon de l'enveloppe limite de l'atome, en Angstrom
    datoms - dictionnaire avec l'info des lignes ATOM des fichiers PDB (sortie de la fonction ATOMparser())
    nomChaine - nom de la chaine pour laquelle on applique la fonction (clé du dico datoms)
    """
    X, Y, Z = ListesCoord(datoms, nomChaine, modeleNum)[0:3]
    CV = []
    
    for res in xrange(len(X)):
        CV.append(CVi(rc, X[res], Y[res], Z[res], datoms,  nomChaine, modeleNum))
    
    return CV

        
def CV_res(CV_atoms, ResNum, NumRes_liste):
    """
    CV d'un residu : Moyenne des CVi des atomes d'un residu. 
    
    CV_atoms - CV de chaque atome de la proteine
    ResNum - Numero du residu pour lequel on veut le CV. 
    NumRes_liste - la liste qui donne a quel aa appartient latome. 
    """
    cvi_res = 0
    no_atom_res = 0 #Nombre total d'atomes dans le residu
    i = 0 #indice
    while i < len(NumRes_liste) and NumRes_liste[i] <= ResNum :
        if NumRes_liste[i] == ResNum:
            cvi_res += CV_atoms[i]
            no_atom_res += 1
        i += 1
    return cvi_res/float(no_atom_res)
 
def CV_res_tous(CV_atoms, NumRes_liste):
    """
    CV de tous les residus d'une proteine
    
    CV_atoms - CV de chaque atome de la proteine
    NumRes_liste - la liste qui donne a quel aa appartient latome. 
    """
    CV_Res = []
    for res in sorted(list(set(NumRes_liste))):
        CV_Res.append(CV_res(CV_atoms, res, NumRes_liste))
    return CV_Res
    
def Niveau_enfouissement(x):
    return 0


if __name__ == "__main__":
    #f = open("/users/Etu9/3301269/Bureau/1cfc.pdb", "r")
    #f = open("/users/Etu9/3301269/Bureau/2bbm.pdb", "r")
    #f = open("/users/Etu9/3301269/Bureau/1cll.pdb", "r")
    
    f = open("3pdz.pdb", "r")
    experience, modeles, resolution = methodexp(f)
    d_3PDZ = DBREFparser(f)
    datoms_3PDZ = ATOMparser(modeles, d_3PDZ.keys(), f)
    
    f2 = open("1fcf_aliSeq.pdb", "r")
    nom_chaines_1FCF = ["Ali1FCF_A"]
    Num_modeles_1FCF = "Not applicable"
    datoms_1FCF = ATOMparser(Num_modeles_1FCF, nom_chaines_1FCF, f2)
    
    #Selection des residus alignes dans la figure de l'enonce
    sel3PDZ = range(21,25) + [26] + range(28,52) + range(53,69)
    sel1FCF = range(159, 164) + range(165, 179) + range(184, 210)
    
    print "La valeur de RMSD de l'alignement de 3PDZ et 1FCF est de", RMSD (datoms_3PDZ, datoms_1FCF, sel3PDZ, sel1FCF, 1, 1, "3PDZ_A", "Ali1FCF_A"), "Angstroms."
    # commande align sur pymol    
    #9 Angstrom c'est beaucoup. Alors que semblables sur pymol : similarité de sequence est faible entre les 2
    
    Contact_Matrix = CarteContacts(datoms_3PDZ, datoms_3PDZ, range(1,97), range(1,97), 1, 1, "3PDZ_A", "3PDZ_A")
    pylab.pcolor(Contact_Matrix)
    plt.colorbar()
    plt.title("Carte de contact de 3PDZ vs 3PDZ")
    plt.xlabel("Residus")
    plt.ylabel("Residus")
    
    CV_atom = CV(20, datoms_3PDZ, "3PDZ_A", 1)
    NumRes_list = ListesCoord(datoms_3PDZ, "3PDZ_A", 1)[3]
    print CV_res(CV_atom, 96, NumRes_list)
    CV_AAs = CV_res_tous(CV_atom, NumRes_list)
    
    f.close()
    f2.close()

##B
#(4)

#Non superposés. RMS =    1.958. Les deux résidus sérine de l'énoncé sensés etre alignés ne sont pas superposés. L'alignement par séquence n'est pas forcément éfficace.