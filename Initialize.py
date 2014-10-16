# This File contains the initialization functions for  "The Race up Mount Improbable" :
#
#  1. create_monomers
#
# copyright Sara Imari Walker 2011



import math
import numpy as np
from operator import itemgetter
import random
import matplotlib.pyplot as plt

from Polymers import *
import itertools


###############################################################################

def initialize_monomers(sequences):
    """This function initializes the first m ID's with the m-monomer species present in the model"""

    import Parameters

    from Parameters import m, M_N

    Parameters.Nmono = 0

    for i in range(0, m):
        """Loop over monomer species in system"""

        ID = Parameters.tot_species
        
        species = Parameters.monomer_species[i]

        sequences[ID] = Polymer(ID, M_N[i], 1, species, [], i, i)
        
        Parameters.Ntot += M_N[i]
        Parameters.tot_species += 1

    
###############################################################################  
def update_propensities(sequences, catalysts):
    """Initializes reaction propensities"""

    import Parameters

    from Parameters import m, kp, kh, kc 
    from Parameters import Ap_p, Ap_h

    """Compute and Update Reaction Propensities"""
    
    Ap_p = 0
    Ap_h = 0
    Ap_cp = 0
    Ap_ch = 0

    for ID1, seq1 in sequences.items():
        """Loop over 3' ends of all sequences"""
        ID1_e3 = seq1.three_prime
       
        for ID2, seq2 in sequences.items():
            """Loop over 5' ends of all sequences"""

            ID2_e5 = seq2.five_prime
            
            """ 3' of seq 1 binds with 5' of seq 2"""
            if ID1 == ID2 and sequences[ID1].tot > 1:
                Ap_p += kp[ID1_e3][ID2_e5]*sequences[ID1].tot*(sequences[ID1].tot - 1)
            else:
                Ap_p += kp[ID1_e3][ID2_e5]*sequences[ID1].tot*sequences[ID2].tot


            bond = str(ID2_e5) + str(ID1_e3)
        
            if bond in catalysts.keys():

                for IDc in catalysts[bond]:

                    if ID1 == ID2 and ID1 == IDc and sequences[ID1].tot > 2:
                        """Both substrates are the catalyst"""
                        Ap_cp += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*(sequences[IDc].tot - 1)*(sequences[IDc].tot - 2)*sequences[IDc].tot

                    elif ID1 == ID2 and ID1 != IDc and sequences[ID1].tot > 1:
                        """Both substrates are the same, but differ from the catalyst"""
                        Ap_cp += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID1].tot*(sequences[ID1].tot - 1)*sequences[IDc].tot

                    elif ID1 != ID2 and ID1 == IDc and sequences[IDc] > 1:
                        """Substrate 1 is the same as the catalyst, but differs from substrate 2"""
                        Ap_cp += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID2].tot*(sequences[IDc].tot - 1)*sequences[IDc].tot

                    elif ID1 != ID2 and IDc == ID2 and sequences[IDc] > 1:
                        """Substrate 2 is the same as the catalyst, but differs from substrate 1"""
                        Ap_cp += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID1].tot*(sequences[ID2].tot - 1)*sequences[IDc].tot

                    elif ID1 != ID2 and ID2 != IDc != ID1:
                        """Catalyst is different from both substrates"""
                        Ap_cp += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID1].tot*sequences[ID2].tot*sequences[IDc].tot


                    

        """
           for c in catalysts[bond]:

              print 'Id = ' + str(c)
              print 'seq = ' + str(sequences[c].sequence)
              print 'F = ' + str(sequences[c].F)
              print 'Func = ' + str(sequences[c].F_motifs)
              print 'bonds = ' + str(sequences[c].B_domains)
        """


        """Calculate hydrolysis rate per polymer"""

        for i in range(seq1.nb):

            e3 = seq1.bond_struct[i][1]
            e5 = seq1.bond_struct[i][0]
            
            Ap_h += kh[e3][e5]*seq1.tot  #THIS CAN BE TIGHTENED UP BY CACULATING THE TOTAL RATE OF BONDS BREAKING AS SUM OVER BONDS THEN MULTIPLYING BY CONC each timestep

            bond = str(e5) + str(e3)

            
            if Parameters.cat_reverse == True and bond in catalysts.keys():

                for IDc in catalysts[bond]:

                    if IDc == ID1 and sequences[ID1].tot > 1:
                        #Substrate is the catalyst
                        Ap_ch += kc[e3][e5]*kh[e3][e5]*(seq1.tot - 1)

                    else:
                        #Catalyst is different from substrate
                        Ap_ch += kc[e3][e5]*kh[e3][e5]*seq1.tot
           
            
   
    Atot = Ap_p + Ap_h + Ap_cp + Ap_ch

    Parameters.Ap_p = Ap_p
    Parameters.Ap_h = Ap_h
    Parameters.Atot = Atot

    Parameters.Ap_cp = Ap_cp
    Parameters.Ap_cp = Ap_ch

    #print 'Ap_p is ' + str(Ap_p)
    #print 'Ap_h is ' + str(Ap_h)
    #print Atot

    #exit()
