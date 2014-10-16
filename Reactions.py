from Parameters import *
from Polymers import *
from decimal import *
#import matplotlib.pyplot as plt
#from pylab import *

import math, random

#######################################################################################################
def polymerization(sequences, catalysts):
    """Polymerization reaction, can link any 5' end in system with any 3' end"""

    import Parameters
    from Parameters import m, kp

    #print 'Binding polymers ...'
   
    dice_roll = random.random()*Parameters.Ap_p

    checkpoint = 0.0

    Found = False

    for ID1, seq1 in sequences.items():
        """Loop over 3' ends"""

        if Found == True:
            break

        ID1_e3 = seq1.three_prime
        
        for ID2, seq2 in sequences.items():
            """Loop over 5' ends"""

            ID2_e5 = seq2.five_prime

            if ID1 == ID2 and sequences[ID1].tot > 1:
                checkpoint += kp[ID1_e3][ID2_e5]*sequences[ID1].tot*(sequences[ID1].tot - 1)
            else:
                checkpoint += kp[ID1_e3][ID2_e5]*sequences[ID1].tot*sequences[ID2].tot

            bond = str(ID2_e5) + str(ID1_e3)
        
            if bond in catalysts.keys():
                """Add propensities for catalyzed assembly"""
                for IDc in catalysts[bond]:

                    if ID1 == ID2 and ID1 == IDc and sequences[ID1].tot > 2:
                        checkpoint += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*(sequences[IDc].tot - 1)*(sequences[IDc].tot - 2)*sequences[IDc].tot
                        
                    elif ID1 == ID2 and ID1 != IDc and sequences[ID1].tot > 1:
                        checkpoint += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID1].tot*(sequences[ID1].tot - 1)*sequences[IDc].tot

                    elif ID1 != ID2 and ID1 == IDc and sequences[IDc] > 1:
                        checkpoint += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID2].tot*(sequences[IDc].tot - 1)*sequences[IDc].tot

                    elif ID1 != ID2 and IDc == ID2 and sequences[IDc] > 1:
                        checkpoint += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID1].tot*(sequences[ID2].tot - 1)*sequences[IDc].tot

                    elif ID1 != ID2 and ID2 != IDc != ID1:
                        checkpoint += kc[ID1_e3][ID2_e5]*kp[ID1_e3][ID2_e5]*sequences[ID1].tot*sequences[ID2].tot*sequences[IDc].tot

            if checkpoint > dice_roll:
                """Matches 3' end of seq ID1 with 5' end of seq ID2"""
                e3 = ID1
                e5 = ID2

                seq1.tot -= 1
                seq2.tot -= 1
                Parameters.Ntot -= 2
                        
                Found = True
                break
                    
   
    new_seq = sequences[e5].sequence + sequences[e3].sequence
    
    found_new = True

    for new_ID, seq in sequences.items():
        """Checks to see if new polymer is in extant species pool"""

        if len(new_seq) == sequences[new_ID].length and new_seq == sequences[new_ID].sequence:

            found_new = False

            seq.tot += 1

            Parameters.Ntot += 1
            
            break


    if found_new == True:
        """If new seq is not in pool, add new sequence to list of extant species"""

        new_ID = Parameters.tot_species

        ID_e3 = sequences[e5].three_prime  #3' end of new polymer is the polymer that reacted its 5' end
        ID_e5 = sequences[e3].five_prime   #5' end of new polymer is the polymer that reacted its 3' end

        sequences[new_ID] = Polymer(new_ID, 1, len(new_seq), new_seq, [], ID_e3, ID_e5)

        #print sequences[e5].bond_struct
        #print sequences[e3].bond_struct

        sequences[new_ID].bond_struct = sequences[e5].bond_struct + [[sequences[e5].five_prime, sequences[e3].three_prime]] + sequences[e3].bond_struct

        Parameters.Ntot += 1
        Parameters.tot_species += 1

        #Check functionality
        for f in range(len(cat_motifs)):
            if cat_motifs[f][0] in new_seq:
                """If motif is in new sequence, add functional tag"""
                sequences[new_ID].F = True
                sequences[new_ID].F_motifs.append(cat_motifs[f][0])
                sequences[new_ID].B_domains.append(cat_motifs[f][1])

                if cat_motifs[f][1] in catalysts.keys():
                    catalysts[cat_motifs[f][1]].append(new_ID)
                else:
                    catalysts[cat_motifs[f][1]] = [new_ID]

               
    Parameters.kp_events += 1
 
#######################################################################################################

def degradation(sequences):
    """Hydrolyses one-bond in the sequence"""

    import Parameters
    from Parameters import m, kh, monomer_species

    #print 'Hydrolyzing bond ...'

    dice_roll = random.random()*Parameters.Ap_h

    checkpoint = 0.0

    Found = False

    for k, seq in sequences.items():

        if Found == True:
            break

        for b in range(seq.nb):

            e3 = seq.bond_struct[b][1]
            e5 = seq.bond_struct[b][0]
            
            checkpoint += kh[e3][e5]*seq.tot    #makes sure seq >= 0 for reaction to occur

            bond = str(e5) + str(e3)

            if Parameters.cat_reverse == True and bond in catalysts.keys():

                for IDc in catalysts[bond]:

                    if IDc == seq.ID and seq.tot > 1:
                        #Substrate is the catalyst
                        checkpoint += kc[e3][e5]*kh[e3][e5]*(seq.tot - 1)

                    else:
                        #Catalyst is different from substrate
                        checkpoint += kc[e3][e5]*kh[e3][e5]*seq.tot

            if dice_roll < checkpoint:

                ID = k
                Found = True

                seq.tot -= 1
                Parameters.Ntot -= 1
   
                break

    """Generate two new sequences at hydrolyzed bond"""
    if sequences[ID].nb == 0:
        bonds1 = []
        bonds2 = []
    else:
        if b == 0:
            bonds1 = []
            bonds2 = sequences[ID].bond_struct[b + 1:]
        else:
            bonds1 = sequences[ID].bond_struct[:b]
            bonds2 = sequences[ID].bond_struct[b + 1:]

    seq1 = sequences[ID].sequence[:b + 1]
    seq2 = sequences[ID].sequence[b + 1:]

  
    found1 = True
    found2 = True

    for k, seq in sequences.items():
        """Check to see if new polymers are in extant species pool"""

        if found1 == False and found2 == False:
            break

        if len(seq1) == sequences[k].length and seq1 == sequences[k].sequence:

            found1 = False

            sequences[k].tot += 1

            Parameters.repeat_seq += 1

            Parameters.Ntot += 1

            new_ID1 = k


        if len(seq2) == sequences[k].length and seq2 == sequences[k].sequence:

            found2 = False

            sequences[k].tot += 1

            Parameters.repeat_seq += 1

            Parameters.Ntot += 1

            new_ID2 = k


    if found1 == True:
        """If new seq 1 not in pool, add to list of extant polymers"""

        new_ID1 = Parameters.tot_species

        for ms in range(m):

            if seq1[:1] == monomer_species[ms]:
                ID_e3 = ms

            if seq1[len(seq1) - 1:] == monomer_species[ms]:
                ID_e5 = ms

        sequences[new_ID1] = Polymer(new_ID1, 1, len(seq1), seq1, [], ID_e3, ID_e5)

        sequences[new_ID1].bond_struct = bonds1

        #Check functionality
        if sequences[ID].F == True:
            """Sequence fragments can only be functional if parent was functional"""
            
            for f in range(len(cat_motifs)):
                if cat_motifs[f][0] in  sequences[new_ID1].sequence:
                    """If motif is in new sequence, add functional tag"""
                    sequences[new_ID1].F = True
                    sequences[new_ID1].F_motifs.append(cat_motifs[f][0])
                    sequences[new_ID1].B_domains.append(cat_motifs[f][1])

                    if cat_motifs[f][1] in catalysts.keys():
                        catalysts[cat_motifs[f][1]].append(new_ID1)
                    else:
                        catalysts[cat_motifs[f][1]] = [new_ID1]



        Parameters.Ntot += 1
        Parameters.tot_species += 1
        

    if found2 == True:
        """If new seq 1 not in pool, add to list of extant polymers"""

        new_ID2 = Parameters.tot_species

        for ms in range(m):

            if seq2[:1] == monomer_species[ms]:
                ID_e3 = ms

            if seq2[len(seq2) - 1:] == monomer_species[ms]:
                ID_e5 = ms

        sequences[new_ID2] = Polymer(new_ID2, 1, len(seq2), seq2, [], ID_e3, ID_e5)

        sequences[new_ID2].bond_struct = bonds2

        #Check functionality
        if sequences[ID].F == True:
            """Sequence fragments can only be functional if parent was functional"""
            
            for f in range(len(cat_motifs)):
                if cat_motifs[f][0] in sequences[new_ID2].sequence:
                    """If motif is in new sequence, add functional tag"""
                    sequences[new_ID2].F = True
                    sequences[new_ID2].F_motifs.append(cat_motifs[f][0])
                    sequences[new_ID2].B_domains.append(cat_motifs[f][1])

                    if cat_motifs[f][1] in catalysts.keys():
                        catalysts[cat_motifs[f][1]].append(new_ID2)
                    else:
                        catalysts[cat_motifs[f][1]] = [new_ID2]


        Parameters.Ntot += 1
        Parameters.tot_species += 1

    Parameters.kh_events += 1
    
