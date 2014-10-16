# This File contains all the class objects for simulating "The Race up Mount Improbable"
# copyright Sara Imari Walker 2011


#import pp
import sys
from time import clock
from Parameters import *
from Polymers import *
from Initialize import *
from Reactions import *
from Output import *

import numpy as np
import random
import math


### polymer structure =  3'.......5'

def main():

    for exp in range(Start_run, Start_run + Total_runs):
    
        start_time = clock()

        """Initialize random number generator for each run"""
        Parameters.run_seed = 100*exp + run_seed      
        random.seed(Parameters.run_seed)

        """Write random number generator and seed to file"""
        file_name = ('%s/%i_seed.txt' % (dirname, exp))
        file = open(file_name, 'a')
        file.write(str(random.getstate()))
        file.close()

        file_name = ('%s/%i_Seeds.txt' % (dirname, exp))
        file = open(file_name, 'a')
        s = str(exp) + '    ' +  str(Parameters.run_seed)
        file.write(s)
        file.close()


        """Clear variables for each new experimental run"""
        sequences.clear()
        initialize_monomers(sequences)

        """Initialize Random Functional motifs"""
        #for i in range(N_T):
        #    #first index is ID, second is sequence, third substrate bond
        #    Parameters.Fun_motifs.append([i, generate_random_sequence(monomer_species, Fun_len)])

        """Initialize Reaction Propensities and Reaction Clock"""
        update_propensities(sequences, catalysts)

        tau = 0
        t   = 0


        """Begin main time evolution loop"""
        while tau < tau_max:

            """Choose Reaction"""

            random.seed()
            dice_roll = random.random()*Parameters.Atot
        
            if(dice_roll < Parameters.Ap_p + Parameters.Ap_cp):
                polymerization(sequences, catalysts)

            elif(dice_roll < Parameters.Ap_p + Parameters.Ap_cp + Parameters.Ap_h):
                degradation(sequences)



            update_propensities(sequences, catalysts)


            #############################################################
            """Check Mass Conservation - If violated evolution loop will exit"""

       
            mass = 0
    
            for ID, v in sequences.items():

                mass += sequences[ID].length*sequences[ID].tot

        
            if sum(M_N) != 0 and mass != sum(M_N):

                print 'Conservation of mass violated, exiting simulation ... '
                break
        

            #############################################################

        
            """ 3. Adjust Time """

            if Parameters.Ntot == 1 and Parameters.h == 0:
                print 'All mass has converged on a single polymer ... '

                break
            

            if( t % t_frequent == 0):
                output_data(exp, t, tau, sequences)

            random.seed()
            dice_roll = random.random()
            tau -= math.log(dice_roll)/Parameters.Atot;

            print tau

            t += 1



        filename = ('%s/%i_Final_Statistics.txt' % (Parameters.dirname, exp))
        file = open(filename, 'w')
            
        for ID, seq in sequences.items():

            #print seq.sequence, seq.tot, seq.bond_struct

            s = str(seq.sequence) + str("   ") + str(seq.tot) + str("   ") + str(seq.bond_struct)

            file.write(s)
            file.write('\n')

        file.close()


        final_data(exp, mass)
            
                
        print 'Saving final data for run number ' + str(exp) + ' ...'
        runtime = clock() - start_time
        
        #output_data(exp, t, tau, sequences, monomers)

    

  
    
# End Function Main

if __name__=="__main__":
    main()
    
