#Python Modules
from math import sqrt, pi, pow
import os
import numpy as np

#####################################################################################################
run_seed = 350
Start_run = 1
Total_runs = 2            #Number of experimental runs - used for cluster computing. For testing purposes set Total_runs = 1

"""System Parameters & Output Specifications"""
t_0 = 0                     #initial time
tau_max = 100             #max sim time

t_frequent = 1
t_infrequent = 10                #total number of Gillespie steps between infrequent outputs

output_all = 1              #If = 0 no polymer concentrations output, if = 1 polymer concentrations output for all polymers at each time step


"""User defined Monomer Species & Statistics"""
monomer_species = ['0', '1']

M_N = [500, 500]                 #Initial number of each monomeric species at start of simulation!

m = len(monomer_species)


"""Microscopic Reaction Rates"""

p = 0.0001                 #assembly/polymerization rate 
h = 0.1                    #degradation rate


"""Matrix of differential chemical rates for monomer attachment if chemistry is NOT the same for all monomer species"""

kp = p*np.ones((m,m))

#kp = [[p]*m for i in range(m)]          #Polymerization rates, kp[0][1] is rate of attachment of 3' of monomer_species[0] to 5' of monomer_species[1]
#kh = [[h]*m for i in range(m)]
kh = h*np.ones((m,m))

kc = np.zeros((m,m))    #catalytic rate enhancement

kc[0][0] = 1000*kp[0][0]

cat_motifs = [['000', '00']] #list of catalytic motifs indexed by bond catalyzed s.t. [motif, bond] where bond = 5'3'

cat_reverse = False #If TRUE catalysts operate on forward and reverse reaction

print cat_motifs

#exit()

#####################################################################################################
### Do not modify parameters below this line!
#####################################################################################################

"""Global variables not to be changed (program will modify)"""

sequences = {}          #Dictionary of all sequences present in the system

catalysts = {}           #Dictionary of all functional sequences, keys are bonds which are catalyzed by functional seq in format '5'3'', values are catalyst IDs for the seq that catalyze the bond key

#length_list = [[], [0,1]]   #NOT CURRENTLY IN CODE ADD THIS for more efficient searches

Npoly = 0   		# Total number of polymers in system (summed over lattice sites)
Ntot = 0                # Total number of molecules (monomers + polymers) in the system
tot_species = 0         # Total number of species that have ever existed in sim

#NBB = [[0]*m for i in range(m)]    # Total number of each distinct type of bond in the system


Atot = 0                # Total global propensity for all events
Ap_p = 0                # Global propensity for uncatalyzed assembly     
Ap_h = 0                # Global propensity for uncatalyzed hydrolysis
                       

Ap_cp = 0               # Global propensity for catalyzed assembly
Ap_ch = 0               # Global propensity for catalyzed hydrolysis


kp_events = 0           # Tracks number of assembly events
kh_events = 0           # Tracks number of degredation events

repeat_seq = 0
null_event = 0


mass = sum(M_N)

dirname = ('data/kp%.4f_kh%.2f_N%d' % (p, h, mass))

if not os.path.exists(dirname):
    os.makedirs(dirname)
