# This File contains all the class objects for simulating "The Race up Mount Improbable"
# copyright Sara Imari Walker 2011

import numpy as np
import Parameters

###################################################################################################

class Polymer(object):
    """A Polymeric Sequence"""
    
    def __init__(self, ID, tot, length, sequence, nbb, three_prime, five_prime):
        print("A new polymer has been born!")
        self.ID = ID
       	self.tot = tot
	self.length = length
        self.sequence = sequence
        self.nb = length - 1
        self.nbb = []              #identifies how many of each bond variety there are in polymer i.e. 
        self.three_prime = three_prime    #ID number identifying the leftmost residue of polymer string
        self.five_prime = five_prime      #ID number identifying the rightmost residue of polymer string
        
        self.bond_struct = []      #identifies the two monomers paired in a bond, i.e. [5', 3'] = [0, 1] 

        
        self.F = False           #if TRUE, sequence contains functional moeity
        self.F_motifs = []       #list of functional motifs in sequence
        self.B_domains = [ ]     #list of binding domains for each functional motif

        
 
    def talk(self):
        print("Hi. I'm a polymer.", "\n")



###################################################################################################
        
class Monomer(object):
    """A Monomer Species"""
    
    def __init__(self, tot_count, species):
        #print("A new polymer has been born!")
       	self.tot = tot_count
	self.species = species
