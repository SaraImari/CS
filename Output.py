import math
import numpy as np
import matplotlib.pyplot as plt


def output_data(exp, t, tau, sequences):

    import Parameters

    from Parameters import t_infrequent


    tot_species(exp, tau)

    avg_length(exp, tau, sequences)

    #extant_species(exp, sequences, tau)

    #diversity(exp, sequences, tau)
    
    mass_fraction(exp, tau, sequences)

    #output_monomers(exp, tau, sequences, monomers)

    

    if Parameters.output_all == 1:
        output_concentrations(exp, tau, sequences)

    #if t % t_infrequent == 0 and Parameters.kr != 0:
    #    sequence_distribution(exp, tau, sequences, monomers)

    #    if Parameters.Fit != 2:
    #        """If the landscape is not flat - output distribution plot of k_values"""
    #        population_landscape(exp, tau, sequences, monomers)

###############################################################################################################################################################
###################### Calculate Ouput ########################################################################################################################
###############################################################################################################################################################

def tot_species(exp, t):

    import Parameters
    
    file_name = ('%s/%i_tot_species.dat' % (Parameters.dirname, exp))
	
    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
            
    s = str(t) + '      ' + str(Parameters.tot_species)
    file.write(s)
    file.write('\n')
    file.close()

#####################################################################################
def avg_length(exp, t, sequences):
    """Calculates and Prints Length distribution to file"""

    import Parameters

    avg_length = 0.0

    tot = 0
    
    for ID, seq in sequences.items():

        if seq.tot != 0:

            avg_length += float(seq.tot)*float(seq.length)

            tot += seq.tot


    avg_length /= tot

    #print 'avg length ' + str(avg_length)
    
    file_name = ('%s/%i_avg_length.dat' % (Parameters.dirname, exp))
	
    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
            
    s = str(t) + '      ' + str(avg_length)
    file.write(s)
    file.write('\n')
    file.close()


################################################################################
def extant_species(exp, sequences, t):

    import Parameters
    
    num_extant = 0

    for ID, v in sequences.items():

        if sequences[ID].tot_count != 0:

            num_extant += 1

            
    """Write Extant Species Count to File"""

    file_name = ('%s/%i_extant_species.dat' % (Parameters.dirname, exp))
	
    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
            
    s = str(t) + '      ' + str(num_extant)
    file.write(s)
    file.write('\n')
    file.close()

################################################################### 
   
def diversity(exp, sequences, t):
    """Calculate Shannon Diversity of Extant Population"""

    import Parameters

    I = 0

    for ID, v in sequences.items():

        if sequences[ID].tot_count != 0:

            p_i = float(sequences[ID].tot_count)/Parameters.Npoly

            I -= p_i*math.log(p_i, 2)

    file_name = ('%s/%i_extant_diversity.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
        
    s = str(t) + "        " + str(I)
    file.write(s)
    file.write('\n')

    file.close()


################################################################################
def mass_fraction(exp, t, sequences):
    """Calculates and prints average fitness to file"""

    import Parameters

    Nmono = 0
    Npoly = 0


    for ID, seq in sequences.items():

        if seq.length == 1:

            Nmono += seq.tot

        else:

            Npoly += seq.tot*seq.length

    filename = ('%s/%i_mono_population.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')

    s = str(t) + '       ' + str(Nmono)
    
    file.write(s)
    file.write('\n')
    file.close()

    filename = ('%s/%i_poly_population.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')

    
    s = str(t) + '       ' + str(Npoly)
    
    file.write(s)
    file.write('\n')
    file.close()


###################################################################
def output_concentrations(exp, t, sequences):
    """Print time-dependent polymer concentrations"""

    import Parameters

    for ID, v in sequences.items():

        if sequences[ID].length <= 3:

            filename = ('%s/%i_sequence_%s.dat' % (Parameters.dirname, exp, sequences[ID].sequence))
	
            if(t == 0):
                file = open(filename, 'w')
            else:
                file = open(filename, 'a')
            s = str(t) + '      ' + str(sequences[ID].tot)
            file.write(s)
            file.write('\n')
            file.close()

###################################################################


def output_monomers(exp, t, sequences, monomers):
    """Print time-dependent polymer concentrations"""

    import Parameters

    for (k, v) in monomers.items():

        filename = ('%s/%i_monomer_%s.dat' % (Parameters.dirname, exp, k))
	
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
            
        s = str(t) + '      ' + str(monomers[k].tot_count)
        file.write(s)
        file.write('\n')
        file.close()

###################################################################

def hamming_populations(exp, t, sequences):
    """Outputs total population size of polymers with given nA/nB ratio"""

    import Parameters

    from Parameters import master_sequence

    N_bins = Parameters.R_L + 1

    bin_abundance = [0]*(N_bins + 1)

    avg_hamming = 0.0
    sig_hamming = 0.0

    for ID, k in sequences.items():

        hamming_distance = 0

        for i in range(Parameters.R_L):

            if sequences[ID].sequence[i] != master_sequence[i]:

                hamming_distance += 1
                

        avg_hamming += hamming_distance*sequences[ID].tot_count
        
        bin_abundance[hamming_distance] += sequences[ID].tot_count



    avg_hamming /= Parameters.Npoly


    for ID, k in sequences.items():

        hamming_distance = 0

        for i in range(Parameters.R_L):

            if sequences[ID].sequence[i] != master_sequence[i]:

                hamming_distance += 1
    
        sig_hamming   += sequences[ID].tot_count*pow(hamming_distance - avg_hamming, 2) 



    sig_hamming = math.sqrt(sig_hamming/(Parameters.Npoly - 1))
    
    
    for j in range(Parameters.R_L):

        filename = ('%s/%i_hamming_population_%i.dat' % (Parameters.dirname, exp, j))

        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')

    
        s = str(t) + '       ' + str(bin_abundance[j])
    
        file.write(s)
        file.write('\n')
        file.close()


    filename = ('%s/%i_avg_hamming.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')

    
    s = str(t) + '       ' + str(avg_hamming)
    
    file.write(s)
    file.write('\n')
    file.close()


    filename = ('%s/%i_sig_hamming.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(filename, 'w')
    else:
        file = open(filename, 'a')

    
    s = str(t) + '       ' + str(sig_hamming)
    
    file.write(s)
    file.write('\n')
    file.close()


##############################################################
# Infrequent Calculations
#############################################################

def population_landscape(exp, t, sequences, monomers):
    """Prints to file the population profile at time t to file"""

    import Parameters

    from Parameters import N_bins

    bin_res = float(Parameters.kr/N_bins)
    
    bin_k = [0]*(N_bins + 1)
    bin_abundance = [0]*(N_bins + 1)


    for i in range(len(bin_k)):

        bin_k[i] = i*bin_res
    
    for ID, v in sequences.items():

        bin_num = int(math.floor((N_bins - 1)*(sequences[ID].k/Parameters.kr)))

        bin_abundance[bin_num] += sequences[ID].tot_count

    fig = plt.figure()   
    ax = fig.add_subplot(111)
    ax.set_ylabel('Abundance')
    ax.set_xlabel('k_value')
    ax.set_xticks(bin_k)
    width = 0.1
   
    rects1 = ax.bar(bin_k, bin_abundance , width, color='cyan', align = 'center')
    
    #plt.setp(plt.gca(), 'yticklabels', [])
    #plt.setp(plt.gca(), 'xticklabels', [])

    file_name = ('%s/%i_population_landscape_%2.f.png' % (Parameters.dirname, exp, t))

    plt.savefig(file_name, format = 'png')

    
################################################################################

def sequence_distribution(exp, t, sequences, monomers):
    """Prints to file the population profile at time t to file"""

    import Parameters

    N_bins = Parameters.R_L + 1
    
    bin_nA = [0]*(N_bins + 1)
    bin_abundance = [0]*(N_bins + 1)

    for i in range(Parameters.R_L + 1):
        """Initialize array to populations"""

        ##############################################
        ### Use binomial coefficient to determine how many sequences have given A/B ratio with: n_A = i, n_B = R_L - i
        ##############################################

        bin_nA[i] = i #number of A monomers in sequence


    for ID, k in sequences.items():

        i = sequences[ID].sequence.count('A')
        bin_abundance[i] += sequences[ID].tot_count
  
 
    fig = plt.figure()   
    ax = fig.add_subplot(111)
    ax.set_ylabel('Abundance')
    ax.set_xlabel('Number of A-monomers in sequence')
    ax.set_xticks(bin_nA)
    width = 0.25
   
    rects1 = ax.bar(bin_nA, bin_abundance , width, color='cyan', align = 'center')
    
    #plt.setp(plt.gca(), 'yticklabels', [])
    #plt.setp(plt.gca(), 'xticklabels', [])

    file_name = ('%s/%i_sequence_distribution_%2.f.png' % (Parameters.dirname, exp, t))

    plt.savefig(file_name, format = 'png')


################################################################################


def final_data(exp, mass):
    """Output final run statistics"""
    import Parameters

    """Write Polymer list to file"""
    filename = ('%s/%i_run_statistics.txt' % (Parameters.dirname, exp))

    
    file = open(filename, 'w')
    """
    file.write("Run parameters \n")
    s = 'Total mass = ' + str(mass) + '\n'
    file.write(s)
    s = 'kr = ' + str(Parameters.kr) + '\n'
    file.write(s)
    s = 'kh = ' + str(Parameters.kh) + '\n'
    file.write(s)
    s = 'km = ' + str(Parameters.km) + '\n'
    file.write(s)
    s = 'kc = ' + str(Parameters.kc) + '\n'
    file.write(s)
    s = 'ks = ' + str(Parameters.ks) + '\n'
    file.write(s)
    s = 'mu = ' + str(Parameters.mu) + '\n'
    file.write(s)
    s = 'Run seed = ' + str(Parameters.run_seed) + '\n'
    file.write(s)
    s = 'Number of Steps = ' + str(Parameters.tot_step) + '\n' + '\n'
    
    """
    file.write("Information about kMC run including run statistics: \n \n")

    file.write("Number of polymerization events = ")
    s = str(Parameters.kp_events)
    file.write(s)
    file.write('\n')

    file.write("Number of degradation events =  ")
    s = str(Parameters.kh_events)
    file.write(s)
    file.write('\n')

    file.close()
