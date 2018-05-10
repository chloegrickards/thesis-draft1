import numpy as np
import numpy.matlib

#Input: intrinsic growth rate r

def geneSystemInputs(r):

    #number of alleles, include gene drive allele
    #This system: 000, 001, 010, 100, 011, 101, 110, 111, and 111g (gene drive)
    num_alleles = 9

    #number of genotypes, calculated from the number of alleles
    # Source = http://scienceprimer.com/hardy-weinberg-equilibrium-calculator
    num_genotypes = int(num_alleles * (num_alleles + 1) / 2)

    #Initial population vector
    init_genotype_pop = np.zeros(num_genotypes)

    # A genotype dictionary that matches alleles to the corresponding genotype
    # Important for building the transition matrix and other functions
    genotype_dict = np.matlib.zeros((num_alleles, num_alleles))

    #Iterating through allele combinations
    counter1 = 0
    counter2 = 0
    for i in range(0, num_alleles):
        for j in range(counter2, num_alleles):

            #Genotype dictionary inputs
            #End result: able to find genotype number from two allele numbers, and vice versa
            genotype_dict[i, j] = counter1
            genotype_dict[j, i] = counter1

            #Initial population inputs
            # Filled in with 10 snails per wt genotype, and a seed of 10 gene drive snails
            if i < num_alleles-1 and j < num_alleles-1:
                init_genotype_pop[counter1] = 10
            else:
                init_genotype_pop[counter1] = 0.000000001

            counter1 = counter1 + 1
        counter2 = counter2 + 1

    #Gene drive seed
    init_genotype_pop[num_genotypes - 1] = 10  # GG

    # Fitness cost for each allele
    # Here, assumes that each receptor comes with a fitness cost of 0.1*r
    # The gene drive comes at an additional cost of 0.1*r
    # For example, the total reproducing power of a 011 allele = r - 0.1*r - 0.1*r
    allele0_fit = 0 #000
    allele1_fit = 0.1 * r #001
    allele2_fit = 0.1 * r #010
    allele3_fit = 0.1 * r #100
    allele4_fit = 0.2 * r #011
    allele5_fit = 0.2 * r #101
    allele6_fit = 0.2 * r #110
    allele7_fit = 0.3 * r #111
    allele_g_fit = 0.4 * r #111 + gene drive
    allele_fitness = [allele0_fit, allele1_fit, allele2_fit, allele3_fit, allele4_fit,
                      allele5_fit, allele6_fit, allele7_fit, allele_g_fit]

    #Does selfing occur?
    selfing = True

    systemInputs = [num_alleles, genotype_dict, init_genotype_pop, allele_fitness, selfing]
    return systemInputs

