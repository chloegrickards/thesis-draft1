#genesystem

import numpy as np
import numpy.matlib

#Information needed: number of alleles, intrinsic growth rate,
# fitness cost for each allele, disease resistance for each allele

#testing everything in 3 alleles, identical to Richard's system

def geneSystemInputs(r):

    #number of alleles, include gene drive
    num_alleles = 9

    #number of genotypes, calculated from the number of alleles
    # uses that equation from that link
    num_genotypes = int(num_alleles * (num_alleles + 1) / 2)

    #Initial population, filled in. Here, every genotype has 10 snails, except for the seed of 5 gene drive snails
    init_genotype_pop = np.zeros(num_genotypes)
    for i in range(0, num_genotypes - 1):
        init_genotype_pop[i] = 10
    init_genotype_pop[num_genotypes - 1] = 5  # GG

    #A genotype dictionary that matches alleles to the corresponding genotype
    # Important for building the transition matrix
    genotype_dict = np.matlib.zeros((num_alleles, num_alleles))
    counter1 = 0
    counter2 = 0
    for i in range(0, num_alleles):
        for j in range(counter2, num_alleles):

            genotype_dict[i, j] = counter1
            genotype_dict[j, i] = counter1

            if i < num_alleles-1 and j < num_alleles-1:
                init_genotype_pop[counter1] = 10
            else:
                init_genotype_pop[counter1] = 0.000000001

            counter1 = counter1 + 1
        counter2 = counter2 + 1

    init_genotype_pop[num_genotypes - 1] = 10  # GG



    # fitness costs for each allele
    allele0_fit = 0 #000
    # insert fitness costs for other, non-gene drive alleles, here
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

    selfing = True

    systemInputs = [num_alleles, genotype_dict, init_genotype_pop, allele_fitness, selfing]
    return systemInputs

