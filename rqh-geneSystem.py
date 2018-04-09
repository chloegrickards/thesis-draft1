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


#transitionmatrix

#info from gene system
#at the end of my transition matrix (or maybe after I calculate transition matrix, and in my actual code,
# remember to convert genotypes --> phenotypes by using codominant alleles
# contains some slight errors, may require debug!
#eventually tailor this to individual use

import numpy as np
import numpy.matlib
import math
from rqh_geneSystem import geneSystemInputs

#arguments could also contain gamma
def myTransitionMatrix(inputs, previous_population, r, mu, beta, number_snail_infections, g, sig, inbr, k):

    num_alleles = inputs[0]
    genotype_dict = inputs[1]
    #initial_genotypes = inputs[2]
    allele_fitness = inputs[3]
    selfing = inputs[4]


    num_genotypes = int(num_alleles * (num_alleles + 1) / 2)

    reproducing_genotypes = np.zeros(num_genotypes)
    genotype_growth = np.zeros(num_genotypes)
    genotype_counter = 0
    counter1 = 0
    for i in range(0, num_alleles):
        for j in range(counter1, num_alleles):

            genotype_growth[genotype_counter] = r - allele_fitness[i] - allele_fitness[j]

            reproducing_genotypes[genotype_counter] = previous_population[genotype_counter] - mu*previous_population[genotype_counter] - (1-mu)*number_snail_infections[genotype_counter]*previous_population[genotype_counter]

            #will have to bug Richard about this later
            # but I think this is important for population calculations down the line
            if reproducing_genotypes[genotype_counter] == 0:
                reproducing_genotypes[genotype_counter] = 0.00000000001

            genotype_counter = genotype_counter + 1
        counter1 = counter1 + 1

    #find total number of reproducing genotypes
    # that is, find genotype distribution for organisms that are still alive and uninfected
    tot_repr_genotypes = np.sum(reproducing_genotypes)

    #genotype proportions
    p_genotype = np.zeros(num_genotypes)
    for i in range(0, num_genotypes):
        p_genotype[i] = reproducing_genotypes[i]/tot_repr_genotypes


    # outcrossing transition matrix
    #each cross can only have four outcomes, so just make everything based on those outcomes
    counter3 = 0
    counter4 = 0
    tm = np.matlib.zeros((num_genotypes, num_genotypes))
    for i in range(0, num_alleles):
        for j in range(counter3, num_alleles):
            for m in range(0, num_alleles):
                for n in range(counter4, num_alleles):

                    cross1 = int(genotype_dict[i, j])  # genotype 1 in the crossing event
                    cross2 = int(genotype_dict[m, n])  # genotype 2 in the crossing event
                    tm_row = cross1

                    # do 4 calculations for each allele combination
                    rate_x1x2 = (genotype_growth[cross1] + genotype_growth[cross2]) / 2

                    #proportion of the cross2 genotype
                    p2 = p_genotype[cross2]

                    # do 4 calculations for each allele combination
                    # im cross
                    tm_col = int(genotype_dict[i, m])

                    # gene drive condition
                    if i != (num_alleles - 1) and m == (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # in cross
                    tm_col = int(genotype_dict[i, n])

                    # gene drive condition
                    if i != (num_alleles - 1) and n == (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and n != (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # jm cross
                    tm_col = int(genotype_dict[j, m])

                    # gene drive condition
                    if j != (num_alleles - 1) and m == (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    elif j == (num_alleles - 1) and m != (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # jn cross
                    tm_col = int(genotype_dict[j, n])

                    # gene drive condition
                    if j != (num_alleles - 1) and n == (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    elif j == (num_alleles - 1) and n != (num_alleles - 1):
                        # the same column gets the leftovers of the gene drive
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                        # the homozygous gene drive column (last column) gets the most of the gene drive
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                        # g2_counter = g2_counter + 1
                counter4 = counter4 + 1
            # g1_counter = g1_counter + 1
            counter4 = 0
        counter3 = counter3 + 1



    #selfing transition matrix
    stm = np.matlib.zeros((num_genotypes, num_genotypes))
    counter5 = 0
    if selfing == True:
        for i in range(0, num_alleles):
            for j in range(counter5, num_alleles):
                m = i
                n = j

                g_self = int(genotype_dict[i,j])
                stm_row = g_self

                # do 4 calculations for each allele combination
                rate_g_self = genotype_growth[g_self]

                # do 4 calculations for each allele combination
                # im cross
                stm_col = int(genotype_dict[i, m])

                # no gene drive condition necessary -- all homozygous
                stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # in cross
                stm_col = int(genotype_dict[i, n])

                # gene drive condition
                if i != (num_alleles - 1) and n == (num_alleles - 1):
                    # the same column gets the leftovers of the gene drive
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 - g) * (1 / 4) * rate_g_self
                    # the homozygous gene drive column (last column) gets the most of the gene drive
                    stm[stm_row, num_genotypes - 1] = stm[stm_row, num_genotypes - 1] + g * (1 / 4) * rate_g_self
                # you never get the condition for i == (num_alleles - 1) and n != (num_alleles - 1):
                else:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # jm cross
                stm_col = int(genotype_dict[j, m])

                # gene drive condition
                # you never meet the condition for if j != (num_alleles - 1) and m == (num_alleles - 1):
                if j == (num_alleles - 1) and m != (num_alleles - 1):
                    # the same column gets the leftovers of the gene drive
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 - g) * (1 / 4) * rate_g_self
                    # the homozygous gene drive column (last column) gets the most of the gene drive
                    stm[stm_row, num_genotypes - 1] = stm[stm_row, num_genotypes - 1] + g * (1 / 4) * rate_g_self
                else:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # jn cross
                stm_col = int(genotype_dict[j, n])

                # no gene drive condition necessary -- all homozygous
                stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self


            counter5 = counter5 + 1



    tot_tm = (1-sig)*tm + sig*inbr*stm

    # 1-(sAA[i-1]*(mu + (1-mu)*rho))/sAA[i-1] + sAA_AA_r/(sAA[i-1]*s_r)*(s_death*k/(s_death + (k - s_death)*math.exp(-s_r)) - s_death)

    #error is below somewhere

    c_genotype = p_genotype * tot_tm
    genotype_tot_c = np.sum(c_genotype)

    #log_growth = (1 / genotype_tot_c) * (tot_repr_genotypes * k / (tot_repr_genotypes + (k - tot_repr_genotypes) * math.exp(-genotype_tot_c)) - tot_repr_genotypes)


    log_growth = (1 / genotype_tot_c) * (tot_repr_genotypes * k / (tot_repr_genotypes + (k - tot_repr_genotypes) * math.exp(-genotype_tot_c)) - tot_repr_genotypes)
    # print(log_growth)
    # print(tot_repr_genotypes)
    # print(math.exp(-genotype_tot_c))
    #print(genotype_tot_c)

    #idk how to FIX THIS
    #ask RICHARD AHAHHH
    #cue lemongrab scream
    if log_growth < 0 :
        log_growth = 1

    ptm = np.matlib.zeros((num_genotypes, num_genotypes))

    for i in range(0, num_genotypes):
        for j in range(0, num_genotypes):


            if i == j:
                #ptm[i, j] =  reproducing_genotypes[i]/previous_population[i] +  (tot_tm[i, j]*p_genotype[i])/(previous_population[i]*genotype_tot_c)*(tot_repr_genotypes*k/(tot_repr_genotypes + (k - tot_repr_genotypes)*math.exp(-genotype_tot_c)) - tot_repr_genotypes)
                #ptm = reproducing_genotypes[i]/previous_population[i] + p_genotype[i]*tot_tm[i, j]
                ptm[i, j] = reproducing_genotypes[i]/previous_population[i] + (tot_tm[i, j]*p_genotype[i])/(previous_population[i])*log_growth
            else:
                #ptm[i, j] = sAA_AA_r/(sAA[i-1]*s_r)*(s_death*k/(s_death + (k - s_death)*math.exp(-s_r)) - s_death)
                #ptm[i, j] = (tot_tm[i, j]*p_genotype[i])/(previous_population[i]*genotype_tot_c)*(tot_repr_genotypes*k/(tot_repr_genotypes + (k - tot_repr_genotypes)*math.exp(-genotype_tot_c)) - tot_repr_genotypes)
                ptm[i, j] = (tot_tm[i, j]*p_genotype[i])/(previous_population[i])*log_growth


            #print((tot_tm[i, j]*p_genotype[i])/(previous_population[i]))
            #print(reproducing_genotypes[i]/previous_population[i])
            #print((tm[i,j] * p_genotype[i] * log_growth)/previous_population[i])


    curr_population = previous_population * ptm
    #print(curr_population)


    output = [tm, ptm, curr_population]

    return output
