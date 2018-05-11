import numpy as np
import numpy.matlib
import math

# This document contains the supporting functions for the DIY_transitionMatrix

# This function creates the genotype dictionary, which makes a matrix of genotype numbers in the row/column that
# matches the allele numbers
def findGenotypeDictionary(n):

    # A genotype dictionary that matches alleles to the corresponding genotype
    # Important for building the transition matrix and other functions
    genotype_dict = np.zeros((n, n))

    # Iterating through allele combinations
    counter1 = 0
    counter2 = 0
    for i in range(0, n):
        for j in range(counter2, n):

            # Genotype dictionary inputs
            # End result: able to find genotype number from two allele numbers, and vice versa
            genotype_dict[i, j] = counter1
            genotype_dict[j, i] = counter1

            counter1 = counter1 + 1
        counter2 = counter2 + 1

    return genotype_dict

# Finds the initial population, according to Option 2 of the DIY_transitionMatrix program
# Makes uniform wild-type alleles, sets up a gene drive seed, and sets heterozygous gene drive genotypes = 0
# Simulates an initial population of wild-type alleles + introduction of gene drive
def findInitialPopulation(n, wt, gd):

    # Define the number of genotypes
    num_genotypes = int(n * (n + 1) / 2)

    # Vector for initial genotypes
    init_genotype_pop = np.zeros(num_genotypes)

    # Iterating through allele combinations
    counter1 = 0
    counter2 = 0
    for i in range(0, n):
        for j in range(counter2, n):
            if i < n - 1 and j < n - 1:
                # Fill in wild type
                init_genotype_pop[counter1] = wt
            else:
                # Fill in ~0 for heterozygous gene drive genotypes
                init_genotype_pop[counter1] = 0.000000001
            counter1 = counter1 + 1
        counter2 = counter2 + 1

    # Gene drive seed
    init_genotype_pop[num_genotypes - 1] = gd  # GG

    return init_genotype_pop

# Creates the transition matrix function, to be called each generation
# Sets up the transition matrix, creates the outcrossing transition matrix, creates the selfing transition matrix,
# creates the overall transition matrix, and returns the new population distribution
def findTransitionMatrix(num_alleles, r, allele_fitness, previous_population, mu, beta, mu_i, g, gd, selfing, sig, inbr, k):

    ## TRANSITION MATRIX SETUP ##

    #Unpacking and defining some important constants
    genotype_dict = findGenotypeDictionary(num_alleles)
    num_genotypes = int(num_alleles * (num_alleles + 1) / 2)

    #Vector for growth rate of genotypes
    #Vecotr for number of reproducing snails, per genotype
    #Genotype counter is used to track current genotype from allele number
    #Counter1 is used to go through each allele combination
    genotype_growth = np.zeros(num_genotypes)
    reproducing_genotypes = np.zeros(num_genotypes)
    genotype_counter = 0
    counter1 = 0
    for i in range(0, num_alleles):
        for j in range(counter1, num_alleles):
            #Growth rate for each genotype, fitness costs applied
            genotype_growth[genotype_counter] = r - allele_fitness[i] - allele_fitness[j]

            #Reproducing snails per genotype, death rate and number of infections applied
            #Reproducing genotypes = snails that are alive and uninfected
            reproducing_genotypes[genotype_counter] = previous_population[genotype_counter] - mu*previous_population[genotype_counter] - (1-mu)*beta*mu_i*previous_population[genotype_counter]

            genotype_counter = genotype_counter + 1
        counter1 = counter1 + 1
    tot_repr_genotypes = np.sum(reproducing_genotypes)

    #Setting up proportions of genotypes within the total population
    p_genotype = np.zeros(num_genotypes)
    for i in range(0, num_genotypes):
        p_genotype[i] = reproducing_genotypes[i]/tot_repr_genotypes

    ## OUTCROSSING TRANSITION MATRIX ##

    # Transition matrix to contain all probabilities of genotype state transitions, from outcrossing events
    tm = np.matlib.zeros((num_genotypes, num_genotypes))

    # This runs through all of the genotype crossing combinations and subsequent allele combinations
    # Genotype ij crosses with genotype mn
    # counter2 goes through all ij combinations, while counter3 goes through all mn combinations
    counter2 = 0
    counter3 = 0
    for i in range(0, num_alleles):
        for j in range(counter2, num_alleles):
            for m in range(0, num_alleles):
                for n in range(counter3, num_alleles):
                    # Find genotypes 1 and 2 in this particular outcrossing event
                    cross1 = int(genotype_dict[i, j])
                    cross2 = int(genotype_dict[m, n])

                    # Row of transition matrix = original genotype = cross 1
                    tm_row = cross1

                    # Offspring from this crossing = average growth rate of the original genotypes
                    rate_x1x2 = (genotype_growth[cross1] + genotype_growth[cross2]) / 2

                    # Proportion of the cross2 genotype
                    p2 = p_genotype[cross2]

                    # Do a calculation for each allele outcome of this genotype cross.
                    # Allele i and allele m:
                    tm_col = int(genotype_dict[i, m])

                    # Gene Drive Condition
                    # If one allele is the gene drive, than the efficiency rate g determines how many genotypes
                    # go to a homozygous gene drive genotype
                    # - Most form a gene drive genotype at a rate g
                    # - Some stay in the heterozygous gene drive state at a rate (1-g)
                    #
                    # Case 1 -- m is the gene drive allele
                    # Case 2 - i is the gene drive allele
                    # Case 3 - the allele is homozygous gene drive, or does not have a gene drive.
                    #
                    # A transition matrix cell comes from the existing value, plus 1/4 of the cross2 population, with
                    # growth rate and gene drive applied.  The 1/4 comes from this particular cross being 1/4 of the outcomes
                    if i != (num_alleles - 1) and m == (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # Allele i and allele n:
                    tm_col = int(genotype_dict[i, n])

                    if i != (num_alleles - 1) and m == (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # Allele j and allele m:
                    tm_col = int(genotype_dict[j, m])

                    if i != (num_alleles - 1) and m == (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # Allele j and allele n:
                    tm_col = int(genotype_dict[j, n])

                    if i != (num_alleles - 1) and m == (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1) and gd == True:
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2
                counter3 = counter3 + 1
            counter3 = 0
        counter2 = counter2 + 1


    ## SELFING TRANSITION MATRIX ##

    # Transition matrix to contain all probabilities of genotype state transitions, from self-reproduction events
    stm = np.matlib.zeros((num_genotypes, num_genotypes))

    # If selfing exists in this system, this runs through all of the genotype selfing combinations
    # (i.e. all homozygous crossings, again)
    # Genotype ij crosses with genotype ij
    # counter4 goes through all ij combinations
    counter4 = 0
    if selfing == True:
        for i in range(0, num_alleles):
            for j in range(counter4, num_alleles):
                # Sets up crossing indexes, copying the same crossing method as in the outcrossing transition matrix
                m = i
                n = j

                #Genotype index for this particular cross
                g_self = int(genotype_dict[i,j])

                # Row of transition matrix = current
                stm_row = g_self

                # No need to take an average - the average is equal to the original
                rate_g_self = genotype_growth[g_self]

                # Because of heterozygous alleles, we still need to do all four allele events, with simplifications
                # Allele i and allele m:
                stm_col = int(genotype_dict[i, m])

                # No gene drive condition necessary -- all homozygous
                stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # Allele i and allele n:
                stm_col = int(genotype_dict[i, n])

                # Gene Drive Condition applies
                # Because of the counter system, only n could be the gene drive allele
                if i != (num_alleles - 1) and n == (num_alleles - 1) and gd == True:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 - g) * (1 / 4) * rate_g_self
                    stm[stm_row, num_genotypes - 1] = stm[stm_row, num_genotypes - 1] + g * (1 / 4) * rate_g_self
                else:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # Allele j and allele m:
                stm_col = int(genotype_dict[j, m])

                # Gene Drive Condition applies
                # Because of the counter system, only j could be the gene drive allele
                if j == (num_alleles - 1) and m != (num_alleles - 1) and gd == True:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 - g) * (1 / 4) * rate_g_self
                    stm[stm_row, num_genotypes - 1] = stm[stm_row, num_genotypes - 1] + g * (1 / 4) * rate_g_self
                else:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # Allele j and allele n:
                stm_col = int(genotype_dict[j, n])

                # No gene drive condition necessary -- all homozygous
                stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self
            counter4 = counter4 + 1


    ## TOTAL TRANSITION MATRIX CALCULATIONS ##

    # If selfing exists, total transition matrix is the sum of the outcrossing and selfing transition matrices
    #       depending on selfing rate and inbreeding cost
    if selfing == True:
        tot_tm = (1 - sig) * tm + sig * inbr * stm
    else:
        tot_tm = tm

    # Current genotype proportions from the previous proportions times the probability of genotype transitions
    c_genotype = p_genotype * tot_tm

    # Total number of snails, after reproduction, infection, and mortality
    genotype_tot_c = np.sum(c_genotype)

    # Logistic growth expression, to be applied on the whole snail population
    log_growth = (1 / genotype_tot_c) * (tot_repr_genotypes * k / (tot_repr_genotypes + (k - tot_repr_genotypes) * math.exp(-genotype_tot_c)) - tot_repr_genotypes)

    # Population transition matrix
    ptm = np.matlib.zeros((num_genotypes, num_genotypes))

    # Fills in the population transition matrix
    # Homozygous genotypes require slightly different math, since this cell in the transition matrix represents
    # no change, so this accounts for the existing genotypes too.
    for i in range(0, num_genotypes):
        for j in range(0, num_genotypes):
            if i == j:
                ptm[i, j] = reproducing_genotypes[i] / previous_population[i] + (tot_tm[i, j] * p_genotype[i]) / (
                previous_population[i]) * log_growth
            else:
                ptm[i, j] = (tot_tm[i, j] * p_genotype[i]) / (previous_population[i]) * log_growth

    # Find the current population distribution
    curr_population = previous_population * ptm

    return curr_population