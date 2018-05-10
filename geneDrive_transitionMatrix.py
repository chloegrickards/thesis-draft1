
import numpy as np
import numpy.matlib
import math

def myTransitionMatrix(inputs, previous_population, r, mu, number_snail_infections, g, sig, inbr, k):

    ## TRANSITION MATRIX SETUP ##

    #Unpacking and defining some important constants
    num_alleles = inputs[0]
    genotype_dict = inputs[1]
    allele_fitness = inputs[3]
    selfing = inputs[4]
    num_genotypes = int(num_alleles * (num_alleles + 1) / 2)

    #Vector for growth rate of genotypes
    #Vecotr for number of reproducing snails, per genotype
    #Genotype counter is used to track current genotype from allele number
    #Counter1 is used to go through each allele combination
    genotype_growth = np.zeros(num_genotypes)
    reproducing_genotypes = np.zeros(num_genotypes)
    genotype_infections = np.zeros(num_genotypes)
    genotype_counter = 0
    counter1 = 0
    for i in range(0, num_alleles):
        for j in range(counter1, num_alleles):
            #Growth rate for each genotype, fitness costs applied
            genotype_growth[genotype_counter] = r - allele_fitness[i] - allele_fitness[j]

            #Reproducing snails per genotype, death rate and number of infections applied
            #Reproducing genotypes = snails that are alive and uninfected
            reproducing_genotypes[genotype_counter] = previous_population[genotype_counter] - mu*previous_population[genotype_counter] - (1-mu)*number_snail_infections[genotype_counter]

            # Infected snails = alive and infected
            genotype_infections[genotype_counter] = (1 - mu) * number_snail_infections[genotype_counter]

            # Accounts for condition of extreme infection and removes double-counting miracidia
            if genotype_infections[genotype_counter] > reproducing_genotypes[genotype_counter]:
                reproducing_genotypes[genotype_counter] = 0
                genotype_infections[genotype_counter] = previous_population[genotype_counter] - mu*previous_population[genotype_counter]

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
                    if i != (num_alleles - 1) and m == (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # Allele i and allele n:
                    tm_col = int(genotype_dict[i, n])

                    if i != (num_alleles - 1) and m == (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # Allele j and allele m:
                    tm_col = int(genotype_dict[j, m])

                    if i != (num_alleles - 1) and m == (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    else:
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 / 4) * rate_x1x2 * p2

                    # Allele j and allele n:
                    tm_col = int(genotype_dict[j, n])

                    if i != (num_alleles - 1) and m == (num_alleles - 1):
                        tm[tm_row, num_genotypes - 1] = tm[tm_row, num_genotypes - 1] + g * (1 / 4) * rate_x1x2 * p2
                        tm[tm_row, tm_col] = tm[tm_row, tm_col] + (1 - g) * (1 / 4) * rate_x1x2 * p2
                    elif i == (num_alleles - 1) and m != (num_alleles - 1):
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

    # This runs through all of the genotype selfing combinations (i.e. all homozygous crossings, again)
    # Genotype ij crosses with genotype ij
    # counter4 goes through all ij combinations
    # If statement for selfing appears, because selfing can be toggled in the gene system, as needed by the user
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
                if i != (num_alleles - 1) and n == (num_alleles - 1):
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 - g) * (1 / 4) * rate_g_self
                    stm[stm_row, num_genotypes - 1] = stm[stm_row, num_genotypes - 1] + g * (1 / 4) * rate_g_self
                else:
                    stm[stm_row, stm_col] = stm[stm_row, stm_col] + (1 / 4) * rate_g_self

                # Allele j and allele m:
                stm_col = int(genotype_dict[j, m])

                # Gene Drive Condition applies
                # Because of the counter system, only j could be the gene drive allele
                if j == (num_alleles - 1) and m != (num_alleles - 1):
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

    # Total transition matrix is the sum of the outcrossing and selfing transition matrices, with the selfing rate applied
    tot_tm = (1-sig)*tm + sig*inbr*stm

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
                ptm[i, j] = reproducing_genotypes[i]/previous_population[i] + (tot_tm[i, j]*p_genotype[i])/(previous_population[i])*log_growth
            else:
                ptm[i, j] = (tot_tm[i, j]*p_genotype[i])/(previous_population[i])*log_growth

    # Find the current population distribution
    curr_population = previous_population * ptm

    # Return: (1) population after reproduction (2) adjusted infections for each genotype
    output = curr_population, genotype_infections

    return output
