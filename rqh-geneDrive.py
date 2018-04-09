##ASSUMPTIONS
# Here resistance assumed to be a base rate that increases with every time step
# But that method doesn't work out so great with the transition matrix I made. Sad.
# every schisto has an equal chance of infecting a snail - shows no preference
#schisto are evenly dispersed across all genotypes, where they face two outcomes: establish infection or die



### To Do
# clean up this code lol wow this is a mess
# adjust infection matrix with contact rate c and force of infection

import numpy as np
from rqh_geneSystem import geneSystemInputs
import numpy.matlib
import scipy
import pylab as pl

from rqh_transitionMatrix import myTransitionMatrix

## CONSTANTS ##

num_gens = 50

r = 40

#natural death rate
mu = 0.2

#infection rate
beta = 0.4

#death due to infection, if applicable to your system
#not applicable to my system
mu_i = 0.5

#gene drive efficiency
g = 0.9

#resistance loss per time step
b = 0

sig = 0.25  # selfing rate

inbr = .7  # cost of inbreeding

k = 200 #carrying capacity

#c = 0.000006  #something like the probability of contact
#in Sokolow 2015 and Feng 2004
c = 0.01

#Likelihood of survival from sporocyst to adult worm
h = 0.5

#Egg production by adult worms
q = 100

#from the gene system function
inputs = geneSystemInputs(r)
num_alleles = inputs[0]
genotype_dict = inputs[1]
initial_genotypes = inputs[2]
allele_fitness = inputs[3]
selfing = inputs[4]

num_receptors = 3
num_genotypes = int(num_alleles * (num_alleles + 1) / 2)
num_strains = 8

receptor_network = np.array([[0, 0, 0],
                            [0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0],
                            [0, 1, 1],
                            [1, 0, 1],
                            [1, 1, 0],
                            [1, 1, 1]])

mucin_network = np.array([[0, 0, 0],
                            [0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0],
                            [0, 1, 1],
                            [1, 0, 1],
                            [1, 1, 0],
                            [1, 1, 1]])

schisto_distribution = np.full(8, 1/8)
miracidia_abundance = 10000

schisto_infection_matrix = np.matlib.zeros((num_gens+1, num_strains))

genotype = np.zeros((num_genotypes, num_gens+1))
for i in range(0, num_genotypes):
    genotype[i,0] = initial_genotypes[i]


for i in range(1, num_gens+1):


    if i == 1 :
        previous_population = initial_genotypes
    else:
        #this is confusing syntax, but it just moves the generations forward
        previous_population = next_generation

    #determine abundance of receptors
    # I need to keep track of the number of infections between snail strain and schisto strain
    # I need a receptor network for all of my genotypes
    # I can make that using something like below or by modifying the gamma thing
    # Or just assume all receptors are dominant
    # Or make a brand new transition matrix that has activation alleles
    # I need to make a call backed up by research probably
    # Webster -- resistance gene dominant trait and I think resistance is mostly dominant
    receptor_abundance_network = np.array(np.zeros(np.shape(receptor_network)))
    receptor_infection_network = np.array(np.zeros(np.shape(receptor_network)))
    genotype_infection_network = np.zeros((num_genotypes, num_receptors))
    counter1 = 0
    for j in range(0, num_alleles):
        for k in range(counter1, num_alleles):
            receptor_phenotype = [0, 0, 0]
            genotype_counter = int(genotype_dict[j,k])
            genotype_abundance = previous_population[genotype_counter]

            if j == 8:
                j = 7
            if k == 8:
                k = 7

            for l in range(0, num_receptors):

                if receptor_network[j,l] == receptor_network[k,l] and receptor_network[j,l] == 0:
                    receptor_phenotype[l] = 0
                else:
                    #expressing dominance!
                    receptor_phenotype[l] = 1
            #match new receptor to genotype
            #store receptor phenotype in a matrix
            #lol note to self but pycharm doesn't like elifs
            if receptor_phenotype == [0, 0, 0]:
                rrow = 0
            if receptor_phenotype == [0, 0, 1]:
                rrow = 1
            if receptor_phenotype == [0, 1, 0]:
                rrow = 2
            if receptor_phenotype == [1, 0, 0]:
                rrow = 3
            if receptor_phenotype == [0, 1, 1]:
                rrow = 4
            if receptor_phenotype == [1, 0, 1]:
                rrow = 5
            if receptor_phenotype == [1, 1, 0]:
                rrow = 6
            if receptor_phenotype == [1, 1, 1]:
                rrow = 7

            for m in range(0, num_receptors):
                receptor_abundance_network[rrow, m] = receptor_abundance_network[rrow, m] + genotype_abundance*receptor_phenotype[m]
                if receptor_phenotype[m] == 0:
                    receptor_phenotype[m] = 1
                elif receptor_phenotype[m] == 1:
                    receptor_phenotype[m] = 0
                receptor_infection_network[rrow,m] = receptor_infection_network[rrow, m] + genotype_abundance*receptor_phenotype[m]
                genotype_infection_network[genotype_counter,m] = genotype_infection_network[genotype_counter, m] + genotype_abundance*receptor_phenotype[m]

        counter1 = counter1 + 1

    susceptible_genotypes = np.zeros((num_genotypes, num_receptors))
    susceptible_genotypes[:] = genotype_infection_network

    schisto_challenge = np.zeros(num_receptors)

    num_snail_infections = np.zeros(num_genotypes)
    num_schisto_infections = np.zeros(num_strains)
    num_schisto_left = np.full(8, miracidia_abundance/8)

    #challenging and infecting each snail and schisto
    # adjust with contact rate and (1-mu)
    for j in range(0, num_strains):
        # every schisto has an equal chance of infecting a snail - shows no preference
        # schisto are evenly dispersed across all genotypes, where they face two outcomes: establish infection or die
        for k in range(0, num_receptors):
            schisto_challenge[k] = (miracidia_abundance*schisto_distribution[j]/num_genotypes)*mucin_network[j,k]*c

        strain_snail_infections = np.matmul(susceptible_genotypes, schisto_challenge)

        #adjusting for population constraints
        for k in range(0, num_genotypes):
            if strain_snail_infections[k] > previous_population[k]:
                strain_snail_infections[k] = previous_population[k]

            num_snail_infections[k] = num_snail_infections[k] + strain_snail_infections[k]

            if num_snail_infections[k] > previous_population[k]:
                num_snail_infections[k] = previous_population[k]

        num_schisto_infections[j] = np.sum(strain_snail_infections)

        num_schisto_left[j] = num_schisto_left[j] - num_schisto_infections[j]

    schisto_infection_matrix[i, :] = num_schisto_infections

    rate_snail_infections = num_snail_infections/previous_population

    #transitionMatrices = myTransitionMatrix(inputs, previous_population, r, mu, beta, num_snail_infections, g, sig, inbr, k)
    transitionMatrices = myTransitionMatrix(inputs, previous_population, r, mu, beta, rate_snail_infections, g, sig, inbr, k)

    current_population = transitionMatrices[2]

    # calculate cohort sizes in each generation i
    next_generation = np.zeros(num_genotypes)
    for j in range(0, num_genotypes):
        next_generation[j] = current_population[0,j]
        genotype[j, i] = current_population[0, j]

    #recalculate schisto distribution
    # 1/12 turnaround, from the amount of successful infections
    # 11/12 = what's outside the snails
    previous_abundance = miracidia_abundance
    for j in range(0, num_strains):
        if i > 12:
            schisto_distribution[j] = (previous_abundance*schisto_distribution[j] - num_schisto_infections[j] + h*q*schisto_infection_matrix[i, j])
        else:
            schisto_distribution[j] = (previous_abundance*schisto_distribution[j] - num_schisto_infections[j])

    miracidia_abundance = np.sum(schisto_distribution)

    for j in range(0, num_strains):
        schisto_distribution[j] = schisto_distribution[j]/miracidia_abundance

#Results I want:
# important snail strains (gene drive, 0, 1, 2, 3 receptors)
# schisto distribution
# number snail infections
# number schisto infections
t = scipy.linspace(0, num_gens, num_gens + 1)
#
pl.figure(1)
g0 = scipy.transpose(genotype[0,:]) #genotype 0
g6 = scipy.transpose(genotype[6,:]) #genotype 0
g7 = scipy.transpose(genotype[7,:]) #genotype 7
g16 = scipy.transpose(genotype[16,:]) #genotype 0

g20 = scipy.transpose(genotype[20,:]) #genotype 20
g26 = scipy.transpose(genotype[26,:]) #genotype 0

g28 = scipy.transpose(genotype[28,:]) #genotype 7

g30 = scipy.transpose(genotype[30,:]) #genotype 30
g36 = scipy.transpose(genotype[36,:]) #genotype 0

g38 = scipy.transpose(genotype[38,:]) #genotype 7

g40 = scipy.transpose(genotype[40,:]) #genotype 40
g42 = scipy.transpose(genotype[num_genotypes-3,:]) #genotype 42
g43 = scipy.transpose(genotype[num_genotypes-2,:]) #genotype 43
g44 = scipy.transpose(genotype[num_genotypes-1,:]) #genotype 44

print(num_genotypes)

pl.plot(t, g0, label='000 & 000')
pl.plot(t, g7, label='g7')
pl.plot(t, g20, label = 'g20')
pl.plot(t, g28, label = 'g28')
pl.plot(t, g30, label = 'g30')
pl.plot(t, g36, label = 'g36')
pl.plot(t, g38, label = 'g38')
pl.plot(t, g40, label = 'g40')
pl.plot(t, g42, label = "g42")
pl.plot(t, g43, label = 'g43')
pl.plot(t, g44, label='111g & 111g')





pl.grid()
pl.xlabel('Generation')
pl.ylabel('Population %')
pl.title('Multi-Receptor Gene Drive System')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))  # sets ticks for x-axis


#TO DO
# figure out important receptors
# put in a c for the schisto success

pl.figure(2)
for i in range(0, num_genotypes):
    X = scipy.transpose(genotype[i,:])

    pl.plot(t, X, label = i)

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.ylabel('Population [# B. glabrata]')
pl.title('Multi-Receptor Gene Drive System')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))  # sets ticks for x-axis


pl.show()
    #challenge the whole snail population with each strain of schisto
    # but once a schisto infects, it infects, so why would it equally infect every snail?
    # maybe divide the schisto by 8 again
    # or 9??? for the gene drive?
    # shit how do I separate infected gene drive snails from infected other snails
    # lol wait none of the gene drive snails end up infected because dominance exists so it doesn't matter yayy
    # I think instead, just change so that we track infections to each genotype instead of infections to each phenotype

    # g = g - b*i
    #
    # if i == 1 :
    #     previous_population = initial_genotypes
    # else:
    #     #this is confusing syntax, but it just moves the generations forward
    #     previous_population = current_population
    #
    # #if time, package these parameters in a better way
    # #transitionMatrices = myTransitionMatrix(r, mu, beta, mu_i, num_alleles, allele_fitness, allele_resistance, previous_population, selfing, g)
    #
    # #current_population = transitionMatrices[2]
    #
    # #current_population = previous_population * selfingTransitionMatrix





    # 8 snail strains
# 6 gene drive alleles
# 7 (8?) schisto strains
# calculate fitness cost for each strain (do at beginning -- this stays the same bc it's based off of number of receptors)
# calculate average fitness cost for B strains (do at each time step -- this depends on previous snail strain distribution)
# calculate infection protection for each strain (do at each time step -- this depends on previous schisto strain distribution)
# keep track of successful schisto infections (do at each time step -- this depends on previous schisto strain and previous snail strain distributions)
# calculate average infection protection for B strains (do at each time step -- this depends on infection protection for each strain)
# Go through the gene drive calculation (do at each time step -- depends on average infeciton protection and average fitness cost)
# Update the strain distribution from gene drive outcomes
# adjust schisto distribution (do at each time step -- depends on successful schisto infections. Remember 1/12 turnover)
