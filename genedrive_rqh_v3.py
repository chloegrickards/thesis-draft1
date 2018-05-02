## ASSUMPTIONS ##
# every schisto has an equal chance of infecting a snail - shows no preference
# schisto are evenly dispersed across all genotypes, where they face two outcomes: establish infection or die
# Resistance is dominant -- snails that are heterozygous for a receptor inherit that receptor. [Source = Webster]
# Not putting in a resistance efficacy term, that's theoretically accounted for with the matching thing
# Assuming density dependence cercarial production. Mutliple miracidia can infect a snail, and worm burden per
# strain is a result of infection intensity in snails, per strain.
# Beta is straight up just calibrated to fit the expected prevalence level. I wish there was a way around it, but
# there's precedence for this in the literature

#Sources cited for parameters:

# Woolhouse 1991  -- Snail life span = 4-12 weeks
# Pointier et al. 1991 -- Snail growth rate, etc
# Costa 2004 -- egg tables for r and sig
# Jarne 1993 - selfing rate

# To Do
# Make lots of plots!
# - plot gene drive at 0% efficiency
# - take picture from even initial conditions
# - capture system at equilibrium, adjust initial conditions (is this even right?)
# - turn on gene drive in the system
# - turn on/off selfing in the system
# - revisit and actually plot unique genotypes, even if they're wrapped up in the sameish thing

# Results I want:
# Worm burden in humans

import numpy as np
from rqh_geneSystem import geneSystemInputs
import scipy
import pylab as pl

from rqh_transitionMatrix_v2 import myTransitionMatrix

## MODELING TOGGLES ##

Stoggle = 1 #force of infection
selftoggle = 1 #selfing toggle
gdtoggle = 1 #gene drive toggle

## PARAMETERS ##
# all parameters are based on literature values, except for beta

# Intrinsic growth rate
# From Costa 2004
# In my system, total population = 370
r = 54

# Natural death rate, per generation
mu = 0.2

# Gene drive efficiency
# Source = Hammond et al. 2016
g = 0.95

# Selfing rate
# Estimated from Jarne 1993
sig = 0.15

# Cost of inbreeding
# From Costa 2004
# Eggs/snail outcrossing = 54
# Eggs/snail inbreeding = 8.1
# therefore, cost = 0.85
inbr = .85  # cost of inbreeding

# Per capita infection rate of snails by a female worm
# Feng says 0.0004, Sokolow says 0.000004
#in Sokolow 2015 and Feng 2004
# I've calibrated it to my expected behavior
beta = 0.0000008

# Worm fecundity
# eggs/snail generation
# difficult to estimate, but it's in the ballpark of Cheever 1968
f = 250

#Dominance coefficient
# Resistant genes are assumed dominant (gamma = 1)
# Source = Webster 2004
gamma = 1

# Human population
H = 1000

# Infection probability of snail to man
# lam for lambda
# Converted from day to snail generation
# Source = Sokolow 2015
lam = 0.0005*365/4

# Estimate from schistosome life span = 3.3 years,
# 3-4 snail generations in a year
# Source = Woolhouse 1991
mu_w = 1/12

## SYSTEM INPUTS and SETUP ##
# this takes inputs from a user-defined gene drive system and sets up for iterations over mutliple generations

#Information from the rqh_geneSystem function:
# Number of alleles
# Genotype-allele dictionary
# Initial genotype population distribution
# Allele fitness costs
# Selfing (True/False)
inputs = geneSystemInputs(r)
num_alleles = inputs[0]
genotype_dict = inputs[1]
initial_population = inputs[2]
allele_fitness = inputs[3]
selfing = inputs[4]

# Inputs specific to this receptor-matching system:
# Number of receptors in snails
# Number of genotypes, from number of alleles
# Number of schisto strains
# Number of snail generations
num_receptors = 3
num_genotypes = int(num_alleles * (num_alleles + 1) / 2)
num_strains = 8
num_gens = 50

## INFECTION CYCLE SETUP ##

# Force of infection
S = np.zeros((num_strains, num_gens+1))
S[:,0] = beta*f*0.5*(lam*(0.015*370))*H*Stoggle


# Mean worm burden
W = np.zeros((num_strains, num_gens+1))
W[:,0] = lam*(0.015*370)

## SNAIL POPULATION SETUP ##
# Many of these vectors will be useful for plotting, later:

# Total initial snail population
total_population0 = np.sum(initial_population)
#Population carrying capacity
k_pop = total_population0

# Genotype abundance (# of snails) and distribution (% of snails), per generation
genotype_abundance = np.zeros((num_genotypes, num_gens+1))
genotype_percent = np.zeros((num_genotypes, num_gens+1))
for i in range(0, num_genotypes):
    genotype_abundance[i,0] = initial_population[i]
    genotype_percent[i,0] = initial_population[i]/total_population0

# Infected snails, per genotype, per generation
infected_snails = np.zeros((num_genotypes, num_gens+1))

# Number of infected and susceptible snails, for use after adjusting for multiple miracidial infections:
num_infections = np.zeros(num_gens + 1)
num_infections[0] = 0.015*total_population0
num_susceptible = np.zeros(num_gens + 1)
num_susceptible[0] = total_population0 - num_infections[0]

## INFECTION CYCLE SETUP ##

# Initial infection prevalence
# This comes from the initial miracidia infection prevalence in snails
# 1-2% of snails show miracidia infection (Source = Sokolow 2015)
initial_prevalence = 0.015*total_population0

# Force of infection
# Initial value of S is a modification of the original equation
S = np.zeros((num_strains, num_gens+1))
S[:,0] = beta*f*0.5*(lam*(initial_prevalence))*H*Stoggle

# Mean worm burden
# Initial value of W is a modification of the original equation
W = np.zeros((num_strains, num_gens+1))
W[:,0] = lam*(initial_prevalence)

## SCHISTOSOME POPULATION SETUP ##
#  Many of these vectors will be useful for plotting, later:

# Input initial miracidia populations here
# Here, value is similar to that of force of infection, but without snail infection probability, beta
miracidia_abundance0 = np.full(8, f*0.5*W[0,0]*H)
total_miracidia0 = np.sum(miracidia_abundance0)

# Miracidia abundance (# of miracidia) and distribution (% of miracidia), per generation
# Miracidia abundance does not account for mortality from egg to miracidia
miracidia_abundance = np.zeros((num_strains, num_gens+1))
miracidia_percent = np.zeros((num_strains, num_gens+1))

for i in range(0, num_strains):
    miracidia_abundance[i,0] = miracidia_abundance0[i]
    miracidia_percent[i,0] = miracidia_abundance0[i]/total_miracidia0

# Vector for number of successfully infectious schistosomes
infectious_schisto = np.zeros((num_strains, num_gens+1))

## RECEPTOR AND MUCIN NETWORKS ##

# Receptor network, for snail alleles
receptor_network = np.array([[0, 0, 0],
                            [0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0],
                            [0, 1, 1],
                            [1, 0, 1],
                            [1, 1, 0],
                            [1, 1, 1]])

# Mucin network, for schisto phenotypes
# Identical to the receptor network, repeated here for clarity
mucin_network = np.array([[0, 0, 0],
                            [0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0],
                            [0, 1, 1],
                            [1, 0, 1],
                            [1, 1, 0],
                            [1, 1, 1]])

## TOGGLE EFFECTS ##

if selftoggle == 0:
    selfing = False

if gdtoggle == 0:
    initial_population[num_genotypes-1] = 0

## SIMULATION ##

#Iterate through each generation
for i in range(1, num_gens+1):

    # Determine resistance phenotype from receptor genotypes
    # Starts with the blank phenotype to be filled in, and for each genotype, determines phenotype abundance
    # End result = network of phenotype abundances, summed across all genotypes
    # This network is useful in calculating overall infections
    susceptible_network = np.zeros((num_genotypes, num_receptors))
    counter1 = 0
    for j in range(0, num_alleles):
        for k in range(counter1, num_alleles):

            #starts with a blank phenoytpe array
            receptor_phenotype = [0, 0, 0]

            # finds current genotype
            genotype_counter = int(genotype_dict[j,k])

            # finds allele numbers
            allele1 = j
            allele2 = k

            # counts 111g as 111
            if allele1 == 8:
                allele1 = 7
            if allele2 == 8:
                allele2 = 7

            # if alleles both = 0, then phenotype = 0
            # if alleles are heterozygous, phenotype = gamma*1
            # if alleles are both = 1, it's 1
            # End result = receptor phenotype for a given genotype
            for l in range(0, num_receptors):
                if receptor_network[allele1,l] == receptor_network[allele2,l] and receptor_network[allele1,l] == 0:
                    receptor_phenotype[l] = 0
                elif receptor_network[allele1,l] == receptor_network[allele2,l]:
                    receptor_phenotype[l] = gamma*1
                else:
                    receptor_phenotype[l] = 1

            # Builds an overall array of the receptor phenotype abundance in a population of snails
            # Row entry = previous entry + abundance of current phenotype
            # Flips the receptor phenotype from presence of receptor (which is easier to input)
            #     to susceptibility (which is easier to calculate infections from)
            for m in range(0, num_receptors):
                if receptor_phenotype[m] == 0:
                    receptor_phenotype[m] = 1
                else:
                    receptor_phenotype[m] = 1 - receptor_phenotype[m]
                susceptible_network[genotype_counter,m] = susceptible_network[genotype_counter, m] + genotype_abundance[genotype_counter, i-1]*receptor_phenotype[m]

        counter1 = counter1 + 1

    # Empty vector, holds the population of miracidia that make up the challenge to the snails
    miracidia_challenge = np.zeros((num_strains, num_receptors))

    # Calculate the infectious miracidia challenge
    for j in range(0, num_strains):
        miracidia_challenge[j, :] = S[j,i-1] * mucin_network[j,:]

    # Infection Network Matrix
    # Calculates the number of infections, per genotype of snail, and per strain of schisto
    # Infections are determined by compatibility
    INM = np.matmul(susceptible_network, miracidia_challenge.T)

    # Find the number of infected snails, per genotype
    for j in range(0, num_genotypes):
        infected_snails[j, i] = np.sum(INM[j, :])

    # Find the number of infectious (successful) schisto, per strain
    for j in range(0, num_strains):
        infectious_schisto[j, i] = np.sum(INM[:,j])

    # Find snail infections for each genotype for this generation
    gen_snail_infections = infected_snails[:,i]

    # Run snail population, after infection, through transition matrix
    tm_out = myTransitionMatrix(inputs, genotype_abundance[:, i-1], r, mu, gen_snail_infections, g, sig, inbr, k_pop)

    # Find new snail population, after outcrossing and selfing
    current_population = tm_out[0]
    total_population = np.sum(current_population)

    # Rewrite infected snails, adjusting for multiple miracidial infections
    infected_snails[:, i] = tm_out[1]
    num_infections[i] = np.sum(infected_snails[:, i])

    # Total susceptible snails
    num_susceptible[i] = total_population - num_infections[i]

    # calculate genotype values sizes for this generation
    #Moves the population abundance backward, in order to calculate the next iteration
    for j in range(0, num_genotypes):
        genotype_abundance[j, i] = current_population[0, j]
        genotype_percent[j, i] = current_population[0, j]/total_population


    # update mean parasite burden in humans and force of infection for snails
    for j in range(0, num_strains):
        W[j,i] = lam*infectious_schisto[j,i] - mu_w * W[j, i-1]

        if W[j, i] < 0:
            W[j, i] = 0

        S[j,i] = beta*f*0.5*W[j,i]*H*Stoggle

        #Miracidia abundance, NOT accounting for miracidia mortality
        miracidia_abundance[j, i] = f*0.5*W[j,i]*H

    total_miracidia = np.sum(miracidia_abundance[:,i])
    for j in range(0, num_strains):
        miracidia_percent[j, i] = miracidia_abundance[j, i]/total_miracidia

print(S)

## PLOTTING ##

t = scipy.linspace(0, num_gens, num_gens + 1)


pl.figure(1)

# for i in range(0, num_genotypes):
#     y = scipy.transpose(genotype_abundance[i, :])
#     pl.plot(t, y, label = i)

g0 = scipy.transpose(genotype_abundance[0,:]) #genotype 0
g7 = scipy.transpose(genotype_abundance[7,:]) #genotype 7
g20 = scipy.transpose(genotype_abundance[20,:]) #genotype 20
g28 = scipy.transpose(genotype_abundance[28,:]) #genotype 7
g30 = scipy.transpose(genotype_abundance[30,:]) #genotype 30
g36 = scipy.transpose(genotype_abundance[36,:]) #genotype 0
g38 = scipy.transpose(genotype_abundance[38,:]) #genotype 7
g40 = scipy.transpose(genotype_abundance[40,:]) #genotype 40
g42 = scipy.transpose(genotype_abundance[num_genotypes-3,:]) #genotype 42
g43 = scipy.transpose(genotype_abundance[num_genotypes-2,:]) #genotype 43
g44 = scipy.transpose(genotype_abundance[num_genotypes-1,:]) #genotype 44

pl.plot(t, g0, label='000 & 000')
pl.plot(t, g7, label='000 & 111')
pl.plot(t, g28, label = '100 & 111')
pl.plot(t, g30, label = '011 & 011')
pl.plot(t, g36, label = '101 & 110')
pl.plot(t, g38, label = '101 & 111g')
pl.plot(t, g40, label = '110 & 111')
pl.plot(t, g42, label = "111 & 111")
pl.plot(t, g43, label = '111 & 111g')
pl.plot(t, g44, label='111g & 111g')

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.ylabel('Population [# B. glabrata]')
pl.title('Biomphalaria glabrata Genotype Distribution')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))  # sets ticks for x-axis

pl.figure(2)

strain0 = miracidia_percent[0, :]
strain3 = miracidia_percent[3, :]
strain6 = miracidia_percent[6, :]
strain7 = miracidia_percent[7, :]

pl.plot(t, strain0, label='0 antigens (000)')
pl.plot(t, strain3, label='1 antigen (001, 101, 100)')
pl.plot(t, strain6, label='2 antigens (011, 101, 110)')
pl.plot(t, strain7, label='3 antigens (111)')

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.ylabel('Population [% S. mansoni strain]')
pl.title('Schistosoma mansoni Strain Distribution')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))  # sets ticks for x-axis

pl.figure(3)

pl.plot(t, num_infections, label = "I")
pl.plot(t, num_susceptible, label = "S")

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.ylabel('Population [# B. glabrata]')
pl.title('Susceptible and Infected B. glabrata')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))  # sets ticks for x-axis

pl.show()
