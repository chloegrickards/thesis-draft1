## ASSUMPTIONS ##
# Resistance is dominant -- snails that are heterozygous for a receptor inherit that receptor. [Source = Webster]
# Receptors are completely effective at clearing miracidiae, so long as they are expressed (i.e. no resistance efficacy term)
# Density dependence cercarial production. Mutliple miracidia can infect a snail, and worm burden per
# strain is a result of infection intensity in snails, per strain.
# Beta is calibrated to fit the expected prevalence level. I made this decision based on literature precedence

# Sources cited for parameters:

# Woolhouse 1991  -- Snail life span = 4-12 weeks
# Pointier et al. 1991 -- Snail growth rate, etc
# Costa 2004 -- egg tables for r and sig
# Jarne 1993 - selfing rate

import numpy as np
from rqh_geneSystem import geneSystemInputs
import scipy
import pylab as pl

from rqh_transitionMatrix_v2 import myTransitionMatrix

## MODELING TOGGLES ##

Stoggle = 1  # force of infection (scale between 0 and 1), normally 1
selftoggle = 0.15  # selfing toggle, toggles selfing rate sigma (scale b/n 0 and 1), normally 0.15
gdtoggle = 1  # gene drive toggle (on/off)
gamtoggle = 1  # dominance toggle (scale b/n 0 and 1), normally 1

## PARAMETERS ##
# all parameters are based on literature values, except for beta

# Intrinsic growth rate
# Source = Costa 2004
r = 54

# Natural death rate, per generation
# Source = Woolhouse 1991
mu = 0.2

# Gene drive efficiency
# Source = Hammond et al. 2016
g = 0.95

# Selfing rate
# normally 0.15
# Source = Jarne 1993
sig = selftoggle

# Cost of inbreeding
# From Costa 2004
# Eggs/snail outcrossing = 54
# Eggs/snail inbreeding = 8.1
# therefore, cost = 0.85
inbr = .85  # cost of inbreeding

# Per capita infection rate of snails by a female worm
# Source 1 = Feng 2004 (beta = 0.0004)
# Source 2 = Sokolow 2015 (beta = 0.000004)
# Here, beta is calibrated to expected infection prevalence
beta = 0.0000008

# Worm fecundity
# eggs/snail generation
# Source = Cheever 1968
f = 250

# Dominance coefficient
# Resistant genes are assumed dominant (gamma = 1)
# Source = Webster 2001
gamma = gamtoggle

# Human population
H = 1000

# Infection probability of snail to man
# lam for lambda
# Converted from day to snail generation
# Source = Sokolow 2015
lam = 0.0005 * 365 / 4

# Estimate from schistosome life span = 3.3 years,
# 3-4 snail generations in a year
# Source = Woolhouse 1991
mu_w = 1 / 12

## SYSTEM INPUTS and SETUP ##
# this takes inputs from a user-defined gene drive system and sets the foundation for iterations over multiple generations

# User inputs from the rqh_geneSystem function, setting up the general genetic system:
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

# User inputs specific to this receptor-matching system:
# Number of receptors in snails
# Number of genotypes, from number of alleles
# Number of schisto strains
# Number of snail generations
num_receptors = 3
num_genotypes = int(num_alleles * (num_alleles + 1) / 2)
num_strains = 8
num_gens = 50

## SNAIL POPULATION SETUP ##
# Many of these vectors will be useful for plotting, later:

# Effect of gene drive toggle
if gdtoggle == 0:
    initial_population[num_genotypes - 1] = 0.000000001

# Total initial snail population
total_population0 = np.sum(initial_population)

# Population carrying capacity
k_pop = total_population0

# Genotype abundance (# of snails) and distribution (% of snails), per generation
genotype_abundance = np.zeros((num_genotypes, num_gens + 1))
genotype_percent = np.zeros((num_genotypes, num_gens + 1))
for i in range(0, num_genotypes):
    genotype_abundance[i, 0] = initial_population[i]
    genotype_percent[i, 0] = initial_population[i] / total_population0

# Infected snails, per genotype, per generation
infected_snails = np.zeros((num_genotypes, num_gens + 1))

# Number of infected and susceptible snails, for use after adjusting for multiple miracidial infections:
# Note: infected_snails contains the total number of miracidia that infect a snail during a generation,
#       while num_infections constains the number of infected snails, period. num_infections does no count the number
#       of miracidia in a snail, just whether or not the snail is infected.
num_infections = np.zeros(num_gens + 1)
num_infections[0] = 0.02 * total_population0
percent_infections = np.zeros(num_gens + 1)
percent_infections[0] = num_infections[0] / total_population0

num_susceptible = np.zeros(num_gens + 1)
num_susceptible[0] = total_population0 - num_infections[0]
percent_susceptible = np.zeros(num_gens + 1)
percent_susceptible[0] = num_susceptible[0] / total_population0

## INFECTION CYCLE SETUP ##

# Initial infection prevalence
# This comes from the initial miracidia infection prevalence in snails
# 1-2% of snails show miracidia infection
# Source = Sokolow 2015
initial_prevalence = 0.02 * total_population0

# Force of infection
# Initial value of S is a modification of the original equation
S = np.zeros((num_strains, num_gens + 1))
S[:, 0] = beta * f * 0.5 * (lam * (initial_prevalence)) * H * Stoggle

# Mean worm burden
# Initial value of W is a modification of the original equation
W = np.zeros((num_strains, num_gens + 1))
W[:, 0] = lam * (initial_prevalence)

# For plotting W:
W_avg = np.zeros(num_gens + 1)
W_avg[0] = np.sum(W[:, 0]) / num_strains

# Vector for number of successfully infectious schistosomes
infectious_schisto = np.zeros((num_strains, num_gens + 1))

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
# Identical to the receptor network, repeated here for clarity and later manipulation
mucin_network = np.array([[0, 0, 0],
                          [0, 0, 1],
                          [0, 1, 0],
                          [1, 0, 0],
                          [0, 1, 1],
                          [1, 0, 1],
                          [1, 1, 0],
                          [1, 1, 1]])

# Effect of selfing toggle
if selftoggle == 0:
    inputs[4] = False

## SIMULATION ##

# Iterate through each generation
for i in range(1, num_gens + 1):

    # Determine resistance phenotype from receptor genotypes
    # Starts with the blank phenotype to be filled in, and for each genotype, determines phenotype abundance
    # End result = network of phenotype abundances, summed across all genotypes
    # This network is useful in calculating overall infections
    susceptible_network = np.zeros((num_genotypes, num_receptors))
    counter1 = 0
    for j in range(0, num_alleles):
        for k in range(counter1, num_alleles):

            # starts with a blank phenotype array
            receptor_phenotype = [0, 0, 0]

            # finds current genotype
            genotype_counter = int(genotype_dict[j, k])

            # finds allele numbers
            allele1 = j
            allele2 = k

            # 111g gene drive genotypes express 111 phenotypes. These conditions express that fact.
            if allele1 == 8:
                allele1 = 7
            if allele2 == 8:
                allele2 = 7

            # if alleles both = 0, then phenotype = 0
            # if alleles are heterozygous, phenotype = gamma
            # if alleles are both = 1, it's 1
            # End result = receptor phenotype for a given genotype
            for l in range(0, num_receptors):
                if receptor_network[allele1, l] == receptor_network[allele2, l] and receptor_network[allele1, l] == 0:
                    receptor_phenotype[l] = 0
                elif receptor_network[allele1, l] != receptor_network[allele2, l]:
                    receptor_phenotype[l] = gamma
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
                susceptible_network[genotype_counter, m] = susceptible_network[genotype_counter, m] + \
                                                           genotype_abundance[genotype_counter, i - 1] * \
                                                           receptor_phenotype[m]

        counter1 = counter1 + 1

    # Empty vector, holds the population of miracidia that make up the challenge to the snails
    miracidia_challenge = np.zeros((num_strains, num_receptors))

    # Calculate the infectious miracidia challenge
    for j in range(0, num_strains):
        miracidia_challenge[j, :] = S[j, i - 1] * mucin_network[j, :]

    # Infection Network Matrix
    # Calculates the number of infections, per genotype of snail, and per strain of schisto
    # Infections are determined by compatibility
    INM = np.matmul(susceptible_network, miracidia_challenge.T)

    # Find the number of infected snails, per genotype
    for j in range(0, num_genotypes):
        infected_snails[j, i] = np.sum(INM[j, :])

    # Find the number of infectious (successful) schisto, per strain
    for j in range(0, num_strains):
        infectious_schisto[j, i] = np.sum(INM[:, j])

    # Find snail infections for each genotype for this generation
    gen_snail_infections = infected_snails[:, i]

    # Run snail population, after infection, through transition matrix
    tm_out = myTransitionMatrix(inputs, genotype_abundance[:, i - 1], r, mu, gen_snail_infections, g, sig, inbr, k_pop)

    # Find new snail population, after outcrossing and selfing
    current_population = tm_out[0]
    total_population = np.sum(current_population)

    # Rewrite infected snails, adjusting for multiple miracidial infections
    infected_snails[:, i] = tm_out[1]
    num_infections[i] = np.sum(infected_snails[:, i])
    percent_infections[i] = num_infections[i] / total_population

    # Total susceptible snails
    num_susceptible[i] = total_population - num_infections[i]
    percent_susceptible[i] = num_susceptible[i] / total_population

    # Calculate genotype values sizes for this generation
    # Moves the population abundance backward, in order to calculate the next iteration
    for j in range(0, num_genotypes):
        genotype_abundance[j, i] = current_population[0, j]
        genotype_percent[j, i] = current_population[0, j] / total_population

    # update mean parasite burden in humans and force of infection for snails
    for j in range(0, num_strains):
        W[j, i] = lam * infectious_schisto[j, i] - mu_w * W[j, i - 1]

        # forbids negative worms
        if W[j, i] < 0:
            W[j, i] = 0

        S[j, i] = beta * f * 0.5 * W[j, i] * H * Stoggle

    # Calculates overall mean worm burden
    W_avg[i] = np.sum(W[:, i]) / num_strains

## PLOTTING ##

t = scipy.linspace(0, num_gens, num_gens + 1)

pl.figure(1)

twopercent = np.full(num_gens + 1, 0.02)

pl.plot(t, percent_infections, label="I")
pl.plot(t, percent_susceptible, label="S")
pl.plot(t, twopercent, color="red")

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.ylabel('Population [% B. glabrata]')
pl.title('Susceptible and Infected B. glabrata')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))

pl.figure(2)

pl.plot(t, W_avg)

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.ylabel('Number of worms per capita')
pl.title('Average Worm Burden in Humans')
pl.xticks(scipy.linspace(0, num_gens, 11))

pl.figure(3)

p0 = scipy.transpose(genotype_percent[0, :])  # genotype 0
p1 = scipy.transpose(genotype_percent[1, :])  # genotype 1
p4 = scipy.transpose(genotype_percent[4, :])  # genotype 4
p7 = scipy.transpose(genotype_percent[7, :])  # genotype 7
p9 = scipy.transpose(genotype_percent[9, :])  # genotype 9
p10 = scipy.transpose(genotype_percent[10, :])  # genotype 10
p12 = scipy.transpose(genotype_percent[12, :])  # genotype 12
p15 = scipy.transpose(genotype_percent[15, :])  # genotype 15
p25 = scipy.transpose(genotype_percent[25, :])  # genotype 25
p30 = scipy.transpose(genotype_percent[30, :])  # genotype 30
p33 = scipy.transpose(genotype_percent[33, :])  # genotype 33
p36 = scipy.transpose(genotype_percent[36, :])  # genotype 36
p42 = scipy.transpose(genotype_percent[num_genotypes - 3, :])  # genotype 42
p43 = scipy.transpose(genotype_percent[num_genotypes - 2, :])  # genotype 43
p44 = scipy.transpose(genotype_percent[num_genotypes - 1, :])  # genotype 44

pl.plot(t, p0, label='000 & 000')
pl.plot(t, p1, label='000 & 001')
pl.plot(t, p4, label='000 & 011')
pl.plot(t, p7, label='000 & 111')
pl.plot(t, p9, label='001 & 001')
pl.plot(t, p10, label='001 & 010')
pl.plot(t, p12, label='001 & 011')
pl.plot(t, p15, label='001 & 111')
pl.plot(t, p25, label='011 & 100')
pl.plot(t, p30, label='011 & 011')
pl.plot(t, p33, '--', label='011 & 111')
pl.plot(t, p36, '--', label='101 & 110')
pl.plot(t, p42, '--', label="111 & 111")
pl.plot(t, p43, 'm--', label='111 & 111g')
pl.plot(t, p44, 'r--', label='111g & 111g')

pl.grid()
pl.xlabel('Time [# B. glabrata generations]')
pl.xticks(scipy.linspace(0, num_gens, 11))

pl.ylabel('Population [% B. glabrata]')
pl.title('Biomphalaria glabrata Genotype Distribution')

ax = pl.subplot(111)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)

pl.show()