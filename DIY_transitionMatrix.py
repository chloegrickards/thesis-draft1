import numpy as np
import usefulFunctions as fn
import scipy
import pylab as pl

## NOTE TO THE USER ##

# Hi! This program allows the user to track changes in genotype over time.  The program takes any number of alleles
# and finds classic Hardy-Weinberg proportions -- not only under normal conditions, but also in  in the context
# gene drive and self-reproduction.  The user needs to input information into the system to make it work.
# Necessary information = number of alleles, intrinsic growth rate, allele-associate fitness costs, population
#   carrying capacity, natural mortality rate, presence of self-reproduction, and presence of gene drive
# If selfing exists, the selfing rate and inbreeding cost become necessary
# If the gene drive exists, the gene drive efficiency becomes necessary
# Optional information = infection rate, infection-induced mortality rate

# This code is currently filled in with example parameters, but feel free to replace the values!


## NECESSARY USER INPUTS ##

# Number of alleles, including the gene drive (if any)
# This program assumes that the gene drive (if any) is the LAST allele
num_alleles = 3

# Intrinsic growth rate
# In other words, the number of offspring from one mating event, no costs applied
# The actual number of offspring is calculated as an average of the two alleles
# True growth rate = (r - allele_fitness_1 + r - allele_fitness_2)/2
r = 25

# Fitness costs associated with each allele
# Ideally a function of r
# The gene drive fitness cost (if any) goes LAST
allele_fitness = [0, 0.2*r, 0.3*r]

# Population carrying capacity
# Max number of individuals in a population
# Used to determine logistic growth
k = 100

# Natural mortality rate
# Percent population that dies each generation
# Value must be b/n 0 and 1
mu = 0.25

# Does self-reproduction occur? True or False
selfing = True

# Is there a gene drive? True or False
gene_drive = True

## CONDITIONALLY NECESSARY USER INPUTS ##

# Selfing rate
# Proportion of population that uses selfing to reproduce
# Value must be b/n 0 and 1
sig = 0.2

# Selfing cost
# Removal rate due to inbreeding
# Value must be b/n 0 and 1
inbr = 0.85

# Gene Drive Efficacy
# CRISPR/Cas9 is 90-99% effective, and on average 98% effective
# Source = Hammond 2016
g = 0.98

## OPTIONAL USER INPUTS ##

# Infection rate
# If no infections, set = 0
# Value must be b/n 0 and 1
beta = 0.3

# Infection-induced mortality rate
# Percent of infected population that dies each generation
# If no infections, set = 0
# Value must be b/n 0 and 1
mu_i = 0.5

# Initial Population Distribution
# Listed as the abundance (# of organisms) for each genotype
# 3 options exist:
#
# 1. Uniform distribution for each genotype:
# initial_population = np.full(# genotypes, abundance)
#
# 2. Uniform wild-type alleles, plus a gene drive seed
wt_abundance = 30
gd_seed = 10
initial_population = fn.findInitialPopulation(num_alleles, wt_abundance, gd_seed)
#
# 3. Input your own values
# Note: The array must be written in the correct allele order -- think of iterating through the
# genotypes as in a nested for loop. For Allele 1, iterate through all the combinations of Allele 2. Then step
# Allele 1 forward and iterate through unique combinations of Allele 2.
# For example in a system of three alleles (A, B, and C), write the genotype abundances in this order:
# AA, AB, AC, BB, BC, CC
# No double counting of alleles. AB = BA and is written as AB.
#
# initial_population = [...]

# Number of generations (i.e. time steps)
num_gens = 100


## SETUP ##

# Establish the number of genotypes from the number of alleles
num_genotypes = int(num_alleles * (num_alleles + 1) / 2)

# Create genotype distribution vector
# Units = abundance of individuals
genotype_distribution = np.zeros((num_gens+1, num_genotypes))
genotype_distribution[0, :] = initial_population

# Create genotype percent vector
# Units = % of total population
total_population0 = np.sum(initial_population)
genotype_percent = np.zeros((num_gens+1, num_genotypes))
genotype_percent[0, :] = initial_population/total_population0


## SIMULATION ##

# Calculate changes in genotype after each genertaion
# May include other functions that happen in a generation
for i in range(1, num_gens + 1):
    genotype_distribution[i, :] = fn.findTransitionMatrix(num_alleles, r, allele_fitness, genotype_distribution[i-1, :], mu, beta, mu_i, g, gene_drive, selfing, sig, inbr, k)
    total_population = np.sum(genotype_distribution[i, :])
    genotype_percent[i, :] = genotype_distribution[i, :]/ total_population

## PLOT ##

# Sets up time domain
t = scipy.linspace(0, num_gens, num_gens + 1)

# Plots every genotype in the system
pl.figure(1)

# Iterate through, transpose, and plot
for i in range(0, num_genotypes):
    y = scipy.transpose(genotype_percent[:, i])
    pl.plot(t, y, label = 'Genotype ' + str(i))

# Note that calling the following shows which alleles match up to each genotype
# Allele 1 = row, allele 2 = column
genotype_dict = fn.findGenotypeDictionary(num_alleles)
print(genotype_dict)

# Plot it
pl.grid()
pl.xlabel('Time [# generations]')
pl.ylabel('Population [# organisms]')
pl.title('Genotype Distribution')
pl.legend(loc='best')
pl.xticks(scipy.linspace(0, num_gens, 11))  # sets ticks for x-axis

pl.show()
