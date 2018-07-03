# thesis-public

This repository contains the code relevant to Chloe Rickard's senior honors thesis (Engineering Resistance: A CRISPR/Cas9-Mediated Gene Drive for Schistosomiasis Control).

UPDATE: This repository code has been taken down in order to conserve it for future publications. For more information, please contact rickards.chloe@gmail.com

geneSystem, geneDrive_transitionMatrix, and geneDrive_sim all work together to create the figures found in the thesis.  The model is modified from geneDrive_sim, which calls upon geneSystem and geneDrive_transitionMatrix as it runs.

geneSystem contains basic information about the system modeled: number of alleles, number of genotypes, initial population distribution, allele fitness costs, and selfing.

geneDrive_transitionMatrix is called during every time step in the geneDrive_sim.  The genotypes in each generation experience outcrossing, selfing, and logistic growth to produce the next generation's genotypes.  This program forms the basis for the DIY_transitionMatrix program.

geneDrive_sim carries out the simulation specific to the honors thesis, in the context of 8 3-receptor alleles + 1 gene drive allele, compatibility-based infection, and a gene drive. 

DIY_transitionMatrix is a user-friendly program that takes any number of alleles to create a transition matrix tailored to the user's desired conditions.  The DIY_transitionMatrix allows the user to toggle gene drive and selfing, and input initial genotype distribution and allele-associated fitness costs.  We hope this program will be helpful to computational biologists who hope to track genotype changes over time, especially in the presence of a gene drive.

The DIY_transitionMatrix script depends on some functions in the usefulFunctions script. Use them together!
