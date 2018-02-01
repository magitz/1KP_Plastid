# 1KP_Plastid
Plastid phylogenomics data and scripts for Gitzendanner et al. 2018 AJB paper:

Gitzendanner, M. A., P. S. Soltis, G. K.-S. Wong, B. R. Ruhfel, D. E. Soltis. 2018. Plastid phylogenomic analysis of green plants: a billion years of evolutionary history. _American Journal of Botany_ in press.

This github repository contains the data, as well as the scripts used for the Gitzendanner et al. 2018 

## Data directory

Within the data directory there is a file called **cpDNA.trimal.gappyout.1Kaa.160715.phy**. This is the aligned, concatenated data file that was analyzed for the phylogenetic analyses reported in the paper.

There is also a file called **Gitzendanner_etal_2018_GreenPlantTree.tree** which is the ExaML bootstrap tree obtained in the analyses descibed in the paper.

The **individual_genes** folder has the individual gene alignments. 

## Scripts directory

Within the scripts directory are the scripts that were used to process the 1KP transcriptome data. It should be emphasized that these scripts are provided to document what was done, not to make a generally applicable pipeline for generic datasets. Many steps in the scripts make assumptions about file and scaffold naming that will be specific to the 1KP project and data organization and cluster used in this paper. Many of the scripts will be helpful as starting points for others, but are not intended to be a readily deployable, general purpose pipeline.
