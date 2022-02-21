# Differential Expression Analysis of the projection neurons in the olfactory bulb.

We find genetic differences of projection neurons using dimensionality reduction techniques and unsupervised machine learning algorithms. The results of the algorithm with the experimental verification is published in https://www.jneurosci.org/content/41/30/6449.abstract 


# The main functions and the data
"OB_DEA.mat" is the main function that performs the Differential Expression Analysis on the data saved in the file "l2_neurons_olfactory.loom" (Zeisel et al., 2018)

The external functions used in "OB_DEA.mat":

- selectGenes.m - searches for the top of overdispersed genes 
- tsne.m  - performs tSNE dimensionality reduction on the top of overdispersed genes (Laurens van der Maaten, 2010)
- HDBSCAN.m - performs clustering of the data in tsne space (Jordan Sorokin, 2017)
- signifGenes.m - finds the differentially expressed genes between the clusters 

The top of overdispersed genes used in the paper is saved in "top_overd_Genes_for_all.mat"
The tsne output used in the paper is saved in "tsne_final.mat" 



