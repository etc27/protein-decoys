# protein-decoys
## Description
Selective degradation of proteins is essential for the proper functioning of the cell. The ubiquitin-proteasome pathway of protein degradation consitutes the major method by which proteins are degraded intracellularly. In this pathway, ubiquitination results in the selective degradation of proteins, and the selectivity is determined by E3 ubiquitin ligases. Another layer of regulation potentially occurs by "endogenous decoys," which have high homology to the substrate recognition domains of E3 ubiquitin ligases, but cannot recruit the ubiquitination machinary. In this way, endogenous decoys may act as inhibitors of ubiquitination. However, the extent and specificity of endogenous decoys is relatively uncharacterized. This project aims to identify putative endogenous decoys for several species. After identification of the endogenous decoys, a heatmap will be generated to visualize the similarity between the endogenous decoys and their E3 ligase "targets."
## Pipeline
1. Obtain target protein sequences as fasta files from UniProt and the coordinates of the domain of interest for each protein
2. Remove F-box or U-box domain using remove_FBDdomain.R or remove_UBXdomain.R to obtain a list of all substrate recognition domains
3. Run edited sequences through customized microProtein Prediction Program (miP3) to obtain a fasta file of all decoys. The microProtein Prediction Program (miP3) is a free software (Copyright (c) 2014 - Niek de Klein, Enrico Magnani, Seung Y. Rhee) that has been edited for use in this workflow.
4. Use evals script to obtain a matrix of evals for each decoy vs. target protein
5. Generate heatmap using make_heatmap.R
## Results
The results from this workflow can be visualized in the RShiny app [here](http://etcorcoran.shinyapps.io/webpage_heatmap).
