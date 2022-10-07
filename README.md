# Trajectory Analysis and Cluster Labeling of Single-Cell mRNA Sequencing Data
<h4 style="text-align:right">
Josh Bishop, Anjana Dissanayaka, Vishal Manickman, Nina Moorman, Jay Wroe
</h4>
## Introduction/Background 
### Topic
Single-cell mRNA sequencing (scRNA-seq) has generated petabytes of data detailing the transcriptomesI of billions of cells, generating insights ranging from biomarkersII of disease states to the discovery of new cell-types. Diverse methods for generating these data cause high variability in the number & similarity of cells sequenced, and mRNAs sequenced per cell[6,11,12]. Developing robust methods for sc-RNAseq analysis is critical to standardizing analytical results and expediting scientific discoveries[6,11,12].

### Dataset & Features
Our primary dataset is a 27297x10396 array[5]. Rows indicate genes, and columns indicate mouse stem cells collected over a 7-day differentiation protocol. Entries are unique observations of mRNAs encoding gene R in cell C. Two authors have identified 5 cell-types at 3 stages of maturity in these data[3,5].
Genes with few observations are statistically unreliable, and dead cells are a source of confounding noise. Regressing gene expression against background processes (e.g. mitosis) further reduces noise. After cleaning, typical scRNAseq analyses focus on 100s-1000s of genes in 1000-5000 cells[6].

## Problem definition & Motivation

![Proposal Figure](proposal_figure.PNG)

Robust methods for interpreting transcriptomicI responses to arbitrary stimuli can yield insights into disease state evolution and aid in the development of novel therapies. Diffusion PseudotimeIII analysis (DPT) is a method for ordering cells along a continuous process and is robust to batch effects on sampling density and sequencing depth, making it a flexible tool for classifying cells[3].
Our goals are: (1) identify distinct populations of cells differentiating from a common origin, and (2) determine the physiological role of those populations.

### Methods:
#### Unsupervised Subproject - Graph Inference
We will obtain a graph of nearest-neighbor cells in transcriptomic space from the unsupervised classifier sklearn.neighbors.NearestNeighbors, then compute a transition matrix representing the probability of each cell transitioning into another using previously described methods[5]. scanpy.tl.dpt can then derive a pseudotemporal ordering of immature cells differentiating into groups of mature cells and the mean positions in the transcriptomic space of mature cells.
#### Supervised Subproject - Graph Annotation
We will determine the physiological role of the differentiated cells by mapping their transcriptomes to representative transcriptomes of labeled cell-types using scanpy.tl.ingest, which takes in a dimensionality-reduced representation of the reference data generated with scanpy.tl.umap.

### Potential results and Discussion
#### Research Questions
* RQ1. Can we order cells along an arbitrary branching process using scRNA-seq data? 

* RQ2. Can we determine the identity of cells near terminal points in the above process using scRNA-seq data?

#### Quantitative Metrics 
_Unsupervised_ - We will use a clustering dispersion metric (ex. Calinski-Harabasz Index) to evaluate the quality of our KNN clusters. We will use kendall rank correlation[3] to compare pseudotemporal ordering of cells to the collection time of those cells. Good correlation will indicate that the transition matrix orders cells according to primarily time-driven processes.

_Supervised_ - Given the ground truths, we can employ matching-based (ex. F-measure), entropy-based (ex. Normalized mutual information), or pairwise measures (ex. Jaccard coefficient) to evaluate supervised classification.

All methods will use a Euclidean distance metric in transcriptomic space.


## Proposed Timeline & Contribution Table

![Proposal Timeline](proposal_timeline.png)


## Definition of Terms
Transcriptome - “A transcriptome is the full range of messenger RNA, or mRNA, molecules expressed by an organism. The term ‘transcriptome’ can also be used to describe the array of mRNA transcripts produced in a particular cell or tissue type.”
Biomarker - A molecule whose presence is indicative of some process, e.g. the development of a disease state, or the metabolization of certain toxins.
Pseudotime - An abstracted notion of time or another independent variable associated with a continuous (time-dependent) process (e.g. the evolution of a system from state A to state B), which can be used to describe and order intermediate positions in that process.

## References
1. Aizarani, N., Saviano, A., Sagar _et al_. A human liver cell atlas reveals heterogeneity and epithelial progenitors. _Nature_ 572, 199–204 (2019). https://doi.org/10.1038/s41586-019-1373-2 
2. Anthony Gitter. Single-cell RNA-seq pseudotime estimation algorithms. 2018. https://doi.org/10.5281/zenodo.1297422 (GitHub:https://github.com/agitter/single-cell-pseudotime) 
3. Haghverdi, L., Büttner, M., Wolf, F. _et al_. Diffusion pseudotime robustly reconstructs lineage branching. Nat Methods 13, 845–848 (2016). https://doi.org/10.1038/nmeth.3971 
4. He Z, Peng C, Li T and Li J. Cell Differentiation Trajectory in Liver Cirrhosis Predicts Hepatocellular Carcinoma Prognosis and Reveals Potential Biomarkers for Progression of Liver Cirrhosis to Hepatocellular Carcinoma. Front Genet 13:858905 (2022). https://doi.org/10.3389/fgene.2022.858905 
5. Klein AM, Mazutis L, Akartuna I, Tallapragada N _et al_. Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells. Cell 161(5):1187-1201 (2015). https://doi.org/10.1016/j.cell.2015.04.044
6. Luecken, M. D., & Theis, F. J. Current best practices in single‐cell RNA‐SEQ Analysis: A tutorial. Mol Syst Biol, 15(6) (2019). https://doi.org/10.15252/msb.20188746 
7. MacParland, S.A., Liu, J.C., Ma, XZ. _et al_. Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nat Commun 9, 4383 (2018). https://doi.org/10.1038/s41467-018-06318-7 
8. Mu, T., Xu, L., Zhong, Y. _et al_. Embryonic liver developmental trajectory revealed by single-cell RNA sequencing in the Foxa2eGFP mouse. Commun Biol 3, 642 (2020). https://doi.org/10.1038/s42003-020-01364-8 
9. Pedregosa, F. _et al_. Scikit-Learn: Machine Learning in Python. Scikit-Learn 12, 2825–2830 (2011). https://dl.acm.org/doi/10.5555/1953048.2078195 
10. Wolf, F., Angerer, P. & Theis, F. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0 
11. Zappia L, Phipson B, Oshlack A. Exploring the single-cell RNA-seq analysis landscape with the scRNA-tools database. PLoS Comput Biol 14(6): e1006245 (2018). https://doi.org/10.1371/journal.pcbi.1006245 
12. Zappia, L., Theis, F.J. Over 1000 tools reveal trends in the single-cell RNA-seq analysis landscape. Genome Biol 22, 301 (2021). https://doi.org/10.1186/s13059-021-02519-4
