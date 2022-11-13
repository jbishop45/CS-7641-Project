# Trajectory Analysis of Single-Cell mRNA Sequencing Data 
**TODO:** MAKE GITHUB PAGE PUBLIC AND SHARE LINK, ENSURE REPO IS PRIVATE BEFORE SUBMISSION
**Authors:** Josh Bishop, Anjana Dissanayaka, Vishal Manickman, Nina Moorman, Jay Wroe

**Instructor:** Prof. Mahdi Roozbahani

**Submitted:** November 13, 2022

## Introduction/Background 
### Topic
Single-cell mRNA sequencing (scRNA-seq) has generated petabytes of data detailing the transcriptomes<sup>I</sup> of billions of cells, generating insights ranging from biomarkers<sup>II</sup> of disease states to the discovery of new cell-types. Diverse methods for generating these data cause high variability in the number & similarity of cells sequenced, and mRNAs sequenced per cell<sup>[6,11,12]</sup>. Developing robust methods for sc-RNAseq analysis is critical to standardizing analytical results and expediting scientific discoveries<sup>[6,11,12]</sup>.

### Dataset & Features
Our primary dataset is a 27297x10396 array<sup>[5]</sup>. Rows indicate genes, and columns indicate mouse stem cells collected over a 7-day differentiation protocol. Entries are unique observations of mRNAs encoding gene R in cell C. Two scientific journal articles have identified 5 cell types at 3 stages of maturity in these data<sup>[3,5]</sup>.
Genes with few observations are statistically unreliable, and dead cells are a source of confounding noise. Regressing gene expression against background processes (e.g. cell division) further reduces noise. After cleaning, typical scRNA-seq analyses focus on 100s-1000s of genes in 1000-5000 cells<sup>[6]</sup>.

## Problem Definition & Motivation
Robust methods for interpreting transcriptomic<sup>I</sup> responses to arbitrary stimuli can yield insights into disease state evolution and aid in the development of novel therapies. Diffusion Pseudotime<sup>III</sup> analysis (DPT) is a method for ordering cells along a continuous process and is robust to batch effects on sampling density and sequencing depth, making it a flexible tool for classifying cells<sup>[3]</sup>.
Our goals are: (1) identify distinct populations of cells differentiating from a common origin, and (2) determine the physiological role of those populations.

![Proposal Figure](figures/proposal_figure.png)

### Potential results and Discussion
#### Research Questions
* **RQ1.** Can we order cells along an arbitrary branching process using scRNA-seq data? 

* **RQ2.** Can we determine the identity of cells near terminal points in the above process using scRNA-seq data?

### Methods:
#### Data Cleaning
Harvesting cells or tissues and creating single-cell suspensions inevitably damages the cell membranes of some fraction of the cells present in the sample, destroying their ability to maintain homeostatic conditions and leading to rapid deterioration of cell surface features and internal environment. Any biological assay aiming to analyze bulk or single cell characteristics must include a cleaning step to identify and exclude these dead and damaged cells.
In scRNA-seq analysis, this is typically done by studying two key parameters of an initial read-count dataset: count depth and unique gene counts. Count depth describes the total number of RNA molecules in a particular cell that were sequenced and successfully mapped to a gene, and can be computed as the sum of the observation counts of all genes in the dataset for a given cell. Unique gene counts indicate how many different genes were observed one or more times within each cell, and can be computed as the number of nonzero column entries for each cell. These parameters are selected for cleaning scRNA-seq data because mRNA molecules in the cytoplasm of a cell are extremely unstable following membrane damage and are quickly destroyed. This leads to a stark reduction in the number of mappable mRNAs present in a cell, which will typically result in low count depth. Cells with low transcriptomic activity, or a large amount of internal compartmentation can deviate from this trend, however.
Unique gene counts can be used to detect these outliers, because the primary source of mRNAs that _are_ shielded inside of organelles within the cell are mitochondria. mRNAs in mitochondria are derived from a limited subset of genes compared to mRNAs found elsewhere in the cell cytoplasm. In a healthy and metabolically active cell, the abundance of mitochondrial mRNA would not be disproportionately reflected in the number of genes that are sequenced, meaning that the number of unique gene counts would remain high. However, in the case of damaged cells, the metabolic activity would likely decrease, which would cause the ratio of mitochondrial mRNA to cytoplasmic mRNA to increase and consequently reduce the number of unique genes per cell. 
We plotted histograms of count depth and unique genes per cell, and in each case observed a bimodal distribution for each of these features, in keeping with the underlying biological theory (**Fig. 1 A,B**). This bimodality is also visible in the two-shelf structure of (**Fig. 1 C**), which shows the count depths of the unfiltered cells in descending order. A plot of count depth versus unique genes showed a single narrow distribution with a weakly exponential association (**Fig. 1 D**), indicating that for our particular dataset, high-quality filtering could be accomplished with only one of these features, but we implemented both to maintain the generalizability of our platform. In keeping with standard practices for the field we visually estimated minimum thresholds of count-depth (**Fig. 1E**) and unique genes (**Fig. 1F**) for the boundary between the “low” and “typical” populations, but we intend to implement a two-component, two-dimensional gaussian mixture model to objectively discriminate these populations in our final report.

![Data Cleaning](figures/Data_Cleaning.png)

Our data cleaning process involved first removing genes which were not observed in any cells from the count matrix, which reduced our dataset from ~27300 features to ~16200 features, removing cells with a count depth below 300, which reduced the dataset from ~10400 to ~8800 observations, and then removing cells with fewer than 500 unique genes, which further reduced the dataset to exactly 8100 observations. 

#### Data Preprocessing
We take care of missing features in our dataset by removing genes with low count values in the data cleaning process. The original dataset includes a reference list of all genes known to be expressed in mice, and it is normal for many of them to be absent in actual sequencing data. With >16000 out of 27300 possible features remaining in the dataset, confident identification of cells' position and trajectory in transcriptomic space is still very achievable. 
The first step in our data preprocessing is gene count normalization. We rely on the ``sc.pp.normalize_total`` function in order to convert the raw gene count for each cell into counts out of 10,000. By doing this, the effect of cell-to-cell variations in count depth can be removed from measures of cell-to-cell variance in the relative levels of expression of individual genes. Following this normalization, we reduce the feature-space of our dataset by removing low-variance features, which hold little information about the differences between cells in our dataset which we leverage to cluster cells and predict their trajectory in gene expression space. We isolate the most variable genes using ``sc.pp.highly_variable_genes``. We determine the min and max mean input parameters as well as the min dispersion coefficient from prior work</sup>[13]</sup>. This reduced the features space of our sequencing data from 16216 to 2242 features.

![Principle Components Variance](figures/Principle_Components_Variance.png)
# TODO: NEAREST-NEIGHBOR GRAPH

#### Unsupervised Subproject - Graph Inference
# TODO: update
We obtain a graph of nearest-neighbor cells in transcriptomic<sup>I</sup> space from the unsupervised classifier <code>sklearn.neighbors.NearestNeighbors</code>, then compute a transition matrix representing the probability of each cell transitioning into another using previously described methods<sup>[5]</sup> (see figure above). <code>scanpy.tl.dpt</code> can then derive a pseudotemporal ordering of immature cells differentiating into groups of mature cells and the mean positions in the transcriptomic<sup>I</sup> space of mature cells.
#### Supervised Subproject - Graph Annotation
# TODO: update
We determine the physiological role of the differentiated cells by mapping their transcriptomes<sup>I</sup> to representative transcriptomes<sup>I</sup> of labeled cell-types using <code>scanpy.tl.ingest</code>, which takes in a dimensionality-reduced representation of the reference data generated with <code>scanpy.tl.umap</code>.

### Results and Discussion
#### Data Cleaning
# TODO: update
#### Data Pre-Processing
# TODO: update
#### Supervised/Unsupervised Method & Results
# TODO: at least 1 supervised/unsupervised method implemented w/ results and metrics

## Proposed Timeline & Contribution Table
# TODO: update
![Proposal Timeline](figures/proposal_timeline.png)
![Gantt Chart](figures/gantt_chart.png)


## Definition of Terms
1. Transcriptome - “A transcriptome is the full range of messenger RNA, or mRNA, molecules expressed by an organism. The term ‘transcriptome’ can also be used to describe the array of mRNA transcripts produced in a particular cell or tissue type.”
2. Biomarker - A molecule whose presence is indicative of some process, e.g. the development of a disease state, or the metabolization of certain toxins.
3. Pseudotime - An abstracted notion of time or another independent variable associated with a continuous (time-dependent) process (e.g. the evolution of a system from state A to state B), which can be used to describe and order intermediate positions in that process.

## References
- Aizarani, N., Saviano, A., Sagar _et al_. A human liver cell atlas reveals heterogeneity and epithelial progenitors. _Nature_ 572, 199–204 (2019). https://doi.org/10.1038/s41586-019-1373-2 
- Anthony Gitter. Single-cell RNA-seq pseudotime estimation algorithms. 2018. https://doi.org/10.5281/zenodo.1297422 (GitHub:https://github.com/agitter/single-cell-pseudotime) 
- Haghverdi, L., Büttner, M., Wolf, F. _et al_. Diffusion pseudotime robustly reconstructs lineage branching. _Nat Methods_ 13, 845–848 (2016). https://doi.org/10.1038/nmeth.3971 
- He Z, Peng C, Li T and Li J. Cell Differentiation Trajectory in Liver Cirrhosis Predicts Hepatocellular Carcinoma Prognosis and Reveals Potential Biomarkers for Progression of Liver Cirrhosis to Hepatocellular Carcinoma. _Front Genet_ 13:858905 (2022). https://doi.org/10.3389/fgene.2022.858905 
- Klein AM, Mazutis L, Akartuna I, Tallapragada N _et al_. Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells. _Cell_ 161(5):1187-1201 (2015). https://doi.org/10.1016/j.cell.2015.04.044
- Luecken, M. D., & Theis, F. J. Current best practices in single‐cell RNA‐SEQ Analysis: A tutorial. _Mol Syst Biol_, 15(6) (2019). https://doi.org/10.15252/msb.20188746 
- MacParland, S.A., Liu, J.C., Ma, XZ. _et al_. Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. _Nat Commun_ 9, 4383 (2018). https://doi.org/10.1038/s41467-018-06318-7 
- Mu, T., Xu, L., Zhong, Y. _et al_. Embryonic liver developmental trajectory revealed by single-cell RNA sequencing in the Foxa2eGFP mouse. _Commun Biol_ 3, 642 (2020). https://doi.org/10.1038/s42003-020-01364-8 
- Pedregosa, F. _et al_. Scikit-Learn: Machine Learning in Python. _Scikit-Learn_ 12, 2825–2830 (2011). https://dl.acm.org/doi/10.5555/1953048.2078195 
- Wolf, F., Angerer, P. & Theis, F. SCANPY: large-scale single-cell gene expression data analysis. _Genome Biol_ 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0 
- Zappia L, Phipson B, Oshlack A. Exploring the single-cell RNA-seq analysis landscape with the scRNA-tools database. _PLoS Comput Biol_ 14(6): e1006245 (2018). https://doi.org/10.1371/journal.pcbi.1006245 
- Zappia, L., Theis, F.J. Over 1000 tools reveal trends in the single-cell RNA-seq analysis landscape. _Genome Biol_ 22, 301 (2021). https://doi.org/10.1186/s13059-021-02519-4
