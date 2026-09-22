# Single-cell mapping of the Acute Myeloid Leukemia Transcriptome

# 1. Executive Summary and Clinical Rationale

# The Clinical Mission

The objective of this project is to compare the single-cell transcriptome of AML cells in patients with Healthy Bone Marrow cells and develop a deeper and insightful picture of what differentiates the two. This project provides an end-to-end high-throughput computational pipeline using Scanpy(python) to map the transcriptomic landscape of AML.

# The Single-Cell advantage over Bulk Transcriptomics

While Bulk-RNA sequencing also provides important insight into the nature of a diseased tissue through its transriptome, it offers an overall snapshot into the environment of a tumour microenvironment but does not display the intricate and deeper inner workings of the cells involved. The individial cell phenotype is exceptionally important to consider if we wish to put our finger onto a specific cog in this vast sea of moving complex parts and say with confidence, this is the problem, and therefore this is how we fix it.



# 2. Data ingestion, Directory Architecture, and Matrix Transposition

# File Structuring and Automated Tracking

The raw data was downloaded from GEO which contained compressed digital expressoin matrices of the patients the samples were derived from. These samples were extracted using 'glob', 'os', and 'gzip, to retrieve, un-zip, and add cell metadata on the fly using the file naming conventnions for subsequent analysis.

# Folder Architecture

scrna-cancer-project/
├── data/
│   ├── GSM3587923_BM1.dem.txt.gz      <-- Raw compressed healthy matrix
│   ├── GSM3587932_AML1.dem.txt.gz     <-- Raw compressed clinical tumor matrix
│   └── combined_raw.h5ad              <-- Consolidated pipeline master file
└── notebooks/
    └── 01_quality_control.ipynb       <-- Active working Jupyter workspace

# The Matrix Transposition Restraint

To be able to use Scanpy and its respective function calls, the Anndata object being used has to have a specific orientation, in that the rows need to be individual cells and columns need to be genes. But the Anndata object retrieved was the other way around, therefore a tranpositon step necessary:

with g.zip.open(filepath, "rt") as f:
    adata_sample = sc.read_text(f).T

Without the ".T" step the data would be interpreted the other way around and all downstream analysis would be illegible

# Automated Metadata Engineering via Lexical Parsing

Since the file itself only contains raw numbers related to gene count and cell number, there is no metadata about the cells and genes themselves. The file's extensions ('dem.tzt.gz') were removed and using the file's string-names metadata was dynamically constructed and added to the cell metadata table ('adata.obs'). 

adata_sample.obs["sample"] = sample_name

if "BM" in sample_name:
    adata_sample["condition"] = "Normal"
else:
    adata_sample["condition"] = "AML"

This ensures that when all 40 individual files are merged together onto one sheet using 'ad.concat(join="outer")', every single cell can tracked back to its orignial pateint cohort and disease

# 3. Qualty Control and Hard Threshold Filtering Ledger

# The Computational cleaning rationale

Raw single-cell sequencing outputs are often contaminated with cells laced with damage due to perhaps processed in cell isolation, microfluidic capture, and library preperation. Hence, some filtering needs to take place to keep only single that are metabolically stable and able to give us the best results from our downstream analysis.

[RAW UNFILTERED DATA MATRIX] -> 41,090 Cells
                     │
                     ├──> Drops cells with < 500 unique genes  (Purges Ambient RNA)
                     ├──> Drops cells with > 5,000 unique genes (Purges Multiplets)
                     └──> Drops cells with > 5% MT- RNA        (Purges Lysed Cells)
                     │
       [CLEANED STRUCTURAL MATRIX] -> 40,201 Cells


# 1) Minimum Gene Threshold ('min_genes = 500')
    Technical Artifact: Ambient RNA / Emptry Droplets
    Low-level mechanism: An oil droplet is supposed to capture a cell with its content and a barcoded bead. Often times the droplet will capture the bead but not the cell but due to nearby cell suspensions it will pick other cells' RNA leading to an almost ghost like cell RNA capture.

# 2) Maximum Gene Threshold ('max_genes = 5,000')
    Technical Artifact: Doublets/Multiplets
    Low-level mechanism: Due to standard physical crowding it is possible that a single oil droplet picks up more than one cell, perhaps even more than two; this can give signifacntly inflated numbers again hurting our downstream analysis

# 3) Mitochondrial Fraction Ceiling ('pct_counts_mt < 5')
    Biological flaw: Dying / Lysed cells
    Low-level mechanism: Physical stress during clinical biopsy extraction and processing can cause the cell membrane to lyse and rupture. Most of the floating mRNA will this way exit the cell leaving behind the heavy mitochondria with mitochondrial RNA sitting inside the dead cell carcass giving you high mitochondrial reads. Cutting them off at 5% accounts for these damaged lysed cells which would again affect your analysis.


# 4 Machine Learning Doublet Detection (Scrublet)

# Limitations of Hard Thresholding

Even though in the previous hard thresholding filtering took place, there is one potential gap left that needs to be accounted for. The max gene filter can clear out the obvious big doublets or multiplets but there can be instances where a doublet gene reading does amount to a reasonable reading instead of a very high one, essentially masking the presence of a doublet. If left as is they can distort cell-state trajectoties and generate artificial intermediate cell types.

# The Scrublet Algorithmic Framework

A 3-step process is underwent by Scrublet to evaluate the 40,201 surviving cells"

# 1) Doublet Simulation: 
    The engine will randomly average out 2 cells' gene count making artificial doublets
# 2) K-Nearest Neighbors (KNN) Graph Embedding: 
    The real cell observations and the newly generated artificial doublet calculations are plotted together onto a high-dimensional coordinated space. A local neighborhood graph is then constructed which calculates the local geometric neighbors for each single data point.\
# 3) Doublet Scoring & Thresholding: 
    Every cell is then basically audited based on its neighborhood density; cells with a higher density of fake doublets around them get a higher doublet score closer to 1.


SCRUBLET KNN NEIGHBORHOOD AUDIT
               
               [Real Cell] -> Surrounded by clean cells = Low Score (0.02) -> KEEP
               [Real Cell] -> Surrounded by fake doublets = High Score (0.85) -> DROP

# Execution Results & Filtering Verdict

The leukemia dataset was then evaluated by the algorithm using 'expected_doublet_rate=0.06' (6.00%), derived directly from the 10x Genomics microfluidic loading manual for a ~10,000 cell target.

The optimal vertical threshold separation line was mathematically computed to be at Doublet Score = 0.79.

Verdict: The algorithm isolated and flagged exactly 3 hidden homotypic doublets sitting on the dangerous side
Action: Using localized local negation operators ('~') to purge the flagged indices, what was left behind was an absolute pristine standalone matrix:

Slicing out the machine-learning flagged doublets cleanly
adata = adata[~adata.obs["predicted_doublet"]].copy()

The final post-doublet filtration baseline was locked at 40,198 cells

# 5 Normalization, Log-Transformation, & Highly Variable Gene (HVG) Selection

# Erasing Sequencing Depth Bias via Total Count Normalization

Due to stochastic (random) vairations in microfluidic droplet processing and machine sequencing efficiency some cells are read deeper than others giving them a higher gene count simply because they were read more. This way whats being counted is absolute abundance instead of relative abundance.

To standardize the baseline, the pipeline executes the count scaling using the formula:

Normalized Expression: Raw Gene Count / All Raw Counts in Cell) * 1000

Target sum was set to 10,000 'target_sum=1e4' and in doing so every cell's transcript profile is mathematically summed up to an identical total depth of 10,000 molecules; the matrix shifts from absolute raw tallies into relative cell-specific proportions.

# Compressing Linear Scale Variance via Log Transformation ('log1p')

A problem that can arise is that common housekeeping genes could have numbers in the thousands while rare oncogenic transcripts could have a value less than 1.0. In this case the common housekeeping gene would completely trump the rare gene and give skewed results. The pipeline solves this by applying a natural logarithmic transformation:

Y = ln(x+1)

This is executed via 'sc.pp.log1p(adata)', this transformation compresses extreme outliers while ensuring that a raw count of '0' maps perfectly back to a log-transformed value of '0.0'. This consequently eradicates the liner gap making sure all genes are evaluated on a balanced, geometrically equivalent scale.

# Feature Selection: Highly Variable Gene (HVG) Extraction

About 85% of the genes are common housekeeping genes that dilute the signal in the matrix and waste processing power.

The pipeline uses the **Seurat-flavor dispersion algorithm** ('sc.pp.highly_variable_genes') leaving out the uninformative and flat houskeeping genes by calculating every gene's mean expressoin against its variance (dispersion) across all 40,198 cells:

* The Cutoffs: Parameters are tightly bound at 'min_mean=0.0125', 'max_mean=3', and 'min_disp=0.5'.
* The Result: The algorithm flagged and isolated exactly 3,068 Highly Variable Genes that display significant, fluctuating, biological swings acriss the dataset.
* The Action: The uninformative grey genes were dropped, saving a focused workbook to 'processed_hvg.h5ad'.


# 6 High-Dimensional Scaling, Mainfold Mapping, and Graph Clustering

# Preserving Feature Weight usin Z-score Unit variance Scaling

Post log-normalization, individual genes themselves still have different absolute baseline scales, some genes will fluctuate between 5 and 10 while rare oncogenic genes will fluctuate between 0.1 and 0.5. If plotted using these numbers on geometric distance calculators it would distort cell-to-cell distance vectors.

To establish geometric parity, z-score standardization takes place via 'sc.pp.scale(adata, max_value=10)':
* Centering: Every gene's mean expression across the 40,198 cells becomes 0
* Scaling: Evey gene's variance is scaled to a standard unit of 1
* Ceiling: Extreme outlier values beyond 10 standard deviations are truncared at the value 10, this prevents one outlier value from warping the trajectory of downstream linear embeddings

# Compressing the coordinate space (PCA & UMAP)

1. PCA (Prinicpial Component Analysis): Linear data reduction (by 'sc.tl.pca') of the 3,068 highly variable genes leads to the finalization of 40 prinicpial components or PCs. This captures the core structural axes of biological varaitions essentially silencing all the background noise.
2. Neighborhood Graph Construction: A cell-to-cell connectivity graphs ('sc.pp.neighbors') maps every cell to its 10 closest neighbors, in the 40-dimensional PCA space.
3. UMAP Projection: Non-linear manifold approximatoin reduces the complexity of the 40 PCAs into a 2-dimensional graph with the axes of UMAP1 and UMAP2. Since closeness on a UMAP plot indicates transriptomic similarity, cells with similar molecular states naturally condense into clear clusters.

# Community detection via Leiden Algorithm

To put borders around these clusters, we used sc.tl.leiden for a clustering algorithm that uses the neighborhood graph like an interconnected social network that essentially jots lines down areas where connectivity lines begin to get weaker.

LEIDEN RESOLUTION TUNING IMPACT
                 
      Resolution = 0.1 -> Coarse View -> 3 Clustered Continents
      Resolution = 0.5 -> Optimal View -> 21 Distinct Biological Lineages
      Resolution = 1.5 -> Hyper-Split -> 50+ Artificial Micro-Clusters

The industry-standard resolution of 0.5 was used to cluster 21 distinct biological lineages on the graph, making differential gene annotation easier for more downstream analysis.


# 7 Differential Expression (DE) and Quantiative Compositional Analysis 

# The Biological Translation via Raw Matrix T-Testing

Because the multi-dimensional mapping step required a Z-score scaled matrix, the true magnitude of individual gene expression for every gene is warped in the primary workspace. To determine medically valid, un-warped fold changes, the pipeline utilizes the 'use_raw=True' parameter within the differential expression framework ('sc.tl.rank_gene_groups').

This commands a 2-sided T-test to bypass the scaled Z-scores and evaluate the normalized, un-scaled transcript counts inside 'adata.raw'. By running a comparative trial isolating each cluster against the other 20 groups combined, background housekeeping signatures are mathematically subtracted leaving behind the unique specific genetic fingerprints of each cellular community.

# Cell Atlas Compositional Insights (AML vs. Normal Control)

Executing relative abundance cross-tabulation ('pd.crosstab') across the 40,198 cell matrix exposes a profound narrative of competitive microenvironmental subversion and immune hijacking within the clinical cohorts:

| Lineage Profile Identity | AML Cohort Fraction (%) | Normal Control Fraction (%) | Clinical / Biological Infiltration Status |

| **Healthy Ribosomal Progenitors** | 1.19% | 47.41% | **Niche Displacement:** Normal developmental engines are aggressively crowded out. 
| **Antigen Presenting Myeloid Cells (HLA-DR+)** | 10.71% | 0.02% | **Immune Hijacking:** Massive tumor-induced recruitment of immunosuppressive shields. 
| **Naive/Memory T-Cells (IL7R+)** | 19.34% | 10.33% | **Reactive Infiltration:** Host immune cells swarm the niche but are blinded by the tumor. 
| **Pathological Myeloblasts (MPO+)** | 3.97% | 0.13% | **Malignant Core:** The primary mutated diagnostic leukemia blast footprint.

# Core strategic Inferences:

1. **Malignant Proportional Dominance**: The complete collapse of the Healthy Risbosomal Progenitor Pool mathematically proves that the leukemia clone doesn't merely coexist; it actively alters bone marrow signaling to displace normal hematopoiesis.
   
2. **Microenvironmental co-option**: The massive expansion of Myeloid cells and T-cell recutiment reveals that tumor is hijacking the host's immune system to use these hyper-expanded myeloid networks as structures to act as barriers between them and the cytotoxic clearance cells; so basically local immune evasion and shielding.

