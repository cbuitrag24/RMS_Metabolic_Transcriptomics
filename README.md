# RMS_Metabolic_Transcriptomics
Reproducible transcriptomic analysis of metabolic gene expression and myogenic differentiation in rhabdomyosarcoma using GSE108022 and GSE28511.
## RNA-seq data source

The RNA-seq count matrix was obtained directly from the NCBI Gene Expression Omnibus (GEO), accession GSE108022.

GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108022

Original supplementary files:
- GSE108022_counts_matrix.txt.gz
- GSE108022_metadata.txt.gz

The original count matrix contains 106 samples. For the primary differential expression analysis, only samples from platform GPL16791 were selected (64 rhabdomyosarcoma tumours and five normal skeletal-muscle samples).

The R scripts used to prepare and analyse these data are provided in the scripts directory.
 Copiar texto
