# Getting started

## To use openomics in a project

```{code-block} python
import openomics
```

## Import GTEx Tissue-specific gene expression dataset (directly from URL)

We start with tissue-specific gene expression data set from The Genotype-Tissue Expression (GTEx) project. First, we
load the GTEx database from the path https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/, which automatically
download and parse the gene_median_tpm files to create gene expression matrix. The matrix can be accessed at `gtex.data`
.

```{code-block} python
import pandas as pd
from openomics.database import GTEx

gtex = GTEx(path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/")

gtex_transcripts = gtex.data

gene_id = pd.Index(gtex_transcripts["gene_id"].unique())
gene_name = pd.Index(gtex_transcripts["gene_name"].unique())

gtex_transcripts.head()
```
