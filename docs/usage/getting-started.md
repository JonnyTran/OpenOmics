# Getting started

Welcome! This tutorial highlights the OpenOmics API’s core features; for in-depth details and conceptual design, see the links within, or the documentation index which has links to use cases, and API reference sections.

## Loading a single-omics dataframe

Suppose you have your own -omics dataset(s) and you'd like to load them as dataframes. Depending on whether your data is
stored locally, retrieved via an online database, or was already a dataframe, you can instantiate a data object from one
of the class in the `openomics.multiomics` submodule.

````{tab} From local file
If the dataset is a local file in a tabular format, OpenOmics can help you load them to Pandas dataframes. As an example, we downloaded the `TGCA_LUAD` dataset from [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD).

```{code-block} python
from openomics.multiomics import  MessengerRNA

# Load each expression dataframe
mrna = MessengerRNA("TCGA_LUAD/LUAD__geneExp.txt",
                    transpose=True,
                    usecols="GeneSymbol|TCGA",
                    gene_index="GeneSymbol")
```

One thing to pay attention is that the raw data file given is column-oriented where columns corresponds to samples, so we have use the argument `transpose=True` to convert to row-oriented.

> MessengerRNA (576, 20472)

````

````{tab} From multiple files (glob)
If your dataset is large, it may be broken up into multiple files with a similar file name prefix. Assuming all the files have the same tabular format, OpenOmics can load all files and contruct

```python
from openomics.multiomics import Expression # a superclass of MessengerRNA

rnaseq = Expression("TCGA_LUAD/LUAD__*",
                    transpose=True,
                    usecols="GeneSymbol|TCGA",
                    gene_index="GeneSymbol")
```

````

````{tab} From DataFrame

```python
import pandas as pd
import numpy as np

# A random dataframe of microRNA gene_id's.
df = pd.DataFrame(data={"ENSG00000194717": np.random.rand(5),
                        "ENSG00000198973": np.random.rand(5),
                        "ENSG00000198974": np.random.rand(5),
                        "ENSG00000198975": np.random.rand(5),
                        "ENSG00000198976": np.random.rand(5),
                        "ENSG00000198982": np.random.rand(5),
                        "ENSG00000198983": np.random.rand(5)},
                  index=range(5))
mirna = MicroRNA(df, transpose=False, sample_level="sample_id")
mirna.expressions
```

````

## Creating a multi-omics dataset

````{tab} From local files
Assuming you already have your own dataset in local files in a tabular format, OpenOmics can help you load them to Pandas dataframes. As an example, we downloaded the `TGCA_LUAD` dataset from [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD).

```{code-block} python
from openomics.multiomics import MessengerRNA, MicroRNA, LncRNA, SomaticMutation, Protein

# Load each expression dataframe
mRNA = MessengerRNA("TCGA_LUAD/LUAD__geneExp.txt",
                    transpose=True,
                    usecols="GeneSymbol|TCGA",
                    gene_index="GeneSymbol")
miRNA = MicroRNA("TCGA_LUAD/LUAD__miRNAExp__RPM.txt",
                 transpose=True,
                 usecols="GeneSymbol|TCGA",
                 gene_index="GeneSymbol")
lncRNA = LncRNA("TCGA_LUAD/TCGA-rnaexpr.tsv",
                transpose=True,
                usecols="Gene_ID|TCGA",
                gene_index="Gene_ID")
som = SomaticMutation("TCGA_LUAD/LUAD__somaticMutation_geneLevel.txt",
                      transpose=True,
                      usecols="GeneSymbol|TCGA",
                      gene_index="gene_name")
pro = Protein("TCGA_LUAD/protein_RPPA.txt",
              transpose=True,
              usecols="GeneSymbol|TCGA",
              gene_index="GeneSymbol")

# Create an integrated MultiOmics dataset
luad_data = MultiOmics(cohort_name="LUAD", omics_data=[mRNA, mRNA, lncRNA, som, pro])

luad_data.build_samples()
```

One thing to pay attention is that the raw data file given is column-oriented where columns corresponds to samples, so we have use the argument `transpose=True` to convert to row-oriented. The `luad_data` is a :class:`MultiOmics` object builds the samples list from all the samples given in each -omics data.

> MessengerRNA (576, 20472)
> MicroRNA (494, 1870)
> LncRNA (546, 12727)
> SomaticMutation (587, 21070)
> Protein (364, 154)

````

````{tab} From URL
We start with tissue-specific gene expression data set from The Genotype-Tissue Expression (GTEx) project. First, we
load the GTEx database from the path [https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/), which automatically
download and parse the gene_median_tpm files to create gene expression matrix. The matrix can be accessed at `gtex.data`

```python
import pandas as pd
from openomics.database import GTEx

gtex = GTEx(path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/")

gtex_transcripts = gtex.data

gtex_transcripts.head()
```

````

````{tab} From DataFrame

```python
import pandas as pd
import numpy as np

# A random dataframe of microRNA gene_id's.
df = pd.DataFrame(data={"ENSG00000194717": np.random.rand(5),
                        "ENSG00000198973": np.random.rand(5),
                        "ENSG00000198974": np.random.rand(5),
                        "ENSG00000198975": np.random.rand(5),
                        "ENSG00000198976": np.random.rand(5),
                        "ENSG00000198982": np.random.rand(5),
                        "ENSG00000198983": np.random.rand(5)},
                  index=range(5))
df
```

````

## Adding clinical data as sample attributes

```python
from openomics import ClinicalData

clinical = ClinicalData("TCGA_LUAD/nationwidechildrens.org_clinical_patient_luad.txt", patient_index=)

luad_data.add_clinical_data(clinical)

```

##

## Import an external database



## Annotate your expression dataset with attributes

```{tab} Genomics attributes
.
```

```{tab} Sequences
.
```

```{tab} Interactions
.
```
