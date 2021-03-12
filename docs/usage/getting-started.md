# Getting started

Welcome! This tutorial highlights the OpenOmics API’s core features; for in-depth details and conceptual design, see the links within, or the documentation index which has links to use cases, and API reference sections.

## Loading a single-omics dataframe

Suppose you have a single-omics dataset and would like to load them as a dataframe.

As an example, we use the `TGCA` Lung Adenocarcinoma dataset from [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD). Data tables are tab-delimited and have the following format:

| GeneSymbol      | EntrezID  | TCGA-05-4244-01A-01R-1107-07 | TCGA-05-4249-01A-01R-1107-07 | ...  |
| --------------- | --------- | ---------------------------- | ---------------------------- | ---- |
| ENSG00000223972 | 100133144 | 10.8123                      | 3.7927                       | ...  |
| ⋮ | ⋮ | ⋮ | ⋮ |

Depending on whether your data table is stored locally as a single file, splitted into multiple files, or was already a dataframe, you can load it using the class {class}`openomics.transcriptomics.Expression` or any of its subclasses.

````{tab} From a single file
If the dataset is a local file in a tabular format, OpenOmics can help you load them to Pandas dataframe.

```{code-block} python
from openomics.multiomics import MessengerRNA

mrna = MessengerRNA(
    data="https://raw.githubusercontent.com/BioMeCIS-Lab/OpenOmics/master/tests/data/TCGA_LUAD/LUAD__geneExp.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA", # A regex that matches all column name with either "GeneSymbol" or "TCGA substring
    gene_index="GeneSymbol", # This column contains the gene index
    )
```

One thing to pay attention is that the raw data file given is column-oriented where columns corresponds to samples, so we have use the argument `transpose=True` to convert to row-oriented.
> MessengerRNA (576, 20472)
````

````{tab} From multiple files (glob)
If your dataset is large, it may be broken up into multiple files with a similar file name prefix/suffix. Assuming all the files have similar tabular format, OpenOmics can load all files and contruct an integrated data table using the memory-efficient Dask dataframe.

```python
from openomics.multiomics import MessengerRNA

mrna = MessengerRNA("TCGA_LUAD/LUAD__*",
                    transpose=True,
                    usecols="GeneSymbol|TCGA",
                    gene_index="GeneSymbol")
```

> INFO: Files matched: ['LUAD__miRNAExp__RPM.txt', 'LUAD__protein_RPPA.txt', 'LUAD__geneExp.txt']
````

````{tab} From DataFrame
If your workflow already produced a dataframe, you can encapsulate it directly with {class}`openomics.transcriptomics.Expression`.

```python
import pandas as pd
import numpy as np
from openomics.multiomics import MessengerRNA

# A random dataframe of microRNA gene_id's.
df = pd.DataFrame(data={"ENSG00000194717": np.random.rand(5),
                        "ENSG00000198973": np.random.rand(5),
                        "ENSG00000198974": np.random.rand(5),
                        "ENSG00000198975": np.random.rand(5),
                        "ENSG00000198976": np.random.rand(5),
                        "ENSG00000198982": np.random.rand(5),
                        "ENSG00000198983": np.random.rand(5)},
                  index=range(5))
mrna = MessengerRNA(df, transpose=False, sample_level="sample_id")
```
````

To access the

## Creating a multi-omics dataset

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
The `luad_data` is a {class}`MultiOmics` object builds the samples list from all the samples given in each -omics data.

```none
MessengerRNA (576, 20472)
MicroRNA (494, 1870)
LncRNA (546, 12727)
SomaticMutation (587, 21070)
Protein (364, 154)
```



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
