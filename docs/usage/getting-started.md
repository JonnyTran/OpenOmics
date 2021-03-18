# Getting started

Welcome! This tutorial highlights the OpenOmics API’s core features; for in-depth details and conceptual guides, see the links within, or the documentation index which has links to use cases, and API reference sections.

## Loading a single-omics dataframe

Suppose you have a single-omics dataset and would like to load them as a dataframe.

As an example, we use the `TGCA` Lung Adenocarcinoma dataset from [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD). Data tables are tab-delimited and have the following format:

| GeneSymbol | EntrezID  | TCGA-05-4244-01A-01R-1107-07 | TCGA-05-4249-01A-01R-1107-07 | ...  |
| ---------- | --------- | ---------------------------- | ---------------------------- | ---- |
| A1BG       | 100133144 | 10.8123                      | 3.7927                       | ...  |
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

mrna = MessengerRNA("TCGA_LUAD/LUAD__*", # Files must be stored locally
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
---
To access the {class}`DataFrame`, simply use {obj}`mrna.expressions`:
```python
print(mrna.expressions)
```
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>GeneSymbol</th>
      <th>A1BG</th>
      <th>A1BG-AS1</th>
      <th>A1CF</th>
      <th>A2M</th>
    </tr>
    <tr>
      <th>sample_index</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244-01A-01R-1107-07</th>
      <td>26.0302</td>
      <td>36.7711</td>
      <td>0.000</td>
      <td>9844.7858</td>
    </tr>
    <tr>
      <th>TCGA-05-4249-01A-01R-1107-07</th>
      <td>120.1349</td>
      <td>132.1439</td>
      <td>0.322</td>
      <td>25712.6617</td>
    </tr>
  </tbody>
</table>
</div>

<br/>

## Creating a multi-omics dataset

With multiple single-omics, each with different sets of genes and samples, you can use the {class}`openomics.MultiOmics` to integrate them.

```{code-block} python
from openomics.multiomics import MessengerRNA, MicroRNA, LncRNA, SomaticMutation, Protein

path = "https://raw.githubusercontent.com/BioMeCIS-Lab/OpenOmics/master/tests/data/TCGA_LUAD/"

# Load each expression dataframe
mRNA = MessengerRNA(path+"LUAD__geneExp.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")
miRNA = MicroRNA(path+"LUAD__miRNAExp__RPM.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")
lncRNA = LncRNA(path+"TCGA-rnaexpr.tsv",
    transpose=True,
    usecols="Gene_ID|TCGA",
    gene_index="Gene_ID")
som = SomaticMutation(path+"LUAD__somaticMutation_geneLevel.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")
pro = Protein(path+"protein_RPPA.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")

# Create an integrated MultiOmics dataset
luad_data = MultiOmics(cohort_name="LUAD", omics_data=[mRNA, mRNA, lncRNA, som, pro])
# You can also add individual -omics one at a time `luad_data.add_omic(mRNA)`

luad_data.build_samples()
```
The `luad_data` is a {class}`MultiOmics` object builds the samples list from all the samples given in each -omics data.

> MessengerRNA (576, 20472)
> MicroRNA (494, 1870)
> LncRNA (546, 12727)
> SomaticMutation (587, 21070)
> Protein (364, 154)

To access individual -omics data within `luad_data`, such as the {obj}`mRNA`, simply use the `.` accessor with the class name {class}`MessengerRNA`:
```python
luad_data.MessengerRNA
# or
luad_data.data["MessengerRNA"]
```

<br/>

## Adding clinical data as sample attributes

When sample attributes are provided for the study cohort, load it as a data table with the {class}`openomics.clinical.ClinicalData`, then add it to the {class}`openomics.multiomics.MultiOmics` dataset to enable querying for subsets of samples across the multi-omics.

```python
from openomics import ClinicalData

clinical = ClinicalData(
    "https://raw.githubusercontent.com/BioMeCIS-Lab/OpenOmics/master/tests/data/TCGA_LUAD/nationwidechildrens.org_clinical_patient_luad.txt",
    patient_index="bcr_patient_barcode")

luad_data.add_clinical_data(clinical)

luad_data.clinical.patient
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>bcr_patient_uuid</th>
      <th>form_completion_date</th>
      <th>histologic_diagnosis</th>
      <th>prospective_collection</th>
      <th>retrospective_collection</th>
    </tr>
    <tr>
      <th>bcr_patient_barcode</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244</th>
      <td>34040b83-7e8a-4264-a551-b16621843e28</td>
      <td>2010-7-22</td>
      <td>Lung Adenocarcinoma</td>
      <td>NO</td>
      <td>YES</td>
    </tr>
    <tr>
      <th>TCGA-05-4245</th>
      <td>03d09c05-49ab-4ba6-a8d7-e7ccf71fafd2</td>
      <td>2010-7-22</td>
      <td>Lung Adenocarcinoma</td>
      <td>NO</td>
      <td>YES</td>
    </tr>
    <tr>
      <th>TCGA-05-4249</th>
      <td>4addf05f-3668-4b3f-a17f-c0227329ca52</td>
      <td>2010-7-22</td>
      <td>Lung Adenocarcinoma</td>
      <td>NO</td>
      <td>YES</td>
    </tr>
  </tbody>
</table>
</div>

Note that in the clinical data table, `bcr_patient_barcode` is the column with `TCGA-XX-XXXX` patient IDs, which matches
that of the `sample_index` index column in the `mrna.expressions` dataframe.

````{note}
In our `TCGA_LUAD` example, mismatches in the `bcr_patient_barcode` sample index of clinical dataframe may happen because the `sample_index` in `mRNA` may have a longer form `TCGA-XX-XXXX-XXX-XXX-XXXX-XX` that contain the samples number and aliquot ID's. To make them match, you can modify the index strings on-the-fly using the [Pandas's extensible API](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.str.slice.html):
```python
mRNA.expressions.index = mRNA.expressions.index.str.slice(0, 12) # Selects only the first 12 characters
```
````

<br/>

## Import an external database

Next, we may want to annotate the genes list in our RNA-seq expression dataset with genomics annotation. To do so, we'd need to download annotations from the [GENCODE database](https://www.gencodegenes.org/), preprocess annotation files into a dataframe, and then match them with the genes in our dataset.

OpenOmics provides a simple, hassle-free API to download the GENCODE annotation files via FTP with these steps:
1. First, provide the base `path` of the FTP download server - usually found in the direct download link on GENCODE's website. Most of the time, selecting the right base `path` allows you to specify the specific species, genome assembly, and database version for your study.
2. Secondly, use the `file_resources` dict parameter to select the data files and the file paths required to construct the annotation dataframe. For each entry in the `file_resources`, the key is the alias of the file required, and the value is the filename with the FTP base `path`.

   For example, the entry `{"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz"}` indicates the GENCODE class to preprocess a `.gtf` file with the alias `"long_noncoding_RNAs.gtf"`, downloaded from the FTP path `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.long_noncoding_RNAs.gtf.gz`

   To see which file alias keys are required to construct a dataframe, refer to the docstring in {class}`openomics.database.sequence.GENCODE`.

```python
from openomics.database import GENCODE

gencode = GENCODE(
    path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
    file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                    "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                    "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz", # lncRNA sequences
                    "transcripts.fa": "gencode.v32.transcripts.fa.gz" # mRNA sequences
                    },
    npartitions=0, # if > 1, then use Dask partition the dataframe and leverage out-of-core multiprocessing
)
```
To access the attributes constructed from the combination of annotations `long_noncoding_RNAs.gtf` and `
basic.annotation.gtf`, use:

```python
gencode.data
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_id</th>
      <th>gene_name</th>
      <th>index</th>
      <th>seqname</th>
      <th>source</th>
      <th>feature</th>
      <th>start</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000243485</td>
      <td>MIR1302-2HG</td>
      <td>0</td>
      <td>chr1</td>
      <td>HAVANA</td>
      <td>gene</td>
      <td>29554</td>
      <td>31109</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000243485</td>
      <td>MIR1302-2HG</td>
      <td>1</td>
      <td>chr1</td>
      <td>HAVANA</td>
      <td>transcript</td>
      <td>29554</td>
      <td>31097</td>
    </tr>
  </tbody>
</table>
</div>


<br/>

## Annotate your expression dataset with attributes
With the annotation database, you can perform a join operation to add gene attributes to your {class}`openomics.transcriptomics.Expression` dataset. To annotate attributes for the `gene_id` list `mRNA.expression`, you must first select the corresponding column in `gencode.data` with matching `gene_id` keys. The following are code snippets for a variety of database types.

````{tab} Genomics attributes
```python
luad_data.MessengerRNA.annotate_attributes(gencode,
    index="gene_id",
    columns=['gene_name', 'start', 'end', 'strand'] # Add these columns to the .annotations dataframe
    )
```

````

````{tab} Sequences
```python
luad_data.MessengerRNA.annotate_sequences(gencode,
    index="gene_name",
    agg_sequences="all", # Collect all sequences with the gene_name into a list
    )
```
````

````{tab} Disease Associations
```python
from openomics.database.disease import DisGeNet
disgenet = DisGeNet(path="https://www.disgenet.org/static/disgenet_ap1/files/downloads/", curated=True)

luad_data.MessengerRNA.annotate_diseases(disgenet, index="gene_name")
```
````

---
To view the resulting annotations dataframe, use:
```python
luad_data.MessengerRNA.annotations
```


For more detailed guide, refer to the [annotation interfaces API](../modules/openomics.annotate.md).
