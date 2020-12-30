=====
Usage
=====

To use openOmics in a project::

    import openomics

## Importing the openomics library & modules


```python
from openomics import MultiOmics, MessengerRNA, MicroRNA, LncRNA
```

## Import GTEx Tissue-specific gene expression dataset (directly from URL)


```python
from openomics.database import GTEx

gtex = GTEx(path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/")

gtex_transcripts = gtex.data

gene_id = pd.Index(gtex_transcripts["gene_id"].unique())
gene_name = pd.Index(gtex_transcripts["gene_name"].unique())
```

## Import GENCODE human release 32

```python
from openomics.database import GENCODE

gencode = GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                  file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                  "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                                  "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                  "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                  remove_version_num=True)
```

## Load single omics expressions for MessengerRNA, MicroRNA, LncRNA
```python
# Gene Expression
messengerRNA_id = gtex_transcripts_gene_id & pd.Index(gencode.df[gencode.df["gene_type"] == "protein_coding"]["gene_id"].unique())

messengerRNA = MessengerRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(messengerRNA_id)],
                            transpose=True, gene_index="gene_name", usecols=None, npartitions=4)

# MicroRNA expression
microRNA_id = pd.Index(ensembl.df[ensembl.df["gene_biotype"] == "miRNA"]["gene_id"].unique()) & gtex_transcripts_gene_id

microRNA = MicroRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(microRNA_id)],
                    gene_index="gene_id", transpose=True, usecols=None, )

# LncRNA expression
lncRNA_id = pd.Index(gencode.df[gencode.df["gene_type"] == "lncRNA"]["gene_id"].unique()) & gtex_transcripts_gene_id
lncRNA = LncRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(lncRNA_id)],
                gene_index="gene_id", transpose=True, usecols=None, )
```

## Create a MultiOmics dataset
```python
gtex_data = MultiOmics(cohort_name="GTEx Tissue Avg Expressions")

gtex_data.add_omic(messengerRNA)
gtex_data.add_omic(microRNA)
gtex_data.add_omic(lncRNA)

gtex_data.build_samples()
```
