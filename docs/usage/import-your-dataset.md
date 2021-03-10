# Importing your own multi-omics dataset

## Load single omics expressions for MessengerRNA, MicroRNA, LncRNA

We instantiate the MessengerRNA, MicroRNA and LncRNA -omics expression data from `gtex.data`. Since the gene expression
were not seperated by RNA type, we use GENCODE and Ensembl gene annotations to filter the list of mRNA, miRNA, and
lncRNAs.

```{code-block} python
from openomics import MessengerRNA, MicroRNA, LncRNA

# Gene Expression
messengerRNA_id = gtex_transcripts_gene_id & pd.Index(gencode.data[gencode.data["gene_type"] == "protein_coding"]["gene_id"].unique())

messengerRNA = MessengerRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(messengerRNA_id)],
                           transpose=True, gene_index="gene_name", usecols=None, npartitions=4)

# MicroRNA expression
microRNA_id = pd.Index(ensembl.data[ensembl.data["gene_biotype"] == "miRNA"]["gene_id"].unique()) & gtex_transcripts_gene_id

microRNA = MicroRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(microRNA_id)],
                   gene_index="gene_id", transpose=True, usecols=None, )

# LncRNA expression
lncRNA_id = pd.Index(gencode.data[gencode.data["gene_type"] == "lncRNA"]["gene_id"].unique()) & gtex_transcripts_gene_id
lncRNA = LncRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(lncRNA_id)],
               gene_index="gene_id", transpose=True, usecols=None, )
```

## Create a MultiOmics dataset

Now, we create a MultiOmics dataset object by combining the messengerRNA, microRNA, and lncRNA.

```{code-block} python
   from openomics import MultiOmics

   gtex_data = MultiOmics(cohort_name="GTEx Tissue Avg Expressions")

   gtex_data.add_omic(messengerRNA)
   gtex_data.add_omic(microRNA)
   gtex_data.add_omic(lncRNA)

   gtex_data.build_samples()
```
