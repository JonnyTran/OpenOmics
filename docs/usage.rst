=====
Usage
=====

To use openomics in a project
-----------------------------

.. code:: python

   import openomics


Import GTEx Tissue-specific gene expression dataset (directly from URL)
-----------------------------------------------------------------------
We start with tissue-specific gene expression data set from The Genotype-Tissue Expression (GTEx) project. First, we load the GTEx database from the path https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/, which automatically download and parse the gene_median_tpm files to create gene expression matrix. The matrix can be accessed at `gtex.data`.


.. code:: python

   from openomics.database import GTEx

   gtex = GTEx(path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/")

   gtex_transcripts = gtex.data

   gene_id = pd.Index(gtex_transcripts["gene_id"].unique())
   gene_name = pd.Index(gtex_transcripts["gene_name"].unique())

   gtex_transcripts.head()

Import GENCODE human release 32
-------------------------------
Next, we can annotate the genes in our GTEx expression dataset with genomics annotation from GENCODE. In this example, we use the URL path prefix "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/" which specifies the species and release version. We also pass a dictionary `file_resources`, with key-value pairs where the key is name of file and value is the suffix of the file download URL.

For example, file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz"} will download file located at `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.long_noncoding_RNAs.gtf.gz` to process the `long_noncoding_RNAs.gtf` file.

Here, we loads both "long_noncoding_RNAs.gtf" and "basic.annotation.gtf" which builds a dataframe of combined annotations for both lncRNAs and mRNAs. You can specify different annotation files options from GENCODE by modifying the `file_resources` dict argument.

.. code:: python

   from openomics.database import GENCODE, EnsemblGenes

   gencode = ENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                     file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                     "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                                     "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                     "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                     remove_version_num=True)

   # We also loads Ensembl genes to get list of miRNA gene IDs
   ensembl = EnsemblGenes(biomart='hsapiens_gene_ensembl', npartitions=8, )

Load single omics expressions for MessengerRNA, MicroRNA, LncRNA
----------------------------------------------------------------
We instantiate the MessengerRNA, MicroRNA and LncRNA -omics expression data from `gtex.data`. Since the gene expression were not seperated by RNA type, we use GENCODE and Ensembl gene annotations to filter the list of mRNA, miRNA, and lncRNAs.

.. code:: python

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

Create a MultiOmics dataset
---------------------------
Now, we create a MultiOmics dataset object by combining the messengerRNA, microRNA, and lncRNA.

.. code:: python

   from openomics import MultiOmics

   gtex_data = MultiOmics(cohort_name="GTEx Tissue Avg Expressions")

   gtex_data.add_omic(messengerRNA)
   gtex_data.add_omic(microRNA)
   gtex_data.add_omic(lncRNA)

   gtex_data.build_samples()
