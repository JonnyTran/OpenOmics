# External annotation databases

## Import GENCODE human release 32

Next, we can annotate the genes in our GTEx expression dataset with genomics annotation from GENCODE. In this example,
we use the URL path prefix "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/" which specifies the
species and release version. We also pass a dictionary `file_resources`, with key-value pairs where the key is name of
file and value is the suffix of the file download URL.

For example, file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz"} will download file
located at <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.long_noncoding_RNAs.gtf.gz>
to process the `long_noncoding_RNAs.gtf` file.

Here, we loads both "long_noncoding_RNAs.gtf" and "basic.annotation.gtf" which builds a dataframe of combined
annotations for both lncRNAs and mRNAs. You can specify different annotation files options from GENCODE by modifying
the `file_resources` dict argument.

```{code-block} python
from openomics.database import GENCODE, EnsemblGenes

gencode = GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                 file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                 "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                                 "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                 "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                 remove_version_num=True)

# We also loads Ensembl genes to get list of miRNA gene IDs
ensembl = EnsemblGenes(biomart='hsapiens_gene_ensembl', npartitions=8, )
```

## Setting the cache download directory
The package `astropy` is used to automatically cache downloaded files. It defaults to saving the files at
`~/.astropy/cache/`, where the cached content is retrieved given the matching URL. To change the path for the cache download file, run:

```python
import openomics

openomics.set_cache_dir(path="PATH/OF/YOUR/CHOICE/")
```

```{note}
Note that this setting doesn't persist across different programming sessions. Ideally, the cache dir should be in one location to minimize automatic FTP downloads, which may cause unnecessary stress on the database server.
```

