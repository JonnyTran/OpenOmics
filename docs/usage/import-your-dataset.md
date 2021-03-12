# Loading a multi-omics dataset

Suppose you have your own -omics dataset(s) and you'd like to load them. One of OpenOmics's primary goal is to
encapsulate the data import process with one line of code along with a few parameters. Given any processed single-omic
dataset, the library loads the data as a tabular structure where rows correspond to observation samples and columns
correspond to measurements of different biomolecules.

Import TCGA LUAD data included in tests dataset (preprocessed from TCGA-Assembler). It is located at [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD).

```{code-block} python
folder_path = "tests/data/TCGA_LUAD/"
```

Load the multiomics: Gene Expression, MicroRNA expression lncRNA expression, Copy Number Variation, Somatic Mutation, DNA Methylation, and Protein Expression data

```{code-block} python
from openomics import MessengerRNA, MicroRNA, LncRNA, SomaticMutation, Protein

# Load each expression dataframe
mRNA = MessengerRNA(data=folder_path + "LUAD__geneExp.txt",
                    transpose=True, usecols="GeneSymbol|TCGA", gene_index="GeneSymbol", gene_level="gene_name")
miRNA = MicroRNA(data=folder_path + "LUAD__miRNAExp__RPM.txt",
                 transpose=True, usecols="GeneSymbol|TCGA", gene_index="GeneSymbol", gene_level="transcript_name")
lncRNA = LncRNA(data=folder_path + "TCGA-rnaexpr.tsv",
                transpose=True, usecols="Gene_ID|TCGA", gene_index="Gene_ID", gene_level="gene_id")
som = SomaticMutation(data=folder_path + "LUAD__somaticMutation_geneLevel.txt",
                      transpose=True, usecols="GeneSymbol|TCGA", gene_index="gene_name")
pro = Protein(data=folder_path + "protein_RPPA.txt",
              transpose=True, usecols="GeneSymbol|TCGA", gene_index="GeneSymbol", gene_level="protein_name")

# Create an integrated MultiOmics dataset
luad_data = MultiOmics(cohort_name="LUAD")
luad_data.add_clinical_data(
    clinical=folder_path + "nationwidechildrens.org_clinical_patient_luad.txt")

luad_data.add_omic(mRNA)
luad_data.add_omic(miRNA)
luad_data.add_omic(lncRNA)
luad_data.add_omic(som)
luad_data.add_omic(pro)

luad_data.build_samples()
```

Each data is stored as a Pandas DataFrame. Below are all the data imported for TCGA LUAD. For each, the first number represents the number of samples, the second number is the number of features.

    PATIENTS (522, 5)
    SAMPLES (1160, 6)
    DRUGS (461, 4)
    MessengerRNA (576, 20472)
    SomaticMutation (587, 21070)
    MicroRNA (494, 1870)
    LncRNA (546, 12727)
    Protein (364, 154)


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

## Accessing clinical data
Each multi-omics and clinical data can be accessed through luad_data.data[], like:

```{code-block} python
luad_data.data["PATIENTS"]
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gender</th>
      <th>race</th>
      <th>histologic_subtype</th>
      <th>pathologic_stage</th>
    </tr>
    <tr>
      <th>bcr_patient_barcode</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-05-4245</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4249</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4250</th>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4382</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>522 rows × 4 columns</p>
</div>


```{code-block} python
luad_data.data["MessengerRNA"]
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>gene_name</th>
      <th>A1BG</th>
      <th>A1BG-AS1</th>
      <th>A1CF</th>
      <th>A2M</th>
      <th>A2ML1</th>
      <th>A4GALT</th>
      <th>A4GNT</th>
      <th>AAAS</th>
      <th>AACS</th>
      <th>AACSP1</th>
      <th>...</th>
      <th>ZXDA</th>
      <th>ZXDB</th>
      <th>ZXDC</th>
      <th>ZYG11A</th>
      <th>ZYG11B</th>
      <th>ZYX</th>
      <th>ZZEF1</th>
      <th>ZZZ3</th>
      <th>psiTPTE22</th>
      <th>tAKR</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244-01A</th>
      <td>4.756500</td>
      <td>5.239211</td>
      <td>0.000000</td>
      <td>13.265291</td>
      <td>0.431997</td>
      <td>7.043317</td>
      <td>1.033652</td>
      <td>9.348765</td>
      <td>9.652057</td>
      <td>0.763921</td>
      <td>...</td>
      <td>5.350285</td>
      <td>8.197321</td>
      <td>9.907260</td>
      <td>0.763921</td>
      <td>10.088859</td>
      <td>11.471139</td>
      <td>9.768648</td>
      <td>9.170597</td>
      <td>2.932118</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4249-01A</th>
      <td>6.920471</td>
      <td>7.056843</td>
      <td>0.402722</td>
      <td>14.650247</td>
      <td>1.383939</td>
      <td>9.178805</td>
      <td>0.717123</td>
      <td>9.241537</td>
      <td>9.967223</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.980428</td>
      <td>8.950001</td>
      <td>10.204971</td>
      <td>4.411650</td>
      <td>9.622978</td>
      <td>11.199826</td>
      <td>10.153700</td>
      <td>9.433116</td>
      <td>7.499637</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4250-01A</th>
      <td>5.696542</td>
      <td>6.136327</td>
      <td>0.000000</td>
      <td>14.048541</td>
      <td>0.000000</td>
      <td>8.481646</td>
      <td>0.996244</td>
      <td>9.203535</td>
      <td>9.560412</td>
      <td>0.733962</td>
      <td>...</td>
      <td>5.931168</td>
      <td>8.517334</td>
      <td>9.722642</td>
      <td>4.782796</td>
      <td>8.895339</td>
      <td>12.408981</td>
      <td>10.194168</td>
      <td>9.060342</td>
      <td>2.867956</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4382-01A</th>
      <td>7.198727</td>
      <td>6.809804</td>
      <td>0.000000</td>
      <td>14.509730</td>
      <td>2.532591</td>
      <td>9.117559</td>
      <td>1.657045</td>
      <td>9.251035</td>
      <td>10.078124</td>
      <td>1.860883</td>
      <td>...</td>
      <td>5.373036</td>
      <td>8.441914</td>
      <td>9.888267</td>
      <td>6.041142</td>
      <td>9.828389</td>
      <td>12.725186</td>
      <td>10.192589</td>
      <td>9.376841</td>
      <td>5.177029</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
<p>576 rows × 20472 columns</p>
</div>

## To match samples accross different multi-omics, use
```{code-block} python
luad_data.match_samples(modalities=["MicroRNA", "MessengerRNA"])
```

    Index(['TCGA-05-4384-01A', 'TCGA-05-4390-01A', 'TCGA-05-4396-01A',
           'TCGA-05-4405-01A', 'TCGA-05-4410-01A', 'TCGA-05-4415-01A',
           'TCGA-05-4417-01A', 'TCGA-05-4424-01A', 'TCGA-05-4425-01A',
           'TCGA-05-4427-01A',
           ...
           'TCGA-NJ-A4YG-01A', 'TCGA-NJ-A4YI-01A', 'TCGA-NJ-A4YP-01A',
           'TCGA-NJ-A4YQ-01A', 'TCGA-NJ-A55A-01A', 'TCGA-NJ-A55O-01A',
           'TCGA-NJ-A55R-01A', 'TCGA-NJ-A7XG-01A', 'TCGA-O1-A52J-01A',
           'TCGA-S2-AA1A-01A'],
          dtype='object', length=465)

