<p align="center">
  <img src="https://github.com/BioMeCIS-Lab/OpenOmics/raw/master/openomics_web/assets/openomics_logo.png" max-height="200">
</p>

[![PyPI version](https://badge.fury.io/py/openomics.svg)](https://badge.fury.io/py/openomics)
[![Documentation Status](https://readthedocs.org/projects/openomics/badge/?version=latest)](https://openomics.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/125549505.svg)](https://zenodo.org/badge/latestdoi/125549505)
[![OpenOmics](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics/branch/master/graph/badge.svg?token=WAN3PJwM42)](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics)

This Python package provide a series of tools to integrate and query the genomics, transcriptomics, proteomics, and
clinical data (aka multi-omics data). With scalable data-frame manipulation tools, OpenOmics facilitates the common
coding tasks when preparing data for RNA-seq bioinformatics analysis.

Documentation ([Latest](https://openomics.readthedocs.io/en/latest/)
| [Stable](https://openomics.readthedocs.io/en/stable/))
| [OpenOmics at a glance](https://openomics.readthedocs.io/en/stable/usage.html)

## Features
OpenOmics assist in integration of heterogeneous multi-omics bioinformatics data. The library provides a Python API as well as an interactive Dash web interface.
It features support for:
- Genomics, Transcriptomics, Proteomics, and Clinical data.
- Harmonization with 20+ popular annotation, interaction, disease-association databases.

OpenOmics also has an efficient data pipeline that bridges the popular data manipulation Pandas library and Dask distributed processing to address the following use cases:

- Providing a standard pipeline for dataset indexing, table joining and querying, which are transparent and customizable
  for end-users.
- Providing Efficient disk storage for large multi-omics dataset with Parquet data structures.
- Integrating various data types including interactions and sequence data, then exporting to NetworkX graphs or data
  generators for down-stream machine learning.
- Accessible by both developers and scientists with a Python API that works seamlessly with external Galaxy tool
  interface or the built-in Dash web interface (WIP).

<br/>

## Installation via pip:

```
$ pip install openomics
```

<br/>

# How to use OpenOmics:


## Importing the openomics library

```python
from openomics import MultiOmics
```

Import TCGA LUAD data included in tests dataset (preprocessed from TCGA-Assembler). It is located at [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD).

```python
folder_path = "tests/data/TCGA_LUAD/"
```

Load the multiomics: Gene Expression, MicroRNA expression lncRNA expression, Copy Number Variation, Somatic Mutation, DNA Methylation, and Protein Expression data

```python
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

## Annotate LncRNAs with GENCODE genomic annotations

```python
# Import GENCODE database (from URL)
from openomics.database import GENCODE

gencode = GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                  file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                  "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                                  "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                  "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                  remove_version_num=True,
                  npartitions=5)

# Annotate LncRNAs with GENCODE by gene_id
luad_data.LncRNA.annotate_attributes(gencode, on="gene_id",
                                     columns=['feature', 'start', 'end', 'strand', 'tag', 'havana_gene'])

luad_data.LncRNA.annotations.info()
```
    <class 'pandas.core.frame.DataFrame'>
    Index: 13729 entries, ENSG00000082929 to ENSG00000284600
    Data columns (total 6 columns):
     #   Column       Non-Null Count  Dtype
    ---  ------       --------------  -----
     0   feature      13729 non-null  object
     1   start        13729 non-null  object
     2   end          13729 non-null  object
     3   strand       13729 non-null  object
     4   tag          13729 non-null  object
     5   havana_gene  13729 non-null  object
    dtypes: object(6)
    memory usage: 1.4+ MB

Each multi-omics and clinical data can be accessed through luad_data.data[], like:


```python
luad_data.data["PATIENTS"]
```

<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>bcr_patient_barcode</th>
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
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244</th>
      <td>TCGA-05-4244</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-05-4245</th>
      <td>TCGA-05-4245</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4249</th>
      <td>TCGA-05-4249</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4250</th>
      <td>TCGA-05-4250</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4382</th>
      <td>TCGA-05-4382</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>522 rows × 5 columns</p>
</div>




```python
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


```python
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


## To prepare the data for classification


```python
# This function selects only patients with patholotic stages "Stage I" and "Stage II"
X_multiomics, y = luad_data.load_dataframe(modalities=["MessengerRNA", "MicroRNA", "LncRNA"], target=['pathologic_stage'],
                                     pathologic_stages=['Stage I', 'Stage II'])
print(X_multiomics['MessengerRNA'].shape, X_multiomics['MicroRNA'].shape, X_multiomics['LncRNA'].shape, y.shape)
```

    (336, 20472) (336, 1870) (336, 12727) (336, 1)


```python
y
```


<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>pathologic_stage</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4390-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4405-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4410-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4417-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4424-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4427-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4433-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-5423-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-5425-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-5428-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-5715-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-38-4631-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-38-7271-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-38-A44F-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-2655-11A</th>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>336 rows × 1 columns</p>
</div>



## Log2 transform the mRNA, microRNA, and lncRNA expression values


```python
def expression_val_transform(x):
    return np.log2(x+1)
X_multiomics['MessengerRNA'] = X_multiomics['MessengerRNA'].applymap(expression_val_transform)
X_multiomics['MicroRNA'] = X_multiomics['MicroRNA'].applymap(expression_val_transform)
# X_multiomics['LncRNA'] = X_multiomics['LncRNA'].applymap(expression_val_transform)
```

## Classification of Cancer Stage


```python
from sklearn import preprocessing
from sklearn import metrics
from sklearn.svm import SVC, LinearSVC
import sklearn.linear_model
from sklearn.model_selection import train_test_split

```


```python
binarizer = preprocessing.LabelEncoder()
binarizer.fit(y)
binarizer.transform(y)
```


    array([0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1,
           0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1,
           0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0,
           0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0,
           0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1,
           1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
           1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
           0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0,
           1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0])




```python
for omic in ["MessengerRNA", "MicroRNA"]:
    print(omic)
    scaler = sklearn.preprocessing.StandardScaler(copy=True, with_mean=True, with_std=False)
    scaler.fit(X_multiomics[omic])

    X_train, X_test, Y_train, Y_test = \
        train_test_split(X_multiomics[omic], y, test_size=0.3, random_state=np.random.randint(0, 10000), stratify=y)
    print(X_train.shape, X_test.shape)


    X_train = scaler.transform(X_train)

    model = LinearSVC(C=1e-2, penalty='l1', class_weight='balanced', dual=False, multi_class="ovr")
#     model = sklearn.linear_model.LogisticRegression(C=1e-0, penalty='l1', fit_intercept=False, class_weight="balanced")
#     model = SVC(C=1e0, kernel="rbf", class_weight="balanced", decision_function_shape="ovo")

    model.fit(X=X_train, y=Y_train)
    print("NONZERO", len(np.nonzero(model.coef_)[0]))
    print("Training accuracy", metrics.accuracy_score(model.predict(X_train), Y_train))
    print(metrics.classification_report(y_pred=model.predict(X_test), y_true=Y_test))

```

    MessengerRNA
    (254, 20472) (109, 20472)
    NONZERO 0
    Training accuracy 0.6929133858267716
                 precision    recall  f1-score   support

        Stage I       0.69      1.00      0.82        75
       Stage II       0.00      0.00      0.00        34

    avg / total       0.47      0.69      0.56       109

    MicroRNA
    (254, 1870) (109, 1870)
    NONZERO 0
    Training accuracy 0.6929133858267716
                 precision    recall  f1-score   support

        Stage I       0.69      1.00      0.82        75
       Stage II       0.00      0.00      0.00        34

    avg / total       0.47      0.69      0.56       109


Credits
-------

This package was created with Cookiecutter_ and the `pyOpenSci/cookiecutter-pyopensci`_ project template, based off `audreyr/cookiecutter-pypackage`_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`pyOpenSci/cookiecutter-pyopensci`: https://github.com/pyOpenSci/cookiecutter-pyopensci
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
