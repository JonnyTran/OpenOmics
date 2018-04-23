# TCGAMultiOmics

This Python package provide a series of tool to download, integrate, and query multi-omics TCGA data. The TCGA data is downloaded from TCGA-Assembler.

Load all multi-omics TCGA data from a given folder_path with the following folder structure:

    folder_path/
        clinical/
            genome.wustl.edu_biospecimen_sample.txt (optional)
            nationwidechildrens.org_clinical_drug.txt
            nationwidechildrens.org_clinical_patient.txt
        gene_exp/
            geneExp.txt
        mirna/
            miRNAExp__RPM.txt
        cnv/
            copyNumber.txt
        protein_rppa/
            protein_RPPA.txt
        somatic/
            somaticMutation_geneLevel.txt
        lncrna/
            TCGA-rnaexpr.tsv
            

            

# Install with "pip install git+https://github.com/JonnyTran/TCGAMultiOmics"

# Importing the TCGAMultiOmics library


```python
from TCGAMultiOmics.multiomics import MultiOmicsData

import pandas as pd
import numpy as np
```

# Import TCGA LUAD data downloaded from TCGA-Assembler


```python
folder_path ="./data/tcga-assembler/LUAD/"
```




    './data/tcga-assembler/LUAD/'




```python
# Load all modalities: Gene Expression, MicroRNA expression lncRNA expression, Copy Number Variation, Somatic Mutation, DNA Methylation, and Protein Expression data
luad_data = MultiOmicsData(cancer_type="LUAD", folder_path=folder_path,
                           modalities=["GE", "MIR", "LNC", "CNV", "SNP", "PRO"])

```

    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/pandas/core/frame.py:3027: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      return super(DataFrame, self).rename(**kwargs)
    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/TCGAMultiOmics/genomic.py:61: RuntimeWarning: invalid value encountered in log2
      return np.log2(x + 1)
    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/TCGAMultiOmics/genomic.py:61: RuntimeWarning: divide by zero encountered in log2
      return np.log2(x + 1)


    PATIENTS (522, 5)
    DRUGS (461, 4)
    GE (576, 20472)
    SNP (587, 21070)
    MIR (494, 1870)
    LNC (546, 12727)
    CNV (1107, 22516)
    PRO (364, 154)
    SAMPLES (1160, 6)


# Each multi-omics and clinical data can be accessed through luad_data.data[""]


```python
luad_data.data["PATIENTS"]
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
    <tr>
      <th>TCGA-05-4384</th>
      <td>TCGA-05-4384</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4389</th>
      <td>TCGA-05-4389</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4390</th>
      <td>TCGA-05-4390</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4395</th>
      <td>TCGA-05-4395</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4396</th>
      <td>TCGA-05-4396</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4397</th>
      <td>TCGA-05-4397</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4398</th>
      <td>TCGA-05-4398</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4402</th>
      <td>TCGA-05-4402</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-05-4403</th>
      <td>TCGA-05-4403</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4405</th>
      <td>TCGA-05-4405</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4410</th>
      <td>TCGA-05-4410</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4415</th>
      <td>TCGA-05-4415</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4417</th>
      <td>TCGA-05-4417</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4418</th>
      <td>TCGA-05-4418</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4420</th>
      <td>TCGA-05-4420</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4422</th>
      <td>TCGA-05-4422</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4424</th>
      <td>TCGA-05-4424</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4425</th>
      <td>TCGA-05-4425</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-05-4426</th>
      <td>TCGA-05-4426</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4427</th>
      <td>TCGA-05-4427</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4430</th>
      <td>TCGA-05-4430</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4432</th>
      <td>TCGA-05-4432</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4433</th>
      <td>TCGA-05-4433</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4434</th>
      <td>TCGA-05-4434</td>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-05-5420</th>
      <td>TCGA-05-5420</td>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SV</th>
      <td>TCGA-MP-A4SV</td>
      <td>MALE</td>
      <td>[Unknown]</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SW</th>
      <td>TCGA-MP-A4SW</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SY</th>
      <td>TCGA-MP-A4SY</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T2</th>
      <td>TCGA-MP-A4T2</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T4</th>
      <td>TCGA-MP-A4T4</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T6</th>
      <td>TCGA-MP-A4T6</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T7</th>
      <td>TCGA-MP-A4T7</td>
      <td>FEMALE</td>
      <td>[Unknown]</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T8</th>
      <td>TCGA-MP-A4T8</td>
      <td>MALE</td>
      <td>[Unknown]</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T9</th>
      <td>TCGA-MP-A4T9</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TA</th>
      <td>TCGA-MP-A4TA</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TC</th>
      <td>TCGA-MP-A4TC</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TD</th>
      <td>TCGA-MP-A4TD</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TE</th>
      <td>TCGA-MP-A4TE</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TF</th>
      <td>TCGA-MP-A4TF</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TH</th>
      <td>TCGA-MP-A4TH</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TI</th>
      <td>TCGA-MP-A4TI</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TJ</th>
      <td>TCGA-MP-A4TJ</td>
      <td>FEMALE</td>
      <td>[Unknown]</td>
      <td>Lung Signet Ring Adenocarcinoma</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TK</th>
      <td>TCGA-MP-A4TK</td>
      <td>FEMALE</td>
      <td>[Unknown]</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A5C7</th>
      <td>TCGA-MP-A5C7</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YF</th>
      <td>TCGA-NJ-A4YF</td>
      <td>FEMALE</td>
      <td>BLACK OR AFRICAN AMERICAN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YG</th>
      <td>TCGA-NJ-A4YG</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YI</th>
      <td>TCGA-NJ-A4YI</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Papillary Adenocarcinoma</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YP</th>
      <td>TCGA-NJ-A4YP</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Papillary Adenocarcinoma</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YQ</th>
      <td>TCGA-NJ-A4YQ</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55A</th>
      <td>TCGA-NJ-A55A</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55O</th>
      <td>TCGA-NJ-A55O</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Mucinous (Colloid) Carcinoma</td>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55R</th>
      <td>TCGA-NJ-A55R</td>
      <td>MALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A7XG</th>
      <td>TCGA-NJ-A7XG</td>
      <td>MALE</td>
      <td>BLACK OR AFRICAN AMERICAN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-O1-A52J</th>
      <td>TCGA-O1-A52J</td>
      <td>FEMALE</td>
      <td>WHITE</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-S2-AA1A</th>
      <td>TCGA-S2-AA1A</td>
      <td>FEMALE</td>
      <td>BLACK OR AFRICAN AMERICAN</td>
      <td>Lung Bronchioloalveolar Carcinoma Mucinous</td>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>522 rows × 5 columns</p>
</div>




```python
luad_data.data["GE"]
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
    <tr>
      <th>TCGA-05-4384-01A</th>
      <td>7.004132</td>
      <td>7.126555</td>
      <td>0.000000</td>
      <td>15.560200</td>
      <td>0.000000</td>
      <td>7.264146</td>
      <td>1.502127</td>
      <td>9.381223</td>
      <td>9.473371</td>
      <td>2.221939</td>
      <td>...</td>
      <td>6.340285</td>
      <td>9.140128</td>
      <td>10.368188</td>
      <td>3.160501</td>
      <td>9.607078</td>
      <td>11.706702</td>
      <td>10.763482</td>
      <td>9.500392</td>
      <td>8.566551</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4389-01A</th>
      <td>6.090449</td>
      <td>5.419832</td>
      <td>5.218843</td>
      <td>12.929239</td>
      <td>5.231786</td>
      <td>8.645042</td>
      <td>1.743730</td>
      <td>9.484694</td>
      <td>10.458104</td>
      <td>0.417488</td>
      <td>...</td>
      <td>6.331458</td>
      <td>9.291179</td>
      <td>9.983577</td>
      <td>6.186106</td>
      <td>9.065818</td>
      <td>11.325551</td>
      <td>9.409804</td>
      <td>9.558459</td>
      <td>3.593808</td>
      <td>0.740798</td>
    </tr>
    <tr>
      <th>TCGA-05-4390-01A</th>
      <td>7.372546</td>
      <td>5.777325</td>
      <td>0.000000</td>
      <td>13.788390</td>
      <td>0.712552</td>
      <td>9.363766</td>
      <td>0.712552</td>
      <td>9.472821</td>
      <td>10.614398</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.630493</td>
      <td>9.480580</td>
      <td>9.618252</td>
      <td>3.115050</td>
      <td>9.274341</td>
      <td>12.360097</td>
      <td>9.610022</td>
      <td>8.942920</td>
      <td>2.754460</td>
      <td>0.712552</td>
    </tr>
    <tr>
      <th>TCGA-05-4395-01A</th>
      <td>4.521962</td>
      <td>3.933412</td>
      <td>0.000000</td>
      <td>12.171309</td>
      <td>8.614152</td>
      <td>9.284894</td>
      <td>0.620023</td>
      <td>9.815079</td>
      <td>9.507134</td>
      <td>0.000000</td>
      <td>...</td>
      <td>4.096464</td>
      <td>7.625694</td>
      <td>9.544953</td>
      <td>3.928588</td>
      <td>8.052096</td>
      <td>12.273278</td>
      <td>8.871455</td>
      <td>9.643634</td>
      <td>2.077687</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4396-01A</th>
      <td>4.204571</td>
      <td>4.362323</td>
      <td>4.080760</td>
      <td>14.593816</td>
      <td>0.814099</td>
      <td>7.726325</td>
      <td>0.814099</td>
      <td>9.221360</td>
      <td>9.441835</td>
      <td>0.000000</td>
      <td>...</td>
      <td>6.399111</td>
      <td>9.298077</td>
      <td>10.635574</td>
      <td>5.102599</td>
      <td>9.317065</td>
      <td>10.713854</td>
      <td>10.527088</td>
      <td>9.789539</td>
      <td>4.014623</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4397-01A</th>
      <td>6.995131</td>
      <td>6.579800</td>
      <td>0.000000</td>
      <td>12.310564</td>
      <td>4.337419</td>
      <td>9.721597</td>
      <td>2.655558</td>
      <td>9.780029</td>
      <td>10.648534</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.879426</td>
      <td>8.463813</td>
      <td>10.214013</td>
      <td>7.981822</td>
      <td>9.690311</td>
      <td>11.913838</td>
      <td>8.954124</td>
      <td>9.620376</td>
      <td>5.045552</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4398-01A</th>
      <td>8.017215</td>
      <td>7.889831</td>
      <td>0.283922</td>
      <td>13.647286</td>
      <td>1.204704</td>
      <td>9.307084</td>
      <td>0.521051</td>
      <td>9.424975</td>
      <td>9.577727</td>
      <td>1.204704</td>
      <td>...</td>
      <td>4.268195</td>
      <td>7.855086</td>
      <td>9.927825</td>
      <td>6.909081</td>
      <td>9.453897</td>
      <td>12.566950</td>
      <td>9.216107</td>
      <td>9.324318</td>
      <td>5.475643</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4402-01A</th>
      <td>5.628158</td>
      <td>5.506405</td>
      <td>0.000000</td>
      <td>14.315344</td>
      <td>4.705536</td>
      <td>9.778099</td>
      <td>1.291368</td>
      <td>9.498088</td>
      <td>10.254248</td>
      <td>1.141629</td>
      <td>...</td>
      <td>6.383946</td>
      <td>8.658785</td>
      <td>9.869494</td>
      <td>5.938222</td>
      <td>10.209383</td>
      <td>11.783793</td>
      <td>9.524795</td>
      <td>10.243664</td>
      <td>3.596744</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4403-01A</th>
      <td>6.026079</td>
      <td>5.284414</td>
      <td>1.558806</td>
      <td>14.989447</td>
      <td>0.000000</td>
      <td>8.432349</td>
      <td>1.257011</td>
      <td>9.166945</td>
      <td>10.016462</td>
      <td>0.000000</td>
      <td>...</td>
      <td>6.264848</td>
      <td>9.307769</td>
      <td>9.481655</td>
      <td>0.353888</td>
      <td>8.882484</td>
      <td>12.469388</td>
      <td>9.731860</td>
      <td>9.115087</td>
      <td>5.309005</td>
      <td>0.637842</td>
    </tr>
    <tr>
      <th>TCGA-05-4405-01A</th>
      <td>7.200949</td>
      <td>7.221463</td>
      <td>0.000000</td>
      <td>14.970769</td>
      <td>1.217417</td>
      <td>8.508278</td>
      <td>1.217417</td>
      <td>8.521346</td>
      <td>9.974718</td>
      <td>1.868094</td>
      <td>...</td>
      <td>5.226690</td>
      <td>8.101538</td>
      <td>9.955590</td>
      <td>6.282907</td>
      <td>9.940103</td>
      <td>11.826372</td>
      <td>10.514482</td>
      <td>9.411214</td>
      <td>3.361319</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4410-01A</th>
      <td>6.329190</td>
      <td>6.262367</td>
      <td>0.000000</td>
      <td>14.485855</td>
      <td>2.978324</td>
      <td>8.381841</td>
      <td>0.819259</td>
      <td>8.721967</td>
      <td>9.984749</td>
      <td>0.819259</td>
      <td>...</td>
      <td>5.808537</td>
      <td>8.191229</td>
      <td>10.427467</td>
      <td>6.882307</td>
      <td>10.383556</td>
      <td>11.873941</td>
      <td>10.493222</td>
      <td>9.351250</td>
      <td>4.330064</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4415-01A</th>
      <td>5.041865</td>
      <td>4.655352</td>
      <td>0.000000</td>
      <td>11.906874</td>
      <td>1.222372</td>
      <td>8.106781</td>
      <td>1.222372</td>
      <td>10.101100</td>
      <td>9.253453</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.751916</td>
      <td>8.415507</td>
      <td>9.872589</td>
      <td>7.825136</td>
      <td>9.397318</td>
      <td>12.788650</td>
      <td>9.884171</td>
      <td>9.580991</td>
      <td>3.650259</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4417-01A</th>
      <td>6.431387</td>
      <td>6.285079</td>
      <td>0.000000</td>
      <td>15.059303</td>
      <td>2.062433</td>
      <td>9.181066</td>
      <td>0.000000</td>
      <td>9.029690</td>
      <td>10.252259</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.997800</td>
      <td>8.868415</td>
      <td>9.951699</td>
      <td>4.791080</td>
      <td>9.678367</td>
      <td>11.785149</td>
      <td>10.105222</td>
      <td>8.802194</td>
      <td>4.855825</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4418-01A</th>
      <td>7.341244</td>
      <td>6.534743</td>
      <td>0.414027</td>
      <td>12.530094</td>
      <td>4.837247</td>
      <td>8.691617</td>
      <td>0.735262</td>
      <td>10.304520</td>
      <td>10.068605</td>
      <td>1.582171</td>
      <td>...</td>
      <td>4.803305</td>
      <td>7.986392</td>
      <td>9.464089</td>
      <td>5.399619</td>
      <td>8.935386</td>
      <td>12.145546</td>
      <td>9.710621</td>
      <td>9.289660</td>
      <td>7.032920</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4420-01A</th>
      <td>8.662965</td>
      <td>6.721158</td>
      <td>0.000000</td>
      <td>13.582962</td>
      <td>1.633943</td>
      <td>8.320556</td>
      <td>0.000000</td>
      <td>9.309954</td>
      <td>10.436395</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.693470</td>
      <td>8.208319</td>
      <td>9.987547</td>
      <td>6.803512</td>
      <td>9.836823</td>
      <td>10.527368</td>
      <td>8.912698</td>
      <td>10.179033</td>
      <td>2.380508</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4422-01A</th>
      <td>7.023161</td>
      <td>6.166229</td>
      <td>5.786223</td>
      <td>13.367342</td>
      <td>3.296898</td>
      <td>7.206795</td>
      <td>1.646577</td>
      <td>9.570829</td>
      <td>9.823888</td>
      <td>0.000000</td>
      <td>...</td>
      <td>6.245758</td>
      <td>8.693800</td>
      <td>10.493342</td>
      <td>3.054137</td>
      <td>8.732517</td>
      <td>11.023314</td>
      <td>10.642012</td>
      <td>8.952358</td>
      <td>4.098782</td>
      <td>2.015819</td>
    </tr>
    <tr>
      <th>TCGA-05-4424-01A</th>
      <td>8.016638</td>
      <td>7.252499</td>
      <td>0.000000</td>
      <td>15.789319</td>
      <td>1.539333</td>
      <td>8.842524</td>
      <td>1.280897</td>
      <td>8.858908</td>
      <td>10.491365</td>
      <td>0.000000</td>
      <td>...</td>
      <td>6.366431</td>
      <td>9.576732</td>
      <td>9.714064</td>
      <td>2.748397</td>
      <td>9.740832</td>
      <td>11.950419</td>
      <td>10.515524</td>
      <td>9.554956</td>
      <td>9.618447</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4425-01A</th>
      <td>5.588340</td>
      <td>5.643484</td>
      <td>0.000000</td>
      <td>13.022930</td>
      <td>1.083520</td>
      <td>9.166382</td>
      <td>0.000000</td>
      <td>9.126494</td>
      <td>10.201919</td>
      <td>0.641176</td>
      <td>...</td>
      <td>5.759201</td>
      <td>8.211729</td>
      <td>10.037278</td>
      <td>1.421587</td>
      <td>9.489774</td>
      <td>11.729891</td>
      <td>8.864264</td>
      <td>9.310666</td>
      <td>1.695281</td>
      <td>0.641176</td>
    </tr>
    <tr>
      <th>TCGA-05-4426-01A</th>
      <td>8.270363</td>
      <td>7.314247</td>
      <td>0.000000</td>
      <td>12.838338</td>
      <td>0.000000</td>
      <td>7.592988</td>
      <td>0.374511</td>
      <td>9.306595</td>
      <td>9.391180</td>
      <td>1.311503</td>
      <td>...</td>
      <td>6.130573</td>
      <td>8.547571</td>
      <td>10.124224</td>
      <td>0.917775</td>
      <td>10.049525</td>
      <td>11.743545</td>
      <td>9.308620</td>
      <td>9.793686</td>
      <td>4.423040</td>
      <td>0.374511</td>
    </tr>
    <tr>
      <th>TCGA-05-4427-01A</th>
      <td>6.607009</td>
      <td>7.421382</td>
      <td>0.000000</td>
      <td>13.580522</td>
      <td>0.589284</td>
      <td>8.212671</td>
      <td>0.589284</td>
      <td>8.424516</td>
      <td>9.796845</td>
      <td>0.000000</td>
      <td>...</td>
      <td>7.187879</td>
      <td>10.058108</td>
      <td>10.813127</td>
      <td>7.370516</td>
      <td>10.453903</td>
      <td>11.610626</td>
      <td>10.197287</td>
      <td>10.059473</td>
      <td>3.181532</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4430-01A</th>
      <td>8.661710</td>
      <td>7.073164</td>
      <td>0.000000</td>
      <td>15.226886</td>
      <td>0.399445</td>
      <td>9.502002</td>
      <td>1.375679</td>
      <td>8.971878</td>
      <td>9.522796</td>
      <td>0.968570</td>
      <td>...</td>
      <td>5.899471</td>
      <td>9.308734</td>
      <td>9.711907</td>
      <td>4.354572</td>
      <td>9.507071</td>
      <td>12.263314</td>
      <td>9.602505</td>
      <td>9.186902</td>
      <td>4.055716</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4432-01A</th>
      <td>7.098543</td>
      <td>6.144934</td>
      <td>0.899408</td>
      <td>13.877613</td>
      <td>0.657091</td>
      <td>8.732290</td>
      <td>0.657091</td>
      <td>9.446189</td>
      <td>10.531749</td>
      <td>1.106817</td>
      <td>...</td>
      <td>6.348250</td>
      <td>8.979012</td>
      <td>9.966395</td>
      <td>5.937408</td>
      <td>9.561920</td>
      <td>11.111787</td>
      <td>10.716945</td>
      <td>9.177473</td>
      <td>3.680977</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4433-01A</th>
      <td>4.856333</td>
      <td>4.918763</td>
      <td>3.211713</td>
      <td>15.086977</td>
      <td>2.574998</td>
      <td>7.765063</td>
      <td>5.374897</td>
      <td>8.912864</td>
      <td>8.973419</td>
      <td>0.869003</td>
      <td>...</td>
      <td>6.154776</td>
      <td>9.260677</td>
      <td>9.895181</td>
      <td>0.869003</td>
      <td>9.049646</td>
      <td>12.543218</td>
      <td>10.657363</td>
      <td>9.270363</td>
      <td>3.076901</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4434-01A</th>
      <td>6.676981</td>
      <td>7.281937</td>
      <td>0.294077</td>
      <td>14.448851</td>
      <td>9.456478</td>
      <td>9.846292</td>
      <td>0.929337</td>
      <td>9.371375</td>
      <td>9.402555</td>
      <td>2.465139</td>
      <td>...</td>
      <td>5.412805</td>
      <td>8.161518</td>
      <td>10.125686</td>
      <td>1.601601</td>
      <td>9.564241</td>
      <td>12.197218</td>
      <td>10.253564</td>
      <td>9.043825</td>
      <td>4.547388</td>
      <td>0.294077</td>
    </tr>
    <tr>
      <th>TCGA-05-5420-01A</th>
      <td>8.252258</td>
      <td>6.990833</td>
      <td>7.555048</td>
      <td>13.604667</td>
      <td>0.000000</td>
      <td>7.667290</td>
      <td>1.480265</td>
      <td>9.552725</td>
      <td>9.685670</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.367231</td>
      <td>8.473418</td>
      <td>8.477046</td>
      <td>5.749703</td>
      <td>9.228266</td>
      <td>11.280407</td>
      <td>8.193660</td>
      <td>9.507316</td>
      <td>4.932430</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-5423-01A</th>
      <td>9.473126</td>
      <td>7.560004</td>
      <td>0.000000</td>
      <td>13.436623</td>
      <td>0.000000</td>
      <td>6.111280</td>
      <td>2.615463</td>
      <td>8.918219</td>
      <td>8.911650</td>
      <td>0.792939</td>
      <td>...</td>
      <td>5.949673</td>
      <td>8.287044</td>
      <td>9.368666</td>
      <td>5.708245</td>
      <td>10.119402</td>
      <td>11.419291</td>
      <td>8.785593</td>
      <td>9.491214</td>
      <td>2.924746</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>TCGA-MN-A4N5-01A</th>
      <td>8.468247</td>
      <td>8.372505</td>
      <td>0.000000</td>
      <td>13.907189</td>
      <td>0.837863</td>
      <td>7.381526</td>
      <td>1.364460</td>
      <td>9.090363</td>
      <td>9.104877</td>
      <td>0.837863</td>
      <td>...</td>
      <td>5.107043</td>
      <td>8.269035</td>
      <td>10.004198</td>
      <td>5.716300</td>
      <td>9.125363</td>
      <td>12.234197</td>
      <td>10.087473</td>
      <td>9.455699</td>
      <td>3.329856</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SV-01A</th>
      <td>6.034878</td>
      <td>6.305482</td>
      <td>0.000000</td>
      <td>13.538944</td>
      <td>0.000000</td>
      <td>8.308098</td>
      <td>0.891497</td>
      <td>9.748049</td>
      <td>9.813945</td>
      <td>2.890058</td>
      <td>...</td>
      <td>5.755264</td>
      <td>8.512773</td>
      <td>10.454617</td>
      <td>7.648595</td>
      <td>9.845127</td>
      <td>11.939382</td>
      <td>10.320684</td>
      <td>9.814630</td>
      <td>6.899982</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SW-01A</th>
      <td>7.186085</td>
      <td>7.212903</td>
      <td>0.000000</td>
      <td>13.569303</td>
      <td>0.541416</td>
      <td>10.503647</td>
      <td>1.496462</td>
      <td>9.560521</td>
      <td>9.851814</td>
      <td>1.712332</td>
      <td>...</td>
      <td>4.455965</td>
      <td>7.300568</td>
      <td>9.811441</td>
      <td>6.820459</td>
      <td>10.271963</td>
      <td>11.945639</td>
      <td>10.175846</td>
      <td>8.544947</td>
      <td>4.956735</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SY-01A</th>
      <td>5.902602</td>
      <td>5.651545</td>
      <td>0.000000</td>
      <td>13.472570</td>
      <td>2.526620</td>
      <td>8.159199</td>
      <td>1.662023</td>
      <td>9.448007</td>
      <td>10.175531</td>
      <td>2.819873</td>
      <td>...</td>
      <td>5.155191</td>
      <td>8.123810</td>
      <td>9.938329</td>
      <td>5.510696</td>
      <td>9.963570</td>
      <td>12.695075</td>
      <td>9.577088</td>
      <td>9.692545</td>
      <td>4.657869</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T4-01A</th>
      <td>7.741933</td>
      <td>8.626581</td>
      <td>0.568810</td>
      <td>13.997212</td>
      <td>0.568810</td>
      <td>9.213242</td>
      <td>1.552525</td>
      <td>9.197889</td>
      <td>9.090706</td>
      <td>0.000000</td>
      <td>...</td>
      <td>4.345751</td>
      <td>7.116320</td>
      <td>9.938379</td>
      <td>7.696132</td>
      <td>9.493198</td>
      <td>13.061718</td>
      <td>10.023245</td>
      <td>9.370046</td>
      <td>4.906867</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T6-01A</th>
      <td>8.539588</td>
      <td>9.142275</td>
      <td>0.000000</td>
      <td>12.989049</td>
      <td>0.000000</td>
      <td>10.070757</td>
      <td>1.300944</td>
      <td>10.328567</td>
      <td>9.936575</td>
      <td>3.862025</td>
      <td>...</td>
      <td>6.348948</td>
      <td>8.037081</td>
      <td>10.279710</td>
      <td>8.382245</td>
      <td>9.429923</td>
      <td>10.948718</td>
      <td>9.993102</td>
      <td>8.196756</td>
      <td>5.218417</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T7-01A</th>
      <td>5.691891</td>
      <td>6.372331</td>
      <td>0.000000</td>
      <td>14.336930</td>
      <td>2.435922</td>
      <td>8.912333</td>
      <td>1.216175</td>
      <td>10.033705</td>
      <td>10.234050</td>
      <td>1.866592</td>
      <td>...</td>
      <td>4.954001</td>
      <td>8.007236</td>
      <td>10.188966</td>
      <td>1.680594</td>
      <td>8.658848</td>
      <td>12.298280</td>
      <td>10.572970</td>
      <td>8.507940</td>
      <td>6.011683</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T8-01A</th>
      <td>4.730846</td>
      <td>4.539606</td>
      <td>0.000000</td>
      <td>13.658237</td>
      <td>0.000000</td>
      <td>9.356888</td>
      <td>1.815452</td>
      <td>9.434240</td>
      <td>10.304970</td>
      <td>0.879627</td>
      <td>...</td>
      <td>4.756847</td>
      <td>8.512126</td>
      <td>9.417504</td>
      <td>6.589885</td>
      <td>9.523382</td>
      <td>12.379501</td>
      <td>9.696927</td>
      <td>9.308978</td>
      <td>2.124196</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T9-01A</th>
      <td>5.826668</td>
      <td>6.257608</td>
      <td>0.000000</td>
      <td>14.917946</td>
      <td>1.580724</td>
      <td>9.138845</td>
      <td>1.134024</td>
      <td>9.523886</td>
      <td>9.219125</td>
      <td>1.134024</td>
      <td>...</td>
      <td>5.171759</td>
      <td>7.999315</td>
      <td>10.045968</td>
      <td>0.483571</td>
      <td>9.839477</td>
      <td>12.412521</td>
      <td>10.579947</td>
      <td>9.209452</td>
      <td>6.634698</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TA-01A</th>
      <td>6.470221</td>
      <td>5.506605</td>
      <td>4.261350</td>
      <td>12.965836</td>
      <td>2.444508</td>
      <td>9.826329</td>
      <td>1.145221</td>
      <td>10.228232</td>
      <td>9.424810</td>
      <td>4.634181</td>
      <td>...</td>
      <td>5.399260</td>
      <td>8.270566</td>
      <td>10.121895</td>
      <td>7.523469</td>
      <td>9.295379</td>
      <td>12.597671</td>
      <td>9.460817</td>
      <td>9.340125</td>
      <td>2.644156</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TC-01A</th>
      <td>8.867801</td>
      <td>8.701134</td>
      <td>0.000000</td>
      <td>13.510088</td>
      <td>0.000000</td>
      <td>8.915287</td>
      <td>2.472384</td>
      <td>9.163243</td>
      <td>9.856887</td>
      <td>0.541019</td>
      <td>...</td>
      <td>5.451745</td>
      <td>8.590364</td>
      <td>9.321393</td>
      <td>7.684620</td>
      <td>9.935404</td>
      <td>12.555356</td>
      <td>10.183592</td>
      <td>9.637111</td>
      <td>3.959650</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TD-01A</th>
      <td>6.722376</td>
      <td>6.099615</td>
      <td>0.578118</td>
      <td>15.296500</td>
      <td>2.789562</td>
      <td>8.732496</td>
      <td>2.789562</td>
      <td>9.729125</td>
      <td>10.058735</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.471548</td>
      <td>7.889211</td>
      <td>9.996060</td>
      <td>3.151485</td>
      <td>8.977862</td>
      <td>12.245876</td>
      <td>10.139751</td>
      <td>8.330234</td>
      <td>3.151485</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TE-01A</th>
      <td>9.056166</td>
      <td>7.362117</td>
      <td>0.000000</td>
      <td>12.008847</td>
      <td>2.177663</td>
      <td>8.483123</td>
      <td>0.526570</td>
      <td>9.967226</td>
      <td>11.622904</td>
      <td>0.526570</td>
      <td>...</td>
      <td>5.939048</td>
      <td>8.636321</td>
      <td>9.848850</td>
      <td>6.228988</td>
      <td>8.738063</td>
      <td>10.526508</td>
      <td>9.591494</td>
      <td>8.495503</td>
      <td>1.679244</td>
      <td>0.526570</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TF-01A</th>
      <td>5.236672</td>
      <td>5.808406</td>
      <td>4.180219</td>
      <td>12.005311</td>
      <td>0.860526</td>
      <td>7.987860</td>
      <td>0.000000</td>
      <td>9.355703</td>
      <td>9.350305</td>
      <td>1.785341</td>
      <td>...</td>
      <td>4.962521</td>
      <td>8.542805</td>
      <td>10.723703</td>
      <td>1.395776</td>
      <td>9.270621</td>
      <td>12.221057</td>
      <td>9.838189</td>
      <td>10.231783</td>
      <td>6.541789</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TH-01A</th>
      <td>6.188219</td>
      <td>6.700798</td>
      <td>0.000000</td>
      <td>15.205090</td>
      <td>2.080248</td>
      <td>7.979475</td>
      <td>0.942984</td>
      <td>9.452701</td>
      <td>9.092074</td>
      <td>0.547252</td>
      <td>...</td>
      <td>5.667197</td>
      <td>8.415344</td>
      <td>10.678119</td>
      <td>2.806592</td>
      <td>9.077367</td>
      <td>11.745336</td>
      <td>10.982709</td>
      <td>8.823461</td>
      <td>7.318964</td>
      <td>0.547252</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TI-01A</th>
      <td>6.048127</td>
      <td>6.569730</td>
      <td>0.000000</td>
      <td>14.588614</td>
      <td>0.597031</td>
      <td>9.659606</td>
      <td>1.343522</td>
      <td>9.856707</td>
      <td>9.269331</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.338656</td>
      <td>7.813686</td>
      <td>9.478662</td>
      <td>1.017993</td>
      <td>9.136368</td>
      <td>13.101266</td>
      <td>10.261550</td>
      <td>9.512466</td>
      <td>3.491981</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TJ-01A</th>
      <td>6.723098</td>
      <td>7.150324</td>
      <td>0.000000</td>
      <td>14.316598</td>
      <td>1.457699</td>
      <td>9.302794</td>
      <td>1.856030</td>
      <td>9.573077</td>
      <td>10.115430</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.852214</td>
      <td>8.386399</td>
      <td>10.126744</td>
      <td>2.830560</td>
      <td>9.573077</td>
      <td>11.861678</td>
      <td>10.690459</td>
      <td>9.233316</td>
      <td>5.715026</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TK-01A</th>
      <td>7.504929</td>
      <td>6.940841</td>
      <td>0.000000</td>
      <td>14.347271</td>
      <td>0.000000</td>
      <td>9.844355</td>
      <td>2.159339</td>
      <td>9.528442</td>
      <td>9.032145</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.116032</td>
      <td>8.825593</td>
      <td>9.261342</td>
      <td>5.412944</td>
      <td>9.371163</td>
      <td>12.965339</td>
      <td>9.811333</td>
      <td>9.003194</td>
      <td>4.601631</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-MP-A5C7-01A</th>
      <td>6.883095</td>
      <td>6.969186</td>
      <td>0.000000</td>
      <td>13.408743</td>
      <td>0.871844</td>
      <td>5.678430</td>
      <td>2.852938</td>
      <td>9.392369</td>
      <td>9.828102</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.746920</td>
      <td>7.784797</td>
      <td>10.762959</td>
      <td>4.496910</td>
      <td>9.668056</td>
      <td>9.748210</td>
      <td>10.439946</td>
      <td>9.723637</td>
      <td>9.359936</td>
      <td>0.500802</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YF-01A</th>
      <td>8.631193</td>
      <td>7.477151</td>
      <td>0.000000</td>
      <td>12.800867</td>
      <td>0.000000</td>
      <td>7.209679</td>
      <td>0.684819</td>
      <td>10.112063</td>
      <td>10.893280</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.559874</td>
      <td>8.499347</td>
      <td>10.556954</td>
      <td>6.786021</td>
      <td>8.904270</td>
      <td>11.445203</td>
      <td>10.143399</td>
      <td>9.011754</td>
      <td>2.215741</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YG-01A</th>
      <td>7.580348</td>
      <td>7.746573</td>
      <td>1.107487</td>
      <td>13.523957</td>
      <td>0.000000</td>
      <td>6.820296</td>
      <td>1.958583</td>
      <td>9.664667</td>
      <td>9.266543</td>
      <td>2.987030</td>
      <td>...</td>
      <td>4.872386</td>
      <td>7.971544</td>
      <td>10.288866</td>
      <td>4.001037</td>
      <td>9.382152</td>
      <td>11.871504</td>
      <td>9.327452</td>
      <td>9.028908</td>
      <td>5.310471</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YI-01A</th>
      <td>8.437292</td>
      <td>9.143533</td>
      <td>0.741057</td>
      <td>14.250102</td>
      <td>0.741057</td>
      <td>5.940707</td>
      <td>2.816047</td>
      <td>9.735205</td>
      <td>10.237742</td>
      <td>9.940369</td>
      <td>...</td>
      <td>4.834206</td>
      <td>8.201570</td>
      <td>10.323387</td>
      <td>4.834206</td>
      <td>9.691361</td>
      <td>12.204581</td>
      <td>9.872761</td>
      <td>8.667130</td>
      <td>4.307305</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YP-01A</th>
      <td>6.200320</td>
      <td>6.053750</td>
      <td>0.000000</td>
      <td>13.542930</td>
      <td>3.928617</td>
      <td>8.890499</td>
      <td>0.817378</td>
      <td>9.563367</td>
      <td>9.235119</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.033780</td>
      <td>8.459915</td>
      <td>10.558043</td>
      <td>1.182883</td>
      <td>9.569169</td>
      <td>11.972999</td>
      <td>9.825760</td>
      <td>9.398274</td>
      <td>3.828977</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YQ-01A</th>
      <td>7.381001</td>
      <td>7.959587</td>
      <td>0.702835</td>
      <td>13.793663</td>
      <td>0.000000</td>
      <td>7.262493</td>
      <td>0.702835</td>
      <td>9.909071</td>
      <td>9.613185</td>
      <td>7.621003</td>
      <td>...</td>
      <td>5.829581</td>
      <td>8.921819</td>
      <td>10.362085</td>
      <td>6.218128</td>
      <td>9.512625</td>
      <td>12.032069</td>
      <td>10.373736</td>
      <td>9.447964</td>
      <td>4.309787</td>
      <td>0.702835</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55A-01A</th>
      <td>6.649694</td>
      <td>6.437491</td>
      <td>0.000000</td>
      <td>15.766314</td>
      <td>0.000000</td>
      <td>8.725331</td>
      <td>0.906275</td>
      <td>9.337701</td>
      <td>9.202952</td>
      <td>3.818953</td>
      <td>...</td>
      <td>5.786609</td>
      <td>8.302469</td>
      <td>10.431879</td>
      <td>2.168899</td>
      <td>9.402497</td>
      <td>12.072282</td>
      <td>10.637460</td>
      <td>9.379959</td>
      <td>6.148450</td>
      <td>0.906275</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55O-01A</th>
      <td>5.426500</td>
      <td>6.029395</td>
      <td>0.000000</td>
      <td>14.451362</td>
      <td>3.285624</td>
      <td>9.893921</td>
      <td>2.170117</td>
      <td>9.341203</td>
      <td>9.412104</td>
      <td>1.459589</td>
      <td>...</td>
      <td>5.386814</td>
      <td>8.366562</td>
      <td>10.114199</td>
      <td>4.209680</td>
      <td>9.506548</td>
      <td>12.031540</td>
      <td>10.382324</td>
      <td>9.195202</td>
      <td>9.069242</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55R-01A</th>
      <td>8.102809</td>
      <td>7.738950</td>
      <td>5.977124</td>
      <td>14.776279</td>
      <td>0.310689</td>
      <td>7.870827</td>
      <td>1.423309</td>
      <td>9.825515</td>
      <td>11.536570</td>
      <td>8.502341</td>
      <td>...</td>
      <td>5.692087</td>
      <td>8.552138</td>
      <td>10.537639</td>
      <td>5.744756</td>
      <td>9.567679</td>
      <td>11.588149</td>
      <td>10.958776</td>
      <td>9.259183</td>
      <td>9.566765</td>
      <td>1.138421</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A7XG-01A</th>
      <td>7.948200</td>
      <td>7.125734</td>
      <td>0.000000</td>
      <td>12.066421</td>
      <td>1.905967</td>
      <td>8.566999</td>
      <td>0.000000</td>
      <td>9.224335</td>
      <td>9.191943</td>
      <td>0.938022</td>
      <td>...</td>
      <td>4.165373</td>
      <td>7.409182</td>
      <td>10.259693</td>
      <td>3.469157</td>
      <td>9.057056</td>
      <td>10.963035</td>
      <td>10.476904</td>
      <td>8.740755</td>
      <td>2.699352</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-O1-A52J-01A</th>
      <td>7.717136</td>
      <td>8.104844</td>
      <td>0.000000</td>
      <td>14.676937</td>
      <td>0.000000</td>
      <td>5.856433</td>
      <td>2.271366</td>
      <td>9.534633</td>
      <td>10.068412</td>
      <td>3.113617</td>
      <td>...</td>
      <td>7.101856</td>
      <td>9.138731</td>
      <td>10.341686</td>
      <td>0.564134</td>
      <td>9.613415</td>
      <td>12.130652</td>
      <td>10.571349</td>
      <td>8.903853</td>
      <td>7.547854</td>
      <td>0.564134</td>
    </tr>
    <tr>
      <th>TCGA-S2-AA1A-01A</th>
      <td>6.852427</td>
      <td>6.367464</td>
      <td>0.727006</td>
      <td>15.698060</td>
      <td>0.000000</td>
      <td>8.028528</td>
      <td>11.615363</td>
      <td>9.819954</td>
      <td>9.728205</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.937297</td>
      <td>8.512893</td>
      <td>10.838901</td>
      <td>2.916859</td>
      <td>9.285612</td>
      <td>12.384620</td>
      <td>10.629313</td>
      <td>8.996451</td>
      <td>8.319333</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
<p>576 rows × 20472 columns</p>
</div>



# To match samples accross different multi-omics, use


```python
luad_data.match_samples(modalities=["MIR", "GE"])
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



# To prepare the data for classification


```python
# This function selects only patients with patholotic stages "Stage I" and "Stage II"
X_multiomics, y = luad_data.load_data(modalities=["GE", "MIR", "LNC"], target=['pathologic_stage'], 
                                     pathologic_stages=['Stage I', 'Stage II'])
print(X_multiomics['GE'].shape, X_multiomics['MIR'].shape, X_multiomics['LNC'].shape, y.shape)
```

    (336, 20472) (336, 1870) (336, 12727) (336, 1)


    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/TCGAMultiOmics/multiomics.py:152: FutureWarning: 
    Passing list-likes to .loc or [] with any missing label will raise
    KeyError in the future, you can use .reindex() as an alternative.
    
    See the documentation here:
    http://pandas.pydata.org/pandas-docs/stable/indexing.html#deprecate-loc-reindex-listlike
      return self.data["SAMPLES"].loc[matched_samples]



```python
y
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
    <tr>
      <th>TCGA-44-2657-11A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-2665-11A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-44-2668-11A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-5644-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-5645-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6144-11A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6145-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6146-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-44-6148-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6775-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6776-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6776-11A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6777-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6777-11A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-6778-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th>TCGA-J2-8192-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-J2-8194-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-J2-A4AD-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-J2-A4AE-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-J2-A4AG-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-L4-A4E5-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-L4-A4E6-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-L9-A443-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-L9-A444-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MN-A4N1-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MN-A4N4-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MN-A4N5-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SV-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SW-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4SY-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4T4-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TA-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TF-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TH-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TI-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TJ-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-MP-A4TK-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-MP-A5C7-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YF-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YP-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A4YQ-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55A-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55O-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-NJ-A55R-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-O1-A52J-01A</th>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>336 rows × 1 columns</p>
</div>



# Log2 transform the mRNA, microRNA, and lncRNA expression values


```python
def expression_val_transform(x):
    return np.log2(x+1)
X_multiomics['GE'] = X_multiomics['GE'].applymap(expression_val_transform)
X_multiomics['MIR'] = X_multiomics['MIR'].applymap(expression_val_transform)
# X_multiomics['LNC'] = X_multiomics['LNC'].applymap(expression_val_transform)
```

# Classification of Cancer Stage


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

    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/sklearn/preprocessing/label.py:95: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)
    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/sklearn/preprocessing/label.py:128: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)





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
for omic in ["GE", "MIR"]:
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

    GE
    (254, 20472) (109, 20472)


    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/sklearn/utils/validation.py:578: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)


    NONZERO 0
    Training accuracy 0.6929133858267716
                 precision    recall  f1-score   support
    
        Stage I       0.69      1.00      0.82        75
       Stage II       0.00      0.00      0.00        34
    
    avg / total       0.47      0.69      0.56       109
    
    MIR
    (254, 1870) (109, 1870)
    NONZERO 0
    Training accuracy 0.6929133858267716
                 precision    recall  f1-score   support
    
        Stage I       0.69      1.00      0.82        75
       Stage II       0.00      0.00      0.00        34
    
    avg / total       0.47      0.69      0.56       109
    


    /Users/jonny/anaconda/envs/py36/lib/python3.6/site-packages/sklearn/metrics/classification.py:1135: UndefinedMetricWarning: Precision and F-score are ill-defined and being set to 0.0 in labels with no predicted samples.
      'precision', 'predicted', average, warn_for)

