# TCGAMultiOmics

This Python package provide a series of tool to integrate and query the genomics, transcriptomics, proteomics, and clinical TCGA data. By providing a series of data manipulation tools, TCGAMultiOmics facilitaties pan-cancer bioinformatics analysis.

Installation via pip (Python >= 3.6.0): 

    pip install git+https://github.com/JonnyTran/TCGAMultiOmics


The TCGA multi-omics data is downloaded by [TCGA-Assembler](http://www.compgenome.org/TCGA-Assembler/).
Load all multi-omics data files according the following folder structure and naming convention:

    folder_path/
        clinical/
            genome.wustl.edu_biospecimen_sample.txt (optional, not used)
            nationwidechildrens.org_clinical_drug.txt
            nationwidechildrens.org_clinical_patient.txt
        gene_exp/
            geneExp.txt
        mirna/
            miRNAExp__RPM.txt
            TargetScan/ (Obtained from http://www.targetscan.org/vert_72/)
                    TargetScan_miR_Family_Info.txt
                    TargetScan_Predicted_Targets_Context_Scores.default_predictions.txt
                    TargetScan_Predicted_Targets_Info.default_predictions.txt
                    TargetScan_Predicted_Targets_Info_default_predictions.tsv
        cnv/
            copyNumber.txt
        protein_rppa/
            protein_RPPA.txt
        somatic/
            somaticMutation_geneLevel.txt
        lncrna/
            TCGA-rnaexpr.tsv (obtained from http://ibl.mdanderson.org/tanric/_design/basic/download.html)
            HGNC_RNA_long_non-coding.txt (obtained from http://ensembl.org)
            

The microRNA and lncRNA data requires additional external databases, e.g. TargetScan, microRNA family, HGNC long non-coding RNA names, etc. Put these files in the same directory as the data files:


                


# How to use TCGAMultiOmics:


## Importing the TCGAMultiOmics library


```python
from TCGAMultiOmics.multiomics import MultiOmicsData
```

## Import TCGA LUAD data downloaded from TCGA-Assembler


```python
folder_path ="./data/tcga-assembler/LUAD/"
```



```python
# Load all modalities: Gene Expression, MicroRNA expression lncRNA expression, Copy Number Variation, Somatic Mutation, DNA Methylation, and Protein Expression data
luad_data = MultiOmicsData(cancer_type="LUAD", folder_path=folder_path,
                           modalities=["GE", "MIR", "LNC", "CNV", "SNP", "PRO"])

```

Each data is stored as a Pandas DataFrame. Below are all the data imported for TCGA LUAD. For each, the first number represents the number of samples, the second number is the number of features.

    PATIENTS (522, 5)
    DRUGS (461, 4)
    GE (576, 20472)
    SNP (587, 21070)
    MIR (494, 1870)
    LNC (546, 12727)
    CNV (1107, 22516)
    PRO (364, 154)
    SAMPLES (1160, 6)


## Each multi-omics and clinical data can be accessed through luad_data.data[""], for instance:


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
luad_data.data["GE"]
```




<div>
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
  </tbody>
</table>
<p>576 rows × 20472 columns</p>
</div>



## To match samples accross different multi-omics, use


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



## To prepare the data for classification


```python
# This function selects only patients with patholotic stages "Stage I" and "Stage II"
X_multiomics, y = luad_data.load_data(modalities=["GE", "MIR", "LNC"], target=['pathologic_stage'], 
                                     pathologic_stages=['Stage I', 'Stage II'])
print(X_multiomics['GE'].shape, X_multiomics['MIR'].shape, X_multiomics['LNC'].shape, y.shape)
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
X_multiomics['GE'] = X_multiomics['GE'].applymap(expression_val_transform)
X_multiomics['MIR'] = X_multiomics['MIR'].applymap(expression_val_transform)
# X_multiomics['LNC'] = X_multiomics['LNC'].applymap(expression_val_transform)
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
    

