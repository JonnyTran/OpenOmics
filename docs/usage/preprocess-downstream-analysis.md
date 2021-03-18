# Preparing data for downstream analyses

## To prepare the data for classification

```python
X_multiomics, y = luad_data.load_data(omics="all", remove_duplicates=True)

print(X_multiomics['MessengerRNA'].shape,
      X_multiomics['MicroRNA'].shape,
      X_multiomics['LncRNA'].shape,
      y.shape)
```

    (338, 20472) (338, 1870) (338, 12727) (338, 1)


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
