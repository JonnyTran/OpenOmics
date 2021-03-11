# Release history

## 0.9.0 (Future)
- Build web-app dashboard interface for importing user multiomics data files

## 0.8.6 (Pending)
- Changed to Github Actions CI from Travis CI
- Revamped openomics.readthedocs.io
- Fixed bugs from pyOpenSci reviews.

## 0.8.5 (2020-02-19)
- Importing GENCODE gtf files using dask now works with gzip compressed files.
- Improved coverage and performance of automated tests.

## 0.8.4 (2020-01-07)
- Enabled the support for Dask dataframes for ExpressionData and Annotation Dataset's. To use this feature, simply use
  the npartitions argument when instantiating an ExpressionData or Annotation/Sequence/Interaction Dataset.

## 0.7.2 (2019-09-01)
- Added compatibility for Python 2.7
- Refactored ClinicalData
- Built working documentations with Sphinx on readthedocs
- Added pytests for MultiOmicsData
- First successful build on Travis CI on Python 3.4-3.7, 2.7

## 0.6.0 (2019-08-31)
- First release on PyPI.
