---
title: 'OpenOmics: A bioinformatics API to integrate multi-omics datasets and interface with public databases.'
tags:
  - Python
  - bioinformatics
  - multiomics
  - data integration
  - big data
authors:
  - name: Nhat C. Tran^[corresponding author]
    orcid: 0000-0002-2575-9633
    affiliation: 1
  - name: Jean X. Gao
    affiliation: 1
affiliations:
  - name: Department of Computer Science and Engineering, The University of Texas at Arlington
    index: 1
date: 25 January 2021
bibliography: paper.bib
---

# Summary

Leveraging large-scale multi-omics data is emerging as the primary approach for systemic research of human diseases and
general biological processes. As data integration and feature engineering are the vital steps in these bioinformatics
projects, there currently lacks a tool for standardized preprocessing of heterogeneous multi-omics and annotation data
within the context of a clinical cohort. OpenOmics is a Python library for integrating heterogeneous multi-omics data
and interfacing with popular public annotation databases, e.g., GENCODE, Ensembl, BioGRID. The library is designed to be
highly flexible to allow the user to parameterize the construction of integrated datasets, interactive to assist complex
data exploratory analyses, and scalable to facilitate working with large datasets on standard machines. In this paper,
we demonstrate the software design choices to support the wide-ranging use cases of OpenOmics with the goal of
maximizing usability and reproducibility of the data integration framework.

# Statement of need

Recent advances in sequencing technology and computational methods have enabled the means to generate large-scale,
high-throughput multi-omics data [@lappalainen2013transcriptome], providing unprecedented research opportunities for
cancer and other diseases. These methods have already been applied to a number of problems within bioinformatics, and
indeed several integrative disease
studies [@zhang2014proteogenomic; @cancer2014comprehensive; @ren2016integration; @hassan2020integration]. In addition to
the genome-wide measurements of different genetic characterizations, the growing public knowledge-base of functional
annotations [@rnacentral2016rnacentral; @derrien2012gencode], experimentally-verified
interactions [@chou2015mirtarbase; @yuan2013npinter; @chou2017mirtarbase; @oughtred2019biogrid], and gene-disease
associations [@huang2018hmdd; @pinero2016disgenet; @chen2012lncrnadisease] also provides the prior-knowledge essential
for system-level analyses. Leveraging these data sources allow for a systematic investigation of disease mechanisms at
multiple molecular and regulatory layers; however, such task remains nontrivial due to the complexity of multi-omics
data.

While researchers have developed several mature tools to access or analyze a particular single omic data
type [@wolf2018scanpy; @stuart2019integrative], the current state of integrative data platforms for multi-omics data is
lacking due to three reasons. First, pipelines for data integration carry out a sequential tasks that does not process
multi-omics datasets holistically. Second, the vast size and heterogeneity of the data poses a challenge on the
necessary data storage and computational processing. And third, implementations of data pipelines are close-ended for
down-stream analysis or not conductive to data exploration use-cases. Additionally, there is currently a need for
increased transparency in the process of multi-omics data integration, and a standardized data preprocessing strategy is
important for the interpretation and exchange of bioinformatic projects. Currently, there exist very few systems that,
on the one hand, supports standardized handling of multi-omics datasets but also allows to query the integrated dataset
within the context of a clinical cohort.

# Related works

There are several existing platforms that aids in the integration of multi-omics data, such as Galaxy, Anduril, MixOmics
and O-Miner. First, Galaxy [@boekel2015multi] and Anduril [@cervera2019anduril] are mature platforms and has an
established workflow framework for genomic and transcriptomic data analysis. Galaxy contains hundreds of
state-of-the-art tools of these core domains for processing and assembling high-throughput sequencing data. Second,
MixOmics [@rohart2017mixomics] is an R library dedicated to the multivariate analysis of biological data sets with a
specific focus on data exploration, dimension reduction and visualisation. Third, O-Miner [@sangaralingam2019multi] is
web tool that provides a pipeline for analysis of both transcriptomic and genomic data starting from raw image files
through in-depth bioinformatics analysis. However, as large-scale multi-omic data analysis demands continue to grow, the
technologies and data analysis needs continually change to adapt with `big data`. For instance, the data manipulation
required for multi-omics integration requires a multitude of complex operations, but the point and click interface given
in existing Galaxy tools can be limiting or not computationally efficient. Although the MixOmics toolkit provides an R
programming interface, it doesn't yet leverage high-performance distributed storage or computing resources. Finally,
while O-Miner can perform end-to-end analysis in an integrated platform, its interim analysis results cannot be exported
elsewhere for down-stream analysis.

![Overall OpenOmics System Architecture, Data Flow, and Use Cases.\label{architecture}](figure.pdf)

# The OpenOmics library

OpenOmics consists of two core modules: multi-omics integration and annotation interface. An overview visualization of
the OpenOmics system architecture is provided in \autoref{architecture}.

## Multi-omics integration

Tabular data are everywhere in bioinformatics. To record expression quantifications, annotations, or variant calls, data
are typically stored in various tabular-like formats, such as BED, GTF, MAF, and VCF, which can be preprocessed and
normalized to row indexed formats. Given any processed single-omic dataset, the library generalizes the data as a
tabular structure where rows correspond to observation samples and columns correspond to measurements of different
biomolecules. The core functionality of the Multi-omics Integration module is to integrate the multiple single-omic
datasets for the overlapping samples. By generating multi-omics data for the same set of samples, our tool can provide
the necessary data structure to develop insights into the flow of biological information across multiple genome,
epigenome, transcriptome, proteome, metabolome and phenome levels. The user can import and integrate the following
supported omic types:

- Genomics: single nucleotide variants (SNV), copy number variation (CNV)
- Epigenomics: DNA methylation
- Transcriptomics: RNA-Seq, miRNA expression, lncRNA expression, microarrays
- Proteomics: reverse phase protein array (RPPA), iTRAQ

After importing each single omics data, OpenOmics stores a Pandas Dataframe that is flexible for a wide range of tabular
operations. For instance, the user is presented with several functions for preprocessing of the expression
quantifications to normalize, filter outliers, or reduce noise.

Within a study cohort, the clinical characteristics are crucial for the study of a disease or biological phenomenon. The
user can characterize the set of samples using the Clinical Data structure, which is comprised of two levels: Patient
and Biospecimen. A Patient can have attribute fields on demographics, clinical diagnosis, disease progression, treatment
responses, and survival outcomes. Typically, multi-omics data observations are captured at the Biospecimen level and
each Patient can have multiple Biospecimens. OpenOmics tracks the ID's of biospecimens and the patient it belongs to, so
the multi-omics data are organized in a hierarchical order to enable aggregated operations.

## Annotation interface

After importing and integrating the multi-omic data, the user can supplement their dataset with various annotation
attributes from public data repositories such as GENCODE, Ensembl, and RNA Central. With just a few operations, the user
can easily download a data repository of choice, select relevant attributes, and efficiently join a variable number of
annotation columns to their genomics, transcriptomics, and proteomics data. The full list of databases and the
availability of annotation attributes is listed in Table 1.

For each public database, the Annotation Interface module provides a series of interfaces to perform specific importing,
preprocessing, and annotation tasks. At the import step, the module can either fetch the database files via a
file-transfer-protocol (ftp) URL or load a locally downloaded file. At this step, the user can specify the species,
genome build, and version of the database by providing a ftp URL of choice. To streamline this process, the module
automatically caches downloaded file to disk, uncompress them, and handle different file extensions, including FASTA,
GTF, VCF, and other tabular formats. Then, at the preprocessing step, the module selects only the relevant attribute
fields specified by the user and perform necessary data cleanings. Finally, the annotation data can be annotated to an
omics dataset by performing a SQL-like join operation on a user-specified index of the biomolecule name or ID. If the
user wishes to import an annotation database not yet included in OpenOmics, they can extend the Annotation Dataset API
to specify their own importing, preprocessing, and annotation tasks in an object-oriented manner.

An innovative feature of our integration module is the ability to cross-reference the gene ID's between different
annotation systems or data sources. When importing a dataset, the user can specify the level of genomic index, such as
at the gene, transcript, protein, or peptide level, and whether it is a gene name or gene ID. Since multiple
single-omics datasets can use different gene nomenclatures, the user is able to convert between the different gene
indexing methods by reindexing the annotation dataframe with a index column of choice. This not only allows the
Annotation Interface to select and join the annotation data to the correct index level, but also allow the user to
customize the selection and aggregation of biological measurements at different levels.

| Data Repository | Annotation Data Available                       | Index     | # entries  |
| --------------- | ----------------------------------------------- | --------- | ---------- |
| GENCODE         | Genomic annotations, primary sequence           | RNAs      | 60660      |
| Ensembl         | Genomic annotations & &                        | Genes     | 232,186    |
| MiRBase         | MicroRNA sequences and annotatinos              | MicroRNAs | 38589      |
| RNA Central     | ncRNA sequence and annotation collection        | ncRNAs    | 14,784,981 |
| NONCODE         | lncRNA sequences and annotations                | LncRNAs   | 173,112    |
| lncrnadb        | lncRNA functional annotations                   | LncRNAs   | 100        |
| Pfam            | Protein family annotation                       | Proteins  | 18,259     |
| Rfam            | RNA family annotations                          | ncRNAs    | 2,600      |
| Gene Ontology   | Functional, cellular, and molecular annotations | Genes     | 44,117     |
| KEGG            | High-level functional pathways                  | Genes     | 22,409     |
| DisGeNet        | gene-disease associations                       | Genes     | 1,134,942  |
| HMDD            | microRNA-disease associations                   | MicroRNAs | 35,547     |
| lncRNAdisease   | lncRNA-disease associations                     | LncRNAs   | 3,000      |
| OMIM            | Ontology of human diseases                      | Diseases  | 25,670     |

Table 1: Public annotation databases and availability of data in the Human genome.

# System design

This section describes the various implementation details behind the scalable processing and efficient data storage, and
the design choices in the development operations.

While the in-memory Pandas dataframes utilized in our data structures are fast, they have size and speed limitations
when the dataset size approaches the system memory limit. When this is an issue, the user can enable out-of-memory
distributed data processing on all OpenOmics operations, implemented by the Dask framework^[https://dask.org/]. When
memory resources is limited, data in a Dask dataframe can be read directly from disk and is only brought into memory
when needed during computations (also called lazy evaluations). When performing data query operations on Dask
dataframes, a task graph containing each operation is built and is only evaluated on command, in a process called lazy
loading.

Operations on Dask dataframes are the same as Pandas dataframes, but can utilize multiple workers and can scale up to
clusters by connecting to a cluster client with minimal configuration. To enable this feature in OpenOmics, the user
simply needs to explicitly enable an option when importing an omics dataset, importing an annotation/interaction
database, or importing a MultiOmics file structure on disk.

## Software requirements

OpenOmics is distributed as a readily installable Python package from the Python Package Index (PyPI) repository. For
users to install OpenOmics in their own Python environment, several software dependencies are automatically downloaded
to reproduce the computing environment.

OpenOmics is compatible with Python 3.6 or higher, and is operational on both Linux and Windows operating systems. The
software requires as little as 4 GB of RAM and 2 CPU cores, and can computationally scale up to large-memory
multi-worker distributed systems such as a compute cluster. To take advantage of increased computational resource,
OpenOmics simply requires one line of code to activate parallel computing functionalities.

## Development operations

We developed OpenOmics following modern software best-practices and package publishing standards. For the version
control of our source-code, we utilized a public GitHub repository which contains two branches, master and develop. The
master branch contains stable and well-tested releases of the package, while the develop branch is used for building new
features or software refactoring. Before each version is released, we utilize Github Actions for continuous integration,
building, and testing for version and dependency compatibility. Our automated test suite covers essential functions of
the package and a reasonable range of inputs and conditions.

# Conclusion

A standardized data preprocessing strategy is essential for the interpretation and exchange of bioinformatics research.
OpenOmics provides researchers with the means to consistently describe the processing and analysis of their experimental
datasets. It equips the user, a bioinformatician, with the ability to preprocess, query, and analyze data with modern
and scalable software technology. As the wide array of tools and methods available in the public domain are largely
isolated, OpenOmics aims toward a uniform framework that can effectively process and analyze multi-omics data in an
end-to-end manner along with biologist-friendly visualization and interpretation.

# Acknowledgements

N/A.

# References
