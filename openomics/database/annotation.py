import os
import re
from os.path import expanduser

import dask.dataframe as dd
import pandas as pd
from bioservices import BioMart
from io import StringIO
from pandas.errors import ParserError

from .base import Database
from ..transforms.agg import get_agg_func

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openomics")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), ".openomics", "databases")


class ProteinAtlas(Database):
    """Loads the  database from  .

        Default path:  .
        Default file_resources: {
            "": "",
            "": "",
            "": "",
        }
        """
    COLUMNS_RENAME_DICT = {
        "Gene": "protein_name",
        "Ensembl": "gene_id",
    }

    def __init__(self, path="https://www.proteinatlas.org/download/", file_resources=None,
                 col_rename=COLUMNS_RENAME_DICT, blocksize=0, verbose=False):
        """
        Args:
            path:
            file_resources:
            col_rename:
            blocksize:
            verbose:
        """
        if file_resources is None:
            file_resources = {}
            file_resources["proteinatlas.tsv.zip"] = "proteinatlas.tsv.zip"

        super().__init__(path, file_resources, col_rename=col_rename, blocksize=blocksize, verbose=verbose)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        if blocksize:
            df = dd.read_table(file_resources["proteinatlas.tsv"],
                               blocksize=None if isinstance(blocksize, bool) else blocksize)
        else:
            df = pd.read_table(file_resources["proteinatlas.tsv"])

        return df

    def get_expressions(self, index="gene_name", type="Tissue RNA"):
        """Returns (NX) expressions from the proteinatlas.tsv table. :param
        index: a column name to index by. If column contain multiple values,
        then aggregate by median values. :param type: one of {"Tissue RNA",
        "Cell RNA", "Blood RNA", "Brain RNA", "RNA - "}. If "RNA - ", then
        select all types of expressions.

        Args:
            index:
            type:

        Returns:
            expressions (pd.DataFrame):
        """
        columns = "|".join([type, index])
        expressions = self.data.filter(regex=columns).groupby(
            index).median()
        return expressions


class RNAcentral(Database):
    """Loads the RNAcentral database from https://rnacentral.org/ .

        Default path: ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/ .
        Default file_resources: {
            "rnacentral_rfam_annotations.tsv": "go_annotations/rnacentral_rfam_annotations.tsv.gz",
            "database_mappings/gencode.tsv": "id_mapping/database_mappings/gencode.tsv",
            "database_mappings/mirbase.tsv": "id_mapping/database_mappings/mirbase.tsv",
            ...
        }
    """
    COLUMNS_RENAME_DICT = {
        'ensembl_gene_id': 'gene_id',
        'external id': 'transcript_id',
        'GO terms': 'go_id',
        'gene symbol': 'gene_id',
    }

    def __init__(self, path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/", file_resources=None,
                 col_rename=COLUMNS_RENAME_DICT, species_id: str = '9606',
                 index_col="RNAcentral id",
                 remove_version_num=True, **kwargs):
        """
        Args:
            path:
            file_resources:
            col_rename:
            species_id:
        """
        self.species_id = species_id
        self.remove_version_num = remove_version_num

        if file_resources is None:
            file_resources = {}
            file_resources["rnacentral_rfam_annotations.tsv.gz"] = "go_annotations/rnacentral_rfam_annotations.tsv.gz"
            file_resources["database_mappings/gencode.tsv"] = "id_mapping/database_mappings/gencode.tsv"
            file_resources["database_mappings/mirbase.tsv"] = "id_mapping/database_mappings/mirbase.tsv"

        super().__init__(path, file_resources, col_rename=col_rename, index_col=index_col, **kwargs)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        # Build transcripts ids by combining `database_mappings/` files
        transcripts_df = []
        for filename in file_resources:
            if "database_mappings" in filename:
                args = dict(low_memory=True, header=None,
                            names=["RNAcentral id", "database", "external id", "species_id", "RNA type", "gene symbol"],
                            dtype={'gene symbol': 'str',
                                   'database': 'category', 'species_id': 'category', 'RNA type': 'category', })

                if blocksize:
                    if filename.endswith('.tsv'):
                        id_mapping: dd.DataFrame = dd.read_table(file_resources[filename],
                                                                 blocksize=None if isinstance(blocksize,
                                                                                              bool) else blocksize,
                                                                 **args)
                    elif filename.endswith('.parquet'):
                        id_mapping: dd.DataFrame = dd.read_parquet(file_resources[filename],
                                                                   blocksize=None if isinstance(blocksize,
                                                                                                bool) else blocksize, )

                    id_mapping = id_mapping.set_index("RNAcentral id", sorted=True)
                else:
                    if filename.endswith('.tsv'):
                        id_mapping = pd.read_table(file_resources[filename], index_col="RNAcentral id", **args)
                    elif filename.endswith('.parquet'):
                        id_mapping = pd.read_parquet(file_resources[filename], index_col="RNAcentral id")

                if self.remove_version_num and 'gene symbol' in id_mapping.columns:
                    id_mapping["gene symbol"] = id_mapping["gene symbol"].str.replace("[.].\d*", "", regex=True)

                transcripts_df.append(id_mapping)

        # Concatenate multiple `database_mappings` files from different databases
        if blocksize:
            transcripts_df = dd.concat(transcripts_df, axis=0, join='outer')
        else:
            transcripts_df = pd.concat(transcripts_df, axis=0, join='outer')

        if self.species_id:
            transcripts_df = transcripts_df.where(transcripts_df["species_id"] == self.species_id)

        # Join go_id and Rfams annotations to each "RNAcentral id" from 'rnacentral_rfam_annotations.tsv'
        args = dict(low_memory=True, names=["RNAcentral id", "GO terms", "Rfams"])
        if self.remove_version_num:
            args['converters'] = {"RNAcentral id": lambda s: re.sub("_(.*)", '', s)}
        if blocksize:
            if 'rnacentral_rfam_annotations.tsv' in file_resources and isinstance(
                file_resources['rnacentral_rfam_annotations.tsv'], str):
                anns: dd.DataFrame = dd.read_table(file_resources["rnacentral_rfam_annotations.tsv"], **args)
            else:
                anns = dd.read_table(file_resources["rnacentral_rfam_annotations.tsv.gz"], compression="gzip", **args)
            anns = anns.set_index("RNAcentral id", sorted=True, npartitions=10)

            # Filter annotations by "RNAcentral id" in `transcripts_df`
            idx = transcripts_df.index.compute()
            anns = anns.loc[anns.index.isin(idx)]

            # Groupby on index
            anns_groupby: dd.DataFrame = anns \
                .groupby(by=lambda idx: idx) \
                .agg({col: get_agg_func('unique', use_dask=True) for col in ["GO terms", 'Rfams']}) \
                .persist()
        else:
            anns = pd.read_table(file_resources["rnacentral_rfam_annotations.tsv"], index_col='RNAcentral id', **args)
            idx = transcripts_df.index.drop_duplicates()
            anns = anns.loc[anns.index.isin(idx)]
            anns_groupby = anns.groupby("RNAcentral id").agg({col: 'unique' for col in ["GO terms", 'Rfams']})

        transcripts_df = transcripts_df.merge(anns_groupby, how='left', left_index=True, right_index=True)

        return transcripts_df


class GTEx(Database):
    """Loads the  database from https://www.gtexportal.org/home/ .

    Default path: "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/" .
    Default file_resources: {
        "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct": "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt": "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct": "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
    }
    """
    COLUMNS_RENAME_DICT = {
        "Name": "gene_id",
        "Description": "gene_name"
    }

    def __init__(self, path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/",
                 file_resources=None, col_rename=None, blocksize=0, verbose=False):
        """
        Args:
            path:
            file_resources:
            col_rename:
            blocksize:
            verbose:
        """
        if file_resources is None:
            file_resources = {}

            file_resources["GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"] = \
                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
            file_resources["GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"] = \
                "https://storage.googleapis.com/gtex_analysis_v8/annotations/" \
                "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
            file_resources["GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"] = \
                "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"

        super().__init__(path, file_resources, col_rename=None, blocksize=blocksize, verbose=verbose)

    def load_dataframe(self, file_resources, blocksize=None):  # type: (dict) -> pd.DataFrame
        """
        Args:
            file_resources:
            blocksize:
        """
        gene_exp_medians = pd.read_csv(
            self.file_resources["GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"],
            sep='\t', header=1, skiprows=1)
        gene_exp_medians["Name"] = gene_exp_medians["Name"].str.replace("[.]\d*", "", regex=True)
        gene_exp_medians = gene_exp_medians.rename(columns=self.COLUMNS_RENAME_DICT)  # Must be done here
        gene_exp_medians.set_index(["gene_id", "gene_name"], inplace=True)

        # # Sample attributes (needed to get tissue type)
        # SampleAttributes = pd.read_table(
        #     self.file_resources["GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"],
        # )
        # SampleAttributes.set_index("SAMPID", inplace=True)
        #
        # # Transcript expression for all samples
        # transcript_exp = pd.read_csv(
        #     self.file_resources["GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct"],
        #     sep='\t', header=1, skiprows=1)
        # print("transcript_exp", transcript_exp.columns)
        # transcript_exp["gene_id"] = transcript_exp["gene_id"].str.replace("[.]\d*", "")
        # transcript_exp["transcript_id"] = transcript_exp["transcript_id"].str.replace("[.]\d*", "")
        # transcript_exp.set_index(["gene_id", "transcript_id"], inplace=True)
        #
        # # Join by sample with tissue type, group expressions by tissue type, and compute medians for each
        # transcript_exp_medians = transcript_exp.T \
        #     .join(SampleAttributes["SMTSD"], how="left") \
        #     .groupby("SMTSD") \
        #     .median()
        #
        # # Reset multilevel index
        # transcript_exp_medians.index.rename(name=None, inplace=True)
        # transcript_exp_medians = transcript_exp_medians.T.set_index(
        #     pd.MultiIndex.from_tuples(tuples=transcript_exp_medians.T.index, names=["gene_id", "transcript_id"]))
        #
        # gene_transcript_exp_medians = pd.concat([gene_exp_medians, transcript_exp_medians], join="inner", copy=True)
        # print("gene_transcript_exp_medians \n", gene_transcript_exp_medians)
        return gene_exp_medians


class NONCODE(Database):
    """Loads the NONCODE database from http://noncode.org .

    Default path: "http://www.noncode.org/datadownload" .
    Default file_resources: {
        "NONCODEv6_human.fa": "NONCODEv6_human.fa.gz",
        "": "",
        "": "",
    }
    """

    def __init__(self, path="http://www.noncode.org/datadownload", file_resources=None, col_rename=None, verbose=False,
                 blocksize=None):
        """
        Args:
            path:
            file_resources:
            col_rename:
            verbose:
            blocksize:
        """
        if file_resources is None:
            file_resources = {}
            file_resources["NONCODEv5_source"] = os.path.join(path, "NONCODEv5_source")
            file_resources["NONCODEv5_Transcript2Gene"] = os.path.join(path, "NONCODEv5_Transcript2Gene")
            file_resources["NONCODEv5_human.func"] = os.path.join(path, "NONCODEv5_human.func")

        super().__init__(path, file_resources, col_rename=col_rename, blocksize=blocksize, verbose=verbose)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        source_df = pd.read_table(file_resources["NONCODEv5_source"], header=None)
        source_df.columns = ["NONCODE Transcript ID", "name type", "Gene ID"]

        transcript2gene_df = pd.read_table(file_resources["NONCODEv5_Transcript2Gene"], header=None)
        transcript2gene_df.columns = ["NONCODE Transcript ID", "NONCODE Gene ID"]

        if blocksize:
            self.noncode_func_df = dd.read_table(file_resources["NONCODEv5_human.func"], header=None,
                                                 blocksize=None if isinstance(blocksize, bool) else blocksize)
        else:
            self.noncode_func_df = pd.read_table(file_resources["NONCODEv5_human.func"], header=None)
        self.noncode_func_df.columns = ["NONCODE Gene ID", "GO terms"]
        self.noncode_func_df.set_index("NONCODE Gene ID", inplace=True)

        # Convert to NONCODE transcript ID for the functional annotation data
        self.noncode_func_df["NONCODE Transcript ID"] = self.noncode_func_df.index.map(
            pd.Series(transcript2gene_df['NONCODE Transcript ID'].values,
                      index=transcript2gene_df['NONCODE Gene ID']).to_dict())

        # Convert NONCODE transcript ID to gene names
        source_gene_names_df = source_df[source_df["name type"] == "NAME"].copy()

        self.noncode_func_df["Gene Name"] = self.noncode_func_df["NONCODE Transcript ID"].map(
            pd.Series(source_gene_names_df['Gene ID'].values,
                      index=source_gene_names_df['NONCODE Transcript ID']).to_dict())


class BioMartManager:
    """
    A base class with functions to query Ensembl Biomarts "https://www.ensembl.org/biomart".
    """
    DTYPES = {
        'entrezgene_id': 'str',
        'gene_biotype': 'category',
        'transcript_biotype': 'category',
        'chromosome_name': 'category',
        'transcript_start': 'int',
        'transcript_end': 'int',
        'transcript_length': 'int',
        'mirbase_id': 'str'}

    def __init__(self, dataset, attributes, host, filename):
        """
        Args:
            dataset:
            attributes:
            host:
            filename:
        """
        pass  # Does not instantiate

    def retrieve_dataset(self, host, dataset, attributes, filename, blocksize=None):
        """
        Args:
            host:
            dataset:
            attributes:
            filename:
            blocksize:
        """
        filename = os.path.join(DEFAULT_CACHE_PATH, f"{filename}.tsv")

        args = dict(
            sep="\t",
            low_memory=True,
            dtype=self.DTYPES,
        )

        if os.path.exists(filename):
            if blocksize:
                df = dd.read_csv(filename, blocksize=None if isinstance(blocksize, bool) else blocksize, **args)
            else:
                df = pd.read_csv(filename, **args)
        else:
            df = self.query_biomart(host=host, dataset=dataset, attributes=attributes,
                                    cache=True, save_filename=filename)
        return df

    def cache_dataset(self, dataset, dataframe, save_filename):
        """
        Args:
            dataset:
            dataframe:
            save_filename:
        """
        if not os.path.exists(DEFAULT_CACHE_PATH):
            os.makedirs(DEFAULT_CACHE_PATH, exist_ok=True)

        if save_filename is None:
            save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(dataset))

        dataframe.to_csv(save_filename, sep="\t", index=False)
        return save_filename

    def query_biomart(self, dataset, attributes, host="www.ensembl.org", cache=True, save_filename=None,
                      blocksize=None):
        """
        Args:
            dataset:
            attributes:
            host:
            cache:
            save_filename:
            blocksize:
        """
        bm = BioMart(host=host)
        bm.new_query()
        bm.add_dataset_to_xml(dataset)
        for at in attributes:
            bm.add_attribute_to_xml(at)
        xml_query = bm.get_xml()

        print("Querying {} from {} with attributes {}...".format(dataset, host, attributes))
        results = bm.query(xml_query)

        try:
            if blocksize:
                df = dd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", low_memory=True,
                                 dtype=self.DTYPES, blocksize=None if isinstance(blocksize, bool) else blocksize)
            else:
                df = pd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", low_memory=True,
                                 dtype=self.DTYPES)
        except Exception as e:
            raise ParserError(f'BioMart Query Result: {results}')

        if cache:
            self.cache_dataset(dataset, df, save_filename)
        return df


class EnsemblGenes(BioMartManager, Database):
    COLUMNS_RENAME_DICT = {'ensembl_gene_id': 'gene_id',
                           'external_gene_name': 'gene_name',
                           'ensembl_transcript_id': 'transcript_id',
                           'external_transcript_name': 'transcript_name',
                           'rfam': 'Rfams'}

    def __init__(self, biomart="hsapiens_gene_ensembl",
                 attributes=None, host="www.ensembl.org", blocksize=None):
        # Do not call super().__init__()
        """
        Args:
            biomart:
            attributes:
            host:
            blocksize:
        """
        if attributes is None:
            attributes = ['ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id',
                          'external_transcript_name',
                          'chromosome_name', 'transcript_start', 'transcript_end', 'transcript_length',
                          'gene_biotype', 'transcript_biotype', ]
        self.filename = "{}.{}".format(biomart, self.__class__.__name__)

        self.biomart = biomart
        self.host = host
        self.data = self.load_data(dataset=biomart, attributes=attributes, host=self.host,
                                   filename=self.filename, blocksize=blocksize)

        self.data = self.data.rename(columns=self.COLUMNS_RENAME_DICT)

    def name(self):
        return f"{super().name()} {self.biomart}"

    def load_data(self, dataset, attributes, host, filename=None, blocksize=None):
        """
        Args:
            dataset:
            attributes:
            host:
            filename:
            blocksize:
        """
        df = self.retrieve_dataset(host, dataset, attributes, filename, blocksize=blocksize)
        return df

class EnsemblGeneSequences(EnsemblGenes):
    def __init__(self, biomart="hsapiens_gene_ensembl",
                 attributes=None, host="www.ensembl.org", blocksize=None):
        """
        Args:
            biomart:
            attributes:
            host:
            blocksize:
        """
        if attributes is None:
            attributes = ['ensembl_gene_id', 'gene_exon_intron', 'gene_flank', 'coding_gene_flank', 'gene_exon',
                          'coding']
        self.filename = "{}.{}".format(biomart, self.__class__.__name__)

        self.biomart = biomart
        self.host = host
        self.df = self.load_data(dataset=biomart, attributes=attributes, host=self.host,
                                 filename=self.filename, blocksize=blocksize)
        self.data = self.data.rename(columns=self.COLUMNS_RENAME_DICT)


class EnsemblTranscriptSequences(EnsemblGenes):
    def __init__(self, biomart="hsapiens_gene_ensembl",
                 attributes=None, host="www.ensembl.org", blocksize=None):
        """
        Args:
            biomart:
            attributes:
            host:
            blocksize:
        """
        if attributes is None:
            attributes = ['ensembl_transcript_id', 'transcript_exon_intron', 'transcript_flank',
                          'coding_transcript_flank',
                          '5utr', '3utr']
        self.filename = "{}.{}".format(biomart, self.__class__.__name__)

        self.biomart = biomart
        self.host = host
        self.df = self.load_data(dataset=biomart, attributes=attributes, host=self.host,
                                 filename=self.filename, blocksize=blocksize)
        self.data = self.data.rename(columns=self.COLUMNS_RENAME_DICT)

class EnsemblSNP(EnsemblGenes):
    def __init__(self, biomart="hsapiens_snp",
                 attributes=None, host="www.ensembl.org", blocksize=None):
        """
        Args:
            biomart:
            attributes:
            host:
            blocksize:
        """
        if attributes is None:
            attributes = ['synonym_name', 'variation_names', 'minor_allele',
                          'associated_variant_risk_allele',
                          'ensembl_gene_stable_id', 'ensembl_transcript_stable_id',
                          'phenotype_name',
                          'chr_name', 'chrom_start', 'chrom_end']
        self.filename = "{}.{}".format(biomart, self.__class__.__name__)

        self.biomart = biomart
        self.host = host
        self.data = self.data.rename(columns=self.COLUMNS_RENAME_DICT)


class EnsemblSomaticVariation(EnsemblGenes):
    def __init__(self, biomart="hsapiens_snp_som",
                 attributes=None, host="www.ensembl.org", blocksize=None):
        """
        Args:
            biomart:
            attributes:
            host:
            blocksize:
        """
        if attributes is None:
            attributes = ['somatic_variation_name', 'somatic_source_name', 'somatic_allele', 'somatic_minor_allele',
                          'somatic_clinical_significance', 'somatic_validated', 'somatic_transcript_location',
                          'somatic_mapweight',
                          'somatic_chromosome_start', 'somatic_chromosome_end']
        self.filename = "{}.{}".format(biomart, self.__class__.__name__)

        self.biomart = biomart
        self.host = host
        self.data = self.data.rename(columns=self.COLUMNS_RENAME_DICT)


class TANRIC(Database):
    def __init__(self, path, file_resources=None, col_rename=None, blocksize=0, verbose=False):
        """
        Args:
            path:
            file_resources:
            col_rename:
            blocksize:
            verbose:
        """
        super().__init__(path, file_resources, col_rename=col_rename, blocksize=blocksize, verbose=verbose)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        pass

    def get_expressions(self, genes_index):
        """Preprocess LNCRNA expression file obtained from TANRIC MDAnderson,
        and replace ENSEMBL gene ID to HUGO gene names (HGNC). This function
        overwrites the GenomicData.process_expression_table() function which
        processes TCGA-Assembler data. TANRIC LNCRNA expression values are log2
        transformed

        Args:
            genes_index:
        """
        df = pd.read_table(self.file_resources["TCGA-LUAD-rnaexpr.tsv"])
        df[genes_index] = df[genes_index].str.replace("[.]\d*", "")  # Removing .# ENGS gene version number at the end
        df = df[~df[genes_index].duplicated(keep='first')]  # Remove duplicate genes

        # Drop NA gene rows
        df.dropna(axis=0, inplace=True)

        # Transpose matrix to patients rows and genes columns
        df.index = df[genes_index]
        df = df.T.iloc[1:, :]

        # Change index string to bcr_sample_barcode standard
        def change_patient_barcode(s):
            if "Normal" in s:
                return s[s.find('TCGA'):] + "-11A"
            elif "Tumor" in s:
                return s[s.find('TCGA'):] + "-01A"
            else:
                return s

        df.index = df.index.map(change_patient_barcode)
        df.index.name = "gene_id"

        return df
