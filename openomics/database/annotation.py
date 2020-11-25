import os
from io import StringIO
from os.path import expanduser

from bioservices import BioMart

import openomics
from openomics.database.base import Dataset
from openomics.utils.df import concat_uniques
from openomics.utils.io import mkdirs

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openomics")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), ".openomics", "databases")

if openomics.__BACKEND__ == "dask":
    import dask.dataframe as dd
else:
    import pandas as pd

class TANRIC(Dataset):
    def __init__(self, path, file_resources=None, col_rename=None, npartitions=0, verbose=False):
        super(TANRIC, self).__init__(path, file_resources, col_rename, npartitions, verbose)

    def load_dataframe(self, file_resources):
        pass

    def get_expressions(self, genes_index):
        """
        Preprocess LNCRNA expression file obtained from TANRIC MDAnderson, and replace ENSEMBL gene ID to HUGO gene names (HGNC). This function overwrites the GenomicData.process_expression_table() function which processes TCGA-Assembler data.
        TANRIC LNCRNA expression values are log2 transformed
        """
        df = pd.read_table(self.file_resources["TCGA-LUAD-rnaexpr.tsv"])
        df[genes_index] = df[genes_index].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
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


class ProteinAtlas(Dataset):
    COLUMNS_RENAME_DICT = {
        "Gene": "protein_name",
        "Ensembl": "gene_id",
    }

    def __init__(self, path="https://www.proteinatlas.org/download/", file_resources=None,
                 col_rename=COLUMNS_RENAME_DICT,
                 npartitions=0, verbose=False):
        if file_resources is None:
            file_resources = {}
            file_resources["proteinatlas.tsv"] = "proteinatlas.tsv.zip"

        super(ProteinAtlas, self).__init__(path, file_resources, col_rename, npartitions, verbose)

    def load_dataframe(self, file_resources):
        df = pd.read_table(file_resources["proteinatlas.tsv"])

        return df

    def get_expressions(self, index="gene_name", type="Tissue RNA"):
        """
        Returns (NX) expressions from the proteinatlas.tsv table.
        Args:
            index: a column name to index by. If column contain multiple values, then aggregate by median values.
            type: one of {"Tissue RNA", "Cell RNA", "Blood RNA", "Brain RNA", "RNA - "}. If "RNA - ", then select all types of expressions.

        Returns:
            expressions (pd.DataFrame):
        """
        columns = "|".join([type, index])
        expressions = self.df.filter(regex=columns).groupby(
            index).median()
        return expressions


class RNAcentral(Dataset):
    COLUMNS_RENAME_DICT = {'ensembl_gene_id': 'gene_id',
                           'gene symbol': 'gene_name',
                           'external id': 'transcript_id',
                           'GO terms': 'go_id'}

    def __init__(self, path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/",
                 file_resources=None, col_rename=COLUMNS_RENAME_DICT, npartitions=0, species=9606, verbose=False):
        self.species = species

        if file_resources is None:
            file_resources = {}
            file_resources["rnacentral_rfam_annotations.tsv"] = "go_annotations/rnacentral_rfam_annotations.tsv.gz"
            file_resources["database_mappings/gencode.tsv"] = "id_mapping/database_mappings/gencode.tsv"
            file_resources["database_mappings/mirbase.tsv"] = "id_mapping/database_mappings/mirbase.tsv"

        super(RNAcentral, self).__init__(path, file_resources, col_rename=col_rename, npartitions=npartitions,
                                         verbose=verbose)

    def load_dataframe(self, file_resources):
        go_terms = pd.read_table(file_resources["rnacentral_rfam_annotations.tsv"],
                                 low_memory=True, header=None, names=["RNAcentral id", "GO terms", "Rfams"])
        go_terms["RNAcentral id"] = go_terms["RNAcentral id"].str.split("_", expand=True, n=2)[0]

        gene_ids = []
        for file in file_resources:
            if "database_mappings" in file:
                id_mapping = pd.read_table(file_resources[file],
                                           low_memory=True, header=None,
                                           names=["RNAcentral id", "database", "external id", "species", "RNA type",
                                                  "gene symbol"])
                gene_ids.append(id_mapping)
        gene_ids = pd.concat(gene_ids, join="inner")

        gene_ids["species"] = gene_ids["species"].astype("O")
        if self.species is not None:
            gene_ids = gene_ids[gene_ids["species"] == self.species]

        lnc_go_terms = go_terms[go_terms["RNAcentral id"].isin(gene_ids["RNAcentral id"])].groupby("RNAcentral id")[
            "GO terms"].apply(lambda x: "|".join(x.unique()))
        lnc_rfams = go_terms[go_terms["RNAcentral id"].isin(gene_ids["RNAcentral id"])].groupby("RNAcentral id")[
            "Rfams"].apply(lambda x: "|".join(x.unique()))

        gene_ids["GO terms"] = gene_ids["RNAcentral id"].map(lnc_go_terms)
        gene_ids["Rfams"] = gene_ids["RNAcentral id"].map(lnc_rfams)
        gene_ids = gene_ids[gene_ids["GO terms"].notnull() | gene_ids["Rfams"].notnull()]

        return gene_ids


class GTEx(Dataset):
    COLUMNS_RENAME_DICT = {
        "Name": "gene_id",
        "Description": "gene_name"
    }

    def __init__(self, path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/",
                 file_resources=None, col_rename=None, npartitions=0, verbose=False):
        if file_resources is None:
            file_resources = {}

            file_resources[
                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"] = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
            file_resources[
                "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"] = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
            file_resources[
                "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct"] = "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"

        super(GTEx, self).__init__(path, file_resources, col_rename=None, npartitions=npartitions, verbose=verbose)

    def load_dataframe(self, file_resources):  # type: (dict) -> pd.DataFrame
        gene_exp_medians = pd.read_csv(
            self.file_resources["GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"],
            sep='\t', header=1, skiprows=1)
        gene_exp_medians["Name"] = gene_exp_medians["Name"].str.replace("[.].*", "")
        gene_exp_medians = gene_exp_medians.rename(columns=self.COLUMNS_RENAME_DICT)  # Must be done here
        gene_exp_medians.set_index(["gene_id", "gene_name"], inplace=True)
        #
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
        # transcript_exp["gene_id"] = transcript_exp["gene_id"].str.replace("[.].*", "")
        # transcript_exp["transcript_id"] = transcript_exp["transcript_id"].str.replace("[.].*", "")
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


class BioMartManager:
    def __init__(self, dataset, attributes, host, filename):
        pass

    def query_biomart(self, dataset, attributes, host="www.ensembl.org", cache=True, save_filename=None):
        bm = BioMart(host=host)
        bm.new_query()
        bm.add_dataset_to_xml(dataset)
        for at in attributes:
            bm.add_attribute_to_xml(at)
        xml_query = bm.get_xml()

        print("Querying {} from {} with attributes {}...".format(dataset, host, attributes))
        results = bm.query(xml_query)
        df = pd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", low_memory=True)

        if cache:
            self.cache_dataset(dataset, df, save_filename)
        return df

    def cache_dataset(self, dataset, dataframe, save_filename):
        if not os.path.exists(DEFAULT_CACHE_PATH):
            mkdirs(DEFAULT_CACHE_PATH)

        if save_filename is None:
            save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(dataset))

        dataframe.to_csv(save_filename, sep="\t", index=False)
        return save_filename

    def retrieve_dataset(self, host, dataset, attributes, filename):
        filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(filename))
        if os.path.exists(filename):
            df = pd.read_csv(filename, sep="\t", low_memory=True)
        else:
            df = self.query_biomart(host=host, dataset=dataset, attributes=attributes,
                                    cache=True, save_filename=filename)
        return df


class EnsemblGenes(BioMartManager, Dataset):
    COLUMNS_RENAME_DICT = {'ensembl_gene_id': 'gene_id',
                           'external_gene_name': 'gene_name',
                           'ensembl_transcript_id': 'transcript_id',
                           'external_transcript_name': 'transcript_name',
                           'rfam': 'Rfams'}

    def __init__(self, dataset="hsapiens_gene_ensembl",
                 attributes=None,
                 host="www.ensembl.org", filename=False):
        if attributes is None:
            attributes = ['ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id',
                          'external_transcript_name',
                          'chromosome_name', 'transcript_start', 'transcript_end', 'transcript_length',
                          'gene_biotype', 'transcript_biotype', ]
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)

        self.df.rename(columns=self.COLUMNS_RENAME_DICT,
                       inplace=True)
        print(self.name(), self.df.columns.tolist())

    def load_dataframe(self, datasets, attributes, host, filename=None, ):
        return self.retrieve_dataset(host, datasets, attributes, filename)

    def get_rename_dict(self, from_index="gene_id", to_index="gene_name"):
        geneid_to_genename = self.df[self.df[to_index].notnull()]\
            .groupby(from_index)[to_index] \
            .apply(concat_uniques).to_dict()
        return geneid_to_genename

    def get_functional_annotations(self, omic, index):
        geneid_to_go = self.df[self.df["go_id"].notnull()]\
            .groupby(index)["go_id"]\
            .apply(lambda x: "|".join(x.unique())).to_dict()
        return geneid_to_go


class EnsemblGeneSequences(EnsemblGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl",
                 attributes=None,
                 host="www.ensembl.org", filename=False):
        if attributes is None:
            attributes = ['ensembl_gene_id', 'gene_exon_intron', 'gene_flank', 'coding_gene_flank', 'gene_exon',
                          'coding']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, filename=self.filename, host=self.host,
                                      attributes=attributes, )
        self.df.rename(columns=self.COLUMNS_RENAME_DICT,
                       inplace=True)


class EnsemblTranscriptSequences(EnsemblGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl",
                 attributes=None,
                 host="www.ensembl.org", filename=False):
        if attributes is None:
            attributes = ['ensembl_transcript_id', 'transcript_exon_intron', 'transcript_flank',
                          'coding_transcript_flank',
                          '5utr', '3utr']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)
        self.df.rename(columns=self.COLUMNS_RENAME_DICT,
                       inplace=True)


class EnsemblSNP(EnsemblGenes):
    def __init__(self, dataset="hsapiens_snp",
                 attributes=None,
                 host="www.ensembl.org", filename=False):
        if attributes is None:
            attributes = ['synonym_name', 'variation_names', 'minor_allele',
                          'associated_variant_risk_allele',
                          'ensembl_gene_stable_id', 'ensembl_transcript_stable_id',
                          'phenotype_name',
                          'chr_name', 'chrom_start', 'chrom_end']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)


class EnsemblSomaticVariation(EnsemblGenes):
    def __init__(self, dataset="hsapiens_snp_som",
                 attributes=None,
                 host="www.ensembl.org", filename=False):
        if attributes is None:
            attributes = ['somatic_variation_name', 'somatic_source_name', 'somatic_allele', 'somatic_minor_allele',
                          'somatic_clinical_significance', 'somatic_validated', 'somatic_transcript_location',
                          'somatic_mapweight',
                          'somatic_chromosome_start', 'somatic_chromosome_end']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)


class NONCODE(Dataset):
    def __init__(self, path, file_resources=None, col_rename=None, verbose=False):
        if file_resources is None:
            file_resources = {}
            file_resources["NONCODEv5_source"] = os.path.join(path, "NONCODEv5_source")
            file_resources["NONCODEv5_Transcript2Gene"] = os.path.join(path, "NONCODEv5_Transcript2Gene")
            file_resources["NONCODEv5_human.func"] = os.path.join(path, "NONCODEv5_human.func")

        super(NONCODE, self).__init__(path, file_resources, col_rename, verbose=verbose)

    def load_dataframe(self, file_resources):
        source_df = dd.read_table(file_resources["NONCODEv5_source"], header=None)
        source_df.columns = ["NONCODE Transcript ID", "name type", "Gene ID"]

        transcript2gene_df = dd.read_table(file_resources["NONCODEv5_Transcript2Gene"], header=None)
        transcript2gene_df.columns = ["NONCODE Transcript ID", "NONCODE Gene ID"]

        self.noncode_func_df = dd.read_table(file_resources["NONCODEv5_human.func"], header=None)
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


