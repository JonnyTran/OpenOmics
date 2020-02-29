import os
from io import StringIO
from os.path import expanduser

import dask.dataframe as dd
import pandas as pd
from Bio import SeqIO
from bioservices import BioMart
from gtfparse import read_gtf

import openomics
from openomics.database.base import Dataset
from openomics.database.sequence import SequenceDataset
from openomics.utils.df import concat_uniques
from openomics.utils.io import mkdirs

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openomics")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), ".openomics", "databases")


class RNAcentral(Dataset):
    COLUMNS_RENAME_DICT = {'ensembl_gene_id': 'gene_id',
                           'gene symbol': 'gene_name',
                           'external id': 'transcript_id',
                           'GO terms': 'go_id'}

    def __init__(self, path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/",
                 file_resources=None, col_rename=COLUMNS_RENAME_DICT, npartitions=0, species=9606):
        self.species = species

        if file_resources is None:
            file_resources = {}
            file_resources["rnacentral_rfam_annotations.tsv"] = "go_annotations/rnacentral_rfam_annotations.tsv.gz"
            file_resources["gencode.tsv"] = "id_mapping/database_mappings/gencode.tsv"

        super(RNAcentral, self).__init__(path, file_resources, col_rename=col_rename, npartitions=npartitions)

    def load_dataframe(self, file_resources):
        go_terms = pd.read_table(file_resources["rnacentral_rfam_annotations.tsv"],
                                 low_memory=True, header=None, names=["RNAcentral id", "GO terms", "Rfams"])
        go_terms["RNAcentral id"] = go_terms["RNAcentral id"].str.split("_", expand=True, n=2)[0]

        gencode_id = pd.read_table(file_resources["gencode.tsv"],
                                   low_memory=True, header=None,
                                   names=["RNAcentral id", "database", "external id", "species", "RNA type",
                                          "gene symbol"])

        gencode_id["species"] = gencode_id["species"].astype("O")
        if self.species is not None:
            gencode_id = gencode_id[gencode_id["species"] == self.species]

        lnc_go_terms = go_terms[go_terms["RNAcentral id"].isin(gencode_id["RNAcentral id"])].groupby("RNAcentral id")[
            "GO terms"].apply(lambda x: "|".join(x.unique()))
        lnc_rfams = go_terms[go_terms["RNAcentral id"].isin(gencode_id["RNAcentral id"])].groupby("RNAcentral id")[
            "Rfams"].apply(lambda x: "|".join(x.unique()))

        gencode_id["GO terms"] = gencode_id["RNAcentral id"].map(lnc_go_terms)
        gencode_id["Rfams"] = gencode_id["RNAcentral id"].map(lnc_rfams)
        gencode_id = gencode_id[gencode_id["GO terms"].notnull() | gencode_id["Rfams"].notnull()]

        return gencode_id


class GTEx(Dataset):
    COLUMNS_RENAME_DICT = {
        "Name": "gene_id",
        "Description": "gene_name"
    }

    def __init__(self, path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/",
                 file_resources=None, col_rename=None, npartitions=0):
        if file_resources is None:
            file_resources = {}

            file_resources[
                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"] = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
            file_resources[
                "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"] = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
            file_resources[
                "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct"] = "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"

        super(GTEx, self).__init__(path, file_resources, col_rename=None, npartitions=npartitions)

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


class GENCODE(SequenceDataset):
    def __init__(self, path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                 file_resources=None, col_rename=None, npartitions=0, agg_sequences="all",
                 replace_U2T=True, remove_version_num=False):
        if file_resources is None:
            file_resources = {"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                              "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                              "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                              "transcripts.fa": "gencode.v32.transcripts.fa.gz"}

        self.remove_version_num = remove_version_num

        super(GENCODE, self).__init__(path=path, file_resources=file_resources, col_rename=col_rename,
                                      npartitions=npartitions,
                                      agg_sequences=agg_sequences, replace_U2T=replace_U2T)

    def load_dataframe(self, file_resources):
        dfs = []
        for gtf_file in file_resources:
            if '.gtf' in gtf_file:
                # Parse lncRNA gtf
                df = read_gtf(file_resources[gtf_file])  # Returns a dask dataframe
                dfs.append(df)
        annotation_df = pd.concat(dfs)

        if self.remove_version_num:
            annotation_df['gene_id'] = annotation_df['gene_id'].str.replace("[.].*", "")
            annotation_df['transcript_id'] = annotation_df['transcript_id'].str.replace("[.].*", "")
        return annotation_df

    def read_fasta(self, fasta_file):
        entries = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            record_dict = {"gene_id": record.id.split("|")[1],
                           "gene_name": record.id.split("|")[5],
                           "transcript_id": record.id.split("|")[0],
                           "transcript_name": record.id.split("|")[4],
                           "transcript_length": record.id.split("|")[6],
                           "transcript_biotype": record.id.split("|")[7],
                           "sequence": str(record.seq),
                           }

            entries.append(record_dict)

        entries_df = pd.DataFrame(entries)
        if self.replace_U2T:
            entries_df["sequence"] = entries_df["sequence"].replace("U", "T")
        if self.remove_version_num:
            entries_df['gene_id'] = entries_df['gene_id'].str.replace("[.].*", "")
            entries_df['transcript_id'] = entries_df['transcript_id'].str.replace("[.].*", "")
        return entries_df

    def get_sequences(self, index, omic, agg_sequences, biotypes=None):
        if agg_sequences == "all":
            agg_func = lambda x: list(x)
        elif agg_sequences == "shortest":
            agg_func = lambda x: min(x, key=len)
        elif agg_sequences == "longest":
            agg_func = lambda x: max(x, key=len)
        else:
            raise Exception("agg_sequences argument must be one of {'all', 'shortest', 'longest'}")

        # Parse lncRNA & mRNA fasta
        if omic == openomics.MessengerRNA.name():
            fasta_file = self.file_resources["transcripts.fa"]
        elif omic == openomics.LncRNA.name():
            fasta_file = self.file_resources["lncRNA_transcripts.fa"]
        else:
            raise Exception("omic argument must be one of {'MessengerRNA', 'LncRNA'}")

        entries_df = self.read_fasta(fasta_file)

        if "gene" in index:
            if biotypes:
                entries_df = entries_df[entries_df["transcript_biotype"].isin(biotypes)]
            else:
                print("INFO: You can pass in a list of transcript biotypes to filter using the argument 'biotypes'.")

            return entries_df.groupby(index)["sequence"].agg(agg_func)
        elif "transcript" in index:
            return entries_df.groupby(index)["sequence"].first()
        else:
            raise Exception(
                "The level argument must be one of {'gene_id', 'transcript_id', or 'gene_name', or 'transcript_name'}")

    def get_rename_dict(self, from_index='gene_id', to_index='gene_name'):
        ensembl_id_to_gene_name = pd.Series(self.df[to_index].values,
                                            index=self.df[from_index]).to_dict()
        return ensembl_id_to_gene_name


class MirBase(Dataset):
    def __init__(self, path="ftp://mirbase.org/pub/mirbase/CURRENT/",
                 sequence="hairpin", species="Homo sapiens", species_id=9606,
                 file_resources=None, col_rename=None,
                 npartitions=0, agg_sequences="longest", replace_U2T=True):
        """

        Args:
            path:
            file_resources:
            col_rename:
            species (int): Species code, e.g., 9606 for human
            agg_sequences (str): {"longest", "shortest", "all"}
                Whether to select the longest, shortest, or a list of all transcript sequences when aggregating transcript sequences by gene_id or gene_name.
            replace_U2T:
            :param npartitions:
        """
        if file_resources is None:
            file_resources = {}
            file_resources["aliases.txt"] = "aliases.txt.gz"
            file_resources["mature.fa"] = "mature.fa.gz"
            file_resources["hairpin.fa"] = "hairpin.fa.gz"
            file_resources["rnacentral.mirbase.tsv"] = \
                "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/mirbase.tsv"

        self.sequence = sequence
        self.agg_sequences = agg_sequences
        self.replace_U2T = replace_U2T
        self.species_id = species_id
        self.species = species
        super(MirBase, self).__init__(path, file_resources, col_rename=col_rename, npartitions=npartitions)

    def load_dataframe(self, file_resources):
        rnacentral_mirbase = pd.read_table(file_resources["rnacentral.mirbase.tsv"], low_memory=True, header=None,
                                           names=["RNAcentral id", "database", "mirbase id", "species", "RNA type",
                                                  "gene name"])
        rnacentral_mirbase = rnacentral_mirbase.set_index("mirbase id")
        rnacentral_mirbase["species"] = rnacentral_mirbase["species"].astype("O")
        if self.species_id is not None:
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species"] == self.species_id]

        mirbase_aliases = pd.read_table(file_resources["aliases.txt"], low_memory=True, header=None,
                                        names=["mirbase id", "gene_name"], dtype="O").set_index("mirbase id")
        mirbase_aliases = mirbase_aliases.join(rnacentral_mirbase, how="inner")

        # Expanding miRNA names in each MirBase Ascension ID
        mirna_names = mirbase_aliases.apply(lambda x: pd.Series(x['gene_name'].split(";")[:-1]),
                                            axis=1).stack().reset_index(
            level=1, drop=True)
        mirna_names.name = "gene_name"
        mirbase_aliases = mirbase_aliases.drop('gene_name', axis=1).join(mirna_names)

        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        return mirbase_aliases

    def get_sequences(self, index="gene_name", omic=None):
        if hasattr(self, "seq_dict"):
            return self.seq_dict

        self.seq_dict = {}
        if self.sequence == "hairpin":
            file = self.file_resources["hairpin.fa"]
        elif self.sequence == "mature":
            file = self.file_resources["mature.fa"]
        else:
            raise Exception("sequence must be either 'hairpin' or 'mature'")

        for record in SeqIO.parse(file, "fasta"):
            gene_name = str(record.name)
            sequence_str = str(record.seq)

            if self.replace_U2T:
                sequence_str = sequence_str.replace("U", "T")

            if self.agg_sequences == "shortest":
                if gene_name not in self.seq_dict:
                    self.seq_dict[gene_name] = sequence_str
                else:
                    if len(self.seq_dict[gene_name]) > len(sequence_str):
                        self.seq_dict[gene_name] = sequence_str
            elif self.agg_sequences == "longest":
                if gene_name not in self.seq_dict:
                    self.seq_dict[gene_name] = sequence_str
                else:
                    if len(self.seq_dict[gene_name]) < len(sequence_str):
                        self.seq_dict[gene_name] = sequence_str
            elif self.agg_sequences == "all":
                if gene_name not in self.seq_dict:
                    self.seq_dict[gene_name] = [sequence_str, ]
                else:
                    self.seq_dict[gene_name].append(sequence_str)
            else:
                self.seq_dict[gene_name] = sequence_str

        return self.seq_dict


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
    # TODO need more fix
    def __init__(self, path, file_resources=None, col_rename=None):
        if file_resources is None:
            file_resources = {}
            file_resources["NONCODEv5_source"] = os.path.join(path, "NONCODEv5_source")
            file_resources["NONCODEv5_Transcript2Gene"] = os.path.join(path, "NONCODEv5_Transcript2Gene")
            file_resources["NONCODEv5_human.func"] = os.path.join(path, "NONCODEv5_human.func")

        super(NONCODE, self).__init__(path, file_resources, col_rename)

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


