import os
from abc import abstractmethod
from collections import defaultdict
from typing import Union, List, Callable

import pandas as pd
import tqdm
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from dask import dataframe as dd

import openomics
from openomics.utils.read_gtf import read_gtf
from .base import Database


# from gtfparse import read_gtf


class SequenceDatabase(Database):
    """Provides a series of methods to extract sequence data from
    SequenceDataset.
    """
    def __init__(self, replace_U2T=False, **kwargs):
        """
        Args:
            replace_U2T:
            **kwargs:
        """
        self.replace_U2T = replace_U2T

        super().__init__(**kwargs)

    @abstractmethod
    def read_fasta(self, fasta_file:str, replace_U2T:bool, npartitions=None):
        """Returns a pandas DataFrame containing the fasta sequence entries.
        With a column named 'sequence'.

        Args:
            fasta_file (str): path to the fasta file, usually as
                self.file_resources[<file_name>]
            replace_U2T (bool):
            npartitions:
        """
        raise NotImplementedError

    @abstractmethod
    def get_sequences(self, index:str, omic:str, agg_sequences:str, **kwargs):
        """Returns a dictionary where keys are 'index' and values are
        sequence(s).

        Args:
            index (str): {"gene_id", "gene_name", "transcript_id",
                "transcript_name"}
            omic (str): {"lncRNA", "microRNA", "messengerRNA"}
            agg_sequences (str): {"all", "shortest", "longest"}
            **kwargs: any additional argument to pass to
                SequenceDataset.get_sequences()
        """
        raise NotImplementedError

    @staticmethod
    def aggregator_fn(agg: Union[str, Callable] = None) -> Callable:
        """Returns a function used aggregate a list of sequences from a groupby
        on a given key.

        Args:
            agg: One of ("all", "shortest", "longest"), default "all". If "all",
                then return a list of sequences.
        """
        if agg == "all":
            agg_func = lambda x: list(x) if not isinstance(x, str) else x
        elif agg == "shortest":
            agg_func = lambda x: min(x, key=len) if isinstance(x, list) else x
        elif agg == "longest":
            agg_func = lambda x: max(x, key=len) if isinstance(x, list) else x
        elif agg == 'first':
            agg_func = lambda x: x[0] if isinstance(x, list) else x
        elif callable(agg):
            return agg
        else:
            raise Exception(
                "agg_sequences argument must be one of {'all', 'shortest', 'longest'}"
            )
        return agg_func

class GENCODE(SequenceDatabase):
    """Loads the GENCODE database from https://www.gencodegenes.org/ .

    Default path: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/ .
    Default file_resources: {
        "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
        "long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
        "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
        "transcripts.fa": "gencode.v32.transcripts.fa.gz",
    }
    """
    def __init__(
        self,
        path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
        file_resources=None,
        col_rename=None,
        npartitions=0,
        replace_U2T=False,
        remove_version_num=False,
    ):
        """
        Args:
            path:
            file_resources:
            col_rename:
            npartitions:
            replace_U2T (bool): Whether to replace nucleotides from U to T on
                the RNA primary sequences.
            remove_version_num (bool): Whether to drop the version number on the
                ensembl ID.
        """
        if file_resources is None:
            file_resources = {
                "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                "long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                "transcripts.fa": "gencode.v32.transcripts.fa.gz",
            }

        self.remove_version_num = remove_version_num

        super().__init__(
            path=path,
            file_resources=file_resources,
            col_rename=col_rename,
            replace_U2T=replace_U2T,
            npartitions=npartitions,
        )

    def load_dataframe(self, file_resources, npartitions=None):
        """
        Args:
            file_resources:
            npartitions:
        """
        dfs = []
        for filename, content in file_resources.items():
            if ".gtf" in filename:
                df = read_gtf(content,
                              npartitions=npartitions,
                              compression="gzip")
                dfs.append(df)

        if npartitions:
            annotation_df = dd.concat(dfs)
        else:
            annotation_df = pd.concat(dfs)

        if self.remove_version_num:
            annotation_df["gene_id"] = annotation_df["gene_id"].str.replace("[.]\d*", "", regex=True)
            annotation_df["transcript_id"] = annotation_df["transcript_id"].str.replace("[.]\d*", "", regex=True)

        return annotation_df

    def read_fasta(self, fasta_file, npartitions=None):
        """
        Args:
            fasta_file:
            replace_U2T:
            npartitions:
        """
        if hasattr(self, '_seq_df_dict') and fasta_file in self._seq_df_dict:
            return self._seq_df_dict[fasta_file]

        entries = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            record_dict = {
                "gene_id": record.id.split("|")[1],
                "gene_name": record.id.split("|")[5],
                "transcript_id": record.id.split("|")[0],
                "transcript_name": record.id.split("|")[4],
                "transcript_length": record.id.split("|")[6],
                "transcript_biotype": record.id.split("|")[7],
                "sequence": str(record.seq),
            }

            entries.append(record_dict)

        seq_df = pd.DataFrame(entries)
        if npartitions:
            seq_df = dd.from_pandas(seq_df)

        if self.replace_U2T:
            seq_df["sequence"] = seq_df["sequence"].str.replace("U", "T")

        if self.remove_version_num:
            seq_df["gene_id"] = seq_df["gene_id"].str.replace("[.]\d*", "", regex=True)
            seq_df["transcript_id"] = seq_df["transcript_id"].str.replace("[.]\d*", "", regex=True)

        # Cache the seq_df
        if not hasattr(self, '_seq_df_dict'):
            self._seq_df_dict = {}
        self._seq_df_dict[fasta_file] = seq_df

        return seq_df

    def get_sequences(self, index: Union[str, List[str]], omic: str, agg_sequences: str, biotypes: List[str] = None):
        """
        Args:
            index (str):
            omic (str):
            agg_sequences (str):
            biotypes (List[str]):
        """
        agg_func = self.aggregator_fn(agg_sequences)

        # Parse lncRNA & mRNA fasta
        if omic == openomics.MessengerRNA.name():
            fasta_file = self.file_resources["transcripts.fa"]
        elif omic == openomics.LncRNA.name():
            fasta_file = self.file_resources["lncRNA_transcripts.fa"]
        else:
            raise Exception("omic argument must be one of {'MessengerRNA', 'LncRNA'}")

        seq_df = self.read_fasta(fasta_file)

        if "gene" in index:
            if biotypes:
                seq_df = seq_df[seq_df["transcript_biotype"].isin(biotypes)]
            else:
                print("INFO: You can pass in a list of transcript biotypes to filter using the argument 'biotypes'.")

            return seq_df.groupby(index)["sequence"].agg(agg_func)

        else:
            return seq_df.groupby(index)["sequence"].first()

    def get_rename_dict(self, from_index="gene_id", to_index="gene_name"):
        """
        Args:
            from_index:
            to_index:
        """
        ensembl_id_to_gene_name = pd.Series(
            self.data[to_index].values, index=self.data[from_index]).to_dict()
        return ensembl_id_to_gene_name


class UniProt(SequenceDatabase):
    COLUMNS_RENAME_DICT = {
        "UniProtKB-AC": 'protein_id',
        "UniProtKB-ID": 'protein_name',
        "Ensembl": "gene_id",
        "Ensembl_TRS": "transcript_id",
        "Ensembl_PRO": "protein_embl_id",
        "NCBI-taxon": "species_id",
        "GeneID (EntrezGene)": "entrezgene_id",
        "GO": "go_id",
    }

    def __init__(self, path="https://ftp.uniprot.org/pub/databases/uniprot/current_release/",
                 species="HUMAN", species_id="9606",
                 file_resources=None, col_rename=COLUMNS_RENAME_DICT, verbose=False,
                 npartitions=None):
        """
        Args:
            path:
            file_resources:
            col_rename:
            verbose:
            npartitions:
        """
        self.species = species
        self.species_id = species_id
        if file_resources is None:
            file_resources = {}
            file_resources['uniprot_sprot.xml'] = os.path.join(path, "knowledgebase/uniprot_sprot.xml.gz")
            file_resources["idmapping_selected.tab"] = os.path.join(path, "knowledgebase/idmapping/",
                                                                    'idmapping_selected.tab.gz')

        if species:
            file_resources['uniprot_sprot.xml'] = os.path.join(path, "knowledgebase/taxonomic_divisions/",
                                                               f'uniprot_sprot_{species.lower()}.xml.gz')
            file_resources["idmapping_selected.tab"] = os.path.join(path, "knowledgebase/idmapping/by_organism/",
                                                                    f'{species}_{species_id}_idmapping_selected.tab.gz')

        super().__init__(path=path, file_resources=file_resources, col_rename=col_rename, verbose=verbose,
                         npartitions=npartitions)

    def load_dataframe(self, file_resources, npartitions=None):
        """
        Args:
            file_resources:
            npartitions:
        """

        options = dict(
            names=['UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100',
                   'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL',
                   'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed'],
            usecols=['UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO',
                     'NCBI-taxon', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO'],
            dtype={'GeneID (EntrezGene)': 'str', 'NCBI-taxon': 'str'})

        if npartitions:
            idmapping: dd.DataFrame = dd.read_table(file_resources["idmapping_selected.tab"], **options)
        else:
            idmapping: pd.DataFrame = pd.read_table(file_resources["idmapping_selected.tab"], **options)

        for col in ['PDB', 'GI', 'GO', 'RefSeq']:
            # Split string to list
            idmapping[col] = idmapping[col].str.split("; ")

        for col in ['Ensembl', 'Ensembl_TRS', 'Ensembl_PRO']:
            # Removing .# ENGS gene version number at the end
            idmapping[col] = idmapping[col].str.replace("[.]\d*", "", regex=True)
            idmapping[col] = idmapping[col].str.split("; ")

            if col == 'Ensembl_PRO':
                # Prepend species_id to ensembl protein ids to match with STRING PPI
                idmapping['protein_external_id'] = idmapping[col]
                idmapping['protein_external_id'] = idmapping[["NCBI-taxon", 'protein_external_id']].apply(
                    lambda x: [".".join([x['NCBI-taxon'], protein_id]) for protein_id in x['protein_external_id']] \
                        if isinstance(x['protein_external_id'], list) else None,
                    axis=1)

        return idmapping

    def get_sequences(self, index: str, omic: str = None, agg: str = "all", **kwargs):
        agg_func = self.aggregator_fn(agg)

        # Parse lncRNA & mRNA fasta
        seq_df = self.read_fasta(self.file_resources["uniprot_sprot.xml"], npartitions=self.npartitions)

        if "protein" in index:
            return seq_df.groupby(index)["sequence"].agg(agg_func)
        else:
            return seq_df.groupby(index)["sequence"].first()

    def read_fasta(self, fasta_file: str, replace_U2T=False, npartitions=None) -> pd.DataFrame:
        if hasattr(self, '_seq_df_dict') and fasta_file in self._seq_df_dict:
            return self._seq_df_dict[fasta_file]

        records = []
        seqfeats = []
        for record in tqdm.tqdm(SeqIO.parse(fasta_file, "uniprot-xml")):
            # Sequence features
            annotations = defaultdict(lambda: None, record.annotations)
            record_dict = {
                'protein_id': record.id,
                "protein_name": record.name,
                'description': record.description,
                'molecule_type': annotations['molecule_type'],
                'gene_name': annotations['gene_name_primary'],
                'created': annotations['created'],
                'ec_id': annotations['type'],
                'subcellular_location': annotations['comment_subcellularlocation_location'],
                'taxonomy': annotations['taxonomy'],
                'keywords': annotations['keywords'],
                'sequence_mass': annotations['sequence_mass'],
                "sequence": str(record.seq),
            }
            records.append(record_dict)

            # Sequence interval features
            _parse_interval = lambda sf: pd.Interval(left=sf.location.start, right=sf.location.end, )

            feature_type_intervals = defaultdict(lambda: [])
            for sf in record.features:
                if isinstance(sf.location.start, ExactPosition) and isinstance(sf.location.end, ExactPosition):
                    feature_type_intervals[sf.type].append(_parse_interval(sf))

            features_dict = {type: pd.IntervalIndex(intervals, name=type) \
                             for type, intervals in feature_type_intervals.items()}
            seqfeats.append({"protein_id": record.id,
                             "protein_name": record.name,
                             **features_dict})

        records_df = pd.DataFrame(records) if not npartitions else dd.from_pandas(records)
        seqfeats_df = pd.DataFrame(seqfeats) if not npartitions else dd.from_pandas(seqfeats)

        # Join new metadata to self.data
        if len(records_df.columns.intersection(self.data.columns)) <= 2:
            exclude_cols = records_df.columns.intersection(self.data.columns)
            self.data = self.data.join(
                records_df.set_index('protein_id').drop(columns=exclude_cols, errors="ignore"), on='protein_id')
            self.seq_feats = seqfeats_df.set_index(['protein_id', 'protein_name'])

        # Cache the seq_df
        if not hasattr(self, '_seq_df_dict'):
            self._seq_df_dict = {}
        self._seq_df_dict[fasta_file] = records_df

        return records_df


class MirBase(SequenceDatabase):
    """Loads the MirBase database from https://mirbase.org .

    Default path: "ftp://mirbase.org/pub/mirbase/CURRENT/" .
    Default file_resources: {
        "aliases.txt": "aliases.txt.gz",
        "mature.fa": "mature.fa.gz",
        "hairpin.fa": "hairpin.fa.gz",
        "rnacentral.mirbase.tsv": "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/mirbase.tsv",
    }
    """

    def __init__(
        self,
        path="ftp://mirbase.org/pub/mirbase/CURRENT/",
        file_resources=None,
        sequence: str = "mature",
        species: str = "Homo sapiens",
        species_id: str = 9606,
        col_rename=None,
        npartitions=0,
        replace_U2T=False,
    ):
        """
        Args:
            path:
            file_resources:
            sequence (str):
            species (str): Species code, e.g., 9606 for human
            species_id (str):
            col_rename:
            npartitions:
            replace_U2T:
        """
        if file_resources is None:
            file_resources = {}
            file_resources["aliases.txt"] = "aliases.txt.gz"
            file_resources["mature.fa"] = "mature.fa.gz"
            file_resources["hairpin.fa"] = "hairpin.fa.gz"
            file_resources[
                "rnacentral.mirbase.tsv"] = "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/mirbase.tsv"

        self.sequence = sequence
        self.species_id = species_id
        self.species = species
        super().__init__(
            path=path,
            file_resources=file_resources,
            col_rename=col_rename,
            npartitions=npartitions,
            replace_U2T=replace_U2T,
        )

    def load_dataframe(self, file_resources, npartitions=None):
        """
        Args:
            file_resources:
            npartitions:
        """
        rnacentral_mirbase = pd.read_table(
            file_resources["rnacentral.mirbase.tsv"],
            low_memory=True,
            header=None,
            names=[
                "RNAcentral id",
                "database",
                "mirbase id",
                "species",
                "RNA type",
                "NA",
            ],
        )

        rnacentral_mirbase = rnacentral_mirbase.set_index("mirbase id")
        rnacentral_mirbase["species"] = rnacentral_mirbase["species"].astype(
            "O")
        if self.species_id is not None:
            rnacentral_mirbase = rnacentral_mirbase[
                rnacentral_mirbase["species"] == self.species_id]

        mirbase_aliases = pd.read_table(
            file_resources["aliases.txt"],
            low_memory=True,
            header=None,
            names=["mirbase id", "gene_name"],
            dtype="O",
        ).set_index("mirbase id")
        mirbase_aliases = mirbase_aliases.join(rnacentral_mirbase, how="inner")

        # Expanding miRNA names in each MirBase Ascension ID
        mirna_names = (mirbase_aliases.apply(
            lambda x: pd.Series(x["gene_name"].split(";")[:-1]),
            axis=1).stack().reset_index(level=1, drop=True))

        mirna_names.name = "gene_name"

        if npartitions:
            mirbase_aliases = dd.from_pandas(mirbase_aliases,
                                             npartitions=npartitions)

        mirbase_aliases = mirbase_aliases.drop("gene_name",
                                               axis=1).join(mirna_names)

        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        return mirbase_aliases

    def read_fasta(self, fasta_file, npartitions=None):
        """
        Args:
            fasta_file:
            npartitions:
        """
        entries = []
        for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
            if i == 0: print(record)
            record_dict = {
                "gene_id": record.id,
                "gene_name": str(record.name).lower(),
                "mirbase_id": record.description.split(" ")[1],
                # "mir_name": record.description.split(" ")[5],
                "species": " ".join(record.description.split(" ")[2:4]),
                "sequence": str(record.seq),
            }

            entries.append(record_dict)

        entries_df = pd.DataFrame(entries)
        if npartitions:
            entries_df = dd.from_pandas(entries_df)

        if self.replace_U2T:
            entries_df["sequence"] = entries_df["sequence"].str.replace("U", "T")

        return entries_df

    def get_sequences(self,
                      index="gene_name",
                      omic=None,
                      agg_sequences="all",
                      **kwargs):
        """
        Args:
            index:
            omic:
            agg_sequences:
            **kwargs:
        """
        if hasattr(self, "_seq_dict"):
            print("Using cached sequences dict")
            return self._seq_dict

        if self.sequence == "hairpin":
            file = self.file_resources["hairpin.fa"]
        elif self.sequence == "mature":
            file = self.file_resources["mature.fa"]
        else:
            raise Exception("sequence must be either 'hairpin' or 'mature'")

        seq_df = self.read_fasta(file)

        self._seq_dict = seq_df.groupby(index)["sequence"].agg(
            self.aggregator_fn(agg_sequences))

        return self._seq_dict
