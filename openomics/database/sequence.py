import os
import re
from abc import abstractmethod
from collections import defaultdict
from typing import Union, List, Callable, Dict, Tuple

import pandas as pd
import tqdm
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from dask import dataframe as dd
from pyfaidx import Fasta
from six.moves import intern

import openomics
from openomics.utils.read_gtf import read_gtf
from .base import Database

SEQUENCE_COL = 'sequence'


class SequenceDatabase(Database):
    """Provides a series of methods to extract sequence data from
    SequenceDataset.
    """

    def __init__(self, **kwargs):
        """
        Args:
            **kwargs:
        """
        super().__init__(**kwargs)

    @abstractmethod
    def load_sequences(self, fasta_file: str, keys: Union[pd.Index, List[str]] = None, blocksize=None, **kwargs):
        """Returns a pandas DataFrame containing the fasta sequence entries.
        With a column named 'sequence'.

        Args:
            fasta_file (str): path to the fasta file, usually as
                self.file_resources[<file_name>]
            keys (pd.Index): list of keys to
            blocksize:
        """
        raise NotImplementedError

    @abstractmethod
    def get_sequences(self, index_name: str, omic: str, agg: str, **kwargs):
        """Returns a dictionary where keys are 'index' and values are
        sequence(s).

        Args:
            index_name (str): {"gene_id", "gene_name", "transcript_id",
                "transcript_name"}
            omic (str): {"lncRNA", "microRNA", "messengerRNA"}
            agg (str): {"all", "shortest", "longest"}
            **kwargs: any additional argument to pass to
                SequenceDataset.get_sequences()
        """
        raise NotImplementedError

    @staticmethod
    def aggregator_fn(agg: Union[str, Callable] = None) -> Callable:
        """Returns a function used aggregate a list of sequences from a groupby
        on a given key.

        Args:
            agg: One of ("all", "shortest", "longest", "first"), default "all". If "all",
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
        blocksize=0,
        replace_U2T=False,
        remove_version_num=False,
    ):
        """
        Args:
            path:
            file_resources:
            col_rename:
            blocksize:
            replace_U2T (bool): Whether to replace nucleotides from U to T on
                the RNA primary sequences.
            remove_version_num (bool): Whether to drop the version number on the
                ensembl ID.
        """
        if file_resources is None:
            file_resources = {
                "basic.annotation.gtf.gz": "gencode.v32.basic.annotation.gtf.gz",
                "long_noncoding_RNAs.gtf.gz": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                "lncRNA_transcripts.fa.gz": "gencode.v32.lncRNA_transcripts.fa.gz",
                "transcripts.fa.gz": "gencode.v32.transcripts.fa.gz",
            }

        self.remove_version_num = remove_version_num

        super().__init__(path=path, file_resources=file_resources, col_rename=col_rename, blocksize=blocksize)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        dfs = []
        if blocksize:
            for filename in file_resources:
                if filename.endswith(".gtf.gz"):
                    df = read_gtf(file_resources[filename], blocksize=blocksize, compression="gzip")
                    dfs.append(df)
            annotation_df = dd.concat(dfs)
        else:
            for filename in file_resources:
                if filename.endswith(".gtf"):
                    df = read_gtf(file_resources[filename])
                    dfs.append(df)

            annotation_df = pd.concat(dfs)

        if self.remove_version_num:
            annotation_df["gene_id"] = annotation_df["gene_id"].str.replace("[.]\d*", "", regex=True)
            annotation_df["transcript_id"] = annotation_df["transcript_id"].str.replace("[.]\d*", "", regex=True)

        return annotation_df

    def load_sequences(self, fasta_file: str, keys: pd.Index = None, blocksize=None, **kwargs):
        """
        Args:
            keys ():
            fasta_file:
            replace_U2T:
            blocksize:
        """
        if hasattr(self, '_seq_df_dict') and fasta_file in self._seq_df_dict:
            return self._seq_df_dict[fasta_file]

        def get_transcript_id(x):
            key = x.split('|')[0]  # transcript_id
            if self.remove_version_num:
                return re.sub("[.]\d*", "", key)
            else:
                return key

        fa = Fasta(fasta_file, key_function=get_transcript_id, as_raw=True)

        entries = []
        for key, record in fa.items():
            if isinstance(keys, (set, list, pd.Index, pd.Series)) and key not in keys: continue

            record_dict = {
                "transcript_id": record.long_name.split("|")[0],
                "gene_id": record.long_name.split("|")[1],
                "gene_name": record.long_name.split("|")[5],
                "transcript_name": record.long_name.split("|")[4],
                "transcript_length": record.long_name.split("|")[6],
                "transcript_biotype": intern(record.long_name.split("|")[7]),
                SEQUENCE_COL: str(record),
            }

            entries.append(record_dict)

        seq_df = pd.DataFrame(entries)
        if blocksize:
            seq_df = dd.from_pandas(seq_df, chunksize=blocksize)

        if self.remove_version_num:
            seq_df["gene_id"] = seq_df["gene_id"].str.replace("[.]\d*", "", regex=True)
            seq_df["transcript_id"] = seq_df["transcript_id"].str.replace("[.]\d*", "", regex=True)

        # Cache the seq_df
        if not hasattr(self, '_seq_df_dict'):
            self._seq_df_dict = {}
        self._seq_df_dict[fasta_file] = seq_df

        return seq_df

    def get_sequences(self, index: Union[str, Tuple[str]], omic: str, agg: str = 'all', biotypes: List[str] = None):
        """
        Args:
            index (str):
            omic (str):
            agg (str):
            biotypes (List[str]):
        """
        agg_func = self.aggregator_fn(agg)

        # Parse lncRNA & mRNA fasta
        if omic == openomics.MessengerRNA.name():
            fasta_file = self.file_resources["transcripts.fa"]
        elif omic == openomics.LncRNA.name():
            fasta_file = self.file_resources["lncRNA_transcripts.fa"]
        else:
            raise Exception("omic argument must be one of {'MessengerRNA', 'LncRNA'}")
        assert isinstance(fasta_file,
                          str), f"Fasta file provided in `file_resources` must be a non-compressed .fa file. Given {fasta_file}"

        seq_df = self.load_sequences(fasta_file)

        if "gene" in index:
            if biotypes:
                seq_df = seq_df[seq_df["transcript_biotype"].isin(biotypes)]
            else:
                print("INFO: You can pass in a list of transcript biotypes to filter using the argument 'biotypes'.")

            return seq_df.groupby(index)[SEQUENCE_COL].agg(agg_func)

        else:
            return seq_df.groupby(index)[SEQUENCE_COL].first()

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

    SPECIES_ID_NAME = {
        '10090': 'MOUSE', '10116': 'RAT', '226900': 'BACCR', '243273': 'MYCGE', '284812': 'SCHPO', '287': 'PSEAI',
        '3702': 'ARATH', '99287': 'SALTY', '44689': 'DICDI', '4577': 'MAIZE', '559292': 'YEAST', '6239': 'CAEEL',
        '7227': 'DROME', '7955': 'DANRE', '83333': 'ECOLI', '9606': 'HUMAN', '9823': 'PIG', }

    SPECIES_ID_TAXONOMIC = {
        'HUMAN': 'human', 'MOUSE': 'rodents', 'RAT': 'rodents', 'BACCR': 'bacteria', 'MYCGE': 'bacteria',
        'SCHPO': 'fungi', 'PSEAI': 'bacteria', 'ARATH': 'plants', 'SALTY': 'bacteria', 'DICDI': 'bacteria',
        'MAIZE': 'plants', 'YEAST': 'fungi', 'CAEEL': 'vertebrates', 'DROME': 'invertebrates', 'DANRE': 'vertebrates',
        'ECOLI': 'bacteria', 'PIG': 'mammals',
    }

    def __init__(self, path="https://ftp.uniprot.org/pub/databases/uniprot/current_release/",
                 species_id: str = "9606",
                 file_resources: Dict[str, str] = None, col_rename=COLUMNS_RENAME_DICT, verbose=False,
                 blocksize=None):
        """
        Loads the UniProt database from https://uniprot.org/ .

        Default path: 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/'
        Default file_resources: {
            file_resources['uniprot_sprot.xml.gz'] = "knowledgebase/complete/uniprot_sprot.xml.gz
            file_resources['uniprot_trembl.xml.gz'] = "knowledgebase/complete/uniprot_trembl.xml.gz
            file_resources["idmapping_selected.tab.gz"] = "knowledgebase/idmapping/idmapping_selected.tab.gz'
        }

        Args:
            path:
            file_resources:
            col_rename:
            verbose:
            blocksize:
        """
        self.species_id = species_id
        self.species = UniProt.SPECIES_ID_NAME[species_id] if isinstance(species_id, str) else None
        self.taxonomic_id = UniProt.SPECIES_ID_TAXONOMIC[self.species] if isinstance(self.species, str) else None

        if file_resources is None:
            file_resources = {}

            file_resources['uniprot_sprot.xml.gz'] = os.path.join(path, "knowledgebase/complete/uniprot_sprot.xml.gz")
            file_resources['uniprot_trembl.xml.gz'] = os.path.join(path, "knowledgebase/complete/uniprot_trembl.xml.gz")
            file_resources["idmapping_selected.tab.gz"] = os.path.join(path, "knowledgebase/idmapping/",
                                                                       'idmapping_selected.tab.gz')

            if self.species:
                file_resources['uniprot_sprot.xml.gz'] = os.path.join(
                    path, "knowledgebase/taxonomic_divisions/", f'uniprot_sprot_{self.taxonomic_id}.xml.gz')
                file_resources['uniprot_trembl.xml.gz'] = os.path.join(
                    path, "knowledgebase/taxonomic_divisions/", f'uniprot_trembl_{self.taxonomic_id}.xml.gz')
                file_resources["idmapping_selected.tab.gz"] = os.path.join(
                    path, "knowledgebase/idmapping/by_organism/",
                    f'{self.species}_{self.species_id}_idmapping_selected.tab.gz')

        if 'proteomes.tsv' not in file_resources:
            file_resources["proteomes.tsv"] = \
                "https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Corganism_id&format=tsv&query=%28%2A%29%20AND%20%28proteome_type%3A1%29"

        super().__init__(path=path, file_resources=file_resources, col_rename=col_rename, verbose=verbose,
                         blocksize=blocksize)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        # Load idmapping_selected
        options = dict(
            names=['UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100',
                   'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL',
                   'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed'],
            usecols=['UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO',
                     'NCBI-taxon', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO'],
            dtype='str')

        if blocksize:
            if "idmapping_selected.tab.gz" not in file_resources:
                idmapping = dd.read_table(file_resources["idmapping_selected.tab"], **options)
            else:
                idmapping: dd.DataFrame = dd.read_table(file_resources["idmapping_selected.tab.gz"], compression="gzip",
                                                        **options)
        else:
            idmapping: pd.DataFrame = pd.read_table(file_resources["idmapping_selected.tab"], index_col='UniProtKB-AC',
                                                    **options)


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
                    lambda x: [".".join([x['NCBI-taxon'], protein_id]) \
                               for protein_id in x['protein_external_id']] \
                        if isinstance(x['protein_external_id'], list) else None,
                    meta=pd.Series([['', '']]),
                    axis=1)

        # Load proteome
        proteomes = pd.read_table(file_resources["proteomes.tsv"],
                                  usecols=['Organism Id', 'Proteome Id'],
                                  dtype={'Organism Id': 'str', 'Proteome Id': 'str'}) \
            .rename(columns={'Organism Id': 'NCBI-taxon', 'Proteome Id': 'proteome_id'}) \
            .dropna()
        idmapping = idmapping.join(proteomes.set_index('NCBI-taxon'), on='NCBI-taxon')

        return idmapping

    def get_sequences(self, index_name: str, omic: str = None, agg: str = "all", **kwargs):
        agg_func = self.aggregator_fn(agg)

        # Parse lncRNA & mRNA fasta
        seq_df = self.load_sequences(self.file_resources["uniprot_sprot.xml"], blocksize=self.blocksize)

        if "uniprot_trembl.xml" in self.file_resources:
            trembl_seq_df = self.load_sequences(self.file_resources["uniprot_trembl.xml"], blocksize=self.blocksize)
            seq_df.update(trembl_seq_df)

        if "gene" in index_name:
            return seq_df.groupby(index_name)[SEQUENCE_COL].agg(agg_func)
        else:
            return seq_df.groupby(index_name)[SEQUENCE_COL].first()

    def load_sequences(self, fasta_xml_file: str, keys=None, blocksize=None) -> pd.DataFrame:
        if hasattr(self, '_seq_df_dict') and fasta_xml_file in self._seq_df_dict:
            return self._seq_df_dict[fasta_xml_file]

        records = []
        seqfeats = []
        for record in tqdm.tqdm(SeqIO.parse(fasta_xml_file, "uniprot-xml")):
            # Sequence features
            annotations = defaultdict(lambda: None, record.annotations)
            record_dict = {
                'protein_id': record.id,
                "protein_name": record.name,
                'gene_name': annotations['gene_name_primary'],
                'description': record.description,
                'molecule_type': annotations['molecule_type'],
                'created': annotations['created'],
                'ec_id': annotations['type'],
                'subcellular_location': annotations['comment_subcellularlocation_location'],
                'taxonomy': annotations['taxonomy'],
                'keywords': annotations['keywords'],
                'sequence_mass': annotations['sequence_mass'],
                SEQUENCE_COL: str(record.seq),
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
                             **features_dict})

        records_df = pd.DataFrame(records) if not blocksize else dd.from_pandas(records, chunksize=blocksize)
        records_df = records_df.set_index(['protein_id'])

        seqfeats_df = pd.DataFrame(seqfeats) if not blocksize else dd.from_pandas(seqfeats, chunksize=blocksize)
        seqfeats_df = seqfeats_df.set_index(['protein_id'])
        seqfeats_df.columns = [f"seq/{col}" for col in seqfeats_df.columns]

        # Join new metadata to self.data
        if SEQUENCE_COL not in self.data.columns:
            exclude_cols = records_df.columns.intersection(self.data.columns)
            self.data = self.data.join(records_df.drop(columns=exclude_cols, errors="ignore"),
                                       on='protein_id', how="left")
        else:
            self.data.update(records_df, overwrite=False)

        # Add new seq features
        if len(seqfeats_df.columns.difference(self.data.columns)):
            self.data = self.data.join(seqfeats_df.drop(columns=seqfeats_df.columns.intersection(self.data.columns)),
                                       on='protein_id', how="left")
        # Fillna seq features
        if len(seqfeats_df.columns.intersection(self.data.columns)):
            self.data.update(seqfeats_df.filter(seqfeats_df.columns.intersection(self.data.columns), axis='columns'),
                             overwrite=False)

        # Cache the seq_df
        if not hasattr(self, '_seq_df_dict'):
            self._seq_df_dict = {}
        self._seq_df_dict[fasta_xml_file] = records_df

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
        species_id: str = '9606',
        sequence: str = "mature",
        col_rename=None,
        blocksize=None,
    ):
        """
        Args:
            path:
            file_resources:
            sequence (str):
            species_id (str): Species code, e.g., 9606 for human
            col_rename:
            blocksize:
        """
        if file_resources is None:
            file_resources = {}
            file_resources["aliases.txt.gz"] = "aliases.txt.gz"
            file_resources["mature.fa.gz"] = "mature.fa.gz"
            file_resources["hairpin.fa.gz"] = "hairpin.fa.gz"
            file_resources["rnacentral.mirbase.tsv"] = "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/" \
                                                       "id_mapping/database_mappings/mirbase.tsv"

        self.sequence = sequence
        self.species_id = species_id
        super().__init__(path=path, file_resources=file_resources, col_rename=col_rename, blocksize=blocksize)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        rnacentral_mirbase = pd.read_table(file_resources["rnacentral.mirbase.tsv"], low_memory=True, header=None,
                                           names=["RNAcentral id", "database", "mirbase id", "species_id", "RNA type",
                                                  "NA", ],
                                           dtype={"species_id": "str"})

        rnacentral_mirbase = rnacentral_mirbase.set_index("mirbase id")
        if self.species_id is not None:
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species_id"] == self.species_id]

        mirbase_aliases = pd.read_table(file_resources["aliases.txt"], low_memory=True, header=None,
                                        names=["mirbase id", "gene_name"], dtype="O", ).set_index("mirbase id")
        mirbase_aliases = mirbase_aliases.join(rnacentral_mirbase, how="inner")

        # Expanding miRNA names in each MirBase Ascension ID
        mirna_names = mirbase_aliases \
            .apply(lambda x: pd.Series(x["gene_name"].split(";")[:-1]), axis=1) \
            .stack() \
            .reset_index(level=1, drop=True)
        mirna_names.name = "gene_name"

        if blocksize:
            mirbase_aliases = dd.from_pandas(mirbase_aliases, chunksize=blocksize)

        mirbase_aliases = mirbase_aliases.drop("gene_name", axis=1).join(mirna_names)

        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        return mirbase_aliases

    def load_sequences(self, fasta_file, keys=None, blocksize=None):
        """
        Args:
            keys ():
            fasta_file:
            blocksize:
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
                SEQUENCE_COL: str(record.seq),
            }

            entries.append(record_dict)

        entries_df = pd.DataFrame(entries)
        if blocksize:
            entries_df = dd.from_pandas(entries_df, chunksize=blocksize)

        return entries_df

    def get_sequences(self,
                      index_name="gene_name",
                      omic=None,
                      agg="all",
                      **kwargs):
        """
        Args:
            index_name:
            omic:
            agg:
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

        seq_df = self.load_sequences(file)

        self._seq_dict = seq_df.groupby(index_name)[SEQUENCE_COL].agg(
            self.aggregator_fn(agg))

        return self._seq_dict
