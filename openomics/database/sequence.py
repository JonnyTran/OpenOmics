import os
import re
from abc import abstractmethod
from collections import defaultdict, OrderedDict
from collections.abc import Iterable
from typing import Union, List, Callable, Dict, Tuple

import numpy as np
import pandas as pd
import tqdm
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from dask import dataframe as dd
from pyfaidx import Fasta
from six.moves import intern

import openomics
from openomics.io.read_gtf import read_gtf
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
    def load_sequences(self, fasta_file: str, index=None, keys: Union[pd.Index, List[str]] = None, blocksize=None) \
        -> Union[pd.DataFrame, Dict]:
        """Returns a pandas DataFrame containing the fasta sequence entries.
        With a column named 'sequence'.

        Args:
            index ():
            fasta_file (str): path to the fasta file, usually as
                self.file_resources[<file_name>]
            keys (pd.Index): list of keys to
            blocksize:
        """
        raise NotImplementedError

    @abstractmethod
    def get_sequences(self, index_name: str, omic: str, agg: str, **kwargs) -> Union[pd.Series, Dict]:
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
        remove_version_num=False,
    ):
        """
        Args:
            path:
            file_resources:
            col_rename:
            blocksize:
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
                if '.gtf' in filename and isinstance(file_resources[filename], str):
                    df = read_gtf(file_resources[filename], blocksize=blocksize,
                                  compression="gzip" if filename.endswith(".gz") else None)
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

    def load_sequences(self, fasta_file: str, index=None, keys: pd.Index = None, blocksize=None):
        """
        Args:
            index ():
            keys ():
            fasta_file:
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
        for key, record in tqdm.tqdm(fa.items(), desc=str(fasta_file)):
            if keys is not None and key not in keys: continue

            attrs = record.long_name.split("|")
            record_dict = {
                "transcript_id": attrs[0],
                "gene_id": attrs[1],
                "gene_name": attrs[5],
                "transcript_name": attrs[4],
                "transcript_length": attrs[6],
                "transcript_biotype": intern(attrs[7]),
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
        if keys is not None:
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

        assert isinstance(fasta_file, str), \
            f"Fasta file provided in `file_resources` must be a non-compressed .fa file. Given {fasta_file}"

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
        # FASTA headers
        "OS": 'species', "OX": 'species_id', 'GN': 'gene_name', 'PE': 'ProteinExistence', 'SV': "version",
        # UniProt XML headers
        "accession": "UniProtKB-AC", "name": "protein_name", "gene": "gene_name",
        "geneLocation": "subcellular_location",
        "keyword": "keywords",
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
                 file_resources: Dict[str, str] = None,
                 species_id: str = "9606", remove_version_num=True,
                 index='UniProtKB-AC', keys=None,
                 col_rename=COLUMNS_RENAME_DICT,
                 **kwargs):
        """
        Loads the UniProt database from https://uniprot.org/ .

        Default path: 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/'
        Default file_resources: {
            file_resources['uniprot_sprot.xml.gz'] = "knowledgebase/complete/uniprot_sprot.xml.gz
            file_resources['uniprot_trembl.xml.gz'] = "knowledgebase/complete/uniprot_trembl.xml.gz
            file_resources["idmapping_selected.tab.gz"] = "knowledgebase/idmapping/idmapping_selected.tab.gz'
            file_resources["proteomes.tsv"] = "https://rest.uniprot.org/proteomes/stream?compressed=true&
                fields=upid%2Corganism%2Corganism_id&format=tsv&query=%28%2A%29%20AND%20%28proteome_type%3A1%29"
            file_resources['speclist.txt'] = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/
                knowledgebase/complete/docs/speclist'
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
        self.remove_version_num = remove_version_num

        if file_resources is None:
            file_resources = {}

            file_resources['uniprot_sprot.xml.gz'] = os.path.join(path, "knowledgebase/complete/uniprot_sprot.xml.gz")
            file_resources['uniprot_trembl.xml.gz'] = os.path.join(path, "knowledgebase/complete/uniprot_trembl.xml.gz")
            file_resources["idmapping_selected.tab.gz"] = os.path.join(path,
                                                                       "knowledgebase/idmapping/idmapping_selected.tab.gz")

            if self.species:
                file_resources['uniprot_sprot.xml.gz'] = os.path.join(
                    path, "knowledgebase/taxonomic_divisions/", f'uniprot_sprot_{self.taxonomic_id}.xml.gz')
                file_resources['uniprot_trembl.xml.gz'] = os.path.join(
                    path, "knowledgebase/taxonomic_divisions/", f'uniprot_trembl_{self.taxonomic_id}.xml.gz')
                file_resources["idmapping_selected.tab.gz"] = os.path.join(
                    path, "knowledgebase/idmapping/by_organism/",
                    f'{self.species}_{self.species_id}_idmapping_selected.tab.gz')

            file_resources["proteomes.tsv"] = \
                "https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Corganism_id&format=tsv&query=%28%2A%29%20AND%20%28proteome_type%3A1%29"

        super().__init__(path=path, file_resources=file_resources, index=index, keys=keys, col_rename=col_rename,
                         **kwargs)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        # Load idmapping_selected.tab
        args = dict(
            names=['UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100',
                   'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL',
                   'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed'],
            usecols=['UniProtKB-AC', 'UniProtKB-ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO',
                     'NCBI-taxon', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO'],
            dtype='str')

        if blocksize:
            if "idmapping_selected.parquet" in file_resources and \
                isinstance(file_resources["idmapping_selected.parquet"], str):
                idmapping = dd.read_parquet(file_resources["idmapping_selected.parquet"])

            elif "idmapping_selected.tab" in file_resources and isinstance(file_resources["idmapping_selected.tab"],
                                                                           str):
                idmapping = dd.read_table(file_resources["idmapping_selected.tab"], blocksize=blocksize, **args)
            else:
                idmapping = dd.read_table(file_resources["idmapping_selected.tab.gz"], compression="gzip", **args, )

            if not idmapping.index.name:
                idmapping = idmapping.set_index('UniProtKB-AC', sorted=False)
        else:
            idmapping: pd.DataFrame = pd.read_table(file_resources["idmapping_selected.tab"], index_col='UniProtKB-AC',
                                                    **args)

        # Filter UniProt accession keys
        if self.keys is not None and idmapping.index.name != None:
            idmapping = idmapping.loc[idmapping.index.isin(self.keys)]

        # Convert string of list elements to a np.array
        list2array = lambda x: np.array(x) if x is not None else x
        args = dict(meta=pd.Series([np.array([''])])) if isinstance(idmapping, dd.DataFrame) else {}
        for col in ['PDB', 'GI', 'GO', 'RefSeq']:
            if col not in idmapping.columns or idmapping[col].head(300).map(type).isin({Iterable}).any(): continue
            # Split string to list
            idmapping[col] = idmapping[col].str.split("; ").map(list2array, **args)

        for col in ['Ensembl', 'Ensembl_TRS', 'Ensembl_PRO']:
            if col not in idmapping.columns or idmapping[col].head(300).map(type).isin({Iterable}).any(): continue
            # Removing .# ENGS gene version number at the end
            if self.remove_version_num:
                idmapping[col] = idmapping[col].str.replace("[.]\d*", "", regex=True)
            idmapping[col] = idmapping[col].str.split("; ").map(list2array, **args)

            if col == 'Ensembl_PRO':
                # Prepend species_id to ensembl protein ids to match with STRING PPI
                idmapping['protein_external_id'] = idmapping[col]
                idmapping['protein_external_id'] = idmapping[["NCBI-taxon", 'protein_external_id']].apply(
                    lambda x: np.array(
                        [".".join([x['NCBI-taxon'], protein_id]) for protein_id in x['protein_external_id']]) \
                        if isinstance(x['protein_external_id'], list) else None,
                    axis=1, **args)

        # Join metadata from uniprot_sprot.parquet
        if 'uniprot_sprot.parquet' in file_resources or 'uniprot_trembl.parquet' in file_resources:
            uniprot = self.load_uniprot_parquet(file_resources, blocksize=blocksize)
            to_join = uniprot[uniprot.columns.difference(idmapping.columns)]
            idmapping = idmapping.join(to_join, how='left')

        # Load proteome.tsv
        if "proteomes.tsv" in file_resources:
            proteomes = pd.read_table(file_resources["proteomes.tsv"],
                                      usecols=['Organism Id', 'Proteome Id'],
                                      dtype={'Organism Id': 'str', 'Proteome Id': 'str'}) \
                .rename(columns={'Organism Id': 'NCBI-taxon', 'Proteome Id': 'proteome_id'}) \
                .dropna().set_index('NCBI-taxon')
            idmapping = idmapping.join(proteomes, on='NCBI-taxon')

        # Load species info from speclist.txt
        if 'speclist.txt' in file_resources:
            speclist = pd.read_fwf(file_resources['speclist.txt'],
                                   names=['species_code', 'Taxon', 'species_id', 'attr'],
                                   comment="==", skipinitialspace=True, skiprows=59, skipfooter=4)
            speclist = speclist.drop(index=speclist.index[~speclist['attr'].str.contains("=")])
            speclist['species_id'] = speclist['species_id'].str.rstrip(":")

            speclist = speclist.fillna(method='ffill')
            speclist = speclist.groupby(speclist.columns[:3].tolist())['attr'] \
                .apply('|'.join) \
                .apply(lambda s: dict(map(str.strip, sub.split('=', 1)) for sub in s.split("|") if '=' in sub)) \
                .apply(pd.Series)

            speclist = speclist.rename(columns={'N': 'Official (scientific) name', 'C': 'Common name', 'S': 'Synonym'}) \
                .reset_index() \
                .set_index('species_id')
            speclist['Taxon'] = speclist['Taxon'].replace(
                {'A': 'archaea', 'B': 'bacteria', 'E': 'eukaryota', 'V': 'viruses', 'O': 'others'})
            speclist.index.name = 'NCBI-taxon'
            idmapping = idmapping.join(speclist, on='NCBI-taxon')

        return idmapping

    def load_uniprot_parquet(self, file_resources: Dict[str, str], blocksize=None) -> Union[dd.DataFrame, pd.DataFrame]:
        dfs = []
        for filename, file_path in file_resources.items():
            if filename not in ['uniprot_sprot.parquet', 'uniprot_trembl.parquet']: continue

            if blocksize:
                df = dd.read_parquet(file_path).rename(columns=UniProt.COLUMNS_RENAME_DICT)
                if self.keys is not None and self.index_col:
                    df = df.loc[df[self.index_col].isin(self.keys)]

                if not df.index.name:
                    try:
                        df = df.set_index(self.index_col, sorted=True)
                    except Exception as e:
                        print(file_path, e)
                        df = df.set_index(self.index_col, sorted=False)
            else:
                df = pd.read_parquet(file_path, index_col=self.index_col).rename(columns=UniProt.COLUMNS_RENAME_DICT)
                if self.keys is not None:
                    df_keys = df.index if df.index.name == self.index_col else df[self.index_col]
                    df = df.loc[df_keys.isin(self.keys)]
            dfs.append(df)

        if dfs:
            dfs = dd.concat(dfs) if blocksize else pd.concat(dfs)

            return dfs
        else:
            return None

    def load_uniprot_xml(self, file_path: str, keys=None, blocksize=None) -> pd.DataFrame:
        records = []
        seqfeats = []
        if isinstance(keys, str):
            index = keys
            keys_set = self.data.index if keys == self.data.index.name else self.data[keys]
        if isinstance(keys, (dd.Index, dd.Series)):
            index = keys.name
            keys_set = keys.compute()
        else:
            index = keys_set = None

        for record in tqdm.tqdm(SeqIO.parse(file_path, "uniprot-xml"), desc=str(file_path)):
            # Sequence features
            annotations = defaultdict(None, record.annotations)
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
            if index is not None:
                if record_dict[keys] not in keys_set: continue

            records.append(record_dict)

            # Sequence interval features
            _parse_interval = lambda sf: pd.Interval(left=sf.location.start, right=sf.location.end, )
            feature_type_intervals = defaultdict(lambda: [])
            for sf in record.features:
                if isinstance(sf.location.start, ExactPosition) and isinstance(sf.location.end, ExactPosition):
                    feature_type_intervals[sf.type].append(_parse_interval(sf))

            features_dict = {type: pd.IntervalIndex(intervals, name=type) \
                             for type, intervals in feature_type_intervals.items()}
            seqfeats.append({"protein_id": record.id, **features_dict})

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

        return records_df

    def load_sequences(self, fasta_file: str, index=None, keys: Union[pd.Index, List[str]] = None, blocksize=None) \
        -> OrderedDict:
        def get_id(s: str):
            if index == 'protein_id':
                return s.split('|')[1]
            elif index == 'protein_name':
                return s.split('|')[2]
            else:
                return s.split('|')[1]

        fa = Fasta(fasta_file, key_function=get_id, as_raw=True, )

        return fa.records

    def get_sequences(self, index: str, omic: str = None, agg: str = "first", **kwargs):
        assert index, '`index` must be either "protein_id" or "protein_name"'

        # Parse lncRNA & mRNA fasta
        seq_df = self.load_sequences(self.file_resources["uniprot_sprot.fasta"], index=index, blocksize=self.blocksize)
        if "uniprot_trembl.fasta" in self.file_resources:
            trembl_seq_df = self.load_sequences(self.file_resources["uniprot_trembl.fasta"], index=index,
                                                blocksize=self.blocksize)
            seq_df.update(trembl_seq_df)

        return seq_df

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

        if 'rnacentral.mirbase.tsv' not in file_resources:
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
        rnacentral_mirbase = pd.read_table(
            file_resources["rnacentral.mirbase.tsv"], low_memory=True, header=None,
            names=["RNAcentral id", "database", "mirbase id", "species_id", "RNA type", "NA", ],
            usecols=["RNAcentral id", "database", "mirbase id", "species_id", "RNA type"],
            dtype={'mirbase id': 'str', "species_id": "str"})

        if isinstance(self.species_id, str):
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species_id"] == self.species_id]
        elif isinstance(self.species_id, list):
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species_id"].isin(set(self.species_id))]

        mirbase_df = pd.read_table(file_resources["aliases.txt"], low_memory=True, header=None,
                                   names=["mirbase id", "gene_name"],
                                   dtype='str', )
        mirbase_df = mirbase_df.merge(rnacentral_mirbase, on='mirbase id', how="left")

        # Expanding miRNA names in each MirBase Ascension ID
        mirbase_df['gene_name'] = mirbase_df['gene_name'].str.rstrip(";").str.split(";")
        mirbase_df = mirbase_df.explode(column='gene_name')

        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        return mirbase_df

    def load_sequences(self, fasta_file, index=None, keys=None, blocksize=None):
        """
        Args:
            index ():
            fasta_file:
            keys ():
            blocksize:
        """
        if hasattr(self, '_seq_df_dict') and fasta_file in self._seq_df_dict:
            return self._seq_df_dict[fasta_file]

        fa = Fasta(fasta_file, read_long_names=True, as_raw=True)

        entries = []
        for key, record in tqdm.tqdm(fa.items(), desc=str(fasta_file)):
            attrs = record.long_name.split(" ")
            record_dict = {
                "gene_id": attrs[0],
                "mirbase_id": attrs[1],
                "gene_name": attrs[-2],
                "species": " ".join(attrs[2:-3]),
                "type": attrs[-1],
                SEQUENCE_COL: str(record),
            }

            entries.append(record_dict)

        df = pd.DataFrame(entries)
        # if blocksize:
        #     df = dd.from_pandas(df, chunksize=blocksize)

        if not hasattr(self, '_seq_df_dict'):
            self._seq_df_dict = {}
        self._seq_df_dict[fasta_file] = df

        return df

    def get_sequences(self,
                      index="gene_name",
                      omic=None,
                      agg="all",
                      **kwargs):
        """
        Args:
            index:
            omic:
            agg:
            **kwargs:
        """
        dfs = []
        for filename in self.file_resources:
            if filename.endswith('.fa'):
                seq_df = self.load_sequences(self.file_resources[filename])
                dfs.append(seq_df)
        seq_df = pd.concat(dfs, axis=0)

        seq_df = seq_df.groupby(index)[SEQUENCE_COL].agg(self.aggregator_fn(agg))

        return seq_df
