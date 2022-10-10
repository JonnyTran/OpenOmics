import os
import re
from abc import abstractmethod
from collections import defaultdict, OrderedDict
from typing import Union, List, Callable, Dict, Tuple, Optional, Iterable

import numpy as np
import pandas as pd
import tqdm
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from dask import dataframe as dd
from logzero import logger
from pyfaidx import Fasta
from six.moves import intern

import openomics
from openomics.io.read_gtf import read_gtf
from .base import Database
from ..transforms.agg import get_agg_func
from ..transforms.df import drop_duplicate_columns

__all__ = ['GENCODE', 'UniProt', 'MirBase', 'RNAcentral']

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
        -> pd.DataFrame:
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
    def get_sequences(self, index: str, omic: str, agg: str, **kwargs) -> Union[pd.Series, Dict]:
        """Returns a dictionary where keys are 'index' and values are
        sequence(s).

        Args:
            index (str): {"gene_id", "gene_name", "transcript_id",
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
            f"Fasta file provided in `file_resources` must be an uncompressed .fa file. Given {fasta_file}."

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
        # idmapping_selected.tab
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
        "accession": "UniProtKB-AC", "name": "protein_name", "gene": "gene_name", "keyword": "keywords",
        "geneLocation": "subcellular_location",

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
                 index_col='UniProtKB-AC', keys=None,
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
        self.species = UniProt.SPECIES_ID_NAME[species_id] if species_id in UniProt.SPECIES_ID_NAME else None
        self.taxonomic_id = UniProt.SPECIES_ID_TAXONOMIC[
            self.species] if self.species in UniProt.SPECIES_ID_TAXONOMIC else None
        self.remove_version_num = remove_version_num

        if file_resources is None:
            file_resources = {}

            file_resources['uniprot_sprot.xml.gz'] = os.path.join(path, "knowledgebase/complete/uniprot_sprot.xml.gz")
            file_resources['uniprot_trembl.xml.gz'] = os.path.join(path, "knowledgebase/complete/uniprot_trembl.xml.gz")
            file_resources["idmapping_selected.tab.gz"] = os.path.join(
                path, "knowledgebase/idmapping/idmapping_selected.tab.gz")

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

        super().__init__(path=path, file_resources=file_resources, index_col=index_col, keys=keys,
                         col_rename=col_rename,
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

            elif "idmapping_selected.tab" in file_resources and \
                isinstance(file_resources["idmapping_selected.tab"], str):
                idmapping = dd.read_table(file_resources["idmapping_selected.tab"], blocksize=blocksize, **args)
            else:
                idmapping = dd.read_table(file_resources["idmapping_selected.tab.gz"], compression="gzip", **args, )

            idmapping: dd.DataFrame
        else:
            if "idmapping_selected.parquet" in file_resources and \
                isinstance(file_resources["idmapping_selected.parquet"], str):
                idmapping = pd.read_parquet(file_resources["idmapping_selected.parquet"])
            else:
                idmapping = pd.read_table(file_resources["idmapping_selected.tab"], index_col=self.index_col, **args)

        # Filter UniProt accession keys
        if self.keys is not None and idmapping.index.name == self.index_col:
            idmapping = idmapping.loc[idmapping.index.isin(self.keys)]
        elif self.keys is not None and idmapping.index.name != self.index_col:
            idmapping = idmapping.loc[idmapping[self.index_col].isin(self.keys)]

        if idmapping.index.name != self.index_col:
            idmapping = idmapping.set_index(self.index_col, sorted=False)
        if not idmapping.known_divisions:
            idmapping.divisions = idmapping.compute_current_divisions()

        # Transform list columns
        if isinstance(idmapping, dd.DataFrame):
            idmapping = idmapping.assign(**self.assign_transforms(idmapping))
        else:
            idmapping = idmapping.assign(**self.assign_transforms(idmapping))

        # Join metadata from uniprot_sprot.parquet
        if any(fn.startswith('uniprot') and fn.endswith('.parquet') for fn in file_resources):
            uniprot_anns = self.load_uniprot_parquet(file_resources, blocksize=blocksize)
            uniprot_anns = uniprot_anns[uniprot_anns.columns.difference(idmapping.columns)]
            uniprot_anns = drop_duplicate_columns(uniprot_anns)
            assert idmapping.index.name == uniprot_anns.index.name, f"{idmapping.index.name} != {uniprot_anns.index.name}"
            idmapping = idmapping.join(uniprot_anns, on=idmapping.index.name, how='left')

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

    def assign_transforms(self, idmapping: pd.DataFrame) -> Dict[str, Union[dd.Series, pd.Series]]:
        # Convert string of list elements to a np.array
        list2array = lambda x: np.array(x) if isinstance(x, Iterable) else x
        assign_fn = {}
        for col in {'PDB', 'GI', 'GO', 'RefSeq'}.intersection(idmapping.columns):
            try:
                # Split string to list
                assign_fn[col] = idmapping[col].str.split("; ").map(list2array)
            except:
                continue

        for col in {'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO'}.intersection(idmapping.columns):
            # Removing .# ENGS gene version number at the end
            try:
                if self.remove_version_num:
                    series = idmapping[col].str.replace("[.]\d*", "", regex=True)
                else:
                    series = idmapping[col]

                assign_fn[col] = series.str.split("; ").map(list2array)

                if col == 'Ensembl_PRO':
                    # Prepend species_id to ensembl protein ids to match with STRING PPI
                    concat = dd.concat([idmapping["NCBI-taxon"], assign_fn[col]]) \
                        if isinstance(idmapping, dd.DataFrame) else \
                        pd.concat([idmapping["NCBI-taxon"], assign_fn[col]])

                    assign_fn['protein_external_id'] = concat.apply(
                        lambda row: np.char.add(row['NCBI-taxon'] + ".", row['Ensembl_PRO']) \
                            if isinstance(row['Ensembl_PRO'], Iterable) else None,
                        axis=1)
            except:
                continue

        return assign_fn

    def load_uniprot_parquet(self, file_resources: Dict[str, str], blocksize=None) -> Union[dd.DataFrame, pd.DataFrame]:
        dfs = []
        for filename, file_path in file_resources.items():
            if not ('uniprot' in filename and filename.endswith('.parquet')): continue
            if blocksize:
                df: dd.DataFrame = dd.read_parquet(file_path) \
                    .rename(columns=UniProt.COLUMNS_RENAME_DICT)
                if df.index.name in UniProt.COLUMNS_RENAME_DICT:
                    df.index = df.index.rename(UniProt.COLUMNS_RENAME_DICT[df.index.name])

                if self.keys is not None:
                    if self.index_col in df.columns:
                        df = df.loc[df[self.index_col].isin(self.keys)]
                    elif df.index.name == self.index_col:
                        df = df.loc[df.index.isin(self.keys)]

                if df.index.size.compute() == 0: continue

                if df.index.name != self.index_col:
                    try:
                        df = df.set_index(self.index_col, sorted=True)
                    except Exception as e:
                        print(file_path, e)
                        df = df.set_index(self.index_col, sorted=False)

                if not df.known_divisions:
                    df.divisions = df.compute_current_divisions()

            else:
                df = pd.read_parquet(file_path).rename(columns=UniProt.COLUMNS_RENAME_DICT).set_index(self.index_col)

                if self.keys is not None:
                    df_keys = df.index if df.index.name == self.index_col else df[self.index_col]
                    df = df.loc[df_keys.isin(self.keys)]
                if df.index.size == 0: continue

            dfs.append(df)

        if dfs:
            dfs = dd.concat(dfs, interleave_partitions=True) if blocksize else pd.concat(dfs)

            return dfs
        else:
            return None

    def load_uniprot_xml(self, file_path: str, keys=None, blocksize=None) -> pd.DataFrame:
        records = []
        seqfeats = []
        if isinstance(keys, str):
            index = keys
            keys_set = self.data.index if keys == self.data.index.name else self.data[keys]
        elif isinstance(keys, (dd.Index, dd.Series)):
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
        path="http://mirbase.org/ftp/CURRENT/",
        file_resources=None,
        species_id: Optional[str] = '9606',
        index_col: str = "mirbase_id",
        col_rename=None,
        **kwargs,
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

        self.species_id = species_id
        super().__init__(path=path, file_resources=file_resources, index_col=index_col, col_rename=col_rename, **kwargs)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        rnacentral_mirbase = pd.read_table(
            file_resources["rnacentral.mirbase.tsv"], low_memory=True, header=None,
            names=["RNAcentral id", "database", "mirbase_id", "species_id", "RNA type", "NA"],
            usecols=["RNAcentral id", "database", "mirbase_id", "species_id", "RNA type"],
            index_col="RNAcentral id",
            dtype={'mirbase_id': 'str', "species_id": "category", 'database': 'category', 'RNA type': 'category'})

        if isinstance(self.species_id, str):
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species_id"] == self.species_id]
        elif isinstance(self.species_id, Iterable):
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species_id"].isin(set(self.species_id))]

        mirbase_df = pd.read_table(file_resources["aliases.txt"], low_memory=True, header=None,
                                   names=["mirbase_id", "mirbase_name"], index_col=self.index_col,
                                   dtype='str', )
        if mirbase_df.index.name == 'mirbase id':
            mirbase_df = mirbase_df.join(rnacentral_mirbase, on=self.index_col, how="left", rsuffix='_rnacentral')
        else:
            mirbase_df = mirbase_df.merge(rnacentral_mirbase, on=self.index_col, how="left")

        # Expanding miRNA names in each MirBase Ascension ID
        mirbase_df['mirbase_name'] = mirbase_df['mirbase_name'].str.rstrip(";").str.split(";")

        seq_dfs = []
        for filename in file_resources:
            if filename.endswith('.fa') or filename.endswith('.fasta'):
                assert isinstance(file_resources[filename], str), f"Must provide a path to an uncompressed .fa file. " \
                                                                  f"Given {file_resources[filename]}"
                df = self.load_sequences(file_resources[filename], index=self.index_col, keys=self.keys)
                seq_dfs.append(df)

        if len(seq_dfs):
            seq_dfs = pd.concat(seq_dfs, axis=0)
            mirbase_df = mirbase_df.join(seq_dfs, how='left', on=self.index_col)
        else:
            logger.info('Missing sequence data because no "hairpin.fa" or "mature.fa" file were given.')

        # mirbase_df = mirbase_df.explode(column='gene_name')
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        return mirbase_df

    def load_sequences(self, fasta_file, index=None, keys=None, blocksize=None):
        """
        Args:
            fasta_file:
            index ():
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
                "species": intern(" ".join(attrs[2:4])),
                "gene_name": attrs[4],
                "type": intern(attrs[5]) if len(attrs) >= 6 else None,
                SEQUENCE_COL: str(record),
            }
            if keys is not None and index:
                if record_dict[index] not in keys:
                    del record_dict
                    continue

            entries.append(record_dict)

        df = pd.DataFrame(entries)
        if index:
            df = df.set_index(index)
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


class RNAcentral(SequenceDatabase):
    """Loads the RNAcentral database from https://rnacentral.org/ .

        Default path: https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/ .
        Default file_resources: {
            "rnacentral_rfam_annotations.tsv": "go_annotations/rnacentral_rfam_annotations.tsv.gz",
            "database_mappings/gencode.tsv": "id_mapping/database_mappings/gencode.tsv",
            "gencode.fasta": "sequences/by-database/gencode.fasta",
            ...
        }
    """
    COLUMNS_RENAME_DICT = {
        'ensembl_gene_id': 'gene_id',
        'external id': 'transcript_id',
        'GO terms': 'go_id',
        'gene symbol': 'gene_id',
    }

    def __init__(self, path="https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/", file_resources=None,
                 col_rename=COLUMNS_RENAME_DICT, species_id: Union[List[str], str, None] = None,
                 index_col="RNAcentral id", keys=None,
                 remove_version_num=True, remove_species_suffix=True, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            col_rename ():
            species_id ():
            index_col ():
            keys ():
            remove_version_num ():
            remove_species_suffix ():
            **kwargs ():
        """
        self.species_id = species_id
        self.remove_version_num = remove_version_num
        self.remove_species_suffix = remove_species_suffix

        if file_resources is None:
            file_resources = {}
            file_resources["rnacentral_rfam_annotations.tsv.gz"] = "go_annotations/rnacentral_rfam_annotations.tsv.gz"
            file_resources["database_mappings/ensembl_gencode.tsv"] = "id_mapping/database_mappings/ensembl_gencode.tsv"
            file_resources["database_mappings/mirbase.tsv"] = "id_mapping/database_mappings/mirbase.tsv"

        super().__init__(path=path, file_resources=file_resources, col_rename=col_rename, index_col=index_col,
                         keys=keys, **kwargs)

    def load_dataframe(self, file_resources, blocksize=None):
        """
        Args:
            file_resources:
            blocksize:
        """
        # Build transcripts ids by combining `database_mappings/` files
        transcripts_df = []
        for filename in (fname for fname in file_resources if "database_mappings" in fname):
            args = dict(low_memory=True, header=None,
                        names=["RNAcentral id", "database", "external id", "species_id", "RNA type", "gene symbol"],
                        dtype={'gene symbol': 'str',
                               'database': 'category', 'species_id': 'category', 'RNA type': 'category', })

            if blocksize:
                if filename.endswith('.tsv'):
                    id_mapping: dd.DataFrame = dd.read_table(
                        file_resources[filename], blocksize=None if isinstance(blocksize, bool) else blocksize, **args)
                elif filename.endswith('.parquet'):
                    id_mapping: dd.DataFrame = dd.read_parquet(
                        file_resources[filename], blocksize=None if isinstance(blocksize, bool) else blocksize, )
                else:
                    id_mapping = None
            else:
                if filename.endswith('.tsv'):
                    id_mapping = pd.read_table(file_resources[filename], **args)
                elif filename.endswith('.parquet'):
                    id_mapping = pd.read_parquet(file_resources[filename])
                else:
                    id_mapping = None

            if id_mapping is None:
                raise Exception("Must provide a file with 'database_mappings/(*).tsv' in file_resources")

            # Filter by species
            if isinstance(self.species_id, str):
                id_mapping = id_mapping.where(id_mapping["species_id"] == self.species_id)
            elif isinstance(self.species_id, Iterable):
                id_mapping = id_mapping.where(id_mapping["species_id"].isin(self.species_id))

            # Filter by index
            if self.keys and id_mapping.index.name == self.index_col:
                id_mapping = id_mapping.loc[id_mapping.index.isin(self.keys)]
            elif self.keys and id_mapping.index.name != self.index_col:
                id_mapping = id_mapping.loc[id_mapping[self.index_col].isin(self.keys)]

            # Add species_id prefix to index values
            if "RNAcentral id" in id_mapping.columns:
                id_mapping["RNAcentral id"] = id_mapping["RNAcentral id"] + "_" + id_mapping["species_id"].astype(str)
            elif "RNAcentral id" == id_mapping.index.name:
                id_mapping.index = id_mapping.index + "_" + id_mapping["species_id"].astype(str)

            # Add sequence column from FASTA file of the same database
            fasta_filename = f"{filename.split('/')[-1].split('.')[0]}.fasta"
            if fasta_filename in file_resources:
                seq_df = self.load_sequences(file_resources[fasta_filename])
                id_mapping = id_mapping.merge(seq_df, how='left',
                                              left_on="RNAcentral id",
                                              left_index=True if id_mapping.index.name == "RNAcentral id" else False,
                                              right_index=True)
            else:
                logger.info(f"{fasta_filename} not provided for `{filename}` so missing sequencing data")

            # Remove species_id prefix from index values if necessary
            if self.remove_version_num and 'gene symbol' in id_mapping.columns:
                id_mapping["gene symbol"] = id_mapping["gene symbol"].str.replace("[.].\d*", "", regex=True)
            if self.remove_species_suffix:
                id_mapping["RNAcentral id"] = id_mapping["RNAcentral id"].str.replace("_(\d*)", '', regex=True)

            # Set index
            args = dict(sorted=True) if blocksize else {}
            id_mapping = id_mapping.set_index(self.index_col, **args)
            if not id_mapping.known_divisions:
                id_mapping.divisions = id_mapping.compute_current_divisions()

            transcripts_df.append(id_mapping)

        # Concatenate multiple `database_mappings` files from different databases
        if blocksize:
            transcripts_df = dd.concat(transcripts_df, axis=0, interleave_partitions=True, join='outer')
        else:
            transcripts_df = pd.concat(transcripts_df, axis=0, join='outer')

        # Join go_id and Rfams annotations to each "RNAcentral id" from 'rnacentral_rfam_annotations.tsv'
        args = dict(low_memory=True, names=["RNAcentral id", "GO terms", "Rfams"])
        if blocksize:
            if 'rnacentral_rfam_annotations.tsv' in file_resources and isinstance(
                file_resources['rnacentral_rfam_annotations.tsv'], str):
                anns = dd.read_table(file_resources["rnacentral_rfam_annotations.tsv"], **args)
            else:
                anns = dd.read_table(file_resources["rnacentral_rfam_annotations.tsv.gz"], compression="gzip", **args)

            # Filter annotations by "RNAcentral id" in `transcripts_df`
            anns = anns.loc[anns["RNAcentral id"].isin(transcripts_df.index.compute())]

            anns = anns.set_index("RNAcentral id", sorted=True)
            if not anns.known_divisions:
                anns.divisions = anns.compute_current_divisions()

            # Groupby on index
            anns_groupby: dd.DataFrame = anns \
                .groupby(by=lambda idx: idx) \
                .agg({col: get_agg_func('unique', use_dask=True) for col in ["GO terms", 'Rfams']})
        else:
            anns = pd.read_table(file_resources["rnacentral_rfam_annotations.tsv"], index_col='RNAcentral id', **args)
            idx = transcripts_df.index.compute() if isinstance(transcripts_df, dd.DataFrame) else transcripts_df.index
            anns = anns.loc[anns.index.isin(set(idx))]
            anns_groupby = anns.groupby("RNAcentral id").agg({col: 'unique' for col in ["GO terms", 'Rfams']})

        transcripts_df = transcripts_df.merge(anns_groupby, how='left', left_index=True, right_index=True)

        return transcripts_df

    def load_sequences(self, fasta_file: str, index=None, keys=None, blocksize=None):
        """
        Args:
            index ():
            fasta_file:
            keys ():
            blocksize:
        """
        if hasattr(self, '_seq_df_dict') and fasta_file in self._seq_df_dict:
            return self._seq_df_dict[fasta_file]

        fa = Fasta(fasta_file, as_raw=True)

        entries = []
        for key, record in tqdm.tqdm(fa.items(), desc=str(fasta_file)):
            id = re.sub("_(\d*)", '', key) if self.remove_species_suffix else key
            if keys is not None and self.index_col == 'RNAcentral id' and id not in keys:
                continue

            if ") " in record.long_name:
                desc = record.long_name.split(") ")[-1].strip()
            elif "\\" in record.long_name:
                desc = record.long_name.split("\\", maxsplit=1)[-1].strip()
            elif "microRNA " in record.long_name:
                desc = record.long_name.split("microRNA ", maxsplit=1)[-1].strip()
            else:
                desc = record.long_name.split(" ", maxsplit=3)[-1]

            record_dict = {
                'RNAcentral id': key,
                'description': desc,
                SEQUENCE_COL: str(record),
            }

            entries.append(record_dict)

        df = pd.DataFrame(entries).set_index("RNAcentral id")

        if not hasattr(self, '_seq_df_dict'):
            self._seq_df_dict = {}
        self._seq_df_dict[fasta_file] = df

        return df

    def get_sequences(self,
                      index="RNAcentral id",
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
            if filename.endswith('.fa') or filename.endswith('.fasta'):
                seq_df = self.load_sequences(self.file_resources[filename])
                dfs.append(seq_df)
        seq_df = pd.concat(dfs, axis=0)

        seq_df = seq_df.groupby(index)[SEQUENCE_COL].agg(self.aggregator_fn(agg))

        return seq_df
