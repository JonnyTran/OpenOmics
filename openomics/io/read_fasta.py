import dask.dataframe as dd
from pyfaidx import Faidx


def read_fasta(filepath, ) -> dd.DataFrame:
    """
    Take in a path of a .fasta file and use Faidx to index it, parse the index names, then build a lazy-loading
    Dask DataFrame with keys on the index and sequence as one of the columns.

    Args:
        filepath ():
    """
    fa = Faidx(filepath)
