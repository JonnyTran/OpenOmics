import gzip
import os
import zipfile
from typing import Tuple, Union, TextIO
from urllib.error import URLError
from logzero import logger

import dask.dataframe as dd
import filetype
import rarfile
import requests
import sqlalchemy as sa
import validators
from astropy.utils import data
from requests.adapters import HTTPAdapter, Retry

import openomics


# @astropy.config.set_temp_cache(openomics.config["cache_dir"])
def get_pkg_data_filename(baseurl: str, filepath: str):
    """Downloads a remote file given the url, then caches it to the user's home
    folder.

    Args:
        baseurl: Url to the download path, excluding the file name
        filepath: The file path to download

    Returns:
        filename (str): A file path on the local file system corresponding to
        the data requested in data_name.
    """
    # Split data url and file name if the user provided a whole path in `file_resources`
    if validators.url(filepath):
        base, filepath = os.path.split(filepath)
        base = base + "/"
    else:
        base, filepath = baseurl, filepath

    try:
        # TODO doesn't yet save files to 'cache_dir' but to astropy's default cache dir
        # logger.debug(f"Fetching file from: {base}{filepath}, saving to {openomics.config['cache_dir']}")

        with data.conf.set_temp("dataurl", base), data.conf.set_temp("remote_timeout", 30):
            return data.get_pkg_data_filename(filepath, package="openomics.database", show_progress=True)

    except (URLError, ValueError) as e:
        raise Exception(f"Unable to download file at {os.path.join(base, filepath)}. "
                        f"Please try manually downloading the files and add path to `file_resources` arg. \n{e}")


def decompress_file(filepath: str, filename: str, file_ext: filetype.Type) \
    -> Tuple[Union[gzip.GzipFile, TextIO], str]:
    """
    Decompress the `data_file` corresponding to its `file_ext`, then remove the `file_ext` suffix from `filename`.
    Args:
        filepath ():
        filename ():
        file_ext ():

    Returns:
        updated_filename, uncompressed_file
    """
    data = filepath
    # This null if-clause is needed when filetype_ext is None, causing the next clauses to fail
    if file_ext is None:
        data = filepath

    elif file_ext.extension == "gz":
        logger.info("Decompressed gzip file at {}".format(filepath))
        data = gzip.open(filepath, "rt")
        filename = filename.replace(".gz", "")

    elif file_ext.extension == "zip":
        logger.info("Decompressed zip file at {}".format(filepath))
        with zipfile.ZipFile(filepath, "r") as zf:
            filename = filename.replace(".zip", "")

            for subfile in zf.infolist():
                # If the file extension matches
                if os.path.splitext(subfile.filename)[-1] == os.path.splitext(filename)[-1]:
                    data = zf.open(subfile.filename, mode="r")

    elif file_ext.extension == "rar":
        logger.info("Decompressed rar file at {}".format(filepath))
        with rarfile.RarFile(filepath, "r") as rf:
            filename = filename.replace(".rar", "")

            for subfile in rf.infolist():
                # If the file extension matches
                if os.path.splitext(subfile.filename)[-1] == os.path.splitext(filename)[-1]:
                    data = rf.open(subfile.filename, mode="r")
    else:
        print(f"WARNING: filepath_ext.extension {file_ext.extension} not supported.")
        data = filepath

    return data, filename

def read_db(path, table, index_col):
    """
    Args:
        path:
        table:
        index_col:
    """
    engine = sa.create_engine(path)
    # conn = engine.connect()
    m = sa.MetaData()
    table = sa.Table(table, m, autoload=True, autoload_with=engine)

    # conn.execute("create table testtable (uid integer Primary Key, datetime NUM)")
    # conn.execute("insert into testtable values (1, '2017-08-03 01:11:31')")
    # print(conn.execute('PRAGMA table_info(testtable)').fetchall())
    # conn.close()

    uid, dt = list(table.columns)
    q = sa.select([dt.cast(sa.types.String)]).select_from(table)

    daskDF = dd.read_sql_table(table, path, index_col=index_col, parse_dates={'datetime': '%Y-%m-%d %H:%M:%S'})
    return daskDF


def retry(num=5):
    """retry connection.

    define max tries num if the backoff_factor is 0.1, then sleep() will
    sleep for [0.1s, 0.2s, 0.4s, ...] between retries. It will also force a
    retry if the status code returned is 500, 502, 503 or 504.

    Args:
        num:
    """
    s = requests.Session()
    retries = Retry(total=num, backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])
    s.mount('http://', HTTPAdapter(max_retries=retries))

    return s
