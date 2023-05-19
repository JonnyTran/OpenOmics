import gzip
import os
import zipfile
from os.path import exists, getsize
from typing import Tuple, Union, TextIO, Optional, Dict, List
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


def decompress_file(filepath: str, filename: str, file_ext: filetype.Type, write_uncompressed=False) \
    -> Tuple[Union[gzip.GzipFile, TextIO], str]:
    """
    Decompress the `filepath` corresponding to its `file_ext` compression type, then return the uncompressed data (or its path) and
    the `filename` without the `file_ext` suffix.

    Args:
        filepath (str): The file path to the data file
        filename (str): The filename of the data file
        file_ext (filetype.Type): The file extension of the data file
        write_uncompressed (bool): Whether to write the uncompressed file to disk

    Returns:
        uncompressed_file (): The uncompressed file path
        updated_filename (str): The filename without the `file_ext` suffix
    """
    data = filepath

    if file_ext is None:
        return data, filename

    elif file_ext.extension == "gz":
        data = gzip.open(filepath, "rt")

    elif file_ext.extension == "zip":
        with zipfile.ZipFile(filepath, "r") as zf:
            for subfile in zf.infolist():
                # Select first file with matching file extension
                subfile_name = os.path.splitext(subfile.filename)[-1]
                if subfile_name == os.path.splitext(filename.replace(".zip", ""))[-1]:
                    data = zf.open(subfile.filename, mode="r")

    elif file_ext.extension == "rar":
        with rarfile.RarFile(filepath, "r") as rf:

            for subfile in rf.infolist():
                # If the file extension matches
                subfile_name = os.path.splitext(subfile.filename)[-1]
                if subfile_name == os.path.splitext(filename.replace(".rar", ""))[-1]:
                    data = rf.open(subfile.filename, mode="r")

    else:
        logger.warn(f"WARNING: filepath_ext.extension {file_ext.extension} not supported.")
        return data, filename

    filename = get_uncompressed_filepath(filename)
    uncompressed_path = get_uncompressed_filepath(filepath)

    if write_uncompressed and not exists(uncompressed_path):
        with open(uncompressed_path, 'w', encoding='utf8') as f_out:
            logger.info(f"Writing uncompressed {filename} file to {uncompressed_path}")
            f_out.write(data.read())

    if exists(uncompressed_path) and getsize(uncompressed_path) > 0:
        data = uncompressed_path

    return data, filename


def get_uncompressed_filepath(filepath: str) -> str:
    """Return the uncompressed filepath by removing the file extension suffix.

    Args:
        filepath (str): File path to the compressed file

    Returns:
        uncompressed_path (str): File path to the uncompressed file
    """
    uncompressed_path = ''
    if filepath.endswith(".gz"):
        uncompressed_path = filepath.removesuffix(".gz")
    elif filepath.endswith(".zip"):
        uncompressed_path = filepath.removesuffix(".zip")
    elif filepath.endswith(".rar"):
        uncompressed_path = filepath.removesuffix(".rar")
    else:
        uncompressed_path = filepath + ".uncompressed"

    if uncompressed_path and filepath != uncompressed_path:
        return uncompressed_path
    else:
        return ''


def select_files_with_ext(file_resources: Dict[str, str], ext: str, contains: Optional[str] = None) -> Dict[str, str]:
    """Return a list of file paths with the specified file extension. Only string values are considered as file paths.

    Args:
        file_resources (dict): A dictionary of file names and their corresponding file paths
        ext (str): The file extension to filter the file names by
        contains (str): If not None, only return file paths that contain the specified string

    Returns:
        file_paths (dict): A dict of file names and corresponding paths with the specified file extension
    """
    subset_file_resources = {}
    for filename, filepath in file_resources.items():
        if not isinstance(filepath, str): continue
        if filename.endswith(ext) and (contains is None or contains in filename):
            subset_file_resources[filename] = filepath

    return subset_file_resources


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
