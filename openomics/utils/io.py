import gzip
import io
import logging
import os
import zipfile
from pathlib import Path
from typing import Tuple, Union, TextIO
from urllib.error import URLError

import dask.dataframe as dd
import openomics
import rarfile
import requests
import sqlalchemy as sa
import validators
from astropy.utils import data
from filetype import Type
from requests.adapters import HTTPAdapter, Retry


# @astropy.config.set_temp_cache(openomics.config["cache_dir"])
def get_pkg_data_filename(baseurl, filepath):
    """Downloads a remote file given the url, then caches it to the user's home
    folder.

    Args:
        baseurl: Url to the download path, excluding the file name
        filepath: The file path to download

    Returns:
        filename (str): A file path on the local file system corresponding to
        the data requested in data_name.
    """
    # Split data url and file name if the user provided a whole path in file_resources
    if validators.url(filepath):
        baseurl, filepath = os.path.split(filepath)
        baseurl = baseurl + "/"

    try:
        logging.debug("Fetching file from: {}{}, saving to {}".format(baseurl, filepath, openomics.config['cache_dir']))

        with data.conf.set_temp("dataurl", baseurl), data.conf.set_temp("remote_timeout", 30):
            return data.get_pkg_data_filename(filepath, package="openomics.database", show_progress=True)

    except (URLError, ValueError) as e:
        raise Exception(f"Unable to download file at {os.path.join(baseurl, filepath)}. "
                        f"Please try manually downloading the files and add path to `file_resources` arg. \n{e}")


def decompress_file(data_file: Path, filename: str, file_ext: Type) -> Tuple[str, Union[gzip.GzipFile, TextIO]]:
    output = data_file
    # This null if-clause is needed when filetype_ext is None, causing the next clauses to fail
    if file_ext is None:
        output = data_file

    elif file_ext.extension == "gz":
        logging.debug("Decompressed gzip file at {}".format(data_file))
        output = gzip.open(data_file, "rt")
        filename = filename.strip(".gz")

    elif file_ext.extension == "zip":
        logging.debug("Decompressed zip file at {}".format(data_file))
        zf = zipfile.ZipFile(data_file, "r")
        filename = filename.strip(".zip")

        for subfile in zf.infolist():
            # If the file extension matches
            if os.path.splitext(subfile.filename)[-1] == os.path.splitext(filename)[-1]:
                output = zf.open(subfile.filename, mode="r")

    elif file_ext.extension == "rar":
        logging.debug("Decompressed rar file at {}".format(data_file))
        rf = rarfile.RarFile(data_file, "r")
        filename = filename.strip(".rar")

        for subfile in rf.infolist():
            # If the file extension matches
            if os.path.splitext(subfile.filename)[-1] == os.path.splitext(filename)[-1]:
                output = rf.open(subfile.filename, mode="r")
    else:
        print(f"WARNING: filepath_ext.extension {file_ext.extension} not supported.")
        output = data_file

    return filename, output

def read_db(path, table, index_col):
    """
    Args:
        path:
        table:
        index_col:
    """
    engine = sa.create_engine(path)
    conn = engine.connect()
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


def get_decompressed_text_gzip(gzip_file):
    """
    Args:
        gzip_file:
    """
    # compressedFile = StringIO()
    # compressedFile.write(gzip_file.read())
    # compressedFile.seek(0)
    return io.TextIOWrapper(gzip_file)
    # decompressedFile = gzip.GzipFile(fileobj=gzip_file, mode='rb')
    # return decompressedFile
