import errno
import os

import requests
from requests.adapters import HTTPAdapter


def mkdirs(outdir):
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass


def retry(num=5):
    """"retry connection.

        define max tries num
        if the backoff_factor is 0.1, then sleep() will sleep for
        [0.1s, 0.2s, 0.4s, ...] between retries.
        It will also force a retry if the status code returned is 500, 502, 503 or 504.

    """
    s = requests.Session()
    retries = Retry(total=num, backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])
    s.mount('http://', HTTPAdapter(max_retries=retries))

    return s