import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command


# Package meta-data.
NAME = 'TCGAMultiOmics'
DESCRIPTION = 'A toolkit to download, integrate, and query multi-omics TCGA cancer data.'
URL = 'https://github.com/JonnyTran/TCGAMultiOmics'
EMAIL = 'nhat.tran@mavs.uta.edu'
AUTHOR = 'Jonny Tran'
REQUIRES_PYTHON = '>=2.7.*'
VERSION = '0.1'

REQUIRED = [
    'numpy', 'pandas', 'networkx', 'dask'
]

here = os.path.abspath(os.path.dirname(__file__))


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPi via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


setup(
    name=NAME,
    version=VERSION,
    packages=find_packages(),
    package_dir={NAME: 'TCGAMultiOmics'},
    url=URL,
    license='',
    python_requires=REQUIRES_PYTHON,
    install_requires=REQUIRED,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
# $ setup.py publish support.
    cmdclass={
        'upload': UploadCommand,
    },
    include_package_data=True

)
