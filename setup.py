import sys
from os import path
from setuptools import setup, find_packages

from omega import __version__

DESCRIPTION = "Omega is a dNdS method to identify positive selection cancer drivers"

# Check the python compatibility
if sys.hexversion < 0x03050000:
    raise RuntimeError('This package requires Python 3.5 or later.')


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    install_requires = f.read().splitlines()


# Get the long description from the README file
with open(path.join(directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name="omega",
    python_requires='>=3.10.12',
    version=__version__,
    packages=find_packages(),
    author='BBGLab (Barcelona Biomedical Genomics Lab)',
    author_email='bbglab@irbbarcelona.org',
    description=DESCRIPTION,
    license="AGPLv3",
    keywords="",
    url="https://github.com/bbglab/omega",
    # download_url="https://github.com/bbglab/omega/get/" + __version__ + ".tar.gz",
    long_description=long_description,
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'omega = omega.main:omega'
        ]
    }
)
