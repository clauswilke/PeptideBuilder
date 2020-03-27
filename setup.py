from setuptools import setup


INSTALL_REQUIRES = [
    'Biopython'
]

TEST_REQUIRES = [
    # testing and coverage
    'pytest', 'coverage', 'pytest-cov',
    # to be able to run `python setup.py checkdocs`
    'collective.checkdocs', 'pygments',
]


with open("README.md", "r") as f:
    long_description = f.read()

with open("PeptideBuilder/__init__.py", "r") as f:
    init = f.readlines()

for line in init:
    if '__version__' in line:
        __version__ = line.split("'")[-2]

setup(name = 'PeptideBuilder', 
    version = __version__, 
    author = 'Matthew Z. Tien', 
    author_email = 'Matthew.Tien89@gmail.com', 
    description = 'Create peptide PDB files with specified geometry',
    long_description=long_description,
    long_description_content_type='ext/markdown',
    url = 'https://github.com/mtien/PeptideBuilder',
    download_url = 'https://github.com/mtien/PeptideBuilder/releases',
    platforms = 'Tested on Mac OS X and Windows 10',
    packages = ['PeptideBuilder'],
    install_requires = INSTALL_REQUIRES,
    extras_require = {
        'test': TEST_REQUIRES + INSTALL_REQUIRES,
    },
)
