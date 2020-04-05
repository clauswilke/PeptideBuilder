# PeptideBuilder: A simple Python library to generate model peptides.

*Matthew Z. Tien, Dariya K. Sydykova, Austin G. Meyer, and Claus O. Wilke*

[![PyPI version](https://badge.fury.io/py/PeptideBuilder.svg)](https://badge.fury.io/py/PeptideBuilder)
![PyPI - Downloads](https://img.shields.io/pypi/dm/PeptideBuilder)
![PyPI - License](https://img.shields.io/pypi/l/PeptideBuilder)
[![Build Status](https://travis-ci.org/clauswilke/PeptideBuilder.svg?branch=master)](https://travis-ci.org/clauswilke/PeptideBuilder)
[![Coverage Status](https://img.shields.io/codecov/c/github/clauswilke/PeptideBuilder/master.svg)](https://codecov.io/github/clauswilke/PeptideBuilder?branch=master)

## Installation

You can install PeptideBuilder with pip:
```
pip install PeptideBuilder
```
PeptideBuilder has one required dependency: [Biopython](https://pypi.org/project/biopython/)


## Examples

For example usage, we encourage you to checkout the scripts in the `examples` folder and in the `tests` folder. The `examples` folder contains two scripts showing typical usage. The script `simpleExample.py` is a brief example script demonstrating basic use of the PeptideBuilder library. The script `evaluation.py` reproduces the results presented in Table 1 of Tien et al. (2013).

The file `test_PeptideBuilder.py` in `tests` contains extensive tests for the various functions provided by this library and may also be useful if you're looking for example usage.

## Misc

The software is provided to you under the MIT license (see file LICENSE).
The most up-to-date version of this software is available at
https://github.com/clauswilke/PeptideBuilder.

To test whether your installation works properly, run `pytest` in the top-level project folder.

## Contributing

Pull requests are welcome on GitHub. However, to be accepted, contributions must:
1. Be styled with [`black`](https://black.readthedocs.io/en/stable/)
2. Be linted with `pylint`
3. Be type-checked with `mypy`
4. Pass the `pytest` unit tests

Thus, before contributing code make sure the following commands exit without errors when run from the root directory of the Peptide Builder project:

- `pytest`
- `black .`
- `mypy PeptideBuilder/`
- `pylint --rcfile=setup.cfg PeptideBuilder/`

**Reference:**
M. Z. Tien, D. K. Sydykova, A. G. Meyer, C. O. Wilke (2013). PeptideBuilder:
A simple Python library to generate model peptides. PeerJ 1:e80.
