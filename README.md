# Overview

The purpose of this code is to enable reproduction
and facilitate extension of the computational
results associated with Ref. [1].

For step-by-step instructions on reproducing the figures in Ref. [1],
see [REPRODUCEME.md](REPRODUCEME.md).

# Installation and Testing

Python can be obtained at [python.org](https://python.org).
The dependencies required can be installed via

```bash
pip install -r requirements.txt
```

To test the installation, navigate to the parent directory and invoke

```bash
python -m doctest src/wells.py REPRODUCEME.md -v > test/results.txt
```

# Data

The raw data can be found in the directory [data/](data/).

# Documentation

The documentation can be found in [doc/](doc/).
A pdf version of the documentation is available [at doc/manual.pdf](doc/manual.pdf).
Other forms of documentation can be built using [sphinx](https://www.sphinx-doc.org).

# Citing This Work

To cite the manuscript, use Ref. [1].
To cite the software and experimental data, use Ref. [2].

.. [1]: DeJaco, R. F.; Majikes, J. M.; Liddle, J. A.; Kearsley, A. J. Binding, Brightness, or Noise? Extracting Temperature-dependent Properties of Dye Bound to DNA. *Under Review*, 2023.

.. [2]: DeJaco, R. F. Software and Data associated with *Temperature-dependent Thermodynamic and Photophysical Properties of SYTO-13 Dye Bound to DNA*, National Institute of Standards and Technology, 2023, https://doi.org/10.18434/mds2-2762.
