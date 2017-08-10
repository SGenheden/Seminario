Seminario
=========

A tool to compute bond and angle force field parameters with the Seminario method.

Description
===========
Seminario is a tool estimate force fields parameters using the method of Seminario:

J. M. Seminario (1996) Calculation of intramolecular force fields from second-derivative tensors. *Internat. J. Quant. Chem.* 60:1271-1277

There are three required inputs:
1. A structure of your model (e.g. in .gro or .pdb format)
2. A Gaussian09 formated checkpoint file with results from a frequency calculation (`freq`-keyword)
3. A list of bonds for which the equilibrium distance and force constant should be estimated

Installation
============
To install this tool, clone or download the git repository. Change to the downloaded directory and install the software with

```
python setup.py install
```

Examples
========

```
seminario_ff -f model.fchk -s model.gro -b :1@N-:3@CU :1@ND1-:3@CU :2@NE2-:3@CU
```
