QCManyBody
==========

<p align="left">
    <picture>
    <img alt="QCManyBody Logo" src="https://github.com/MolSSI/QCManyBody/blob/main/docs/logo.png" height="150px">
    </picture>
</p>

[![GitHub Actions](https://img.shields.io/github/actions/workflow/status/MolSSI/QCManyBody/ci.yml?logo=github)](https://github.com/MolSSI/QCManyBody/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/loriab/QCManyBody/graph/badge.svg?token=E4S0706HJ0)](https://codecov.io/gh/loriab/QCManyBody)
[![Documentation Status](https://img.shields.io/github/actions/workflow/status/MolSSI/QCManyBody/ci.yml?label=docs&logo=readthedocs&logoColor=white)](https://molssi.github.io/QCManyBody/)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/qcmanybody?color=blue&logo=anaconda&logoColor=white)](https://anaconda.org/conda-forge/qcmanybody)
![python](https://img.shields.io/badge/python-3.8+-blue.svg)

QCManyBody is a python package for running quantum chemistry many-body expansions and interaction calculations in a
package-independent way.

## Features

- **Flexible Many-Body Expansions**: Compute interaction energies, gradients, and Hessians at various n-body levels
- **BSSE Corrections**: Support for nocp, cp (counterpoise), and vmfc (Valiron-Mayer function counterpoise) corrections
- **Multi-Level Calculations**: Use different QC methods at different n-body levels for efficiency
- **Parallel Execution**: Built-in support for multiprocessing and MPI for HPC clusters
- **QC Program Integration**: Works with Psi4, NWChem, CFOUR, and other QCEngine-supported programs
- **QCSchema Compliance**: Full integration with QCElemental and QCEngine ecosystems
- **Command-Line Interface**: Run calculations from JSON/YAML input files without Python coding

## Quick Links

- **[Documentation](https://molssi.github.io/QCManyBody/)** - Complete guides and API reference
- **[Getting Started](https://molssi.github.io/QCManyBody/getting_started.html)** - Your first calculation
- **[Examples](qcmanybody/tests/test_examples.py)** - Working code examples
- **[Parallel Execution](https://molssi.github.io/QCManyBody/parallel_execution_guide.html)** - Speed up calculations
- **[How-To Guides](https://molssi.github.io/QCManyBody/how-to-guides.html)** - Common tasks and recipes

## Installation

QCManyBody is available from [PyPI](https://pypi.org/project/qcmanybody) and from
[conda-forge](https://anaconda.org/conda-forge/qcmanybody).

```bash
# Installation from PyPI
pip install qcmanybody

# Installation from conda-forge
conda install -c conda-forge qcmanybody
```

To install the latest development version directly from
[GitHub](https://github.com/MolSSI/QCManyBody), you can use the following command:

```bash
pip install git+https://github.com/MolSSI/QCManyBody.git
```

### Optional Dependencies

```bash
# For parallel execution with MPI (optional)
pip install qcmanybody[mpi]

# For CLI with enhanced formatting (optional)
pip install qcmanybody[cli]

# For QC program execution
pip install qcengine
```

## Python API Quick Start

```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Define molecule with fragments
mol = Molecule(
    symbols=["He", "He", "He"],
    geometry=[[0, 0, 0], [0, 0, 2], [0, 0, 4]],
    fragments=[[0], [1], [2]],
)

# Create input specification
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {"max_nbody": 2, "bsse_type": ["cp"]},
        "specification": {
            "scf/sto-3g": {
                "program": "psi4",
                "model": {"method": "scf", "basis": "sto-3g"},
                "driver": "energy",
            }
        },
    },
)

# Run calculation (with optional parallel execution)
result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
)

# Access results
print(f"Interaction energy: {result.properties.cp_corrected_interaction_energy_through_2_body} Eh")
```

See the [Getting Started Guide](https://molssi.github.io/QCManyBody/getting_started.html) for a complete tutorial.

## Command-Line Interface (CLI)

QCManyBody provides a command-line interface for running many-body calculations without writing Python code.

### Quick Start

Create an input file (`input.json`):

```json
{
  "schema_name": "qcmanybody_cli_input",
  "schema_version": 1,
  "molecule": {
    "source": "inline",
    "inline": {
      "symbols": ["He", "He", "He"],
      "geometry": [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 0.0, 6.0]],
      "fragments": [[0], [1], [2]],
      "units": "angstrom"
    }
  },
  "calculation": {
    "type": "single",
    "single": {
      "driver": "energy",
      "method": "mp2",
      "basis": "cc-pvdz",
      "program": "psi4"
    }
  },
  "bsse": {
    "type": ["cp"]
  },
  "manybody": {
    "max_nbody": 3
  }
}
```

Run the calculation:

```bash
# Validate input
qcmanybody validate input.json

# Preview execution plan
qcmanybody plan input.json

# Run calculation
qcmanybody run input.json -o results.json
```

### Available Commands

- **`qcmanybody run`**: Execute many-body calculations
- **`qcmanybody plan`**: Show execution plan without running
- **`qcmanybody validate`**: Validate input files
- **`qcmanybody convert`**: Convert between JSON and YAML formats

See the [CLI User Guide](docs/cli_guide.md) for comprehensive documentation.

## Documentation

Full documentation is available at [https://molssi.github.io/QCManyBody/](https://molssi.github.io/QCManyBody/)

## Authors

* Asem Alenaizan, [@alenaizan](https://github.com/alenaizan), original Psi4 implementations of vmfc Hessians, multi-level, and embedded point charges
* Lori A. Burns, [@loriab](https://github.com/loriab), ManyBody QCSchema and high-level interface
* Benjamin P. Pritchard, [@bennybp](https://github.com/bennybp), core interface and QCArchive integration
* Daniel G. A. Smith, [@dgasmith](https://github.com/dgasmith), original Psi4 implementations of nocp, cp, and vmfc single-level e/g/H and distributed driver integration

## Citation [![doi](https://img.shields.io/badge/doi-10.1063/5.0231843-5077AB.svg)](https://doi.org/10.1063/5.0231843)

The journal article reference describing QCManyBody is:

L. A. Burns, C. D. Sherrill, B. P. Pritchard,
"QCManyBody: A Flexible Implementation of the Many-Body Expansion",
J. Chem. Phys. 161(15) 152501 (2024).

## Demonstration

An example of the core and high-level interfaces can be found in [test_examples](qcmanybody/tests/test_examples.py) with
directions in [tests/README](qcmanybody/tests/README.md).
