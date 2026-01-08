# Getting Started with QCManyBody

This guide will help you get started with QCManyBody for quantum chemistry many-body expansion calculations.

## Installation

QCManyBody is available from PyPI and conda-forge:

```bash
# Install from PyPI
pip install qcmanybody

# Or install from conda-forge
conda install -c conda-forge qcmanybody
```

For parallel execution support with MPI:

```bash
pip install qcmanybody[mpi]
```

For CLI support with enhanced formatting:

```bash
pip install qcmanybody[cli]
```

Development version from GitHub:

```bash
pip install git+https://github.com/MolSSI/QCManyBody.git
```

## Your First Calculation

Here's a simple example to get you started with a water dimer energy calculation:

```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Define a water dimer with fragment specification
water_dimer = Molecule(
    symbols=["O", "H", "H", "O", "H", "H"],
    geometry=[
        [0.000, 0.000, 0.000],    # O atom (fragment 1)
        [0.757, 0.586, 0.000],    # H atom (fragment 1)
        [-0.757, 0.586, 0.000],   # H atom (fragment 1)
        [3.000, 0.000, 0.000],    # O atom (fragment 2)
        [3.757, 0.586, 0.000],    # H atom (fragment 2)
        [2.243, 0.586, 0.000],    # H atom (fragment 2)
    ],
    fragments=[[0, 1, 2], [3, 4, 5]],  # Define two water molecules
    molecular_charge=0.0,
    molecular_multiplicity=1,
)

# Create input specification for 2-body expansion
manybodyinput = ManyBodyInput(
    molecule=water_dimer,
    specification={
        "driver": "energy",
        "keywords": {
            "bsse_type": ["cp"],  # Counterpoise correction
            "max_nbody": 2,       # Up to 2-body interactions
        },
        "specification": {
            "scf/sto-3g": {
                "program": "psi4",
                "model": {"method": "scf", "basis": "sto-3g"},
                "driver": "energy",
            }
        },
    },
)

# Run the calculation
result = ManyBodyComputer.from_manybodyinput(manybodyinput)

# Access results
print(f"Total energy: {result.return_result}")
print(f"Interaction energy: {result.properties.cp_corrected_interaction_energy_through_2_body}")
```

## Understanding the Results

The `result` object contains various properties:

- `return_result`: The final computed property (energy, gradient, or hessian)
- `properties`: Contains detailed many-body analysis results
  - Interaction energies at different n-body levels
  - Total energies with different corrections (CP, VMFC, noCP)
  - Individual fragment contributions

Example accessing specific properties:

```python
# Get 2-body interaction energy with CP correction
ie_2body = result.properties.cp_corrected_interaction_energy_through_2_body

# Get total energy through 2-body with CP correction
total_2body = result.properties.cp_corrected_total_energy_through_2_body

# Access all available properties
print(result.properties.dict().keys())
```

## Common Calculation Types

### Energy Calculation

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {"max_nbody": 3, "bsse_type": ["cp"]},
        "specification": {"method/basis": {...}},
    },
)
```

### Gradient Calculation

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "gradient",
        "keywords": {"max_nbody": 2, "bsse_type": ["nocp"]},
        "specification": {"method/basis": {...}},
    },
)
```

### Multi-Level Calculations

Use different methods at different n-body levels:

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "levels": {
                1: "ccsd/tz",      # 1-body at CCSD/cc-pVTZ
                2: "mp2/dz",       # 2-body at MP2/cc-pVDZ
            },
            "bsse_type": ["cp"],
        },
        "specification": {
            "ccsd/tz": {
                "program": "psi4",
                "model": {"method": "ccsd", "basis": "cc-pvtz"},
                "driver": "energy",
            },
            "mp2/dz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"},
                "driver": "energy",
            },
        },
    },
)
```

## BSSE Correction Methods

QCManyBody supports three BSSE correction schemes:

- **`nocp`** (No counterpoise): Standard MBE without basis set corrections
- **`cp`** (Counterpoise corrected): Full counterpoise correction
- **`vmfc`** (Valiron-Mayer function counterpoise): Advanced correction scheme

```python
# Request multiple correction schemes
"keywords": {
    "bsse_type": ["nocp", "cp", "vmfc"],
    "max_nbody": 2,
}
```

## Using the Command-Line Interface

QCManyBody also provides a CLI for running calculations from input files:

Create an input file `input.json`:

```json
{
  "molecule": {
    "source": "inline",
    "inline": {
      "symbols": ["He", "He", "He"],
      "geometry": [[0, 0, 0], [0, 0, 2], [0, 0, 4]],
      "fragments": [[0], [1], [2]]
    }
  },
  "calculation": {
    "type": "single",
    "single": {
      "driver": "energy",
      "method": "scf",
      "basis": "sto-3g",
      "program": "psi4"
    }
  },
  "manybody": {
    "max_nbody": 3,
    "bsse_type": ["cp"]
  }
}
```

Run the calculation:

```bash
qcmanybody run input.json -o results.json
```

See the [CLI Guide](cli_guide.md) for more details.

## Parallel Execution

For larger calculations, enable parallel execution:

```python
from qcmanybody import ManyBodyComputer

# Simple parallel execution with 4 workers
result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
)
```

For advanced parallel execution control:

```python
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

config = ExecutorConfig(n_workers=8, timeout_per_task=3600)
executor = MultiprocessingExecutor(config)

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    executor=executor,
)
```

See the [Parallel Execution Guide](parallel_execution_guide.md) for details on MPI and advanced features.

## Next Steps

Now that you've run your first calculation, explore:

- [High-Level Interface](high-level-interface.md) - Complete interface details
- [Core Interface](core-interface.md) - Lower-level control
- [CLI Guide](cli_guide.md) - Using the command-line interface
- [Parallel Execution Guide](parallel_execution_guide.md) - Speed up calculations
- [How-To Guides](how-to-guides.md) - Common tasks and recipes
- [Keywords](keywords.md) - All available options

## Common Issues

### QCEngine Not Found

If you get an error about QCEngine not being available:

```bash
pip install qcengine
```

### Psi4 Not Available

Install Psi4 via conda:

```bash
conda install -c conda-forge psi4
```

### Fragment Specification Errors

Ensure your molecule fragments are correctly specified:
- Fragment indices must be 0-based
- All atoms must be assigned to exactly one fragment
- Fragments should not overlap

```python
# Correct fragment specification for 3 atoms
fragments=[[0], [1], [2]]  # Each atom in its own fragment

# Or for a dimer
fragments=[[0, 1, 2], [3, 4, 5]]  # Two 3-atom fragments
```

### Memory Issues

For large systems, use parallel execution with memory limits:

```python
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

config = ExecutorConfig(
    n_workers=4,
    timeout_per_task=7200,  # 2 hours per task
)
executor = MultiprocessingExecutor(config)
```

## Getting Help

- Documentation: [https://molssi.github.io/QCManyBody/](https://molssi.github.io/QCManyBody/)
- GitHub Issues: [https://github.com/MolSSI/QCManyBody/issues](https://github.com/MolSSI/QCManyBody/issues)
- API Reference: [API Documentation](api.md)
