# High-Level Interface

The high-level interface provides a comprehensive workflow for many-body expansion calculations using [QCSchema](qcschema.md) data structures. It handles the complete calculation pipeline from input specification to final results, with built-in support for running quantum chemistry calculations through [QCEngine](https://github.com/MolSSI/QCEngine).

## Overview

The high-level interface is built around two main components:

1. **`ManyBodyInput`** - QCSchema-based input specification
2. **`ManyBodyComputer`** - Computation engine that processes the input and returns results

The high-level interface uses the [core interface](core-interface.md) under the hood but provides additional conveniences:

- Automatic QCEngine integration for running QC calculations
- QCSchema validation for inputs and outputs
- Built-in support for multiple QC programs (Psi4, NWChem, CFOUR, etc.)
- Parallel execution support
- Comprehensive result packaging

## When to Use the High-Level Interface

Use the high-level interface when:

- You want a complete, turnkey solution for many-body calculations
- You're using QCEngine to run quantum chemistry programs
- You need QCSchema-compliant input/output for integration with other tools
- You want automatic validation of input specifications
- You're running calculations through QCFractal or similar systems

Consider the [core interface](core-interface.md) if you need more control over how individual calculations are executed.

## Quick Start

Here's a basic example of a 2-body energy calculation on a neon trimer:

```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Define molecule with fragments
mol = Molecule(
    symbols=["Ne", "Ne", "Ne"],
    geometry=[[0, 0, 0], [0, 0, 2], [0, 0, 4]],
    fragments=[[0], [1], [2]],
)

# Create input specification
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "max_nbody": 2,
            "bsse_type": ["cp"],
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

# Run calculation
result = ManyBodyComputer.from_manybodyinput(manybodyinput)

# Access results
print(f"Total energy: {result.return_result}")
print(f"Interaction energy: {result.properties.cp_corrected_interaction_energy_through_2_body}")
```

## Detailed Examples

### Basic Energy Calculation

A simple counterpoise-corrected energy calculation:

```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Water dimer
water_dimer = Molecule(
    symbols=["O", "H", "H", "O", "H", "H"],
    geometry=[
        [0.000, 0.000, 0.000],
        [0.757, 0.586, 0.000],
        [-0.757, 0.586, 0.000],
        [3.000, 0.000, 0.000],
        [3.757, 0.586, 0.000],
        [2.243, 0.586, 0.000],
    ],
    fragments=[[0, 1, 2], [3, 4, 5]],
)

manybodyinput = ManyBodyInput(
    molecule=water_dimer,
    specification={
        "driver": "energy",
        "keywords": {
            "max_nbody": 2,
            "bsse_type": ["nocp", "cp", "vmfc"],
        },
        "specification": {
            "mp2/aug-cc-pvdz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "aug-cc-pvdz"},
                "driver": "energy",
            }
        },
    },
)

result = ManyBodyComputer.from_manybodyinput(manybodyinput)

# Compare different BSSE corrections
print(f"No CP: {result.properties.nocp_corrected_interaction_energy_through_2_body}")
print(f"CP:    {result.properties.cp_corrected_interaction_energy_through_2_body}")
print(f"VMFC:  {result.properties.vmfc_corrected_interaction_energy_through_2_body}")
```

### Multi-Level Calculation

Use different methods at different n-body levels for efficiency:

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "levels": {
                1: "ccsd/tz",      # Accurate for monomers
                2: "mp2/dz",       # Faster for dimers
                3: "mp2/dz",       # Same for trimers
            },
            "bsse_type": ["cp"],
        },
        "specification": {
            "ccsd/tz": {
                "program": "nwchem",
                "model": {"method": "ccsd", "basis": "cc-pvtz"},
                "driver": "energy",
                "keywords": {"ccsd__thresh": 1.e-7},
            },
            "mp2/dz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"},
                "driver": "energy",
            },
        },
    },
)

result = ManyBodyComputer.from_manybodyinput(manybodyinput)
```

### Gradient Calculation

Compute many-body gradients for geometry optimization or dynamics:

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "gradient",
        "keywords": {
            "max_nbody": 2,
            "bsse_type": ["cp"],
        },
        "specification": {
            "scf/6-31g": {
                "program": "psi4",
                "model": {"method": "scf", "basis": "6-31g"},
                "driver": "gradient",
            }
        },
    },
)

result = ManyBodyComputer.from_manybodyinput(manybodyinput)

# result.return_result is now a gradient (Natom x 3 array)
print(result.return_result.shape)  # (Natom, 3)
```

### Supersystem Calculation

Include all higher-order terms exactly using the supersystem:

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "levels": {
                1: "mp2/dz",
                2: "mp2/dz",
                "supersystem": "ccsd/tz",  # Full system at high level
            },
            "bsse_type": ["cp"],
        },
        "specification": {
            "mp2/dz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"},
                "driver": "energy",
            },
            "ccsd/tz": {
                "program": "psi4",
                "model": {"method": "ccsd", "basis": "cc-pvtz"},
                "driver": "energy",
            },
        },
    },
)

result = ManyBodyComputer.from_manybodyinput(manybodyinput)
```

### With Parallel Execution

Enable parallel execution for faster calculations:

```python
# Simple parallel execution
result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
)

# Advanced parallel control
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

config = ExecutorConfig(
    n_workers=8,
    timeout_per_task=3600,
    checkpoint_file="checkpoint.json",
)
executor = MultiprocessingExecutor(config)

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    executor=executor,
)
```

See the [Parallel Execution Guide](parallel_execution_guide.md) for details.

## Input Specification

### ManyBodyInput Structure

The `ManyBodyInput` has three main sections:

```python
ManyBodyInput(
    molecule=...,        # QCElemental Molecule object
    specification={      # Calculation specification
        "driver": ...,
        "keywords": {...},
        "specification": {...},
    },
)
```

### Molecule Specification

The molecule must be a `qcelemental.models.Molecule` with fragments defined:

```python
from qcelemental.models import Molecule

mol = Molecule(
    symbols=["He", "He", "He"],
    geometry=[[0, 0, 0], [0, 0, 2], [0, 0, 4]],  # Bohr by default
    fragments=[[0], [1], [2]],                    # 0-based indices
    molecular_charge=0.0,
    molecular_multiplicity=1,
)
```

**Fragment Requirements:**
- Each atom must belong to exactly one fragment
- Fragment indices are 0-based
- Fragments can have different sizes
- Fragment charges/multiplicities can be specified

See [Molecule Input](keywords.md#molecule) for details.

### Driver

Specifies what property to compute:

- `"energy"` - Electronic energy
- `"gradient"` - Energy gradient (forces)
- `"hessian"` - Energy Hessian (second derivatives)

### Keywords

Many-body expansion options (see [Keywords](keywords.md) for full list):

```python
"keywords": {
    "max_nbody": 3,                    # Maximum n-body level
    "bsse_type": ["nocp", "cp"],      # BSSE correction schemes
    "return_total_data": True,         # Include total properties
    "levels": {1: "method1", 2: "method2"},  # Multi-level specification
    "supersystem_ie_only": False,      # Use supersystem for IE only
}
```

### Specification

Defines the quantum chemistry methods to use:

```python
"specification": {
    "label": {                         # Arbitrary label
        "program": "psi4",            # QC program
        "model": {
            "method": "mp2",          # QC method
            "basis": "cc-pvdz",       # Basis set
        },
        "driver": "energy",           # What this spec computes
        "keywords": {...},            # Program-specific keywords
    }
}
```

Multiple specifications can be defined for multi-level calculations.

## Output: ManyBodyResult

The `from_manybodyinput()` method returns a `ManyBodyResult` object:

```python
result = ManyBodyComputer.from_manybodyinput(manybodyinput)
```

### Main Result

```python
result.return_result  # Final computed property (energy, gradient, or hessian)
```

### Properties

Detailed many-body analysis in `result.properties`:

```python
# Interaction energies
result.properties.cp_corrected_interaction_energy_through_2_body
result.properties.nocp_corrected_interaction_energy_through_3_body

# Total energies
result.properties.cp_corrected_total_energy_through_2_body

# Individual n-body contributions
result.properties.cp_corrected_2_body_contribution_to_interaction_energy
```

See [Results](results.md) for complete documentation.

### Component Results

Access individual task results:

```python
result.nbody_components  # Dict of individual calculation results
```

### Provenance

Metadata about the calculation:

```python
result.provenance  # Software versions, timing, etc.
```

## API Reference

### ManyBodyComputer.from_manybodyinput()

::: qcmanybody.computer.ManyBodyComputer.from_manybodyinput
    options:
      show_source: false

### ManyBodyInput

::: qcmanybody.models.ManyBodyInput
    options:
      show_source: false

### ManyBodyResult

::: qcmanybody.models.ManyBodyResult
    options:
      show_source: false

## Common Patterns

### Running Without QCEngine

If you want to handle QC calculations yourself, use `build_tasks=False`:

```python
computer = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    build_tasks=False,
)

# computer.task_list contains the required calculations
for task in computer.task_list.values():
    # Run task yourself
    # Store results in computer.results
    pass

# Analyze results
computer.analyze()
```

See the [Core Interface](core-interface.md) for more control.

### Integration with Psi4

Psi4 uses QCManyBody by inheriting from `ManyBodyComputer`:

```python
from qcmanybody import ManyBodyComputer

class Psi4ManyBodyComputer(ManyBodyComputer):
    def compute_nbody_components(self):
        # Run calculations using Psi4's internal methods
        pass
```

This allows Psi4 to use QCManyBody's validation and structure while controlling execution.

### Checkpointing and Resumption

Save progress for long calculations:

```python
from qcmanybody.parallel import CheckpointManager

checkpoint_mgr = CheckpointManager(checkpoint_file="calculation.chk")

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
    checkpoint_manager=checkpoint_mgr,
)

# If interrupted, resume by running the same command
# Completed tasks will be loaded from checkpoint
```

## Troubleshooting

### Fragment Specification Errors

**Error:** "Not all atoms assigned to fragments"

**Solution:** Ensure all atoms are included in exactly one fragment:

```python
# For 6 atoms, all must be assigned
fragments=[[0, 1, 2], [3, 4, 5]]  # Correct
fragments=[[0, 1], [3, 4, 5]]     # Error: atom 2 missing
```

### QCEngine Not Available

**Error:** "QCEngine is required for high-level interface"

**Solution:** Install QCEngine:

```bash
pip install qcengine
```

### Program Not Found

**Error:** "Program 'psi4' not found"

**Solution:** Install the quantum chemistry program:

```bash
conda install -c conda-forge psi4
```

### Memory Issues

For large systems with many fragments, consider:

- Use multi-level calculations (lower levels for higher n-body)
- Enable parallel execution with appropriate worker limits
- Use `supersystem_ie_only=True` for high-order terms

## Next Steps

- [Core Interface](core-interface.md) - Lower-level control
- [Keywords](keywords.md) - All available options
- [Results](results.md) - Understanding output
- [Parallel Execution Guide](parallel_execution_guide.md) - Speed up calculations
- [How-To Guides](how-to-guides.md) - Common recipes

## Working Example

See `test_highlevel_interface_example()` in [test_examples.py](https://github.com/MolSSI/QCManyBody/blob/main/qcmanybody/tests/test_examples.py) for a complete working example.
