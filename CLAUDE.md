# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QCManyBody is a Python package for running quantum chemistry many-body expansions and interaction calculations in a package-independent way. It provides both a Python API and a command-line interface for computing many-body expansions with various BSSE (Basis Set Superposition Error) treatments.

## Common Commands

### Testing

```bash
# Run all tests
pytest -v qcmanybody/

# Run tests without external QC programs
pytest -v qcmanybody/ -m 'not addon'

# Run specific test categories
pytest -v qcmanybody/ -k "examples"              # Example tests only
pytest -v qcmanybody/ -k "not (3b or 4b)"        # Skip expensive 3/4-body tests

# Run single test file
pytest -v qcmanybody/tests/test_utils.py

# Run tests with coverage
pytest -v qcmanybody/ --cov=qcmanybody --cov-report=html
```

### Code Quality

```bash
# Format code (excludes test files)
black --line-length=120 qcmanybody/

# Sort imports
isort --profile black --line-length=120 qcmanybody/

# Type checking
mypy qcmanybody/

# Pre-commit hooks
pre-commit run --all-files
```

### CLI Commands

```bash
# Validate input file
qcmanybody validate input.json

# Show execution plan
qcmanybody plan input.json

# Run calculation
qcmanybody run input.json -o results.json

# Add verbosity
qcmanybody -v run input.json
qcmanybody -vvv run input.json  # Maximum verbosity

# Convert between formats
qcmanybody convert input.json -o input.yaml
```

### Building/Installing

```bash
# Install in development mode
pip install -e .

# Install with optional dependencies
pip install -e .[cli]      # CLI support (YAML, rich formatting)
pip install -e .[mpi]      # MPI parallel support
pip install -e .[tests]    # Testing dependencies
pip install -e .[dev]      # Development tools
```

### Documentation

```bash
# Build documentation
QCMANYBODY_MAX_NBODY=5 PYTHONPATH=docs/extensions mkdocs build

# Serve documentation locally
mkdocs serve
```

## Architecture

### Core Components

**ManyBodyCore** (`qcmanybody/core.py`):
- Low-level interface for many-body calculations
- Handles BSSE treatment logic (CP, NOCP, VMFC)
- Builds compute task lists via `builder.py`
- Manages multilevel calculations (different QC methods at different n-body levels)
- Returns raw computational results without orchestrating QC execution

**ManyBodyComputer** (`qcmanybody/computer.py`):
- High-level interface that orchestrates full calculations
- Wraps `ManyBodyCore` and handles QC program execution via QCEngine
- Primary entry point: `ManyBodyComputer.from_manybodyinput()`
- Manages `AtomicComputer` objects for individual QC tasks
- Converts results to `ManyBodyResult` with full provenance

**builder.py**:
- Generates lists of computations needed for BSSE treatments
- Key function: `build_nbody_compute_list()`
- Returns dictionary mapping BSSE types to sets of (fragment, basis) index tuples
- Handles supersystem-only calculations

**dependency.py**:
- `NBodyDependencyGraph` class for task dependency management
- Used by parallel execution to determine task ordering and parallelization opportunities

### Data Models

Located in `qcmanybody/models/`:
- **v1/**: Pydantic v1 models (current production)
  - `ManyBodyInput`: Input specification schema
  - `ManyBodyResult`: Output results schema
  - `ManyBodyKeywords`: Configuration options
  - `BsseEnum`: BSSE treatment types (cp, nocp, vmfc)
- **v2/**: Pydantic v2 models (in development, forward compatibility)

Models use QCElemental's `Molecule` and `AtomicResult` types for interoperability.

### Parallel Execution

**qcmanybody/parallel.py** and **qcmanybody/parallel/** module:
- `ParallelManyBodyComputer`: Extension of ManyBodyComputer with parallel support
- Multiple executor backends:
  - `SequentialExecutor`: Reference implementation, no parallelism
  - `MultiprocessingExecutor`: Single-node, multi-core via `multiprocessing`
  - `MPIExecutor`: Multi-node HPC via MPI (requires mpi4py)
- `TaskScheduler`: Priority-based task scheduling
- `CheckpointManager`: Fault tolerance for long calculations
- See `qcmanybody/parallel/README.md` for detailed documentation

Typical usage:
```python
result = ManyBodyComputer.from_manybodyinput(input_model, parallel=True, n_workers=4)
```

### CLI Module

Located in `qcmanybody/cli/`:
- **main.py**: Entry point, argument parsing
- **input_parser.py**: Parses JSON/YAML input files
- **molecule_loader.py**: Loads molecules from various sources (inline, XYZ files, PDB, etc.)
- **converter.py**: Converts between JSON and YAML formats
- **commands/**: Individual command implementations (run, plan, validate, convert)
- **schemas/**: CLI-specific input schemas

CLI input schema uses `qcmanybody_cli_input` schema with structured molecule/calculation/bsse/manybody sections.

### Utilities

**utils.py**:
- `labeler()` / `delabeler()`: Convert between string and tuple labels
- `resize_gradient()` / `resize_hessian()`: Adapt derivatives to full system size
- `collect_vars()`: Gather variables across n-body levels
- `print_nbody_energy()`: Format energy output tables

## Key Concepts

### BSSE Treatments
- **CP (Counterpoise)**: Ghost functions for all fragments in all calculations
- **NOCP (No Counterpoise)**: No ghost functions
- **VMFC (Valiron-Mayer Function Counterpoise)**: Selective ghost functions

### Multilevel Calculations
Different QC methods can be used at different n-body levels:
```python
levels = {
    1: "ccsd(t)/cc-pvtz",  # 1-body at high level
    2: "mp2/cc-pvdz",       # 2-body at medium level
    3: "hf/cc-pvdz"         # 3-body at low level
}
```

### Fragment Indexing
- Fragments are 0-indexed in the molecule specification
- Compute tasks use tuples like `((1, 2), (1, 2, 3))` meaning fragments {1,2} with basis {1,2,3}
- The `labeler()` function converts tuples to strings like `"(1, 2)-(1, 2, 3)"`

## Testing Structure

See `qcmanybody/tests/README.md` for details. Key categories:

1. **Examples** (`test_examples.py`): Easy-to-read demonstrations of core/high-level API
2. **Unit Tests**: Fast, no external dependencies (test_utils.py, test_schema_keywords.py)
3. **Static-Data Regression**: Pre-computed results, requires zstandard (test_core_*.py)
4. **End-to-End**: Actual QC calculations, requires QCEngine + QC program (test_computer_*.py)

Test markers:
- `@pytest.mark.addon`: Requires external software
- `@pytest.mark.qcengine`: Needs QCEngine
- `@pytest.mark.psi4`, `@pytest.mark.nwchem`, etc.: Specific QC program required

## Development Conventions

### Code Style
- Line length: 120 characters
- Black formatter with `--line-length=120`
- isort with black profile
- Test files (test_*.py) excluded from black formatting

### Pydantic Compatibility
The codebase supports both Pydantic v1 and v2:
- Production code uses v1 API (imported as `from pydantic.v1 import ...`)
- v2 models exist in `models/v2/` for forward compatibility
- Comments with `# v2:` mark code for eventual v2 migration

### Versioning
- Uses `setuptools-scm` for automatic versioning from git tags
- Version available via `qcmanybody.__version__`

### Environment Variables
- `QCMANYBODY_EMBEDDING_CHARGES=1`: Enable experimental embedding charges for EE-MBE
- `QCMANYBODY_MAX_NBODY=5`: Set max n-body for documentation generation

## Special Considerations

### Molecule Units
- QCElemental expects bohr by default
- CLI accepts "angstrom" or "bohr" in input files
- Conversion happens in `molecule_loader.py:70` (be careful with unit handling)

### Git Workflow
- Main branch: `main`
- Use feature branches for development
- CI runs on PRs to main
- Multiple test configurations in `.github/workflows/ci.yml` covering Python 3.8-3.13

### Files Pending Removal
The `files_for_removal/` directory contains old development artifacts tagged for cleanup. Ignore this directory.

## External Dependencies

**Required**:
- numpy
- pydantic (v1 or v2)
- qcelemental (>=0.28.0, <0.70.0)

**Optional**:
- qcengine: For high-level interface and QC program execution
- pyyaml: For YAML input file support in CLI
- rich: For enhanced terminal output in CLI
- mpi4py: For MPI parallel execution
- zstandard: For compressed test data

**QC Programs** (at least one needed for actual calculations):
- Psi4 (recommended)
- NWChem
- CFOUR

## Common Patterns

### Creating a Calculation
```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Define molecule with fragments
mol = Molecule(symbols=["He", "He", "He"],
               geometry=[[0,0,0], [0,0,3], [0,0,6]],
               fragments=[[0], [1], [2]],
               units="angstrom")

# Create input
mbin = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {"bsse_type": ["cp"], "max_nbody": 3},
        "specification": {
            "mp2/cc-pvdz": {"program": "psi4",
                           "model": {"method": "mp2", "basis": "cc-pvdz"}}
        }
    }
)

# Execute
result = ManyBodyComputer.from_manybodyinput(mbin)
```

### Working with Results
```python
# Access final result
interaction_energy = result.return_result

# Access component energies
for nbody in [1, 2, 3]:
    energy = result.nbody_energy[nbody]

# Access individual task results
for label, atomic_result in result.properties.calcinfo_natoms.items():
    # atomic_result is AtomicResult from QCElemental
    pass
```
