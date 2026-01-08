# How-To Guides

Practical guides for common tasks with QCManyBody.

## Table of Contents

1. [Set Up a Basic Calculation](#set-up-a-basic-calculation)
2. [Use Parallel Execution](#use-parallel-execution)
3. [Configure BSSE Corrections](#configure-bsse-corrections)
4. [Use the CLI](#use-the-cli)
5. [Integrate with Psi4/NWChem](#integrate-with-psi4nwchem)
6. [Checkpoint and Resume Calculations](#checkpoint-and-resume-calculations)
7. [Debug Calculation Failures](#debug-calculation-failures)
8. [Optimize Performance](#optimize-performance)

---

## Set Up a Basic Calculation

### Goal
Run a simple 2-body energy calculation on a molecular cluster.

### Steps

**1. Define your molecule with fragments:**

```python
from qcelemental.models import Molecule

# Water dimer example
mol = Molecule(
    symbols=["O", "H", "H", "O", "H", "H"],
    geometry=[
        [0.000, 0.000, 0.000],
        [0.757, 0.586, 0.000],
        [-0.757, 0.586, 0.000],
        [3.000, 0.000, 0.000],
        [3.757, 0.586, 0.000],
        [2.243, 0.586, 0.000],
    ],
    fragments=[[0, 1, 2], [3, 4, 5]],  # Two water molecules
)
```

**2. Create input specification:**

```python
from qcmanybody.models import ManyBodyInput

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
```

**3. Run the calculation:**

```python
from qcmanybody import ManyBodyComputer

result = ManyBodyComputer.from_manybodyinput(manybodyinput)

# Access results
print(f"Interaction energy: {result.properties.cp_corrected_interaction_energy_through_2_body} Eh")
```

### Tips
- Start with small basis sets (sto-3g) for testing
- Use `max_nbody=2` for initial calculations
- Always specify fragments correctly

---

## Use Parallel Execution

### Goal
Speed up calculations by running tasks in parallel.

### Basic Parallel Execution

**Enable with simple flag:**

```python
result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
)
```

### Advanced Parallel Control

**Configure executor settings:**

```python
from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig

config = ExecutorConfig(
    n_workers=8,
    timeout_per_task=3600,  # 1 hour per task
    max_retries=2,
)
executor = MultiprocessingExecutor(config)

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    executor=executor,
)
```

### MPI for HPC Clusters

**Use MPI for distributed computing:**

```python
from qcmanybody.parallel import MPIExecutor, ExecutorConfig

config = ExecutorConfig(
    n_workers=None,  # Auto-detect from MPI communicator
    timeout_per_task=7200,
)
executor = MPIExecutor(config)

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    executor=executor,
)
```

**SLURM batch script:**

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00

module load python mpi4py psi4

mpirun -n 64 python my_calculation.py
```

### Performance Tips
- Use 1 worker per physical core (not threads)
- For QC programs, account for their internal parallelism
- Monitor memory usage per worker
- Use checkpointing for long calculations

See [Parallel Execution Guide](parallel_execution_guide.md) for details.

---

## Configure BSSE Corrections

### Goal
Understand and apply different BSSE correction schemes.

### Available Corrections

**1. No Counterpoise (nocp):**
- Standard MBE without basis set corrections
- Fastest but includes BSSE
- Use for large basis sets where BSSE is negligible

```python
"keywords": {
    "bsse_type": ["nocp"],
}
```

**2. Counterpoise (cp):**
- Full counterpoise correction
- Most common choice
- Corrects for BSSE but may overcorrect

```python
"keywords": {
    "bsse_type": ["cp"],
}
```

**3. Valiron-Mayer Function Counterpoise (vmfc):**
- Advanced correction scheme
- More accurate than CP for some systems
- More expensive computationally

```python
"keywords": {
    "bsse_type": ["vmfc"],
}
```

### Compare All Three

**Request multiple corrections in one calculation:**

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "max_nbody": 2,
            "bsse_type": ["nocp", "cp", "vmfc"],  # All three
        },
        "specification": {...},
    },
)

result = ManyBodyComputer.from_manybodyinput(manybodyinput)

# Compare results
print("Interaction Energies:")
print(f"  NoCP: {result.properties.nocp_corrected_interaction_energy_through_2_body}")
print(f"  CP:   {result.properties.cp_corrected_interaction_energy_through_2_body}")
print(f"  VMFC: {result.properties.vmfc_corrected_interaction_energy_through_2_body}")
```

### Recommendations
- **Small basis sets (e.g., 6-31G)**: Use CP or VMFC
- **Large basis sets (e.g., aug-cc-pV5Z)**: NoCP may suffice
- **Publication quality**: Compare CP and VMFC
- **Testing**: Use NoCP for speed

---

## Use the CLI

### Goal
Run calculations from input files without writing Python code.

### Create Input File

**Save as `input.json`:**

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
    "bsse_type": ["cp", "nocp"]
  },
  "output": {
    "format": "json",
    "file": "results.json"
  }
}
```

### Run Calculation

```bash
qcmanybody run input.json -o results.json
```

### With Parallel Execution

```bash
qcmanybody run input.json -o results.json --parallel --n-workers 4
```

### From XYZ File

**Input file with XYZ molecule:**

```json
{
  "molecule": {
    "source": "xyz",
    "file": "water_dimer.xyz",
    "fragments": [[0, 1, 2], [3, 4, 5]]
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
  "manybody": {
    "max_nbody": 2,
    "bsse_type": ["cp"]
  }
}
```

### Plan Before Running

**Preview tasks without executing:**

```bash
qcmanybody plan input.json
```

This shows what calculations will be performed.

### Validate Input

**Check input file for errors:**

```bash
qcmanybody validate input.json
```

See [CLI Guide](cli_guide.md) for complete documentation.

---

## Integrate with Psi4/NWChem

### Goal
Use specific QC programs with appropriate settings.

### Psi4 Integration

**Basic Psi4 calculation:**

```python
"specification": {
    "mp2/cc-pvdz": {
        "program": "psi4",
        "model": {
            "method": "mp2",
            "basis": "cc-pvdz",
        },
        "driver": "energy",
        "keywords": {
            "scf_type": "df",
            "mp2_type": "df",
            "freeze_core": True,
        },
    }
}
```

**Psi4-specific keywords:**

```python
"keywords": {
    "reference": "rhf",
    "scf_type": "df",
    "mp2_type": "df",
    "freeze_core": True,
    "e_convergence": 1e-8,
    "d_convergence": 1e-6,
}
```

### NWChem Integration

**Basic NWChem calculation:**

```python
"specification": {
    "ccsd/cc-pvtz": {
        "program": "nwchem",
        "model": {
            "method": "ccsd",
            "basis": "cc-pvtz",
        },
        "driver": "energy",
        "keywords": {
            "ccsd__thresh": 1.e-7,
            "scf__thresh": 1.e-8,
        },
    }
}
```

**NWChem-specific keywords:**

```python
"keywords": {
    "scf__thresh": 1.e-8,
    "ccsd__thresh": 1.e-7,
    "memory": "2000 mb",
}
```

### Multi-Program Setup

**Use different programs at different levels:**

```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={
        "driver": "energy",
        "keywords": {
            "levels": {
                1: "nwchem-ccsd/tz",  # Accurate for monomers
                2: "psi4-mp2/dz",      # Faster for dimers
            },
            "bsse_type": ["cp"],
        },
        "specification": {
            "nwchem-ccsd/tz": {
                "program": "nwchem",
                "model": {"method": "ccsd", "basis": "cc-pvtz"},
                "driver": "energy",
            },
            "psi4-mp2/dz": {
                "program": "psi4",
                "model": {"method": "mp2", "basis": "cc-pvdz"},
                "driver": "energy",
            },
        },
    },
)
```

### Installation

**Psi4:**
```bash
conda install -c conda-forge psi4
```

**NWChem:**
```bash
conda install -c conda-forge nwchem
```

---

## Checkpoint and Resume Calculations

### Goal
Save progress for long calculations and resume if interrupted.

### Basic Checkpointing

**Enable checkpointing:**

```python
from qcmanybody.parallel import CheckpointManager

checkpoint_mgr = CheckpointManager(checkpoint_file="calculation.chk")

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
    checkpoint_manager=checkpoint_mgr,
)
```

### Resume After Interruption

**Run the same command again:**

```python
# If interrupted, just run again
checkpoint_mgr = CheckpointManager(checkpoint_file="calculation.chk")

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    parallel=True,
    n_workers=4,
    checkpoint_manager=checkpoint_mgr,
)
# Completed tasks loaded from checkpoint automatically
```

### Configure Checkpoint Frequency

**Save every N completed tasks:**

```python
from qcmanybody.parallel import ExecutorConfig, MultiprocessingExecutor

config = ExecutorConfig(
    n_workers=8,
    checkpoint_interval=5,  # Save after every 5 tasks
    checkpoint_file="calculation.chk",
)
executor = MultiprocessingExecutor(config)

result = ManyBodyComputer.from_manybodyinput(
    manybodyinput,
    executor=executor,
)
```

### CLI Checkpointing

```bash
qcmanybody run input.json --checkpoint checkpoint.chk

# If interrupted, resume:
qcmanybody run input.json --checkpoint checkpoint.chk
```

### Best Practices
- Use checkpointing for calculations expected to take >1 hour
- Checkpoint files can be large - ensure adequate disk space
- Keep checkpoint files until calculation completes
- Use unique checkpoint filenames for different calculations

---

## Debug Calculation Failures

### Goal
Diagnose and fix common calculation errors.

### Enable Verbose Logging

```python
import logging

logging.basicConfig(level=logging.DEBUG)

result = ManyBodyComputer.from_manybodyinput(manybodyinput)
```

### Common Errors and Solutions

#### Error: "Program not found"

**Problem:** QC program not installed or not in PATH

**Solution:**
```bash
# Check if program is available
which psi4

# Install if missing
conda install -c conda-forge psi4
```

#### Error: "Fragment specification invalid"

**Problem:** Fragments don't cover all atoms or overlap

**Solution:**
```python
# Check fragment specification
natoms = len(molecule.symbols)
all_atoms = set()
for frag in molecule.fragments:
    all_atoms.update(frag)

print(f"Atoms in fragments: {sorted(all_atoms)}")
print(f"Expected: {list(range(natoms))}")

# Fix: ensure all atoms assigned exactly once
fragments = [[0, 1, 2], [3, 4, 5]]  # For 6 atoms
```

#### Error: "QCEngine task failed"

**Problem:** Individual QC calculation failed

**Solution:**
```python
# Check individual task results
for task_id, task_result in result.nbody_components.items():
    if not task_result.get('success', True):
        print(f"Failed task: {task_id}")
        print(f"Error: {task_result.get('error')}")
```

#### Error: "Memory exceeded"

**Problem:** QC program ran out of memory

**Solution:**
```python
# Reduce worker count to give each task more memory
config = ExecutorConfig(n_workers=2)  # Instead of 8

# Or limit QC program memory
"keywords": {
    "memory": "2 GB",  # Psi4
}
```

#### Error: "Timeout"

**Problem:** Task exceeded time limit

**Solution:**
```python
config = ExecutorConfig(
    timeout_per_task=7200,  # Increase to 2 hours
)
```

### Inspect Failed Calculations

```python
# Get detailed error information
if not result.success:
    print("Calculation failed!")
    print(f"Error: {result.error}")

# Check individual components
for task_id, component in result.nbody_components.items():
    if hasattr(component, 'error') and component.error:
        print(f"Task {task_id} failed: {component.error}")
```

---

## Optimize Performance

### Goal
Make calculations run faster and more efficiently.

### Choose Appropriate N-Body Level

**Don't compute higher levels than needed:**

```python
# For interaction energies, max_nbody=2 often sufficient
"keywords": {
    "max_nbody": 2,  # Instead of 3 or 4
}
```

**Benchmark:**
- 2-body: N(N-1)/2 tasks (for N fragments)
- 3-body: N(N-1)(N-2)/6 tasks
- 4-body: N(N-1)(N-2)(N-3)/24 tasks

For 5 fragments:
- 2-body: 10 tasks
- 3-body: 10 tasks
- 4-body: 5 tasks
- Total through 4-body: ~30 tasks

### Use Multi-Level Calculations

**Cheap methods for high n-body:**

```python
"keywords": {
    "levels": {
        1: "ccsd/tz",      # Expensive, accurate
        2: "mp2/dz",       # Medium cost
        3: "scf/dz",       # Cheap approximation
    },
}
```

### Leverage Supersystem

**For >3 fragments, use supersystem for high-order terms:**

```python
"keywords": {
    "levels": {
        1: "mp2/dz",
        2: "mp2/dz",
        "supersystem": "mp2/dz",  # Captures all >2-body exactly
    },
    "supersystem_ie_only": True,  # Only use for interaction energy
}
```

### Parallel Execution Best Practices

**Optimal worker count:**

```python
import os

# Use all physical cores
n_workers = os.cpu_count() // 2  # Assumes hyperthreading

config = ExecutorConfig(n_workers=n_workers)
```

**Account for QC program parallelism:**

```python
# If Psi4 uses 4 threads per task
n_psi4_threads = 4
n_workers = os.cpu_count() // n_psi4_threads

"keywords": {
    "num_threads": n_psi4_threads,
}
```

### Use Smaller Basis Sets for Testing

**Test with minimal basis:**

```python
# Testing phase
"model": {"method": "scf", "basis": "sto-3g"}

# Production (after testing)
"model": {"method": "ccsd(t)", "basis": "aug-cc-pvqz"}
```

### Enable Density Fitting

**Much faster for MP2 and beyond:**

```python
"keywords": {
    "scf_type": "df",      # Density-fitted SCF
    "mp2_type": "df",      # Density-fitted MP2
}
```

### Reduce BSSE Corrections

**If BSSE is small, skip corrections:**

```python
# Only nocp (fastest)
"keywords": {
    "bsse_type": ["nocp"],
}

# Instead of all three (3x slower)
"keywords": {
    "bsse_type": ["nocp", "cp", "vmfc"],
}
```

### Monitor Resource Usage

```python
import time

start = time.time()
result = ManyBodyComputer.from_manybodyinput(manybodyinput)
elapsed = time.time() - start

print(f"Calculation took {elapsed:.1f} seconds")
print(f"Tasks completed: {len(result.nbody_components)}")
print(f"Time per task: {elapsed / len(result.nbody_components):.1f} s")
```

---

## Additional Resources

- [Getting Started](getting_started.md) - First steps with QCManyBody
- [High-Level Interface](high-level-interface.md) - Complete API documentation
- [CLI Guide](cli_guide.md) - Command-line interface
- [Parallel Execution Guide](parallel_execution_guide.md) - Detailed parallel execution
- [Keywords](keywords.md) - All available options
- [Results](results.md) - Understanding output
