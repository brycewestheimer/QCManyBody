# Results

QCManyBody provides comprehensive analysis results from many-body expansion calculations, including energies, gradients, or hessians at various n-body levels with different BSSE corrections.

## Overview

Results are returned differently depending on the interface:

- **[High-Level Interface](high-level-interface.md)**: Returns a `ManyBodyResult` QCSchema object
- **[Core Interface](core-interface.md)**: Results stored in nested dictionaries

This page focuses on the high-level interface `ManyBodyResult` object, which provides structured access to all computed properties.

## ManyBodyResult Structure

The `ManyBodyResult` object contains:

```python
result = ManyBodyComputer.from_manybodyinput(manybodyinput)

result.return_result     # Main result (energy, gradient, or hessian)
result.properties        # Detailed many-body properties
result.nbody_components  # Individual task results
result.provenance        # Calculation metadata
result.success           # Boolean success flag
result.input_data        # Copy of input specification
```

## Main Result: return_result

The `return_result` field contains the final computed property based on the driver:

### Energy Driver

```python
result.return_result  # Float: total energy in Hartree
```

Example:
```python
manybodyinput = ManyBodyInput(
    molecule=mol,
    specification={"driver": "energy", ...},
)
result = ManyBodyComputer.from_manybodyinput(manybodyinput)
print(f"Total energy: {result.return_result} Eh")
```

### Gradient Driver

```python
result.return_result  # Array shape (Natoms, 3): gradient in Eh/Bohr
```

Example:
```python
result = ManyBodyComputer.from_manybodyinput(
    ManyBodyInput(molecule=mol, specification={"driver": "gradient", ...})
)
gradient = result.return_result  # Shape: (Natoms, 3)
forces = -gradient  # Forces are negative gradient
```

### Hessian Driver

```python
result.return_result  # Array shape (3*Natoms, 3*Natoms): Hessian in Eh/Bohr²
```

## Properties: Detailed Many-Body Analysis

The `properties` field contains all computed many-body properties as a `ManyBodyProperties` object. Property names follow a systematic naming convention:

```
{correction}_{property_type}_through_{nbody}_body
{correction}_{nbody}_body_contribution_to_{property_type}
```

### Property Name Components

**Correction types:**
- `nocp_corrected` - No counterpoise correction
- `cp_corrected` - Counterpoise corrected
- `vmfc_corrected` - Valiron-Mayer function counterpoise

**Property types:**
- `total_energy` - Total energy of the system
- `interaction_energy` - Interaction energy (deviation from monomers)

**N-body levels:**
- `1_body`, `2_body`, `3_body`, etc.
- `through_N_body` - Sum through N-body level

### Common Properties

#### Total Energies

Cumulative energies through each n-body level:

```python
# CP-corrected total energies
result.properties.cp_corrected_total_energy_through_1_body
result.properties.cp_corrected_total_energy_through_2_body
result.properties.cp_corrected_total_energy_through_3_body

# No CP total energies
result.properties.nocp_corrected_total_energy_through_2_body

# VMFC total energies
result.properties.vmfc_corrected_total_energy_through_2_body
```

#### Interaction Energies

Cumulative interaction energies (relative to isolated monomers):

```python
# CP-corrected interaction energies
result.properties.cp_corrected_interaction_energy_through_1_body  # Always 0.0
result.properties.cp_corrected_interaction_energy_through_2_body
result.properties.cp_corrected_interaction_energy_through_3_body

# No CP interaction energies
result.properties.nocp_corrected_interaction_energy_through_2_body

# VMFC interaction energies
result.properties.vmfc_corrected_interaction_energy_through_3_body
```

#### Individual N-Body Contributions

Individual contribution at each n-body level:

```python
# 2-body contribution to interaction energy
result.properties.cp_corrected_2_body_contribution_to_interaction_energy
result.properties.nocp_corrected_2_body_contribution_to_interaction_energy

# 3-body contribution
result.properties.cp_corrected_3_body_contribution_to_interaction_energy
```

### Accessing All Properties

List all available properties:

```python
# Get all property names
property_names = result.properties.dict().keys()

# Print all properties
for name, value in result.properties.dict().items():
    if value is not None:
        print(f"{name}: {value}")
```

Filter by correction type:

```python
# Only CP-corrected properties
cp_properties = {
    name: value
    for name, value in result.properties.dict().items()
    if name.startswith("cp_corrected")
}
```

## Component Results: nbody_components

The `nbody_components` dictionary contains results from individual quantum chemistry calculations:

```python
result.nbody_components  # Dict[str, Dict]
```

Each key is a task identifier, and the value contains:

- Task specification (molecule, method, basis, etc.)
- Individual calculation result (energy, gradient, or hessian)
- Provenance information

Example:

```python
for task_id, task_result in result.nbody_components.items():
    print(f"Task: {task_id}")
    print(f"  Energy: {task_result.get('energy')}")
```

## Provenance Information

Metadata about the calculation:

```python
result.provenance = {
    "creator": "QCManyBody",
    "version": "0.2.0",
    "routine": "qcmanybody.computer",
    ...
}
```

## Complete Example

```python
from qcelemental.models import Molecule
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Water trimer calculation
water_trimer = Molecule(
    symbols=["O", "H", "H"] * 3,
    geometry=[...],  # 9 atoms total
    fragments=[[0, 1, 2], [3, 4, 5], [6, 7, 8]],
)

manybodyinput = ManyBodyInput(
    molecule=water_trimer,
    specification={
        "driver": "energy",
        "keywords": {
            "max_nbody": 3,
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

# Main result
print(f"Final energy: {result.return_result} Eh")

# Compare BSSE corrections
print("\n2-Body Interaction Energies:")
print(f"  No CP: {result.properties.nocp_corrected_interaction_energy_through_2_body} Eh")
print(f"  CP:    {result.properties.cp_corrected_interaction_energy_through_2_body} Eh")
print(f"  VMFC:  {result.properties.vmfc_corrected_interaction_energy_through_2_body} Eh")

# N-body contributions
print("\nCP-Corrected Contributions:")
print(f"  1-body: {result.properties.cp_corrected_1_body_contribution_to_interaction_energy} Eh")
print(f"  2-body: {result.properties.cp_corrected_2_body_contribution_to_interaction_energy} Eh")
print(f"  3-body: {result.properties.cp_corrected_3_body_contribution_to_interaction_energy} Eh")

# Total energies
print("\nCP-Corrected Total Energies:")
print(f"  Through 1-body: {result.properties.cp_corrected_total_energy_through_1_body} Eh")
print(f"  Through 2-body: {result.properties.cp_corrected_total_energy_through_2_body} Eh")
print(f"  Through 3-body: {result.properties.cp_corrected_total_energy_through_3_body} Eh")
```

## Result Units

All results are returned in atomic units:

| Property | Unit |
|----------|------|
| Energy | Hartree (Eh) |
| Gradient | Hartree/Bohr (Eh/a₀) |
| Hessian | Hartree/Bohr² (Eh/a₀²) |
| Geometry | Bohr (a₀) in results |

Note: Input geometries are typically in Angstroms but are converted internally.

### Unit Conversion

Convert to common units:

```python
# Hartree to kcal/mol
HARTREE_TO_KCAL = 627.509474
ie_kcal = result.properties.cp_corrected_interaction_energy_through_2_body * HARTREE_TO_KCAL

# Hartree to kJ/mol
HARTREE_TO_KJ = 2625.499639
ie_kj = result.properties.cp_corrected_interaction_energy_through_2_body * HARTREE_TO_KJ

# Gradients to forces (negative gradient)
forces = -result.return_result  # If driver="gradient"
```

## Exporting Results

### JSON Export

```python
import json

# Export to JSON
result_dict = result.dict()
with open("results.json", "w") as f:
    json.dump(result_dict, f, indent=2, default=str)
```

### Properties Only

```python
# Export only properties
properties_dict = result.properties.dict()
with open("properties.json", "w") as f:
    json.dump(properties_dict, f, indent=2)
```

### Pandas DataFrame

Convert properties to DataFrame for analysis:

```python
import pandas as pd

# Create DataFrame of properties
properties_dict = result.properties.dict()
df = pd.DataFrame([
    {"property": k, "value": v}
    for k, v in properties_dict.items()
    if v is not None
])

# Save to CSV
df.to_csv("results.csv", index=False)
```

## Special Cases

### 1-Body Results

For 1-body calculations:
- Interaction energies are always 0.0 (no interaction with only monomers)
- Total energy equals sum of monomer energies
- All correction schemes (nocp, cp, vmfc) give identical results

```python
# All equal for 1-body
result.properties.cp_corrected_total_energy_through_1_body ==     result.properties.nocp_corrected_total_energy_through_1_body ==     result.properties.vmfc_corrected_total_energy_through_1_body
```

### Supersystem Calculations

When using `supersystem` in levels, results include exact higher-order terms:

```python
"keywords": {
    "levels": {1: "mp2/dz", 2: "mp2/dz", "supersystem": "ccsd/tz"},
}

# Result includes exact (3-body + 4-body + ... + N-body) from supersystem
```

### Multi-Level Results

For multi-level calculations, properties reflect the mixed-level approach:

```python
"keywords": {
    "levels": {1: "ccsd/tz", 2: "mp2/dz"},
}

# 1-body terms computed at CCSD/cc-pVTZ
# 2-body terms computed at MP2/cc-pVDZ
# Total energy combines both levels
```

## Understanding Property Names

### Cumulative vs. Contribution

**Cumulative ("through N-body"):**
```python
# Sum of 1-body + 2-body contributions
result.properties.cp_corrected_interaction_energy_through_2_body
```

**Contribution ("N-body contribution"):**
```python
# Only the 2-body contribution
result.properties.cp_corrected_2_body_contribution_to_interaction_energy
```

### Total vs. Interaction

**Total Energy:**
- Absolute energy of the system
- Never zero

**Interaction Energy:**
- Energy relative to isolated monomers
- Negative = attractive interaction
- Zero for 1-body (no interaction)

```python
# Relationship
total_energy = monomer_energies + interaction_energy
```

## Troubleshooting

### Property Not Found

If a property is `None` or missing:

```python
# Check if property was requested
if result.properties.vmfc_corrected_total_energy_through_3_body is None:
    print("VMFC not computed - check bsse_type includes 'vmfc'")
```

Ensure the BSSE correction was requested:

```python
"keywords": {
    "bsse_type": ["nocp", "cp", "vmfc"],  # Request all three
}
```

### Missing N-Body Levels

Properties only exist up to `max_nbody`:

```python
# If max_nbody=2, these don't exist
result.properties.cp_corrected_3_body_contribution_to_interaction_energy  # None
```

### Gradient Shape Issues

Gradient shape is `(Natoms, 3)`:

```python
gradient = result.return_result
natoms = len(molecule.symbols)
assert gradient.shape == (natoms, 3)

# Access gradient for atom i
grad_i = gradient[i, :]  # [dx, dy, dz]
```

## API Reference

### ManyBodyResult

::: qcmanybody.models.ManyBodyResult
    options:
      show_source: false

### ManyBodyProperties

::: qcmanybody.models.ManyBodyResultProperties
    options:
      show_source: false

## Next Steps

- [High-Level Interface](high-level-interface.md) - How to generate results
- [Keywords](keywords.md) - Control what properties are computed
- [How-To Guides](how-to-guides.md) - Common analysis patterns
- [QCSchema](qcschema.md) - Understanding the data structures
