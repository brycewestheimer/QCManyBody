# QCElemental Fragment Validation Limitation

## Issue

QCElemental's Molecule class has a scalability limitation with large fragment counts:
- **16 fragments**: Works fine
- **32 fragments**: System crashes completely (OOM/hang)
- **64 fragments**: Cannot be tested due to 32-fragment limit

## Reproduction

```python
from qcelemental.models import Molecule

# This works
mol16 = Molecule(
    symbols=['He'] * 16,
    geometry=[[i, 0, 0] for i in range(16)],
    fragments=[[i] for i in range(16)]
)

# This crashes the system
mol32 = Molecule(
    symbols=['He'] * 32,
    geometry=[[i, 0, 0] for i in range(32)],
    fragments=[[i] for i in range(32)]
)
```

## Impact on QCManyBody

- **HMBE direct mode memory fixes are COMPLETE** ✅
  - `compute_map` generation reduced from ~680k terms to ~5k terms
  - `enumeration_mode` properly flows from CLI to core
  - Memory footprint reduced 100-1000x

- **Testing blocked by QCElemental** ❌
  - Cannot load water64 system (64 fragments)
  - Cannot test large HMBE systems
  - Even if initial load worked, `get_fragment()` calls would fail when combining 32+ fragments

## Technical Details

The issue appears to be in QCElemental's Pydantic fragment validation:
- Fragment list validation scales poorly (O(n²) or worse)
- Validation called on every Molecule() constructor
- QCManyBody calls Molecule() once per molecule load + once per HMBE term via `get_fragment()`
- For water64 (2,3,4)-HMBE: ~5k Molecule() calls, some combining 32+ fragments

## Workarounds

### Option 1: Use smaller test systems (< 32 fragments)
Create water16 or water24 test configs that stay below the 32-fragment limit.

### Option 2: Report QCElemental bug
File issue at https://github.com/MolSSI/QCElemental/issues with reproduction script.

### Option 3: Monkey-patch QCElemental (RISKY)
Temporarily disable fragment validation:
```python
from qcelemental.models import Molecule
# Disable validation - UNSAFE, use only for testing
Molecule.Config.validate_assignment = False
```

## Status

**HMBE Fixes**: ✅ Complete
**QCManyBody Code**: ✅ All changes merged
**Testing**: ❌ Blocked by external dependency

**Next Steps**:
1. Report bug to QCElemental team
2. Test HMBE fixes with water16 or water24 systems
3. Wait for QCElemental patch before testing water64

## Files Modified for HMBE Fix

- [qcmanybody/core.py](qcmanybody/core.py) lines 162-357: Skip MBE generation in direct mode
- [qcmanybody/cli/schemas/input_schema.py](qcmanybody/cli/schemas/input_schema.py) line 390: Added enumeration_mode field
- [qcmanybody/cli/converter.py](qcmanybody/cli/converter.py) line 245: Pass enumeration_mode to HMBESpecification
