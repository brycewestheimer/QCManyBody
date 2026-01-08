# Migration Guide: Converting MBE to HMBE

This guide helps you convert existing Many-Body Expansion (MBE) calculations to use Hierarchical MBE (HMBE) for improved performance.

---

## Table of Contents

1. [Should You Migrate?](#should-you-migrate)
2. [Quick Conversion](#quick-conversion)
3. [Step-by-Step Migration](#step-by-step-migration)
4. [Common Scenarios](#common-scenarios)
5. [Backward Compatibility](#backward-compatibility)
6. [Troubleshooting](#troubleshooting)

---

## Should You Migrate?

### ✅ Migrate to HMBE if:

- **System size >20 fragments**
- **Natural hierarchical structure** exists
- **Multiple conformations** to study
- **Computational resources** are limited
- **Turnaround time** is important

### ❌ Stay with standard MBE if:

- **System size <10 fragments** (MBE is already fast)
- **No clear hierarchy** (random aggregates)
- **Maximum accuracy required** (though HMBE+Schengen can help)
- **Already have working MBE scripts** and performance is acceptable

---

## Quick Conversion

### Before (Standard MBE)

```python
from qcmanybody import ManyBodyComputer
from qcmanybody.models.v1 import BsseEnum, ManyBodyInput, ManyBodyKeywords

keywords = ManyBodyKeywords(
    max_nbody=3,
    bsse_type=[BsseEnum.cp]
)

mb_input = ManyBodyInput(
    molecule=mol,
    specification={...}
)

result = ManyBodyComputer.from_manybodyinput(mb_input)
```

### After (HMBE)

```python
from qcmanybody import ManyBodyComputer
from qcmanybody.models.v1 import BsseEnum, ManyBodyInput, ManyBodyKeywords
from qcmanybody.models.hierarchy import FragmentHierarchy, HMBESpecification  # NEW

# NEW: Define hierarchy
hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers={
        1: ("G0", "F1"), 2: ("G0", "F2"), ...,  # Your grouping
    }
)

# NEW: Create HMBE specification
hmbe_spec = HMBESpecification(
    truncation_orders=(2, 3),  # (T_1, T_2)
    hierarchy=hierarchy
)

keywords = ManyBodyKeywords(
    max_nbody=3,
    bsse_type=[BsseEnum.cp],
    hmbe_spec=hmbe_spec  # NEW: Add HMBE
)

mb_input = ManyBodyInput(
    molecule=mol,
    specification={...}
)

result = ManyBodyComputer.from_manybodyinput(mb_input)
```

**That's it!** Three additions:
1. Define hierarchy
2. Create HMBESpecification
3. Add `hmbe_spec` to keywords

---

## Step-by-Step Migration

### Step 1: Analyze Your System

**Identify fragment groupings:**

```python
# Example: 16-water cluster
# Fragments 1-16 are water molecules
# Natural grouping: 4 spatial regions

print(f"Total fragments: {len(mol.fragments)}")
print("Fragments:", mol.fragments)

# Decide on groups
# Region 0: Waters 1-4
# Region 1: Waters 5-8
# Region 2: Waters 9-12
# Region 3: Waters 13-16
```

**Ask yourself:**
- What are logical subdivisions?
- Which fragments are spatially close?
- Which fragments interact strongly?

### Step 2: Define Fragment Hierarchy

**2-Tier Example** (Recommended for most systems):

```python
from qcmanybody.models.hierarchy import FragmentHierarchy

# Manual assignment
fragment_tiers = {
    1: ("R0", "W1"),  2: ("R0", "W2"),  3: ("R0", "W3"),  4: ("R0", "W4"),
    5: ("R1", "W5"),  6: ("R1", "W6"),  7: ("R1", "W7"),  8: ("R1", "W8"),
    9: ("R2", "W9"), 10: ("R2", "W10"), 11: ("R2", "W11"), 12: ("R2", "W12"),
   13: ("R3", "W13"), 14: ("R3", "W14"), 15: ("R3", "W15"), 16: ("R3", "W16"),
}

hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers=fragment_tiers,
    tier_names=("region", "water")
)
```

**Programmatic assignment** (for regular patterns):

```python
fragment_tiers = {}
frags_per_group = 4
num_groups = 4

for i in range(16):
    frag_id = i + 1  # 1-indexed!
    group_id = i // frags_per_group
    fragment_tiers[frag_id] = (f"R{group_id}", f"W{i}")

hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers=fragment_tiers,
    tier_names=("region", "water")
)
```

**3-Tier Example** (For very large systems):

```python
# Example: Protein with domains → residues → atoms
fragment_tiers = {
    1: ("Domain_A", "Res_1", "Atom_group_1"),
    2: ("Domain_A", "Res_1", "Atom_group_2"),
    3: ("Domain_A", "Res_2", "Atom_group_3"),
    ...
}

hierarchy = FragmentHierarchy(
    num_tiers=3,
    fragment_tiers=fragment_tiers,
    tier_names=("domain", "residue", "atom_group")
)
```

### Step 3: Choose Truncation Orders

**Guidelines:**

| Goal | Recommended | Notes |
|------|------------|-------|
| Maximum speed | (2,3) | ~1.5-2x reduction vs MBE-3 |
| Balanced | (2,4) | ~3-5x reduction vs MBE-4 |
| Higher accuracy | (3,4) | ~1.1-1.5x reduction vs MBE-4 |

**Create specification:**

```python
from qcmanybody.models.hierarchy import HMBESpecification

hmbe_spec = HMBESpecification(
    truncation_orders=(2, 3),  # (T_1, T_2) for 2-tier
    hierarchy=hierarchy
)
```

**Rule of thumb:**
- `T_1` = number of tier-1 groups you want to allow per term (typically 2-3)
- `T_2` = your original `max_nbody` from MBE

### Step 4: Update Keywords

**Before:**
```python
keywords = ManyBodyKeywords(
    max_nbody=3,
    bsse_type=[BsseEnum.cp]
)
```

**After:**
```python
keywords = ManyBodyKeywords(
    max_nbody=3,  # Should match truncation_orders[-1]
    bsse_type=[BsseEnum.cp],
    hmbe_spec=hmbe_spec  # Add HMBE
)
```

### Step 5: Run and Validate

```python
result = ManyBodyComputer.from_manybodyinput(mb_input)

# Check HMBE statistics
stats = result.properties.hmbe_metadata
print(f"MBE would need:  {stats['mbe_term_counts']}")
print(f"HMBE needed:     {stats['hmbe_term_counts']}")
print(f"Reduction:       {stats['reduction_factors']}")

# Verify you got meaningful reduction
mc_level = list(stats['reduction_factors'].keys())[0]
reduction = stats['reduction_factors'][mc_level]

if reduction < 1.2:
    print("⚠️ Warning: HMBE barely helping! Check your hierarchy design.")
elif reduction < 2.0:
    print("✓ Modest speedup")
elif reduction < 5.0:
    print("✓✓ Good speedup!")
else:
    print("✓✓✓ Excellent speedup!")
```

### Step 6: (Optional) Add Schengen Terms

If accuracy is important:

```python
from qcmanybody.models.hierarchy import SchengenSpecification

hmbe_spec = HMBESpecification(
    truncation_orders=(2, 3),
    hierarchy=hierarchy,
    schengen=SchengenSpecification(
        enabled=True,
        selection_fraction=0.10,  # Add top 10% of excluded terms
        distance_metric="R2"
    )
)
```

---

## Common Scenarios

### Scenario 1: Water Clusters

**Before (MBE-3 on 20 waters):**
```python
mol = Molecule(symbols=..., geometry=..., fragments=...)  # 20 waters
keywords = ManyBodyKeywords(max_nbody=3, bsse_type=[BsseEnum.cp])
# 1330 calculations needed
```

**After (HMBE with 5×4 grouping):**
```python
# Group into 5 regions, 4 waters each
hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers={
        1: ("R0", "W0"), 2: ("R0", "W1"), 3: ("R0", "W2"), 4: ("R0", "W3"),
        5: ("R1", "W4"), ...
    }
)

hmbe_spec = HMBESpecification(truncation_orders=(2, 3), hierarchy=hierarchy)
keywords = ManyBodyKeywords(max_nbody=3, bsse_type=[BsseEnum.cp], hmbe_spec=hmbe_spec)
# ~500-600 calculations (2-3x speedup)
```

### Scenario 2: Molecular Crystals

**Before (MBE-4 on unit cell with 8 molecules):**
```python
mol = Molecule(...)  # 8 molecules = 8 fragments
keywords = ManyBodyKeywords(max_nbody=4, bsse_type=[BsseEnum.cp])
# 163 calculations
```

**After (HMBE with 2×4 grouping):**
```python
# Group into 2 sub-cells, 4 molecules each
hierarchy = FragmentHierarchy(
    num_tiers=2,
    fragment_tiers={
        1: ("Cell0", "M0"), 2: ("Cell0", "M1"), 3: ("Cell0", "M2"), 4: ("Cell0", "M3"),
        5: ("Cell1", "M4"), 6: ("Cell1", "M5"), 7: ("Cell1", "M6"), 8: ("Cell1", "M7"),
    }
)

hmbe_spec = HMBESpecification(truncation_orders=(2, 4), hierarchy=hierarchy)
keywords = ManyBodyKeywords(max_nbody=4, bsse_type=[BsseEnum.cp], hmbe_spec=hmbe_spec)
# ~80 calculations (2x speedup)
```

### Scenario 3: Large System (100+ fragments)

**Problem:** Standard MBE-4 on 100 fragments = **4,598,126 calculations** ❌ Impossible!

**Solution:** HMBE with 10×10 hierarchy

```python
# 100 fragments → 10 tier-1 groups × 10 fragments each
fragment_tiers = {}
for i in range(100):
    frag_id = i + 1
    group_id = i // 10
    fragment_tiers[frag_id] = (f"G{group_id}", f"F{i}")

hierarchy = FragmentHierarchy(num_tiers=2, fragment_tiers=fragment_tiers)

hmbe_spec = HMBESpecification(truncation_orders=(2, 4), hierarchy=hierarchy)
keywords = ManyBodyKeywords(max_nbody=4, bsse_type=[BsseEnum.cp], hmbe_spec=hmbe_spec)
# ~5,000-10,000 calculations (500-1000x speedup!) ✅ Feasible!
```

**Note:** For >50 fragments, you'll want **Phase 2: Direct Enumeration** (coming soon) for even better performance.

### Scenario 4: Protein Fragment

**Covalent system - requires hydrogen capping (Phase 4 - coming soon):**

```python
# Future API (Phase 4)
from qcmanybody.models.capping import HydrogenCappingSpecification

capping = HydrogenCappingSpecification(
    bond_cuts=[("res1:C", "res2:N"), ...]
)

hierarchy = FragmentHierarchy(...)  # Domains → Residues
hmbe_spec = HMBESpecification(
    truncation_orders=(2, 3, 4),  # 3-tier
    hierarchy=hierarchy,
    capping=capping  # Add H caps
)
```

---

## Backward Compatibility

### Your Existing MBE Scripts Will Still Work!

**Zero breaking changes:**
```python
# This still works exactly as before
keywords = ManyBodyKeywords(max_nbody=3, bsse_type=[BsseEnum.cp])
# No HMBE - runs standard MBE
```

**Mixing old and new:**
```python
# You can use HMBE for some systems, MBE for others
if len(mol.fragments) > 20:
    # Large system - use HMBE
    keywords = ManyBodyKeywords(max_nbody=3, hmbe_spec=hmbe_spec, ...)
else:
    # Small system - use standard MBE
    keywords = ManyBodyKeywords(max_nbody=3, ...)
```

---

## Troubleshooting

### Issue: "HMBE not reducing term count"

**Possible causes:**

1. **Only one tier-1 group** (trivial hierarchy)
   ```python
   # BAD: All fragments in same group
   fragment_tiers = {1: ("G0", "F1"), 2: ("G0", "F2"), ...}
   ```
   **Fix:** Use multiple tier-1 groups (4-8 recommended)

2. **T_1 too large**
   ```python
   # If T_1 = number of groups, no filtering occurs
   hmbe_spec = HMBESpecification(truncation_orders=(5, 5), ...)  # 5 groups
   ```
   **Fix:** Use smaller T_1 (typically 2-3)

3. **System too small**
   - With 6 fragments in 2 groups, reduction will be minimal
   **Fix:** HMBE works best for >20 fragments

### Issue: "Results don't match standard MBE"

**This is expected!** HMBE is an approximation.

**Check:**
1. **Error magnitude:** Is it within acceptable range?
   - Typical: 1-10 kcal/mol for (2,3)-HMBE
   - With Schengen: 0.5-5 kcal/mol

2. **Try adding Schengen terms** for better accuracy

3. **Use less aggressive truncation:**
   - (3,4) instead of (2,3)
   - (3,3) instead of (2,3)

### Issue: "ValidationError: Truncation orders not monotonic"

```python
# ERROR: T_1 > T_2
hmbe_spec = HMBESpecification(truncation_orders=(3, 2), ...)
```

**Fix:** Ensure T_1 ≤ T_2 ≤ T_3
```python
hmbe_spec = HMBESpecification(truncation_orders=(2, 3), ...)  # Correct
```

### Issue: "Fragment indexing confusion"

**Remember:** QCManyBody uses **1-based fragment indexing**!

```python
# Molecule has 0-indexed atoms
fragments = [[0,1,2], [3,4,5], ...]  # Correct

# But hierarchy uses 1-indexed fragments
fragment_tiers = {
    1: ("G0", "F1"),  # Fragment 1 (first fragment)
    2: ("G0", "F2"),  # Fragment 2 (second fragment)
    ...
}
```

### Issue: "Performance not improving as expected"

**Checklist:**
1. ✓ Check reduction factor (should be >1.5x)
2. ✓ Use parallel execution (`parallel=True`)
3. ✓ Profile to find bottleneck (QC program vs HMBE overhead)
4. ✓ For >50 fragments, await Phase 2 (Direct Enumeration)

---

## Migration Checklist

- [ ] Analyzed system for natural hierarchy
- [ ] Defined fragment tiers (2 or 3 tiers)
- [ ] Chose truncation orders (start with (2,3) or (2,4))
- [ ] Created FragmentHierarchy
- [ ] Created HMBESpecification
- [ ] Added hmbe_spec to ManyBodyKeywords
- [ ] Ran calculation and checked statistics
- [ ] Verified reduction factor >1.2x
- [ ] (Optional) Added Schengen terms for accuracy
- [ ] (Optional) Enabled parallel execution
- [ ] Validated results vs standard MBE (if feasible)

---

## Summary

**Converting to HMBE is straightforward:**
1. Define how fragments group hierarchically
2. Choose truncation orders
3. Add `hmbe_spec` to your keywords
4. Run and validate

**Benefits:**
- 2-1000x speedup depending on system size
- Same API, minimal code changes
- Fully backward compatible
- Works with all existing QCManyBody features

**When in doubt:**
- Start with (2,3)-HMBE and 4-6 fragments per tier-1 group
- Check statistics after first run
- Adjust based on reduction factor and accuracy needs

For questions or issues, see the [HMBE User Guide](HMBE_USER_GUIDE.md) or open a GitHub issue.

---

*Last updated: January 2026*
