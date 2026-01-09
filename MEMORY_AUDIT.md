# QCManyBody Memory Audit - HMBE Focus

**Date**: 2026-01-08
**Focus**: Identifying memory-heavy operations, especially for HMBE with large systems

---

## Summary

**Overall Assessment**: ✅ Good - Direct enumeration mode addresses the main memory concern

**Key Findings**:
1. ✅ Direct enumeration (`enumeration_mode="direct"`) is already implemented
2. ⚠️ Filter mode can be memory-intensive for large systems (>50 fragments)
3. ✅ No unbounded caching detected
4. ✅ Minimal deepcopy usage
5. ⚠️ Some opportunities for optimization in compute_map storage

---

## Memory-Heavy Operations Identified

### 1. **builder.py - MBE Term Generation** ⚠️ HIGH IMPACT

**Location**: `qcmanybody/builder.py` lines 94-145

**Issue**: Generates **ALL** MBE terms up to max_nbody using `itertools.combinations`

```python
for x in itertools.combinations(fragment_range, sublevel):
    # Creates C(N, k) combinations
```

**Memory Impact**:
- 16 fragments, MBE-4: ~2,000 terms ✅ Fine
- 64 fragments, MBE-4: ~635,000 terms ⚠️ ~50 MB
- 100 fragments, MBE-4: ~4.6M terms ❌ ~400 MB+

**Solution**: ✅ **Already Implemented!**
- Use `enumeration_mode="direct"` in HMBESpecification
- Bypasses filter mode entirely for large systems
- Direct enumeration only generates HMBE terms (100-1000x fewer)

**Recommendation**:
```python
# For systems with >30 fragments, use:
hmbe_spec = HMBESpecification(
    truncation_orders=(2, 4),
    hierarchy=hierarchy,
    enumeration_mode="direct"  # or "auto" (chooses based on size)
)
```

### 2. **Sub-cluster Generation** ⚠️ MEDIUM IMPACT

**Location**:
- `qcmanybody/hmbe_filter.py:148-160` (generate_all_subclusters)
- `qcmanybody/core.py:712-783` (_validate_hmbe_completeness)

**Issue**: Generates all sub-clusters for validation/Schengen

```python
for sub_size in range(1, n):
    for sub_cluster in combinations(frag_tuple, sub_size):
        # For n=5: generates C(5,1) + C(5,2) + C(5,3) + C(5,4) = 30 sub-clusters
```

**Memory Impact**:
- Typically small (n ≤ 4 for most HMBE terms)
- Worst case (3,5)-HMBE: n=5 → 31 sub-clusters per term
- Total memory: Usually <10 MB even for large systems

**Status**: ✅ Acceptable - Only computed once, not stored long-term

### 3. **Compute Map Storage** ℹ️ LOW IMPACT

**Location**: `qcmanybody/core.py` - `compute_map` property

**Issue**: Stores (fragment_tuple, basis_tuple) pairs in nested dicts

**Structure**:
```python
{
  "modelchem": {
    "bsse_type": {
      nbody_level: Set[(frag_tuple, bas_tuple)]
    }
  }
}
```

**Memory Estimation**:
- Each tuple pair: ~50-100 bytes
- 1,000 terms: ~100 KB
- 10,000 terms: ~1 MB
- 100,000 terms: ~10 MB

**Status**: ✅ Acceptable - Reasonable overhead

**Potential Optimization** (low priority):
- Use frozenset instead of tuple for fragment/basis (more memory-efficient)
- Implement lazy iteration instead of set storage

### 4. **deepcopy Usage** ✅ MINIMAL

**Found**: Only 1 occurrence in `computer.py:91` for keywords dict

```python
return copy.deepcopy(keywords)
```

**Status**: ✅ Minimal impact - keywords dict is small

---

## Operations That Are Memory-Efficient ✅

### 1. **Direct HMBE Enumeration**
- Generates only HMBE terms (not all MBE)
- Memory scales with HMBE terms, not MBE terms
- 100x - 1000x reduction for large systems

### 2. **Iterator-Based Processing**
- `iterate_molecules()` and `iterate_molecules_by_level()` are generators
- Process fragments one at a time, don't store all in memory

### 3. **No Unbounded Caching**
- No `@lru_cache` decorators found
- No global caches that grow with system size

### 4. **Minimal Array Allocations**
- No large pre-allocated numpy arrays
- Arrays created only when needed for specific calculations

---

## Recommendations

### For Current Users (Water-16, <30 fragments)
✅ **No changes needed** - Current implementation is fine

### For Medium Systems (30-100 fragments)
⚠️ **Use direct enumeration or auto mode**:
```python
hmbe_spec = HMBESpecification(
    truncation_orders=(2, 4),
    hierarchy=hierarchy,
    enumeration_mode="auto"  # Automatically chooses "direct" for ≥30 fragments
)
```

### For Large Systems (100+ fragments)
❌ **Must use direct enumeration**:
```python
hmbe_spec = HMBESpecification(
    truncation_orders=(2, 4),
    hierarchy=hierarchy,
    enumeration_mode="direct"  # Required for large systems
)
```

---

## Memory Usage Estimates

### Water-16 System (16 fragments, 4×4 hierarchy)

**Full MBE-4**:
- Terms: 1,820
- Memory: ~200 KB ✅

**(2,4)-HMBE with filter mode**:
- Full MBE generated: 1,820 terms (~200 KB)
- Filtered to: ~500 terms (~50 KB)
- Total memory: ~250 KB ✅

**(2,4)-HMBE with direct mode**:
- Terms generated: ~500 terms
- Memory: ~50 KB ✅
- **Saves**: ~200 KB (not significant for this size)

### Water-64 System (64 fragments, 4×4×4 hierarchy)

**Full MBE-4**:
- Terms: 635,376
- Memory: ~60 MB ⚠️

**(2,3,4)-HMBE with filter mode**:
- Full MBE generated: 635,376 terms (~60 MB)
- Filtered to: ~5,000 terms (~500 KB)
- Total memory: ~60 MB ⚠️ (wastes ~59.5 MB!)

**(2,3,4)-HMBE with direct mode**:
- Terms generated: ~5,000 terms
- Memory: ~500 KB ✅
- **Saves**: ~59.5 MB (120x reduction!)

### Protein System (256 fragments, 4×4×4×4 hierarchy)

**Full MBE-4**:
- Terms: ~178 million
- Memory: ~17 GB ❌ **WILL CRASH**

**(2,2,3,4)-HMBE with filter mode**:
- Full MBE generated: ~178 million terms ❌ **WILL CRASH**
- Would filter to: ~20,000 terms
- **Cannot use filter mode**

**(2,2,3,4)-HMBE with direct mode**:
- Terms generated: ~20,000 terms
- Memory: ~2 MB ✅
- **Only viable approach**

---

## Implementation Notes

### Current "Auto" Mode Threshold
```python
# From core.py line 182
if enumeration_mode == "auto":
    enumeration_mode = "direct" if self.nfragments >= 30 else "filter"
```

**Recommendation**: ✅ Good threshold
- Filter mode fine for <30 fragments
- Direct mode essential for ≥30 fragments

### Parallel Execution Impact

**Memory per worker**:
- Each worker needs its own copy of compute_map
- For 4 workers with 60 MB compute_map: 240 MB total ⚠️
- With direct enumeration (500 KB): 2 MB total ✅

**Recommendation**: Always use direct mode for parallel execution with >30 fragments

---

## Future Optimizations (Low Priority)

### 1. Lazy Compute Map Construction
Instead of storing all (frag, bas) pairs upfront:
- Generate on-demand during iteration
- Save memory at cost of slight CPU overhead

### 2. Compressed Storage for Large Systems
- Use integer indices instead of tuples
- Compress with delta encoding
- Potential 2-4x memory reduction

### 3. Disk-Based Caching for Very Large Systems
- Store compute_map to disk
- Load fragments on-demand
- For >1000 fragment systems

---

## Conclusion

✅ **QCManyBody HMBE implementation is memory-efficient for typical use cases**

Key points:
1. Direct enumeration mode solves the main memory bottleneck
2. For water-16 and similar systems (<30 fragments): no issues
3. For larger systems: use `enumeration_mode="direct"` or `"auto"`
4. No memory leaks or unbounded caching detected
5. Parallel execution is safe with direct enumeration

**Action Items**:
- ✅ Direct enumeration already implemented
- ✅ Auto mode threshold (30 fragments) is appropriate
- ✅ Documentation should emphasize using direct/auto mode for large systems
- 📝 Add memory usage estimates to user guide
