# Water-16 HMBE Test Inputs

**Purpose**: Real-world validation of HMBE implementations with parallel execution

---

## Test Configurations

### 1. `water16_hmbe_23.json` - (2,3)-HMBE
- **Truncation**: Max 2 groups, max 3-body
- **Expected reduction**: ~1.5-2x vs full MBE-3
- **Use case**: Moderate reduction, good accuracy

### 2. `water16_hmbe_24.json` - (2,4)-HMBE
- **Truncation**: Max 2 groups, max 4-body
- **Expected reduction**: ~1.8-2.5x vs full MBE-4
- **Use case**: Higher accuracy with manageable cost

### 3. `water16_hmbe_34.json` - (3,4)-HMBE
- **Truncation**: Max 3 groups, max 4-body
- **Expected reduction**: ~1.2-1.5x vs full MBE-4
- **Use case**: Maximum accuracy, minimal reduction

---

## System Details

**Molecule**: 16 water molecules
**Hierarchy**: 4×4 two-tier (4 groups of 4 waters)
- Group 0: Waters 1-4
- Group 1: Waters 5-8
- Group 2: Waters 9-12
- Group 3: Waters 13-16

**Method**: RI-MP2/aug-cc-pVDZ
**Driver**: Energy
**BSSE**: NOCP (no counterpoise correction)

---

## Running Tests

### Using QCManyBody CLI

```bash
# Run (2,3)-HMBE calculation
qcmanybody run test_inputs/water16_hmbe_23.json

# Run with different parallel settings
qcmanybody run test_inputs/water16_hmbe_24.json --n-workers 8

# Run without parallel execution
qcmanybody run test_inputs/water16_hmbe_34.json --no-parallel
```

### Using Python API

```python
import json
from qcmanybody import ManyBodyComputer
from qcmanybody.models import ManyBodyInput

# Load input file
with open("test_inputs/water16_hmbe_23.json") as f:
    input_data = json.load(f)

# Convert to ManyBodyInput (requires schema conversion)
# Or use CLI instead for these test files

# Run calculation
result = ManyBodyComputer.from_manybodyinput(
    mbin,
    parallel=True,
    n_workers=4
)
```

### Lightweight Testing (Structure Only)

To test HMBE term enumeration **without running QC calculations**:

```bash
# Run structure tests (fast, no QC required)
pytest qcmanybody/tests/test_hmbe_parallel.py -v

# Test specific configurations
pytest qcmanybody/tests/test_hmbe_parallel.py::TestHMBEStructure::test_hmbe_23_term_count -v
```

---

## Expected Reference Results

From `hmbe_development/ref_results/energy_summary_water16.000.txt`:

| Method | Energy (Ha) | ΔE from Full (kcal/mol) |
|--------|-------------|-------------------------|
| Full RI-MP2 | -1220.4842172 | 0.000 (reference) |
| (2,3)-HMBE | -1220.4731682 | +6.933 |
| (2,4)-HMBE | -1220.4821682 | +1.286 |
| (3,4)-HMBE | -1220.4949320 | -6.724 |

**Note**: (3,4)-HMBE slightly overshoots the reference due to including more high-order terms than full MBE-4.

---

## Memory Requirements

### Filter Mode (default for <30 fragments)
- (2,3)-HMBE: ~200 KB memory
- (2,4)-HMBE: ~300 KB memory
- (3,4)-HMBE: ~400 KB memory

### Direct Mode (auto-selected for ≥30 fragments)
- All configurations: ~50-100 KB memory
- Recommended for parallel execution

**To force direct mode** (optional):
```json
"manybody": {
  "hmbe": {
    "enumeration_mode": "direct"
  }
}
```

---

## Parallel Execution Notes

**Safe configurations**:
- ✅ 2-4 workers: Recommended for laptop/workstation
- ✅ 4-8 workers: Recommended for compute node
- ⚠️ >8 workers: Limited by I/O for water-16 size

**Memory per worker**:
- Each worker needs copy of compute map (~200 KB with filter mode)
- Use direct enumeration mode to reduce memory: ~50 KB per worker
- Total memory = base + (n_workers × per_worker)

**Expected speedup**:
- 4 workers: ~3-3.5x speedup (75-85% efficiency)
- 8 workers: ~5-6x speedup (60-75% efficiency)
- Limited by serial task planning and I/O overhead

---

## Validation Checklist

When validating your HMBE implementation:

- [ ] **Term counts**: Check `hmbe_term_counts` vs `mbe_term_counts` in stats
- [ ] **Reduction factors**: Verify expected reduction (see table above)
- [ ] **Completeness**: No errors from `_validate_hmbe_completeness()`
- [ ] **Energy values**: Within ~10 kcal/mol of reference values
- [ ] **Parallel determinism**: Sequential == parallel results
- [ ] **Thread safety**: No race conditions or deadlocks
- [ ] **Memory usage**: Stable across runs, no leaks

---

## Troubleshooting

### "RuntimeError: HMBE completeness violation"
- **Cause**: Bug in sub-cluster addition (Schengen or direct enumeration)
- **Fix**: Report as bug - should never happen in production

### Out of memory errors
- **Solution**: Use `"enumeration_mode": "direct"` in HMBE spec
- For large systems: Always use direct mode with parallel execution

### Different results between sequential and parallel
- **Cause**: Non-deterministic ordering or floating point differences
- **Check**: Verify differences are <1e-10 (round-off errors)
- If larger: Report as thread-safety bug

### Poor parallel speedup
- **Cause**: System too small, I/O bottleneck, or overhead
- **Solution**: Water-16 is small; try larger systems for better speedup
- Expected efficiency: 70-85% for 4-8 workers

---

## Converting Other Input Files

To convert additional mini-app input files:

```bash
python scripts/convert_hmbe_input.py \
    hmbe_development/ref_results/your_input.inp \
    -o test_inputs/your_output.json \
    --bsse nocp \
    --parallel \
    --n-workers 4
```

See `scripts/convert_hmbe_input.py --help` for all options.

---

## Related Files

- `hmbe_development/ref_results/`: Reference energies and mini-app inputs
- `qcmanybody/tests/test_hmbe_parallel.py`: Comprehensive parallel tests
- `docs/HMBE_USER_GUIDE.md`: Full HMBE documentation
- `MEMORY_AUDIT.md`: Memory usage analysis and recommendations
