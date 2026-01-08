# QCManyBody Development TODOs

This file consolidates TODO comments from the codebase, organized by priority and area.

## High Priority

### Data Validation
- [ ] **Enforce point charge sum validation** (computer.py:167)
  - Enforce that point charge sum equals fragment_charges value
  - Impact: Data integrity

- [ ] **Embedding charges validation** (models/v2/many_body.py:150-151, models/v1/manybody_input_pydv1.py:169-170)
  - Embedding charges should sum to fragment charge - enforce?
  - Embedding charges irrelevant to CP (basis sets always present)? - clarify
  - Impact: Correctness of embedding calculations

## Medium Priority

### HMBE Features

- [ ] **Implement FMO-style distance metric** (hmbe_filter.py:287)
  - Add FMO-style scaling with van der Waals radii for Schengen term selection
  - Currently falls back to R2 metric
  - Impact: Better Schengen term selection for systems where vdW radii are important

### Code Quality & Architecture

- [ ] **Consolidate vmfc lists** (builder.py:125)
  - vmfc_compute_list and vmfc_level_list are identical, consolidate
  - Impact: Code simplification, reduced duplication

- [ ] **Multilevel and supersystem_ie_only interaction** (core.py:93)
  - multilevel and supersystem_ie_only=T not allowed together - enforce
  - Impact: Prevent invalid configurations

- [ ] **Supersystem in levels handling** (core.py:94)
  - supersystem in levels is not to be trusted -- nfrag only and skips levels
  - Need to review and fix or document limitations
  - Impact: Correctness

- [ ] **Reduce duplication in higher nbody** (core.py:253)
  - "This is a bit of a hack. Lots of duplication when reaching higher nbody"
  - Impact: Performance and code clarity

- [ ] **Reconsider nbodies_per_mc_level structure** (computer.py:358)
  - Perhaps change nbodies_per_mc_level into dict of lists so that position/label indexing coincides
  - Impact: API consistency

- [ ] **Review non-contiguous nbody levels** (core.py:119)
  - Test and reenable if appropriate - noncontiguous nb might be fine on core computing side
  - Impact: Feature enablement

### Data Structure Improvements

- [ ] **Use Lori's code for level handling** (core.py:91)
  - Consider using existing implementation
  - Impact: Code reuse

- [ ] **Convert levels dict to list of lists** (core.py:92)
  - Handle non-contiguous levels better
  - Impact: Flexibility

- [ ] **Add computed_field decorators** (computer.py:206, 284-291)
  - Document computed fields with proper Pydantic v2 decorators
  - Extensive comments suggest desired documentation for nbodies_per_mc_level
  - Impact: Better API documentation

## Low Priority

### Test Optimizations

- [ ] **Reduce task counts in tests** (test_computer_he4_singlelevel.py:486, 528, 570)
  - Several tests note that task counts could be reduced (65→61, 50→46, 22→18)
  - Impact: Test performance

- [ ] **Reduce task counts in multilevel tests** (test_computer_he4_multilevel.py:646-647, 695-696)
  - Similar optimization opportunities: 22→18, 4→0
  - Impact: Test performance

### Documentation & Cleanup

- [ ] **Review VMFC-CORRECTED TOTAL ENERGY labels** (multiple test files)
  - Many tests have "TODO remove?" on VMFC-CORRECTED TOTAL ENERGY labels
  - Appears to be legacy naming
  - Files: test_computer_he4_singlelevel.py:297, 326, 353, 380
  - Files: test_computer_he4_multilevel.py:396, 404, 421
  - Impact: Code clarity

- [ ] **Clarify supersystem_ie_only behavior** (test_schema_keywords.py:112)
  - When supersystem_ie_only=T, nbodies_per_mc_level and levels isn't really accurate
  - Document or fix
  - Impact: Documentation accuracy

- [ ] **Enable IE for HET4 gradient** (test_computer_het4_gradient.py:88)
  - "IE should be present" - commented out assertion
  - Investigate why and either fix or remove
  - Impact: Test coverage

- [ ] **Document 1-body special cases** (models/v2/many_body.py:540-541, models/v1/manybody_output_pydv1.py:226-227)
  - Note that TOT 1BODY cp=nocp=vmfc
  - Note that summ INTERACTION ENERGY props return 0.0 for max_nbody=1 for completeness
  - Impact: Documentation

- [ ] **Review component_results handling** (computer.py:591-592)
  - Commented out code for component_results with TODO about when/where to include individual outputs
  - Impact: API clarity

- [ ] **Consider extras vs qcvars** (computer.py:586)
  - All besides nbody may be better candidates for extras than qcvars
  - energy/gradient/hessian_body_dict in particular are too simple for qcvars
  - Impact: Data model design

### Minor Items

- [ ] **Remove obsolete code?** (computer.py:45)
  - "can remove?" - review and remove if obsolete
  - Impact: Code cleanup

- [ ] **Simplify levels handling** (computer.py:272)
  - Consider: levels = {plan.max_nbody: method}
  - Impact: Code simplification

- [ ] **Review logic movement** (computer.py:341)
  - "reconsider logic. move this from levels to here?"
  - Impact: Code organization

- [ ] **Review min() logic** (computer.py:346)
  - "once was return min(v, nfragments)"
  - Document why it changed or restore if appropriate
  - Impact: Correctness

## Technical Debt

### Pydantic Migration (Deferred)
- [ ] **Migrate from pydantic.v1 to Pydantic v2**
  - Current codebase uses pydantic.v1 imports throughout
  - Migration requires updating ~45 files
  - Primary focus: qcmanybody/cli/ module
  - See /home/westh/.claude/plans/rippling-sprouting-wind.md for detailed migration plan
  - Impact: Remove technical debt, enable Pydantic v2 features

## Summary Statistics

- **Total TODOs in code**: 52
- **High Priority**: 3 items (data validation)
- **Medium Priority**: 13 items (HMBE features, architecture and code quality)
- **Low Priority**: 20 items (optimizations and documentation)
- **Technical Debt**: 1 major item (Pydantic v2)

## Contributing

When addressing TODOs:
1. Reference this file in commit messages (e.g., "Fixes TODO.md #5")
2. Update this file when TODOs are completed
3. Add new TODOs here rather than inline code comments (except for critical issues)
4. Consider creating GitHub issues for high-priority items
