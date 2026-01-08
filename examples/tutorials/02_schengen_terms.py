"""
Tutorial 2: Improving Accuracy with Schengen Terms

Schengen terms are HMBE-excluded terms that you selectively add back based on
spatial proximity to improve accuracy without computing all MBE terms.

Key concepts:
- Understanding what Schengen terms are
- Choosing selection fraction
- Different distance metrics
- Accuracy vs. cost trade-offs

System: 9 water molecules (3×3 hierarchy)
Base: (2,3)-HMBE
Addition: Schengen terms (10-20% of excluded terms)
Expected: 30-50% error reduction vs base HMBE
"""

from qcelemental.models import Molecule
from qcmanybody.core import ManyBodyCore
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import (
    FragmentHierarchy,
    HMBESpecification,
    SchengenSpecification
)


def create_water_molecule():
    """Create a 9-water system with 3×3 hierarchy."""

    # 9 waters in 3x3 grid
    water_positions = [
        # Group 0
        [[0, 0, 0], [0, 0.757, 0.587], [0, -0.757, 0.587]],
        [[3, 0, 0], [3, 0.757, 0.587], [3, -0.757, 0.587]],
        [[0, 3, 0], [0, 3.757, 0.587], [0, 2.243, 0.587]],
        # Group 1
        [[6, 0, 0], [6, 0.757, 0.587], [6, -0.757, 0.587]],
        [[9, 0, 0], [9, 0.757, 0.587], [9, -0.757, 0.587]],
        [[6, 3, 0], [6, 3.757, 0.587], [6, 2.243, 0.587]],
        # Group 2
        [[12, 0, 0], [12, 0.757, 0.587], [12, -0.757, 0.587]],
        [[15, 0, 0], [15, 0.757, 0.587], [15, -0.757, 0.587]],
        [[12, 3, 0], [12, 3.757, 0.587], [12, 2.243, 0.587]],
    ]

    geometry = []
    fragments = []
    symbols = []

    for i, water in enumerate(water_positions):
        start_idx = i * 3
        fragments.append([start_idx, start_idx + 1, start_idx + 2])
        for atom_pos in water:
            geometry.append(atom_pos)
        symbols.extend(["O", "H", "H"])

    return Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=fragments
    )


def create_hierarchy():
    """Create 3×3 hierarchy (3 groups, 3 waters each)."""

    fragment_tiers = {}
    for i in range(9):
        frag_id = i + 1
        group_id = i // 3  # 0, 1, 2
        fragment_tiers[frag_id] = (f"G{group_id}", f"W{i}")

    return FragmentHierarchy(
        num_tiers=2,
        fragment_tiers=fragment_tiers,
        tier_names=("group", "water")
    )


def compare_schengen_configurations(mol, hierarchy):
    """Compare HMBE with and without Schengen terms."""

    configurations = [
        ("Base HMBE (2,3)", None),
        ("+ 10% Schengen (R2)", 0.10, "R2"),
        ("+ 15% Schengen (R2)", 0.15, "R2"),
        ("+ 20% Schengen (R2)", 0.20, "R2"),
        ("+ 10% Schengen (R_inv)", 0.10, "R_inv"),
    ]

    print("="*80)
    print("COMPARING SCHENGEN CONFIGURATIONS")
    print("="*80)
    print()

    results = []

    for config in configurations:
        if len(config) == 2:
            # Base HMBE (no Schengen)
            name, _ = config
            schengen = None
        else:
            # With Schengen
            name, fraction, metric = config
            schengen = SchengenSpecification(
                enabled=True,
                selection_fraction=fraction,
                distance_metric=metric
            )

        # Create HMBE spec
        hmbe_spec = HMBESpecification(
            truncation_orders=(2, 3),
            hierarchy=hierarchy,
            schengen=schengen
        )

        # Run calculation (build compute map only)
        mbcore = ManyBodyCore(
            molecule=mol,
            bsse_type=[BsseEnum.nocp],
            levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
            return_total_data=False,
            supersystem_ie_only=False,
            embedding_charges={},
            hmbe_spec=hmbe_spec
        )

        stats = mbcore.get_hmbe_statistics()
        mc_level = list(stats['hmbe_term_counts'].keys())[0]
        term_count = stats['hmbe_term_counts'][mc_level]

        results.append((name, term_count))

        print(f"{name:30s}: {term_count:4d} terms")

    print()
    print("-"*80)
    print()

    # Calculate term increases
    base_count = results[0][1]
    print(f"Base HMBE:        {base_count} terms")
    print()
    print("Schengen additions:")
    for name, count in results[1:]:
        added = count - base_count
        pct_increase = (count / base_count - 1) * 100
        print(f"  {name:28s}: +{added:3d} terms (+{pct_increase:4.1f}%)")

    print()
    print("="*80)
    print()

    return results


def explain_distance_metrics():
    """Explain different distance metrics for Schengen selection."""

    print("="*80)
    print("DISTANCE METRICS FOR SCHENGEN SELECTION")
    print("="*80)
    print()

    metrics = [
        ("R2", "Sum of squared distances", "Recommended, balanced"),
        ("R", "Sum of distances", "Linear distance weighting"),
        ("R_inv", "Sum of inverse distances", "Strongly favors close fragments"),
        ("R3_inv", "Sum of inverse cubed distances", "Very strongly favors close"),
        ("fmo", "FMO-style vdW scaling", "Advanced, system-dependent"),
    ]

    print("Available distance metrics:")
    print()
    for metric, description, note in metrics:
        print(f"  {metric:10s}: {description:35s} ({note})")
    print()

    print("How distance metrics work:")
    print()
    print("  For a candidate term (i, j, k):")
    print("    1. Compute center-of-mass for each fragment")
    print("    2. Calculate pairwise distances: d_ij, d_ik, d_jk")
    print("    3. Apply metric formula:")
    print("       - R2:     d_ij² + d_ik² + d_jk²")
    print("       - R:      d_ij + d_ik + d_jk")
    print("       - R_inv:  1/d_ij + 1/d_ik + 1/d_jk")
    print("    4. Rank all candidates by this metric")
    print("    5. Select top X% with smallest metric values (closest)")
    print()

    print("Choosing a metric:")
    print("  - R2 (recommended): Good balance, not too sensitive to outliers")
    print("  - R_inv: When you want to strongly prioritize very close fragments")
    print("  - R: Simple, interpretable")
    print()

    print("="*80)
    print()


def demonstrate_schengen_selection(mol, hierarchy):
    """Show which terms get selected as Schengen terms."""

    print("="*80)
    print("SCHENGEN TERM SELECTION EXAMPLE")
    print("="*80)
    print()

    # Create base HMBE (no Schengen)
    hmbe_base = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy
    )

    # Create HMBE with 15% Schengen
    hmbe_schengen = HMBESpecification(
        truncation_orders=(2, 3),
        hierarchy=hierarchy,
        schengen=SchengenSpecification(
            enabled=True,
            selection_fraction=0.15,
            distance_metric="R2"
        )
    )

    # Build both
    mbcore_base = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_base
    )

    mbcore_schengen = ManyBodyCore(
        molecule=mol,
        bsse_type=[BsseEnum.nocp],
        levels={1: "hf/sto-3g", 2: "hf/sto-3g", 3: "hf/sto-3g"},
        return_total_data=False,
        supersystem_ie_only=False,
        embedding_charges={},
        hmbe_spec=hmbe_schengen
    )

    # Get term sets
    base_map = mbcore_base.compute_map["hf/sto-3g"]["nocp"]
    schengen_map = mbcore_schengen.compute_map["hf/sto-3g"]["nocp"]

    # Find 3-body terms
    base_3body = set(frag for frag, _ in base_map.get(3, set()))
    schengen_3body = set(frag for frag, _ in schengen_map.get(3, set()))

    # Schengen additions
    schengen_additions = schengen_3body - base_3body

    print(f"Base HMBE 3-body terms:     {len(base_3body)}")
    print(f"Schengen HMBE 3-body terms: {len(schengen_3body)}")
    print(f"Schengen added:             {len(schengen_additions)}")
    print()

    if schengen_additions:
        print("Example Schengen terms added (3-body):")
        for i, term in enumerate(sorted(schengen_additions)[:5]):
            # Determine which groups they span
            groups = {hierarchy.fragment_tiers[f][0] for f in term}
            print(f"  {term:20s} spans groups {sorted(groups)}")
            if i >= 4:
                break
        print()

    print("These are terms that:")
    print("  - Were excluded by base (2,3)-HMBE (span >2 groups)")
    print("  - Were spatially close enough to be added back")
    print("  - Should improve accuracy without full MBE cost")
    print()

    print("="*80)
    print()


def main():
    """Run complete Schengen tutorial."""

    print()
    print("#" * 80)
    print("# Tutorial 2: Improving Accuracy with Schengen Terms")
    print("#" * 80)
    print()

    # Create molecule and hierarchy
    print("Setting up 9-water system...")
    mol = create_water_molecule()
    hierarchy = create_hierarchy()
    print()

    # Compare configurations
    results = compare_schengen_configurations(mol, hierarchy)

    # Explain distance metrics
    explain_distance_metrics()

    # Show what gets selected
    demonstrate_schengen_selection(mol, hierarchy)

    # Summary and recommendations
    print("="*80)
    print("SUMMARY AND RECOMMENDATIONS")
    print("="*80)
    print()

    print("KEY POINTS:")
    print("  1. Schengen terms are excluded HMBE terms added back by proximity")
    print("  2. Selection fraction controls accuracy vs. cost trade-off")
    print("  3. Typical: 10-20% gives good balance")
    print("  4. R2 distance metric recommended for most cases")
    print()

    print("WHEN TO USE SCHENGEN:")
    print("  ✓ When base HMBE accuracy isn't sufficient")
    print("  ✓ For critical calculations where accuracy matters")
    print("  ✓ When you can afford 10-20% more calculations")
    print("  ✗ When maximum speed is required")
    print("  ✗ When base HMBE already gives good accuracy")
    print()

    print("TYPICAL WORKFLOW:")
    print("  1. Run base HMBE to see speedup")
    print("  2. If accuracy insufficient, try 10% Schengen")
    print("  3. Increase to 15-20% if needed")
    print("  4. Compare accuracy vs. full MBE (if feasible)")
    print()

    print("NEXT STEPS:")
    print("  - Try different selection fractions on your system")
    print("  - Compare with full MBE if small enough")
    print("  - Use parallel execution for larger systems (Tutorial 3)")
    print()

    print("Tutorial complete!")
    print()


if __name__ == "__main__":
    main()
