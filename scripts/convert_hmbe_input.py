#!/usr/bin/env python3
"""
Convert old HMBE mini-app input files to QCManyBody CLI input format.

This script transforms hierarchical fragment input files (with PF/SF/TF structure)
into the QCManyBody CLI JSON format for HMBE calculations.

Usage:
    python convert_hmbe_input.py input.inp -o output.json
    python convert_hmbe_input.py input.inp  # outputs to stdout
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert old HMBE mini-app input files to QCManyBody CLI format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s water16.000.inp -o water16_qcmb.json
    %(prog)s water64_isomer1_three_tier.inp --parallel --n-workers 8
    %(prog)s input.inp --bsse cp vmfc --schengen --schengen-fraction 0.2
        """,
    )
    parser.add_argument("input_file", type=Path, help="Input file in old HMBE format (.inp)")
    parser.add_argument("-o", "--output", type=Path, help="Output JSON file (default: stdout)")
    parser.add_argument(
        "--bsse",
        nargs="+",
        default=["nocp"],
        choices=["nocp", "cp", "vmfc"],
        help="BSSE correction type(s) (default: nocp)",
    )
    parser.add_argument(
        "--units",
        default="angstrom",
        choices=["angstrom", "bohr"],
        help="Geometry units (default: angstrom)",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        default=None,
        help="Enable parallel execution (overrides input file setting)",
    )
    parser.add_argument(
        "--no-parallel",
        action="store_true",
        help="Disable parallel execution (overrides input file setting)",
    )
    parser.add_argument("--n-workers", type=int, default=4, help="Number of parallel workers (default: 4)")
    parser.add_argument(
        "--executor-type",
        default="multiprocessing",
        choices=["multiprocessing", "concurrent_futures", "dask"],
        help="Parallel executor type (default: multiprocessing)",
    )
    parser.add_argument(
        "--return-total-data",
        action="store_true",
        default=True,
        help="Return total data (default: True)",
    )
    parser.add_argument(
        "--no-return-total-data",
        action="store_true",
        help="Don't return total data",
    )
    parser.add_argument(
        "--schengen",
        action="store_true",
        help="Enable Schengen filtering",
    )
    parser.add_argument(
        "--schengen-fraction",
        type=float,
        default=0.1,
        help="Schengen selection fraction (default: 0.1)",
    )
    parser.add_argument(
        "--schengen-metric",
        default="R2",
        choices=["R1", "R2", "R3"],
        help="Schengen distance metric (default: R2)",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default=None,
        help="Results output file path",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        default=True,
        help="Pretty print output JSON (default: True)",
    )

    return parser.parse_args()


def load_old_format(input_path: Path) -> dict[str, Any]:
    """Load and parse the old HMBE input format."""
    with open(input_path, "r") as f:
        return json.load(f)


def determine_hierarchy_depth(fragments: list[dict]) -> int:
    """
    Determine the depth of the hierarchical fragment structure.
    
    Returns the number of tiers (2 for PF->SF, 3 for PF->SF->TF, etc.)
    """
    if not fragments:
        return 0
    
    # Check first primary fragment
    pf = fragments[0]
    if "fragments" not in pf or not pf["fragments"]:
        return 1
    
    # Check first secondary fragment
    sf = pf["fragments"][0]
    if "fragments" not in sf or not sf["fragments"]:
        return 2
    
    # Has tertiary fragments
    tf = sf["fragments"][0]
    if "fragments" not in tf or not tf["fragments"]:
        return 3
    
    # Could have more levels, but typically 2-3 tiers
    return 4


def get_tier_names(num_tiers: int) -> list[str]:
    """Generate tier names based on the number of tiers."""
    base_names = ["domain", "subdomain", "fragment", "atom"]
    return base_names[:num_tiers]


def flatten_fragments(
    fragments: list[dict], depth: int
) -> tuple[list[str], list[list[float]], list[list[int]], dict[str, list[str]]]:
    """
    Flatten the hierarchical fragment structure into a flat list of atoms.
    
    Returns:
        - symbols: flat list of all atom symbols
        - geometry: list of [x, y, z] coordinates for each atom
        - fragment_indices: list of atom indices for each leaf fragment
        - fragment_tiers: mapping of fragment index to tier path
    """
    symbols: list[str] = []
    geometry: list[list[float]] = []
    fragment_indices: list[list[int]] = []
    fragment_tiers: dict[str, list[str]] = {}
    
    atom_index = 0
    fragment_index = 1  # 1-indexed for QCManyBody
    
    def process_fragment(
        frag: dict, tier_indices: list[int], current_depth: int, target_depth: int
    ) -> None:
        """
        Process a fragment recursively.
        
        Args:
            frag: The fragment dictionary
            tier_indices: List of 0-indexed tier positions [pf_idx, sf_idx, tf_idx, ...]
            current_depth: Current depth in hierarchy (0=PF, 1=SF, 2=TF, ...)
            target_depth: Maximum depth to traverse
        """
        nonlocal atom_index, fragment_index
        
        if current_depth == target_depth or "fragments" not in frag or not frag["fragments"]:
            # This is a leaf fragment - extract atoms
            if "symbols" in frag and "geometry" in frag:
                frag_symbols = frag["symbols"]
                frag_geometry = frag["geometry"]
                
                # Geometry is stored as flat list [x1, y1, z1, x2, y2, z2, ...]
                num_atoms = len(frag_symbols)
                indices = []
                
                for i in range(num_atoms):
                    symbols.append(frag_symbols[i])
                    # Extract xyz for this atom
                    x = frag_geometry[i * 3]
                    y = frag_geometry[i * 3 + 1]
                    z = frag_geometry[i * 3 + 2]
                    geometry.append([x, y, z])
                    indices.append(atom_index)
                    atom_index += 1
                
                fragment_indices.append(indices)
                
                # Build tier path from indices
                # Format: ["G0", "G0_S1"] for 2-tier, ["G0", "G0_S1", "G0_S1_T2"] for 3-tier
                tier_path = build_tier_path(tier_indices)
                fragment_tiers[str(fragment_index)] = tier_path
                fragment_index += 1
        else:
            # Recurse into sub-fragments
            for i, sub_frag in enumerate(frag["fragments"]):
                new_tier_indices = tier_indices + [i]
                process_fragment(sub_frag, new_tier_indices, current_depth + 1, target_depth)
    
    def build_tier_path(tier_indices: list[int]) -> list[str]:
        """
        Build tier path from indices.
        
        For a 2-tier system with tier_indices [0, 2] (PF1, SF3 -> 0-indexed: 0, 2):
            Returns: ["G0", "G0_S2"]
        
        For a 3-tier system with tier_indices [1, 0, 3] (PF2, SF1, TF4 -> 0-indexed: 1, 0, 3):
            Returns: ["G1", "G1_S0", "G1_S0_T3"]
        """
        prefixes = ["G", "S", "T", "Q", "R"]  # Group, Sub, Tertiary, Quaternary, ...
        tier_path = []
        
        for level in range(len(tier_indices)):
            # Build cumulative tier code up to this level
            codes = []
            for i in range(level + 1):
                prefix = prefixes[i] if i < len(prefixes) else f"L{i}"
                codes.append(f"{prefix}{tier_indices[i]}")
            tier_path.append("_".join(codes))
        
        return tier_path
    
    # Process all primary fragments
    for i, pf in enumerate(fragments):
        process_fragment(pf, [i], 0, depth)
    
    return symbols, geometry, fragment_indices, fragment_tiers


def build_fragment_tiers_from_structure(
    fragments: list[dict], num_tiers: int
) -> dict[str, list[str]]:
    """
    Build fragment_tiers mapping from hierarchical structure.
    
    The tier path for each leaf fragment encodes its position in the hierarchy.
    """
    fragment_tiers: dict[str, list[str]] = {}
    fragment_index = 1
    
    def traverse(frag: dict, tier_path: list[tuple[str, int]], current_depth: int) -> None:
        nonlocal fragment_index
        
        if "fragments" not in frag or not frag["fragments"]:
            # Leaf fragment
            # Build tier codes from path
            codes = []
            prefixes = ["G", "S", "T", "Q"]
            for i, (_, idx) in enumerate(tier_path):
                prefix = prefixes[i] if i < len(prefixes) else f"L{i}"
                codes.append(f"{prefix}{idx}")
            
            fragment_tiers[str(fragment_index)] = codes
            fragment_index += 1
        else:
            for i, sub_frag in enumerate(frag["fragments"]):
                new_path = tier_path + [(frag.get("id", f"D{current_depth}"), i)]
                traverse(sub_frag, new_path, current_depth + 1)
    
    for i, pf in enumerate(fragments):
        traverse(pf, [(pf.get("id", f"PF{i+1}"), i)], 0)
    
    return fragment_tiers


def convert_to_qcmanybody_format(
    old_data: dict[str, Any], args: argparse.Namespace
) -> dict[str, Any]:
    """Convert old HMBE format to QCManyBody CLI format."""
    chemical_system = old_data["chemical_system"]
    runtime_params = old_data["runtime_params"]
    
    fragments = chemical_system["fragments"]
    num_tiers = determine_hierarchy_depth(fragments)
    
    # Flatten the hierarchical structure
    symbols, geometry, fragment_indices, fragment_tiers = flatten_fragments(
        fragments, num_tiers
    )
    
    # Get molecular charge and multiplicity from the top level
    # (use first fragment's values or defaults)
    mol_charge = 0
    mol_mult = 1
    if fragments:
        mol_charge = fragments[0].get("molecular_charge", 0)
        mol_mult = fragments[0].get("molecular_multiplicity", 1)
    
    # Determine parallel setting
    if args.no_parallel:
        parallel = False
    elif args.parallel is not None:
        parallel = args.parallel
    else:
        parallel = runtime_params.get("parallel", False)
    
    # Determine return_total_data
    return_total_data = not args.no_return_total_data
    
    # Build the output structure
    output: dict[str, Any] = {
        "schema_name": "qcmanybody_cli_input",
        "schema_version": 1,
        "molecule": {
            "source": "inline",
            "inline": {
                "symbols": symbols,
                "geometry": geometry,
                "fragments": fragment_indices,
                "units": args.units,
                "molecular_charge": float(mol_charge),
                "molecular_multiplicity": mol_mult,
            },
        },
        "calculation": {
            "type": "single",
            "single": {
                "driver": runtime_params.get("driver", "energy"),
                "method": runtime_params["model"]["method"].lower(),
                "basis": runtime_params["model"]["basis"],
                "program": runtime_params.get("program", "psi4"),
            },
        },
        "bsse": {"type": args.bsse},
        "manybody": {
            "max_nbody": max(runtime_params.get("mbe_orders", [2])),
            "return_total_data": return_total_data,
            "supersystem_ie_only": False,
            "hmbe": {
                "truncation_orders": runtime_params.get("mbe_orders", [2]),
                "hierarchy": {
                    "num_tiers": num_tiers,
                    "fragment_tiers": fragment_tiers,
                    "tier_names": get_tier_names(num_tiers),
                },
            },
        },
        "execution": {
            "parallel": parallel,
        },
        "output": {
            "format": "json",
            "pretty": args.pretty,
            "include_timings": True,
        },
    }
    
    # Add keywords if present
    if runtime_params.get("keywords"):
        output["calculation"]["single"]["keywords"] = runtime_params["keywords"]
    
    # Add parallel execution settings if enabled
    if parallel:
        output["execution"]["n_workers"] = args.n_workers
        output["execution"]["executor_type"] = args.executor_type
    
    # Add Schengen settings if enabled
    if args.schengen:
        output["manybody"]["hmbe"]["schengen"] = {
            "enabled": True,
            "selection_fraction": args.schengen_fraction,
            "distance_metric": args.schengen_metric,
        }
    
    # Add output file if specified
    if args.output_file:
        output["output"]["file"] = args.output_file
    elif args.output:
        # Generate default output filename from input
        output["output"]["file"] = str(args.output.stem) + "_results.json"
    
    return output


def main() -> None:
    """Main entry point."""
    args = parse_arguments()
    
    # Validate input file
    if not args.input_file.exists():
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Load old format
    try:
        old_data = load_old_format(args.input_file)
    except json.JSONDecodeError as e:
        print(f"Error: Failed to parse input file as JSON: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to read input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate required fields
    if "chemical_system" not in old_data:
        print("Error: Input file missing 'chemical_system' field", file=sys.stderr)
        sys.exit(1)
    if "runtime_params" not in old_data:
        print("Error: Input file missing 'runtime_params' field", file=sys.stderr)
        sys.exit(1)
    
    # Convert to new format
    try:
        new_data = convert_to_qcmanybody_format(old_data, args)
    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Output result
    indent = 2 if args.pretty else None
    json_output = json.dumps(new_data, indent=indent)
    
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write(json_output)
        print(f"Converted input written to: {args.output}", file=sys.stderr)
    else:
        print(json_output)


if __name__ == "__main__":
    main()
