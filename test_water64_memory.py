#!/usr/bin/env python
"""Test that water64 doesn't crash due to memory in direct enumeration mode."""

import tracemalloc

# Start memory tracking
tracemalloc.start()

print("Loading water64 config...")
import sys
sys.stdout.flush()

from qcmanybody.cli.input_parser import parse_input_file
print("  - imported parse_input_file")
sys.stdout.flush()

from qcmanybody.cli.converter import convert_to_manybody_input  
print("  - imported convert_to_manybody_input")
sys.stdout.flush()

from qcmanybody.computer import ManyBodyComputer
print("  - imported ManyBodyComputer")
sys.stdout.flush()

# Load using CLI parser (same as qcmanybody run command)
print("Parsing input file...")
sys.stdout.flush()
cli_input = parse_input_file("water64_three_tier.json")
print("  - parsed input file")
sys.stdout.flush()

print("Converting to ManyBodyInput...")
sys.stdout.flush()

# Add timing
import time
t0 = time.time()
mb_input = convert_to_manybody_input(cli_input, input_file_path="water64_three_tier.json")
t1 = time.time()

print(f"  - converted to ManyBodyInput in {t1-t0:.1f}s")
sys.stdout.flush()

print(f"Config loaded. Creating ManyBodyComputer...")
print(f"  Fragments: {mb_input.specification.nfragments}")
print(f"  HMBE spec: {mb_input.specification.hmbe_spec.enumeration_mode if mb_input.specification.hmbe_spec else 'not set'}")
sys.stdout.flush()

# Create computer
print("Creating ManyBodyComputer...")
mbc = ManyBodyComputer(mb_input, "", "")

# Check memory before accessing compute_map
current, peak = tracemalloc.get_traced_memory()
print(f"Memory before compute_map: current={current/1024/1024:.1f}MB peak={peak/1024/1024:.1f}MB")

# Access compute_map (this is where it should crash on the old code)
print("Accessing compute_map...")
cmap = mbc.nbody.compute_map

# Check memory after
current, peak = tracemalloc.get_traced_memory()
print(f"Memory after compute_map: current={current/1024/1024:.1f}MB peak={peak/1024/1024:.1f}MB")

# Count terms
total_terms = 0
for mc, compute_dict in cmap.items():
    for bsse_type, nbody_dict in compute_dict.items():
        if bsse_type == "all":
            for nbody, terms in nbody_dict.items():
                print(f"  MC={mc} BSSE={bsse_type} nbody={nbody}: {len(terms)} terms")
                total_terms += len(terms)

print(f"\nTotal terms in compute_map: {total_terms}")
print("SUCCESS - did not crash!")

tracemalloc.stop()
