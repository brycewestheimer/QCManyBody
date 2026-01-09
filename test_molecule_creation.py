#!/usr/bin/env python3
"""Test QCElemental Molecule creation with 64 fragments."""

import json
import time
import sys

print("Loading input file...", file=sys.stderr)
with open("water64_three_tier.json") as f:
    data = json.load(f)

inline = data["molecule"]["inline"]

print(f"Loaded: {len(inline['symbols'])} atoms, {len(inline['fragments'])} fragments", file=sys.stderr)
print(f"Geometry size: {len(inline['geometry'])}", file=sys.stderr)

# Convert geometry units
from qcelemental import constants

geometry = inline["geometry"]
if inline["units"] == "angstrom":
    conversion = constants.conversion_factor("angstrom", "bohr")
    geometry = [[x * conversion for x in coord] for coord in geometry]
    print("Converted geometry to Bohr", file=sys.stderr)

# Try creating Molecule
from qcelemental.models import Molecule

print("About to create Molecule...", file=sys.stderr)
sys.stderr.flush()

start = time.time()
mol = Molecule(
    symbols=inline["symbols"],
    geometry=geometry,
    fragments=inline["fragments"],
    molecular_charge=inline["molecular_charge"],
    molecular_multiplicity=inline["molecular_multiplicity"],
    fragment_charges=inline.get("fragment_charges"),
    fragment_multiplicities=inline.get("fragment_multiplicities"),
)
elapsed = time.time() - start

print(f"Created Molecule in {elapsed:.2f}s", file=sys.stderr)
print(f"Molecule has {len(mol.symbols)} atoms, {len(mol.fragments)} fragments", file=sys.stderr)
