"""
QCManyBody CLI Molecule Loading Utilities

Handles loading molecules from various sources (inline, XYZ, PDB, QCSchema).
"""

import json
import logging
from pathlib import Path
from typing import List, Optional

from qcelemental.models import Molecule

from qcmanybody.cli.schemas.input_schema import (
    HMBEHierarchicalMoleculeSchema,
    MoleculeSchema,
    MoleculeSourceEnum,
)

logger = logging.getLogger(__name__)


class MoleculeLoadError(Exception):
    """Exception raised when molecule loading fails."""

    pass


def load_xyz_file(file_path: str, fragments: Optional[List[List[int]]] = None) -> Molecule:
    """
    Load molecule from XYZ file.

    Parameters
    ----------
    file_path : str
        Path to XYZ file
    fragments : Optional[List[List[int]]]
        Fragment definitions (0-based atom indices). If None, treats as single fragment.

    Returns
    -------
    Molecule
        QCElemental Molecule object

    Raises
    ------
    MoleculeLoadError
        If XYZ file cannot be loaded
    """
    path = Path(file_path)

    if not path.exists():
        raise MoleculeLoadError(f"XYZ file not found: {file_path}")

    try:
        # Read XYZ file
        with open(path, "r") as f:
            lines = f.readlines()

        # Parse XYZ format
        # Line 1: number of atoms
        # Line 2: comment
        # Lines 3+: symbol x y z

        if len(lines) < 3:
            raise MoleculeLoadError(f"XYZ file {file_path} is too short (need at least 3 lines)")

        natoms = int(lines[0].strip())

        if len(lines) < natoms + 2:
            raise MoleculeLoadError(
                f"XYZ file {file_path} claims {natoms} atoms but only has {len(lines) - 2} coordinate lines"
            )

        symbols = []
        geometry = []

        for i in range(natoms):
            line = lines[i + 2].strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 4:
                raise MoleculeLoadError(f"XYZ file {file_path} line {i+3} has too few fields (expected symbol x y z)")

            symbol = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])

            symbols.append(symbol)
            geometry.append([x, y, z])

        # Create fragments if not provided
        if fragments is None:
            fragments = [list(range(natoms))]
            logger.info(f"No fragments specified for XYZ file, treating as single fragment")
        else:
            # Validate fragments
            all_atoms = set()
            for frag in fragments:
                all_atoms.update(frag)
            if len(all_atoms) != natoms:
                raise MoleculeLoadError(
                    f"Fragment specification covers {len(all_atoms)} atoms but XYZ has {natoms} atoms"
                )

        # Create Molecule (XYZ files conventionally use Angstroms, need to convert to Bohr)
        from qcelemental import constants

        conversion = constants.conversion_factor("angstrom", "bohr")
        geometry_bohr = [[x * conversion for x in coord] for coord in geometry]

        mol = Molecule(symbols=symbols, geometry=geometry_bohr, fragments=fragments)
        logger.info(f"Loaded molecule from XYZ file: {natoms} atoms, {len(fragments)} fragments (units: angstrom)")

        return mol

    except (ValueError, IndexError) as e:
        raise MoleculeLoadError(f"Failed to parse XYZ file {file_path}: {e}") from e
    except IOError as e:
        raise MoleculeLoadError(f"Failed to read XYZ file {file_path}: {e}") from e


def load_qcschema_file(file_path: str, fragments: Optional[List[List[int]]] = None) -> Molecule:
    """
    Load molecule from QCSchema JSON file.

    Parameters
    ----------
    file_path : str
        Path to QCSchema JSON file
    fragments : Optional[List[List[int]]]
        Fragment definitions. If None, uses fragments from file or treats as single fragment.

    Returns
    -------
    Molecule
        QCElemental Molecule object

    Raises
    ------
    MoleculeLoadError
        If QCSchema file cannot be loaded
    """
    path = Path(file_path)

    if not path.exists():
        raise MoleculeLoadError(f"QCSchema file not found: {file_path}")

    try:
        with open(path, "r") as f:
            data = json.load(f)

        # Try to create Molecule from QCSchema
        mol = Molecule(**data)

        # Override fragments if provided
        if fragments is not None:
            mol = Molecule(**{**mol.dict(), "fragments": fragments})
            logger.info(f"Overrode fragments from file with user-specified fragments")

        # Ensure fragments are set
        if mol.fragments is None or len(mol.fragments) == 0:
            natoms = len(mol.symbols)
            mol = Molecule(**{**mol.dict(), "fragments": [list(range(natoms))]})
            logger.info(f"No fragments in QCSchema file, treating as single fragment")

        logger.info(f"Loaded molecule from QCSchema file: {len(mol.symbols)} atoms, {len(mol.fragments)} fragments")

        return mol

    except json.JSONDecodeError as e:
        raise MoleculeLoadError(f"Failed to parse QCSchema JSON file {file_path}: {e}") from e
    except Exception as e:
        raise MoleculeLoadError(f"Failed to create Molecule from QCSchema file {file_path}: {e}") from e


def load_pdb_file(file_path: str, fragments: Optional[List[List[int]]] = None) -> Molecule:
    """
    Load molecule from PDB file.

    Parameters
    ----------
    file_path : str
        Path to PDB file
    fragments : Optional[List[List[int]]]
        Fragment definitions. If None, attempts to detect from PDB chains or treats as single fragment.

    Returns
    -------
    Molecule
        QCElemental Molecule object

    Raises
    ------
    MoleculeLoadError
        If PDB file cannot be loaded

    Notes
    -----
    This is a basic PDB parser. For complex PDB files, consider using specialized tools.
    """
    path = Path(file_path)

    if not path.exists():
        raise MoleculeLoadError(f"PDB file not found: {file_path}")

    try:
        # Read PDB file and extract ATOM/HETATM records
        symbols = []
        geometry = []
        chains = []

        with open(path, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # PDB format (simplified):
                    # Columns:
                    #   13-14: element symbol (right-justified)
                    #   31-38: x coordinate
                    #   39-46: y coordinate
                    #   47-54: z coordinate
                    #   22: chain ID

                    # Extract element (columns 77-78, or fallback to atom name)
                    element = line[76:78].strip()
                    if not element:
                        # Fallback: use first 1-2 chars of atom name (columns 13-14)
                        element = line[12:14].strip()
                        # Remove digits
                        element = "".join([c for c in element if not c.isdigit()])

                    # Extract coordinates (angstroms)
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # Extract chain ID
                    chain = line[21:22].strip() or "A"

                    symbols.append(element)
                    geometry.append([x, y, z])
                    chains.append(chain)

        if not symbols:
            raise MoleculeLoadError(f"No atoms found in PDB file {file_path}")

        # Create fragments
        if fragments is None:
            # Try to auto-detect from chains
            unique_chains = sorted(set(chains))
            if len(unique_chains) > 1:
                fragments = []
                for chain in unique_chains:
                    frag_atoms = [i for i, c in enumerate(chains) if c == chain]
                    fragments.append(frag_atoms)
                logger.info(f"Auto-detected {len(fragments)} fragments from PDB chains")
            else:
                fragments = [list(range(len(symbols)))]
                logger.info(f"Single chain in PDB, treating as single fragment")
        else:
            # Validate fragments
            all_atoms = set()
            for frag in fragments:
                all_atoms.update(frag)
            if len(all_atoms) != len(symbols):
                raise MoleculeLoadError(
                    f"Fragment specification covers {len(all_atoms)} atoms but PDB has {len(symbols)} atoms"
                )

        # Create Molecule (PDB coordinates are in angstroms, need to convert to Bohr)
        from qcelemental import constants

        conversion = constants.conversion_factor("angstrom", "bohr")
        geometry_bohr = [[x * conversion for x in coord] for coord in geometry]

        mol = Molecule(symbols=symbols, geometry=geometry_bohr, fragments=fragments)
        logger.info(
            f"Loaded molecule from PDB file: {len(symbols)} atoms, {len(fragments)} fragments (units: angstrom)"
        )

        return mol

    except (ValueError, IndexError) as e:
        raise MoleculeLoadError(f"Failed to parse PDB file {file_path}: {e}") from e
    except IOError as e:
        raise MoleculeLoadError(f"Failed to read PDB file {file_path}: {e}") from e


def create_inline_molecule(mol_schema: MoleculeSchema) -> Molecule:
    """
    Create molecule from inline specification.

    Parameters
    ----------
    mol_schema : MoleculeSchema
        Molecule schema with inline specification

    Returns
    -------
    Molecule
        QCElemental Molecule object

    Raises
    ------
    MoleculeLoadError
        If inline specification is invalid
    """
    if mol_schema.inline is None:
        raise MoleculeLoadError("Inline molecule specification is missing")

    inline = mol_schema.inline

    try:
        # Convert geometry to Bohr if needed (Molecule always stores in Bohr)
        from qcelemental import constants

        geometry = inline.geometry
        if inline.units == "angstrom":
            # Convert from Angstroms to Bohr
            conversion = constants.conversion_factor("angstrom", "bohr")
            geometry = [[x * conversion for x in coord] for coord in geometry]

        mol = Molecule(
            symbols=inline.symbols,
            geometry=geometry,
            fragments=inline.fragments,
            molecular_charge=inline.molecular_charge,
            molecular_multiplicity=inline.molecular_multiplicity,
            fragment_charges=inline.fragment_charges,
            fragment_multiplicities=inline.fragment_multiplicities,
        )

        logger.info(
            f"Created inline molecule: {len(inline.symbols)} atoms, {len(inline.fragments)} fragments (units: {inline.units})"
        )

        return mol

    except Exception as e:
        raise MoleculeLoadError(f"Failed to create molecule from inline specification: {e}") from e


def convert_hmbe_hierarchical(hmbe_schema: HMBEHierarchicalMoleculeSchema) -> dict:
    """
    Convert HMBE hierarchical molecule schema to internal hierarchy format.

    This function converts the CLI schema format to the dictionary format expected
    by qcmanybody.hmbe.hierarchy.FragmentHierarchy.from_dict().

    Parameters
    ----------
    hmbe_schema : HMBEHierarchicalMoleculeSchema
        HMBE molecule specification from CLI input

    Returns
    -------
    dict
        Hierarchical structure dictionary compatible with FragmentHierarchy.from_dict()
        Contains keys: tiers, max_primary_per_nmer, fragments, units

    Notes
    -----
    The returned dictionary will be stored in ManyBodyKeywords.hmbe_hierarchy
    and processed by the HMBE preprocessor to build the flat molecule.
    """
    # Convert Pydantic model to dictionary format expected by FragmentHierarchy
    def convert_fragment(frag_schema):
        """Recursively convert fragment schema to dict."""
        frag_dict = {"id": frag_schema.id}

        # Add elementary fragment data if present (leaf nodes)
        if frag_schema.symbols is not None:
            frag_dict["symbols"] = frag_schema.symbols
        if frag_schema.geometry is not None:
            frag_dict["geometry"] = frag_schema.geometry
        if frag_schema.molecular_charge is not None:
            frag_dict["molecular_charge"] = frag_schema.molecular_charge
        if frag_schema.molecular_multiplicity is not None:
            frag_dict["molecular_multiplicity"] = frag_schema.molecular_multiplicity

        # Recursively convert sub-fragments
        if frag_schema.sub_fragments is not None and len(frag_schema.sub_fragments) > 0:
            frag_dict["sub_fragments"] = [convert_fragment(sub) for sub in frag_schema.sub_fragments]

        return frag_dict

    # Build hierarchical structure dictionary
    hierarchy_dict = {
        "tiers": hmbe_schema.tiers,
        "max_primary_per_nmer": hmbe_schema.max_primary_per_nmer,
        "fragments": [convert_fragment(frag) for frag in hmbe_schema.fragments],
        "units": hmbe_schema.units,
    }

    logger.info(
        f"Converted HMBE hierarchical molecule: {hmbe_schema.tiers} tiers, "
        f"{len(hmbe_schema.fragments)} primary fragments, units: {hmbe_schema.units}"
    )

    return hierarchy_dict


def load_molecule(mol_schema: MoleculeSchema):
    """
    Load molecule according to schema specification.

    This is the main entry point for molecule loading.

    Parameters
    ----------
    mol_schema : MoleculeSchema
        Molecule specification from input file

    Returns
    -------
    Molecule or dict
        For standard molecule sources (inline, xyz, pdb, qcschema): QCElemental Molecule object
        For hmbe_hierarchical: dict with hierarchical structure (will be processed by converter)

    Raises
    ------
    MoleculeLoadError
        If molecule cannot be loaded

    Notes
    -----
    For HMBE hierarchical molecules, this returns a dict rather than a Molecule.
    The converter.py module detects this and stores the hierarchy in ManyBodyKeywords.hmbe_hierarchy.
    The HMBE preprocessor then builds the flat molecule during computation.

    Examples
    --------
    >>> from qcmanybody.cli.schemas.input_schema import MoleculeSchema, MoleculeSourceEnum
    >>> mol_schema = MoleculeSchema(source=MoleculeSourceEnum.xyz, file="water.xyz")
    >>> mol = load_molecule(mol_schema)
    """
    source = mol_schema.source

    if source == MoleculeSourceEnum.inline:
        return create_inline_molecule(mol_schema)

    elif source == MoleculeSourceEnum.xyz:
        if mol_schema.file is None:
            raise MoleculeLoadError("XYZ file path not specified")
        return load_xyz_file(mol_schema.file, mol_schema.fragments)

    elif source == MoleculeSourceEnum.qcschema:
        if mol_schema.file is None:
            raise MoleculeLoadError("QCSchema file path not specified")
        return load_qcschema_file(mol_schema.file, mol_schema.fragments)

    elif source == MoleculeSourceEnum.pdb:
        if mol_schema.file is None:
            raise MoleculeLoadError("PDB file path not specified")
        return load_pdb_file(mol_schema.file, mol_schema.fragments)

    elif source == MoleculeSourceEnum.hmbe_hierarchical:
        if mol_schema.hmbe is None:
            raise MoleculeLoadError("HMBE hierarchy specification not provided")
        # Return hierarchical structure dict (not a Molecule)
        # The converter will detect this and store it in ManyBodyKeywords.hmbe_hierarchy
        return convert_hmbe_hierarchical(mol_schema.hmbe)

    else:
        raise MoleculeLoadError(f"Unsupported molecule source: {source}")
