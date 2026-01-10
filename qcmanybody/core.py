from __future__ import annotations

import logging
import math
import os
import string
from collections import Counter, defaultdict
from typing import Any, Dict, Iterable, Literal, Mapping, Optional, Sequence, Set, Tuple, Union

import numpy as np
from qcelemental.models import Molecule

from qcmanybody.builder import build_nbody_compute_list
from qcmanybody.dependency import NBodyDependencyGraph
from qcmanybody.models.v1 import BsseEnum
from qcmanybody.models.hierarchy import HMBESpecification
from qcmanybody.utils import (
    all_same_shape,
    collect_vars,
    copy_value,
    delabeler,
    find_shape,
    labeler,
    modelchem_labels,
    print_nbody_energy,
    resize_gradient,
    resize_hessian,
    shaped_zero,
    sum_cluster_data,
)

logger = logging.getLogger(__name__)


from dataclasses import dataclass


@dataclass
class RawMoleculeData:
    """Stores molecule data without QCElemental validation.

    This lightweight storage class holds raw molecular data (symbols, geometry, fragments)
    without triggering QCElemental's expensive fragment validation. This enables QCManyBody
    to work with systems having 64+ fragments, which would otherwise crash during Molecule
    construction due to O(n²) validation overhead in QCElemental.

    Attributes
    ----------
    symbols : list
        Atom element symbols (e.g., ['O', 'H', 'H', 'O', 'H', 'H', ...])
    geometry : np.ndarray
        Nx3 array of atomic coordinates in Bohr
    fragments : list
        List of lists, each inner list contains 0-indexed atom indices for a fragment
        (e.g., [[0,1,2], [3,4,5], ...] for two 3-atom fragments)
    molecular_charge : float, optional
        Total molecular charge (default: 0.0)
    molecular_multiplicity : int, optional
        Molecular spin multiplicity (default: 1)
    fix_com : bool, optional
        Fix center of mass (default: True)
    fix_orientation : bool, optional
        Fix molecular orientation (default: True)
    fix_symmetry : str, optional
        Symmetry to impose (default: "c1")

    Notes
    -----
    This class exists to work around QCElemental's scalability limitations with large
    fragment counts (>32). QCManyBody doesn't actually need a pre-validated full Molecule
    object - it only needs raw data to construct small fragment molecules on-demand for
    QC calculations.
    """
    symbols: list
    geometry: np.ndarray
    fragments: list
    molecular_charge: float = 0.0
    molecular_multiplicity: int = 1
    fix_com: bool = True
    fix_orientation: bool = True
    fix_symmetry: str = "c1"

    @property
    def nfragments(self) -> int:
        """Number of fragments in the system."""
        return len(self.fragments)

    @property
    def natoms(self) -> int:
        """Total number of atoms."""
        return len(self.symbols)

    def get_fragment_atom_indices(self, fragment_indices: Sequence[int]) -> list:
        """Get 0-indexed atom indices for a subset of fragments.

        Parameters
        ----------
        fragment_indices : Sequence[int]
            0-indexed fragment indices to extract

        Returns
        -------
        list
            List of 0-indexed atom indices

        Examples
        --------
        >>> raw_mol = RawMoleculeData(
        ...     symbols=['O', 'H', 'H', 'O', 'H', 'H'],
        ...     geometry=np.array([[0,0,0], [1,0,0], [0,1,0], [3,0,0], [4,0,0], [3,1,0]]),
        ...     fragments=[[0,1,2], [3,4,5]]
        ... )
        >>> raw_mol.get_fragment_atom_indices([0])
        [0, 1, 2]
        >>> raw_mol.get_fragment_atom_indices([0, 1])
        [0, 1, 2, 3, 4, 5]
        """
        atoms = []
        for frag_idx in fragment_indices:
            atoms.extend(self.fragments[frag_idx])
        return atoms


def construct_fragment_molecule(
    raw_mol: RawMoleculeData,
    real_fragments: Tuple[int, ...],  # 1-indexed fragment IDs
    basis_fragments: Tuple[int, ...],  # 1-indexed fragment IDs
) -> Molecule:
    """Construct a small Molecule object on-demand for QC calculations.

    This function creates Molecule objects with a SINGLE fragment definition to avoid
    QCElemental's expensive validation overhead that scales poorly with fragment count.
    Each constructed Molecule represents one computational fragment (monomer, dimer, etc.)
    and contains only the atoms needed for that specific calculation.

    Parameters
    ----------
    raw_mol : RawMoleculeData
        Raw molecule data containing symbols, geometry, and fragment definitions
    real_fragments : Tuple[int, ...]
        1-indexed fragment IDs for real (non-ghost) atoms
    basis_fragments : Tuple[int, ...]
        1-indexed fragment IDs for basis set (real + ghost atoms)

    Returns
    -------
    Molecule
        QCElemental Molecule object with single fragment definition, ready for QC calculation

    Notes
    -----
    - Ghost atoms are basis atoms not in real atoms
    - The returned Molecule has only ONE fragment definition (not multiple), avoiding
      QCElemental's O(n²) validation overhead
    - For CP-corrected calculations, basis_fragments includes all supersystem fragments
    - Atom ordering: real atoms first, then ghost atoms

    Examples
    --------
    >>> # Create dimer (1,2) in trimer basis (1,2,3) for CP correction
    >>> raw_mol = RawMoleculeData(symbols=['He']*3, geometry=..., fragments=[[0],[1],[2]])
    >>> mol = construct_fragment_molecule(raw_mol, (1,2), (1,2,3))
    >>> # mol contains 3 atoms: He1, He2 (real), He3 (ghost)
    >>> # Fragment definition: [[0,1]] (single fragment with 2 real atoms)
    """
    # Convert 1-indexed fragment IDs to 0-indexed
    real_frag_indices = [f - 1 for f in real_fragments]
    basis_frag_indices = [f - 1 for f in basis_fragments]

    # Get atom indices (0-indexed)
    real_atom_indices = raw_mol.get_fragment_atom_indices(real_frag_indices)
    basis_atom_indices = raw_mol.get_fragment_atom_indices(basis_frag_indices)

    # Ghost atoms are in basis but not in real
    ghost_atom_indices = [idx for idx in basis_atom_indices if idx not in real_atom_indices]

    # Combined atom list: real atoms first, then ghosts
    all_atom_indices = real_atom_indices + ghost_atom_indices

    # Extract data for these atoms
    symbols = [raw_mol.symbols[i] for i in all_atom_indices]
    geometry = raw_mol.geometry[all_atom_indices]

    # Create real/ghost mask
    real_mask = [True] * len(real_atom_indices) + [False] * len(ghost_atom_indices)

    # Construct Molecule with SINGLE fragment definition
    # This is key: one fragment avoids QCElemental's expensive multi-fragment validation
    # Fragment must include ALL atoms (real + ghost), with `real` mask distinguishing them
    mol = Molecule(
        symbols=symbols,
        geometry=geometry,
        fragments=[[i for i in range(len(all_atom_indices))]],  # Single fragment with all atoms!
        real=real_mask,
        molecular_charge=raw_mol.molecular_charge,
        molecular_multiplicity=raw_mol.molecular_multiplicity,
        fix_com=raw_mol.fix_com,
        fix_orientation=raw_mol.fix_orientation,
        fix_symmetry=raw_mol.fix_symmetry,
    )

    return mol


__all__ = ["ManyBodyCalculator", "ManyBodyCore", "RawMoleculeData", "construct_fragment_molecule"]


class ManyBodyCore:
    def __init__(
        self,
        molecule: Union[Molecule, RawMoleculeData],
        bsse_type: Sequence[BsseEnum],
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        *,
        return_total_data: bool,
        supersystem_ie_only: bool,
        embedding_charges: Mapping[int, Sequence[float]],
        hmbe_spec: Optional[HMBESpecification] = None,
    ):
        self.hmbe_spec = hmbe_spec
        self.embedding_charges = embedding_charges
        if self.embedding_charges:
            if not bool(os.environ.get("QCMANYBODY_EMBEDDING_CHARGES", False)):  # obscure until further validation
                raise ValueError(
                    f"Embedding charges for EE-MBE are still in testing. Set environment variable QCMANYBODY_EMBEDDING_CHARGES=1 to use at your own risk."
                )

        # Support both Molecule (legacy, fails with 32+ fragments) and RawMoleculeData (new, scalable)
        if isinstance(molecule, RawMoleculeData):
            # New path: raw molecule data, no validation overhead
            self.raw_molecule = molecule
            self.molecule = None  # Will be None in this mode
            self._use_raw_molecule = True
            self.nfragments = molecule.nfragments
        elif isinstance(molecule, dict):
            # Legacy path: construct Molecule from dict
            mol = Molecule(**molecule)
            self.molecule = mol
            self.raw_molecule = None
            self._use_raw_molecule = False
            self.nfragments = len(mol.fragments)
        elif isinstance(molecule, Molecule):
            # Legacy path: copy Molecule
            mol = molecule.copy()
            self.molecule = mol
            self.raw_molecule = None
            self._use_raw_molecule = False
            self.nfragments = len(mol.fragments)
        else:
            raise ValueError(
                f"Molecule input type of {type(molecule)} not understood. "
                f"Expected Molecule, RawMoleculeData, or dict."
            )

        self.bsse_type = [BsseEnum(x) for x in bsse_type]
        self.return_total_data = return_total_data
        self.supersystem_ie_only = supersystem_ie_only

        # Cache fragment range set for embedding charges to avoid recreating on every iteration
        # For 5,000 HMBE terms with 64 fragments, this saves ~320KB + eliminates 5k redundant set creations
        if self.embedding_charges:
            self._fragment_range_set = set(range(1, self.nfragments + 1))

        # Cache fragment range tuple for HMBE and supersystem calculations
        # Saves ~3.5KB and eliminates 1-2 redundant tuple creations
        self._fragment_range_tuple = tuple(range(1, self.nfragments + 1))

        self.levels = levels

        # Levels without supersystem
        self.levels_no_ss = {int(k): v for k, v in levels.items() if k != "supersystem"}

        # Just a set of all the modelchems
        self.mc_levels = set(self.levels.values())

        self.max_nbody = max(self.levels_no_ss.keys())

        if len(self.bsse_type) == 0:
            raise ValueError("No BSSE correction specified")

        if BsseEnum.vmfc in self.bsse_type and len(set(self.levels.values())) == 1:
            # For single-modelchem VMFC, NOCP & sometimes CP are produced for free
            if BsseEnum.nocp not in self.bsse_type:
                self.bsse_type.append(BsseEnum.nocp)
            if BsseEnum.cp not in self.bsse_type and self.max_nbody == self.nfragments:
                self.bsse_type.append(BsseEnum.cp)

        self.return_bsse_type = self.bsse_type[0]

        ###############################
        # Build nbodies_per_mc_level
        # TODO - use Lori's code
        # TODO - dict to list of lists to handle non-contiguous levels
        # TODO multilevel and supersystem_ie_only=T not allowed together
        # TODO supersystem in levels is not to be trusted -- nfrag only and skips levels
        max_level = max(self.levels_no_ss.keys())

        if set(range(1, max_level + 1)) != set(self.levels_no_ss.keys()):
            raise ValueError(f"Levels must be contiguous from 1 to {max_level}")

        self.nbodies_per_mc_level: Dict[str, list] = {mc_level: [] for mc_level in self.mc_levels}
        for k, v in self.levels_no_ss.items():
            self.nbodies_per_mc_level[v].append(k)

        # order nbodies_per_mc_level keys (modelchems) by the lowest n-body level covered; any
        #   supersystem key (replaced below) is at the end. Order nbodies within each modelchem.
        #   Reset mc_levels to match.
        self.nbodies_per_mc_level = {
            k: sorted(v)
            for (k, v) in sorted(self.nbodies_per_mc_level.items(), key=lambda item: sorted(item[1] or [1000])[0])
        }
        assert self.mc_levels == set(self.nbodies_per_mc_level.keys())  # remove after some downstream testing
        self.mc_levels = list(self.nbodies_per_mc_level.keys())

        for mc, nbs in self.nbodies_per_mc_level.items():
            if nbs and ((nbs[-1] - nbs[0]) != len(nbs) - 1):
                raise ValueError(
                    f"QCManyBody: N-Body levels must be contiguous within a model chemistry spec ({mc}: {nbs}). Use an alternate spec name to accomplish this input."
                )
                # TODO - test and reenable if appropriate. my guess is that noncontig nb is fine on the core computing side,
                #   but trouble for computer and nbodies_per_mc_level inverting and indexing. Safer to deflect for now since input tweak allows the calc.

        # Supersystem is always at the end
        if "supersystem" in levels:
            ss_mc = levels["supersystem"]
            self.nbodies_per_mc_level[ss_mc].append("supersystem")

        # To be built on the fly
        self.mc_compute_dict = None
        self._dependency_graph = None

        if self.nfragments == 1:
            # Usually we try to "pass-through" edge cases, so a single-fragment mol would return 0 or ordinary energy,
            #   depending on rtd=T/F. But it seems more likely that user just forgot the fragments field, so we don't
            #   want to start a full energy on monsterMol. Reconsider handling in future.
            raise ValueError("""QCManyBody: Molecule fragmentation has not been specified through `fragments` field.""")

        # Validate fragment contiguity and build size/slices dictionaries
        # Different paths for Molecule vs RawMoleculeData
        if self._use_raw_molecule:
            # RawMoleculeData path: validate fragments are contiguous
            all_atoms = []
            for frag in self.raw_molecule.fragments:
                all_atoms.extend(frag)
            if not np.array_equal(sorted(all_atoms), np.arange(len(self.raw_molecule.symbols))):
                raise ValueError("""QCManyBody: non-contiguous fragments could be implemented but aren't at present""")

            # Build size and slices dictionaries. Assumes fragments are contiguous
            self.fragment_size_dict = {}
            self.fragment_slice_dict = {}
            iat = 0
            for ifr in range(1, self.nfragments + 1):
                nat = len(self.raw_molecule.fragments[ifr - 1])
                self.fragment_size_dict[ifr] = nat
                self.fragment_slice_dict[ifr] = slice(iat, iat + nat)
                iat += nat
        else:
            # Legacy Molecule path
            if not np.array_equal(np.concatenate(self.molecule.fragments), np.arange(len(self.molecule.symbols))):
                raise ValueError("""QCManyBody: non-contiguous fragments could be implemented but aren't at present""")

            # Build size and slices dictionaries. Assumes fragments are contiguous
            self.fragment_size_dict = {}
            self.fragment_slice_dict = {}
            iat = 0
            for ifr in range(1, self.nfragments + 1):
                nat = len(self.molecule.fragments[ifr - 1])
                self.fragment_size_dict[ifr] = nat
                self.fragment_slice_dict[ifr] = slice(iat, iat + nat)
                iat += nat

    @property
    def has_supersystem(self) -> bool:
        return "supersystem" in self.levels

    @property
    def compute_map(self) -> Dict[str, Dict[str, Dict[int, Set["FragBasIndex"]]]]:
        if self.mc_compute_dict is not None:
            return self.mc_compute_dict

        # Build the compute lists
        self.mc_compute_dict = {}

        for mc in self.mc_levels:
            nbodies = self.nbodies_per_mc_level[mc]

            # Apply HMBE filtering if specified
            if self.hmbe_spec is not None:
                # Decide which enumeration mode to use
                enumeration_mode = self.hmbe_spec.enumeration_mode
                if enumeration_mode == "auto":
                    # Auto mode: use direct for >=30 fragments, filter otherwise
                    enumeration_mode = "direct" if self.nfragments >= 30 else "filter"

                if enumeration_mode == "direct":
                    # Validate that Schengen is not enabled with direct mode
                    if self.hmbe_spec.schengen and self.hmbe_spec.schengen.enabled:
                        raise ValueError(
                            "Schengen terms with direct enumeration not yet supported. "
                            "Direct enumeration avoids building the full MBE list (~680k terms for 64 fragments) "
                            "but Schengen candidate selection currently requires the full list. "
                            "Use enumeration_mode='filter' for Schengen calculations, "
                            "or disable Schengen for large systems (>= 30 fragments)."
                        )

                    # Direct enumeration: generate only HMBE terms, skip building full MBE list
                    # This is critical for large systems to avoid memory explosion
                    from qcmanybody.hmbe_enumerate import enumerate_hmbe_terms

                    # Get HMBE fragment tuples (only ~5k terms for 64 fragments)
                    hmbe_frag_tuples = enumerate_hmbe_terms(self.hmbe_spec)

                    # Construct (frag, bas) pairs directly from HMBE fragment tuples
                    # without generating all ~680k MBE combinations
                    fragment_range = self._fragment_range_tuple
                    filtered_dict = {
                        "all": {},
                        "cp": {},
                        "nocp": {},
                        "vmfc_compute": {},
                        "vmfc_levels": {},
                    }

                    # Build compute lists for each BSSE type
                    for frag_tuple in hmbe_frag_tuples:
                        nbody = len(frag_tuple)

                        # CP: Use supersystem basis (all fragments)
                        if BsseEnum.cp in self.bsse_type:
                            if nbody not in filtered_dict["cp"]:
                                filtered_dict["cp"][nbody] = set()
                            filtered_dict["cp"][nbody].add((frag_tuple, fragment_range))

                        # NOCP: Use natural basis (only the fragments in the tuple)
                        if BsseEnum.nocp in self.bsse_type:
                            if nbody not in filtered_dict["nocp"]:
                                filtered_dict["nocp"][nbody] = set()
                            filtered_dict["nocp"][nbody].add((frag_tuple, frag_tuple))

                        # VMFC: Like CP but for each combination
                        if BsseEnum.vmfc in self.bsse_type:
                            if nbody not in filtered_dict["vmfc_compute"]:
                                filtered_dict["vmfc_compute"][nbody] = set()
                                filtered_dict["vmfc_levels"][nbody] = set()
                            filtered_dict["vmfc_compute"][nbody].add((frag_tuple, frag_tuple))
                            filtered_dict["vmfc_levels"][nbody].add((frag_tuple, frag_tuple))

                    # Add monomers in monomer basis if total data requested
                    if self.return_total_data and 1 in nbodies:
                        if 1 not in filtered_dict["nocp"]:
                            filtered_dict["nocp"][1] = set()
                        # Only add monomers that are in HMBE terms
                        monomer_frags = {frag for term in hmbe_frag_tuples for frag in term}
                        for ifr in monomer_frags:
                            filtered_dict["nocp"][1].add(((ifr,), (ifr,)))

                    # Combine into "all"
                    all_nbodies = set()
                    for subdict in [filtered_dict["cp"], filtered_dict["nocp"], filtered_dict["vmfc_compute"]]:
                        all_nbodies.update(subdict.keys())
                    
                    for nbody in all_nbodies:
                        filtered_dict["all"][nbody] = set()
                        for subdict_key in ["cp", "nocp", "vmfc_compute"]:
                            if nbody in filtered_dict[subdict_key]:
                                filtered_dict["all"][nbody].update(filtered_dict[subdict_key][nbody])

                    # Add Schengen terms if enabled
                    # Note: Schengen requires full MBE list for candidates, so we build it only if needed
                    if self.hmbe_spec.schengen and self.hmbe_spec.schengen.enabled:
                        logger.warning(
                            "Schengen terms with direct enumeration requires building full MBE list. "
                            "This may consume significant memory for large systems. "
                            "Consider using filter mode or disabling Schengen for large systems."
                        )
                        
                        # Build full MBE list only for Schengen candidate selection
                        base_compute_dict = build_nbody_compute_list(
                            self.bsse_type,
                            self.nfragments,
                            nbodies,
                            self.return_total_data,
                            self.supersystem_ie_only,
                            self.max_nbody,
                        )

                        from qcmanybody.hmbe_filter import (
                            get_schengen_candidates,
                            select_schengen_terms,
                        )

                        for bsse_key in ["cp", "nocp", "vmfc_compute"]:
                            if bsse_key not in ["all", "vmfc_levels"] and filtered_dict[bsse_key]:
                                # Get candidates (terms in MBE but not in base HMBE)
                                candidates = get_schengen_candidates(
                                    base_compute_dict[bsse_key], filtered_dict[bsse_key], self.hmbe_spec
                                )

                                # Select top fraction by distance metric + required sub-clusters
                                schengen_terms, required_subs = select_schengen_terms(
                                    candidates, self.molecule, self.hmbe_spec
                                )

                                # Add BOTH Schengen terms AND their sub-clusters
                                all_terms_to_add = schengen_terms | required_subs

                                for frag_tuple in all_terms_to_add:
                                    nbody = len(frag_tuple)
                                    if nbody not in filtered_dict[bsse_key]:
                                        filtered_dict[bsse_key][nbody] = set()

                                    # Construct (frag, bas) based on BSSE type
                                    if bsse_key == "cp":
                                        filtered_dict[bsse_key][nbody].add((frag_tuple, fragment_range))
                                    elif bsse_key == "nocp":
                                        filtered_dict[bsse_key][nbody].add((frag_tuple, frag_tuple))
                                    elif bsse_key == "vmfc_compute":
                                        filtered_dict[bsse_key][nbody].add((frag_tuple, frag_tuple))

                                logger.info(
                                    f"[{mc}] Added {len(schengen_terms)} Schengen terms + "
                                    f"{len(required_subs)} sub-clusters to {bsse_key}"
                                )

                        # Rebuild "all" after adding Schengen
                        filtered_dict["all"] = {}
                        for nbody in all_nbodies:
                            filtered_dict["all"][nbody] = set()
                            for subdict_key in ["cp", "nocp", "vmfc_compute"]:
                                if nbody in filtered_dict[subdict_key]:
                                    filtered_dict["all"][nbody].update(filtered_dict[subdict_key][nbody])

                    self.mc_compute_dict[mc] = filtered_dict

                else:
                    # Filter mode: generate all MBE terms, then filter
                    base_compute_dict = build_nbody_compute_list(
                        self.bsse_type,
                        self.nfragments,
                        nbodies,
                        self.return_total_data,
                        self.supersystem_ie_only,
                        self.max_nbody,
                    )

                    from qcmanybody.hmbe_filter import (
                        filter_compute_list,
                        get_schengen_candidates,
                        select_schengen_terms,
                    )

                    # Filter to base HMBE terms
                    filtered_dict = {}
                    for bsse_key, compute_list in base_compute_dict.items():
                        filtered_dict[bsse_key] = filter_compute_list(compute_list, self.hmbe_spec)

                    # Add Schengen terms if enabled
                    if self.hmbe_spec.schengen and self.hmbe_spec.schengen.enabled:
                        for bsse_key in filtered_dict:
                            # Get candidates (terms in MBE but not in base HMBE)
                            candidates = get_schengen_candidates(
                                base_compute_dict[bsse_key], filtered_dict[bsse_key], self.hmbe_spec
                            )

                            # Select top fraction by distance metric + required sub-clusters
                            schengen_terms, required_subs = select_schengen_terms(
                                candidates, self.molecule, self.hmbe_spec
                            )

                            # Add BOTH Schengen terms AND their sub-clusters for completeness
                            all_terms_to_add = schengen_terms | required_subs

                            for frag_tuple in all_terms_to_add:
                                nbody = len(frag_tuple)
                                if nbody not in filtered_dict[bsse_key]:
                                    filtered_dict[bsse_key][nbody] = set()

                                # Find matching (frag, bas) pairs from base MBE
                                for frag, bas in base_compute_dict[bsse_key].get(nbody, set()):
                                    if frag == frag_tuple:
                                        filtered_dict[bsse_key][nbody].add((frag, bas))

                            logger.info(
                                f"[{mc}] Added {len(schengen_terms)} Schengen terms + "
                                f"{len(required_subs)} sub-clusters to {bsse_key}"
                            )

                    self.mc_compute_dict[mc] = filtered_dict
            else:
                # Standard MBE (no filtering)
                base_compute_dict = build_nbody_compute_list(
                    self.bsse_type,
                    self.nfragments,
                    nbodies,
                    self.return_total_data,
                    self.supersystem_ie_only,
                    self.max_nbody,
                )
                self.mc_compute_dict[mc] = base_compute_dict

        return self.mc_compute_dict

    @property
    def dependency_graph(self) -> NBodyDependencyGraph:
        """Get the N-body dependency graph for level-ordered iteration.

        Returns
        -------
        NBodyDependencyGraph
            Dependency graph instance for level-by-level fragment iteration
        """
        if self._dependency_graph is None:
            self._dependency_graph = NBodyDependencyGraph(self.compute_map)
        return self._dependency_graph

    def get_hmbe_statistics(self) -> Optional[Dict[str, Any]]:
        """Get statistics about HMBE vs MBE term counts.

        Returns
        -------
        Optional[Dict[str, Any]]
            Dictionary with term counts, reduction factors, and HMBE metadata.
            Returns None if HMBE is not enabled.

            Keys include:
            - mbe_term_counts: Dict[str, int] - term counts for full MBE by modelchem
            - hmbe_term_counts: Dict[str, int] - term counts for HMBE by modelchem
            - reduction_factors: Dict[str, float] - reduction factor (MBE/HMBE) by modelchem
            - truncation_orders: Tuple[int, ...] - HMBE truncation orders
            - num_tiers: int - number of hierarchical tiers
            - schengen_enabled: bool - whether Schengen terms are included
        """
        if self.hmbe_spec is None:
            return None

        # Determine which enumeration mode was/will be used
        enumeration_mode = self.hmbe_spec.enumeration_mode
        actual_mode = enumeration_mode
        if enumeration_mode == "auto":
            actual_mode = "direct" if self.nfragments >= 30 else "filter"

        # Count base MBE terms
        # For large systems in direct mode, avoid building the full MBE list (memory explosion!)
        # Instead, calculate counts using combinatorial formulas
        mbe_counts = {}
        if actual_mode == "direct" and self.nfragments >= 30:
            # Use combinatorial formula: sum of C(n,k) for k in nbodies
            from math import comb
            
            for mc in self.mc_levels:
                nbodies = self.nbodies_per_mc_level[mc]
                total_count = 0
                for nbody in nbodies:
                    if nbody == "supersystem":
                        total_count += 1
                    else:
                        # C(nfragments, nbody) combinations
                        total_count += comb(self.nfragments, nbody)
                mbe_counts[mc] = total_count
        else:
            # For small systems or filter mode, actually build the list
            for mc in self.mc_levels:
                nbodies = self.nbodies_per_mc_level[mc]
                base_dict = build_nbody_compute_list(
                    self.bsse_type,
                    self.nfragments,
                    nbodies,
                    self.return_total_data,
                    self.supersystem_ie_only,
                    self.max_nbody,
                )
                # Count all terms across all BSSE types
                all_terms = set()
                for bsse_dict in base_dict.values():
                    for terms in bsse_dict.values():
                        all_terms.update(terms)
                mbe_counts[mc] = len(all_terms)

        # Count HMBE terms (from compute_map which has filtering applied)
        hmbe_counts = {}
        for mc, compute_dict in self.compute_map.items():
            all_terms = set()
            for bsse_dict in compute_dict.values():
                for terms in bsse_dict.values():
                    all_terms.update(terms)
            hmbe_counts[mc] = len(all_terms)

        return {
            "mbe_term_counts": mbe_counts,
            "hmbe_term_counts": hmbe_counts,
            "reduction_factors": {
                mc: mbe_counts[mc] / hmbe_counts[mc] if hmbe_counts[mc] > 0 else 0.0
                for mc in self.mc_levels
            },
            "truncation_orders": self.hmbe_spec.truncation_orders,
            "num_tiers": self.hmbe_spec.num_tiers,
            "schengen_enabled": (
                self.hmbe_spec.schengen.enabled if self.hmbe_spec.schengen else False
            ),
            "enumeration_mode": enumeration_mode,
            "actual_enumeration_mode": actual_mode,
        }

    def format_calc_plan(self, sset: str = "all") -> Tuple[str, Dict[str, Dict[int, int]]]:
        """Formulate per-modelchem and per-body job count data and summary text.

        Parameters
        ----------
        sset
            Among {"all", "nocp", "cp", "vmfc_compute"}, which data structure to return.

        Returns
        -------
        info
            A text summary with per- model chemistry and per- n-body-level job counts.
            ```
            Model chemistry "c4-ccsd" (§A):         22
                 Number of 1-body computations:     16 (nocp: 0, cp: 0, vmfc_compute: 16)
                 Number of 2-body computations:      6 (nocp: 0, cp: 0, vmfc_compute: 6)

            Model chemistry "c4-mp2" (§B):          28
                 Number of 1-body computations:     12 (nocp: 0, cp: 0, vmfc_compute: 12)
                 Number of 2-body computations:     12 (nocp: 0, cp: 0, vmfc_compute: 12)
                 Number of 3-body computations:      4 (nocp: 0, cp: 0, vmfc_compute: 4)
            ```
        Dict[str, Dict[int, int]]
            Data structure with outer key mc-label, inner key 1-indexed n-body, and value job count.
        """
        # Rearrange compute_list from key nb having values (species) to compute all of that nb
        #   to key nb having values counting that nb.
        compute_list_count = {}
        for mc, compute_dict in self.compute_map.items():
            compute_list_count[mc] = {}
            for sub in compute_dict:  # all, nocp, cp, vmfc
                all_calcs = set().union(*compute_dict[sub].values())
                compute_list_count[mc][sub] = Counter([len(frag) for (frag, _) in all_calcs])

        mc_labels = modelchem_labels(self.nbodies_per_mc_level, presorted=True)
        full_to_ordinal_mc_lbl = {v[0]: v[1] for v in mc_labels.values()}
        info = []
        for mc, counter in compute_list_count.items():
            all_counter = counter["all"]
            mcheader = f'    Model chemistry "{mc}" ({full_to_ordinal_mc_lbl[mc]}):'
            info.append(f"{mcheader:38} {sum(all_counter.values()):6}")
            for nb, count in sorted(all_counter.items()):
                other_counts = [f"{sub}: {counter[sub][nb]}" for sub in ["nocp", "cp", "vmfc_compute"]]
                info.append(f"        Number of {nb}-body computations: {count:6} ({', '.join(other_counts)})")
            info.append("")
        info = "\n".join(info)

        logger.info(info)
        return info, {mc: dsset[sset] for mc, dsset in compute_list_count.items()}

    def resize_gradient(self, grad: np.ndarray, bas: Tuple[int, ...], *, reverse: bool = False) -> np.ndarray:
        return resize_gradient(grad, bas, self.fragment_size_dict, self.fragment_slice_dict, reverse=reverse)

    def resize_hessian(self, hess: np.ndarray, bas: Tuple[int, ...], *, reverse: bool = False) -> np.ndarray:
        return resize_hessian(hess, bas, self.fragment_size_dict, self.fragment_slice_dict, reverse=reverse)

    def iterate_molecules(self) -> Iterable[Tuple[str, str, Molecule]]:
        """Iterate over all the molecules needed for the computation.

        Yields model chemistry, label, and molecule.
        """

        done_molecules = set()

        if self._use_raw_molecule:
            # New path: construct molecules on-demand from RawMoleculeData
            for mc, compute_dict in self.compute_map.items():
                for compute_list in compute_dict["all"].values():
                    for real_atoms, basis_atoms in compute_list:
                        label = labeler(mc, real_atoms, basis_atoms)
                        if label in done_molecules:
                            continue

                        # Construct small Molecule on-demand (single fragment definition, no validation overhead)
                        mol = construct_fragment_molecule(self.raw_molecule, real_atoms, basis_atoms)

                        if self.embedding_charges:
                            # Use cached fragment range set instead of creating new range every iteration
                            embedding_frags = list(self._fragment_range_set - set(basis_atoms))
                            charges = []
                            for ifr in embedding_frags:
                                # Get fragment geometry directly from raw data
                                frag_atoms = self.raw_molecule.get_fragment_atom_indices([ifr - 1])
                                positions = self.raw_molecule.geometry[frag_atoms].tolist()
                                charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[ifr])])
                            mol.extras["embedding_charges"] = charges

                        done_molecules.add(label)
                        yield mc, label, mol

        else:
            # Legacy path: use Molecule.get_fragment() (fails with 32+ fragments)
            for mc, compute_dict in self.compute_map.items():
                # TODO - this is a bit of a hack. Lots of duplication when reaching higher nbody
                for compute_list in compute_dict["all"].values():
                    for real_atoms, basis_atoms in compute_list:
                        label = labeler(mc, real_atoms, basis_atoms)
                        if label in done_molecules:
                            continue

                        ghost_atoms = list(set(basis_atoms) - set(real_atoms))

                        # Shift to zero-indexing
                        real_atoms_0 = [x - 1 for x in real_atoms]
                        ghost_atoms_0 = [x - 1 for x in ghost_atoms]
                        mol = self.molecule.get_fragment(real_atoms_0, ghost_atoms_0, orient=False, group_fragments=False)
                        updates = {"fix_com": True, "fix_orientation": True}
                        if self.molecule.fix_symmetry == "c1":
                            # symmetry in the parent usually irrelevant to symmetry in the fragments, but
                            #   if parent symmetry is cancelled, catch that and pass it along
                            updates["fix_symmetry"] = "c1"
                        mol = mol.copy(update=updates)

                        if self.embedding_charges:
                            # Use cached fragment range set instead of creating new range every iteration
                            embedding_frags = list(self._fragment_range_set - set(basis_atoms))
                            charges = []
                            for ifr in embedding_frags:
                                positions = self.molecule.get_fragment(ifr - 1).geometry.tolist()
                                charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[ifr])])
                            mol.extras["embedding_charges"] = charges

                        done_molecules.add(label)
                        yield mc, label, mol

    def iterate_molecules_by_level(self) -> Iterable[Tuple[int, str, str, Molecule]]:
        """Iterate over molecules needed for computation, grouped by N-body dependency level.

        This method provides level-by-level iteration that respects mathematical dependencies:
        monomers (level 1) → dimers (level 2) → trimers (level 3) → etc.

        This enables safe parallel execution within each level while respecting dependencies
        between levels.

        Yields
        ------
        Tuple[int, str, str, Molecule]
            Tuple of (level, model_chemistry, label, molecule) where:
            - level: N-body dependency level (1, 2, 3, ...)
            - model_chemistry: String identifying the quantum chemistry method
            - label: Fragment label in JSON format
            - molecule: QCElemental Molecule object for this fragment

        Notes
        -----
        This method preserves the exact same molecule set as iterate_molecules(),
        only changing the ordering to respect N-body dependencies. All fragments
        yielded by iterate_molecules() will be yielded by this method exactly once.

        Examples
        --------
        >>> mbc = ManyBodyCore(molecule, ["cp"], {1: "hf", 2: "mp2"})
        >>> for level, mc, label, mol in mbc.iterate_molecules_by_level():
        ...     print(f"Level {level}: {mc} calculation for {label}")
        Level 1: hf calculation for ["hf", [1], [1]]
        Level 1: hf calculation for ["hf", [2], [2]]
        Level 2: mp2 calculation for ["mp2", [1, 2], [1, 2]]
        """
        done_molecules: Set[str] = set()

        # Performance optimization: pre-compute common values
        has_embedding = bool(self.embedding_charges)

        if self._use_raw_molecule:
            # New path: construct molecules on-demand from RawMoleculeData
            # No upfront validation, enables 64+ fragment systems
            for level, fragments_at_level in self.dependency_graph.iterate_molecules_by_level():
                for fragment_dep in fragments_at_level:
                    mc = fragment_dep.mc
                    label = fragment_dep.label

                    if label in done_molecules:
                        continue

                    real_atoms = fragment_dep.real_atoms
                    basis_atoms = fragment_dep.basis_atoms

                    # Construct small Molecule on-demand (single fragment definition, no validation overhead)
                    mol = construct_fragment_molecule(self.raw_molecule, real_atoms, basis_atoms)

                    if has_embedding:
                        # Use cached fragment range set instead of creating new range every iteration
                        embedding_frags = list(self._fragment_range_set - set(basis_atoms))
                        charges = []
                        for ifr in embedding_frags:
                            # Get fragment geometry directly from raw data
                            frag_atoms = self.raw_molecule.get_fragment_atom_indices([ifr - 1])
                            positions = self.raw_molecule.geometry[frag_atoms].tolist()
                            charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[ifr])])
                        mol.extras["embedding_charges"] = charges

                    done_molecules.add(label)
                    yield level, mc, label, mol

        else:
            # Legacy path: use Molecule.get_fragment() (fails with 32+ fragments)
            fix_c1_symmetry = self.molecule.fix_symmetry == "c1"

            # Base updates dict - avoid creating it for each molecule
            base_updates = {"fix_com": True, "fix_orientation": True}
            if fix_c1_symmetry:
                base_updates["fix_symmetry"] = "c1"

            # Use dependency graph to get level-ordered iteration
            for level, fragments_at_level in self.dependency_graph.iterate_molecules_by_level():
                for fragment_dep in fragments_at_level:
                    mc = fragment_dep.mc  # model chemistry
                    label = fragment_dep.label  # fragment label

                    if label in done_molecules:
                        continue

                    # Performance optimization: use cached real_atoms and basis_atoms from FragmentDependency
                    real_atoms = fragment_dep.real_atoms
                    basis_atoms = fragment_dep.basis_atoms

                    # Performance optimization: use set difference for ghost atoms
                    ghost_atoms = list(set(basis_atoms) - set(real_atoms))

                    # Shift to zero-indexing
                    real_atoms_0 = [x - 1 for x in real_atoms]
                    ghost_atoms_0 = [x - 1 for x in ghost_atoms]
                    mol = self.molecule.get_fragment(real_atoms_0, ghost_atoms_0, orient=False, group_fragments=False)

                    # Use pre-computed updates
                    mol = mol.copy(update=base_updates)

                    if has_embedding:
                        # Use cached fragment range set instead of creating new range every iteration
                        embedding_frags = list(self._fragment_range_set - set(basis_atoms))
                        charges = []
                        for ifr in embedding_frags:
                            positions = self.molecule.get_fragment(ifr - 1).geometry.tolist()
                            charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[ifr])])
                        mol.extras["embedding_charges"] = charges

                    done_molecules.add(label)
                    yield level, mc, label, mol

    def _assemble_nbody_components(
        self,
        property_label: str,
        component_results: Dict[str, Union[float, np.ndarray]],
    ) -> Dict[str, Any]:
        """Assembles N-body components for a single derivative level and a single model chemistry level
        into interaction quantities according to requested BSSE treatment(s).
        """

        # When HMBE filtering is active, the MBE job list is no longer the full combinatorial
        # set. The standard inclusion–exclusion coefficients derived for complete MBE coverage
        # therefore overcount missing higher-order terms and can explode (e.g., huge 3-body
        # contributions). For HMBE we instead do a per-fragment Möbius inversion over the
        # actually enumerated clusters, summing contributions bottom-up.
        if self.hmbe_spec is not None:
            return self._assemble_nbody_components_hmbe(property_label, component_results)

        # which level are we assembling?
        delabeled = [delabeler(k) for k in component_results.keys()]
        mc_level_labels = {x[0] for x in delabeled}

        if len(mc_level_labels) != 1:
            raise RuntimeError(f"Multiple model chemistries passed into _assemble_nbody_components: {mc_level_labels}")

        mc_level = mc_level_labels.pop()
        if mc_level not in self.mc_levels:
            raise RuntimeError(f"Model chemistry {mc_level} not found in {self.mc_levels}")

        # get the range of nbodies and the required calculations for this level
        bsse_type = self.bsse_type
        return_bsse_type = self.return_bsse_type
        nbodies = self.nbodies_per_mc_level[mc_level]
        if "supersystem" in nbodies:
            nbodies = list(range(1, self.max_nbody + 1))
            bsse_type = [BsseEnum.nocp]
            return_bsse_type = BsseEnum.nocp

        max_nbody = max(nbodies)
        compute_dict = self.compute_map[mc_level]

        if not all_same_shape(component_results.values()):
            raise ValueError("All values in data dictionary must have the same shape.")

        # Use first data value to determine shape
        first_key = next(iter(component_results.keys()))
        property_shape = find_shape(component_results[first_key])

        # Accumulation dictionaries
        # * {bsse_type}_by_level is filled by sum_cluster_data to contain for NOCP
        #   & CP the summed total energies (or other property) of each nb-body. That is:
        #   * NOCP: {1: 1b@1b,    2: 2b@2b,      ..., max_nbody: max_nbody-b@max_nbody-b} and
        #   * CP:   {1: 1b@nfr-b, 2: 2b@nfr-b,   ..., max_nbody: max_nbody-b@nfr-b}.
        #   VMFC bookkeeping is different. For key 1 it contains the summed 1b total energies.
        #   But for higher keys, it contains each nb-body (non-additive) contribution to the energy.
        #   * VMFC: {1: 1b@1b,    2: 2b contrib, ..., max_nbody: max_nbody-b contrib}
        cp_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_by_level = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}

        # * {bsse_type}_body_dict is usually filled with total energies (or other property).
        #   Multiple model chemistry levels may be involved.
        #   Generally, all consecutive keys between 1 and max_nbody will be present in the body_dict,
        #   but if supersystem_ie_only=T, only 1b and nfr-b are present, or if "supersystem" in levels, ???
        #   * TOT: {1: 1b@1b, 2: 2b tot prop with bsse_type treatment, ..., max_nbody: max_nbody-b tot prop with bsse_type treatment}
        #   If 1b@1b (monomers in monomer basis) aren't available, which can happen when return_total_data=F
        #   and 1b@1b aren't otherwise needed, body_dict contains interaction energies (or other property).
        #   * IE: {1: shaped_zero, 2: 2b interaction prop using bsse_type, ..., max_nbody: max_nbody-b interaction prop using bsse_type}
        #   For both TOT and IE cases, body_dict values are cummulative, not additive. For TOT, total,
        #   interaction, and contribution data in ManyBodyResultProperties can be computed in
        #   collect_vars. For IE, interaction and contribution data can be computed.
        cp_body_dict = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        nocp_body_dict = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}
        vmfc_body_dict = {n: shaped_zero(property_shape) for n in range(1, nbodies[-1] + 1)}

        # Sum up all of the levels
        # * compute_dict[bt][nb] holds all the computations needed to compute nb
        #   *not* all the nb-level computations, so build the latter
        cp_compute_list = {nb: set() for nb in range(1, nbodies[-1] + 1)}
        nocp_compute_list = {nb: set() for nb in range(1, nbodies[-1] + 1)}

        for nb in nbodies:
            # HMBE filtering may omit certain nbody levels from compute_dict
            if nb in compute_dict["cp"]:
                for v in compute_dict["cp"][nb]:
                    if len(v[1]) != 1:
                        cp_compute_list[len(v[0])].add(v)
            if nb in compute_dict["nocp"]:
                for w in compute_dict["nocp"][nb]:
                    nocp_compute_list[len(w[0])].add(w)

        for nb in range(1, nbodies[-1] + 1):
            cp_by_level[nb] = sum_cluster_data(component_results, cp_compute_list[nb], mc_level)
            nocp_by_level[nb] = sum_cluster_data(component_results, nocp_compute_list[nb], mc_level)
            if nb in compute_dict["vmfc_levels"]:
                vmfc_by_level[nb] = sum_cluster_data(
                    component_results, compute_dict["vmfc_levels"][nb], mc_level, vmfc=True, nb=nb
                )

        # Extract data for monomers in monomer basis for CP total data
        if 1 in nbodies and 1 in compute_dict["nocp"]:
            monomers_in_monomer_basis = [v for v in compute_dict["nocp"][1] if len(v[1]) == 1]
            monomer_sum = sum_cluster_data(component_results, set(monomers_in_monomer_basis), mc_level)
        else:
            monomer_sum = shaped_zero(property_shape)

        # Compute cp
        if BsseEnum.cp in bsse_type:
            for nb in range(1, nbodies[-1] + 1):
                if nb == self.nfragments:
                    cp_body_dict[nb] = cp_by_level[nb] - bsse
                    continue

                for k in range(1, nb + 1):
                    take_nk = math.comb(self.nfragments - k - 1, nb - k)
                    sign = (-1) ** (nb - k)
                    cp_body_dict[nb] += take_nk * sign * cp_by_level[k]

                if nb == 1:
                    bsse = cp_body_dict[nb] - monomer_sum
                    cp_body_dict[nb] = copy_value(monomer_sum)
                else:
                    cp_body_dict[nb] -= bsse

        # Compute nocp
        if BsseEnum.nocp in bsse_type:
            for nb in range(1, nbodies[-1] + 1):
                if nb == self.nfragments:
                    nocp_body_dict[nb] = nocp_by_level[nb]
                    continue

                for k in range(1, nb + 1):
                    take_nk = math.comb(self.nfragments - k - 1, nb - k)
                    sign = (-1) ** (nb - k)
                    nocp_body_dict[nb] += take_nk * sign * nocp_by_level[k]

        # Compute vmfc
        if BsseEnum.vmfc in bsse_type:
            for nb in nbodies:
                for k in range(1, nb + 1):
                    vmfc_body_dict[nb] += vmfc_by_level[k]

        # Collect specific and generalized returns
        results = {
            f"cp_{property_label}_body_dict": cp_body_dict,
            f"nocp_{property_label}_body_dict": nocp_body_dict,
            f"vmfc_{property_label}_body_dict": vmfc_body_dict,
        }

        # Overall return body dict & value for this property
        results[f"{property_label}_body_dict"] = results[f"{return_bsse_type.value}_{property_label}_body_dict"]
        results[f"ret_{property_label}"] = copy_value(results[f"{property_label}_body_dict"][max_nbody])

        if not self.return_total_data:
            results[f"ret_{property_label}"] -= results[f"{property_label}_body_dict"][1]

        return results

    def _validate_hmbe_completeness(
        self,
        fragment_tuples: Set[Tuple[int, ...]],
        mc_level: str,
        bsse_type: str,
    ) -> None:
        """Validate that cluster set satisfies completeness requirement for Möbius inversion.

        Möbius inversion is mathematically exact ONLY when every cluster has all its
        proper sub-clusters present. This function enforces that requirement.

        Parameters
        ----------
        fragment_tuples : Set[Tuple[int, ...]]
            All fragment tuples with computed results
        mc_level : str
            Model chemistry level (for error reporting)
        bsse_type : str
            BSSE type ("nocp", "cp") for error reporting

        Raises
        ------
        RuntimeError
            If any cluster is missing required sub-clusters. This indicates a bug
            in HMBE filtering, Schengen sub-cluster addition, or direct enumeration.

        Notes
        -----
        This check should NEVER fail in production. If it does, it's a QCManyBody bug.
        """
        from itertools import combinations

        missing = []

        for frag_tuple in fragment_tuples:
            n = len(frag_tuple)
            if n == 1:
                continue  # Monomers have no sub-clusters

            # Check all proper sub-clusters (size 1 to n-1)
            for sub_size in range(1, n):
                for sub_cluster in combinations(frag_tuple, sub_size):
                    if sub_cluster not in fragment_tuples:
                        missing.append((frag_tuple, sub_cluster))

        if missing:
            error_msg = (
                f"HMBE completeness violation in {mc_level}/{bsse_type}:\n"
                f"Found {len(missing)} clusters missing required sub-clusters.\n"
                f"This violates the mathematical requirement for Möbius inversion.\n"
                f"\n"
                f"First 5 violations:\n"
            )
            for cluster, missing_sub in missing[:5]:
                error_msg += f"  Cluster {cluster} missing sub-cluster {missing_sub}\n"

            error_msg += (
                f"\n"
                f"This is a bug in QCManyBody. Please report this issue.\n"
                f"Likely causes: Schengen sub-cluster addition failed, or HMBE filtering bug."
            )

            raise RuntimeError(error_msg)

        logger.debug(
            f"HMBE completeness validated for {mc_level}/{bsse_type}: "
            f"{len(fragment_tuples)} clusters OK"
        )

    def _assemble_nbody_components_hmbe(
        self,
        property_label: str,
        component_results: Dict[str, Union[float, np.ndarray]],
    ) -> Dict[str, Any]:
        """Assemble N-body results when HMBE truncation is active.

        Uses Möbius inversion over the enumerated fragment tuples. This is
        mathematically EXACT (not approximate) provided the completeness
        requirement is satisfied.

        Mathematical Foundation
        -----------------------
        For cluster C with computed energy E_total[C], its n-body contribution is:

            E_contrib[C] = E_total[C] - Σ_{S ⊂ C, S ≠ C} E_contrib[S]

        This is the Möbius function on the subset lattice. CRITICAL REQUIREMENT:
        All proper subsets S of C must be present (computed) for this to be exact.

        Completeness Enforcement
        ------------------------
        The completeness requirement is ENFORCED by:
        1. Base HMBE filtering: Guarantees completeness by construction
           (if C passes filter, all subsets pass filter due to hierarchical structure)
        2. Schengen sub-cluster addition: When Schengen adds C, automatically adds
           all required sub-clusters
        3. Validation: _validate_hmbe_completeness() catches any bugs

        Why Standard MBE Assembly Fails for HMBE
        -----------------------------------------
        Standard MBE uses inclusion-exclusion coefficients math.comb(N-k-1, n-k)
        that assume ALL C(N,k) combinations exist. When HMBE filters terms, these
        coefficients overcount missing terms, causing exploding higher-order energies.

        Parameters
        ----------
        property_label : str
            Property to assemble ("energy", "gradient", "hessian")
        component_results : Dict[str, Union[float, np.ndarray]]
            Computed results keyed by labeler(mc, frag, bas)

        Returns
        -------
        Dict[str, Any]
            Assembled n-body results with body_dict and ret_<property>
        """

        # Determine model chemistry for this batch of results
        delabeled = [delabeler(k) for k in component_results.keys()]
        mc_level_labels = {x[0] for x in delabeled}

        if len(mc_level_labels) != 1:
            raise RuntimeError(f"Multiple model chemistries passed into _assemble_nbody_components_hmbe: {mc_level_labels}")

        mc_level = mc_level_labels.pop()
        if mc_level not in self.mc_levels:
            raise RuntimeError(f"Model chemistry {mc_level} not found in {self.mc_levels}")

        compute_dict = self.compute_map[mc_level]

        # All values must share shape
        if not all_same_shape(component_results.values()):
            raise ValueError("All values in data dictionary must have the same shape.")

        first_key = next(iter(component_results.keys()))
        property_shape = find_shape(component_results[first_key])

        max_nbody = max([n for n in self.nbodies_per_mc_level[mc_level] if n != "supersystem"], default=0)

        def zero():
            return shaped_zero(property_shape)

        results: Dict[str, Any] = {}

        # Currently only nocp/cp are supported for HMBE. VMFC would require bespoke handling.
        if BsseEnum.vmfc in self.bsse_type:
            raise NotImplementedError("VMFC with HMBE is not supported yet.")

        for bt in self.bsse_type:
            bt_key = bt.value
            contribution_by_level = {n: zero() for n in range(1, max_nbody + 1)}

            # Map fragment tuple -> total property for this bsse type
            energy_by_frag: Dict[Tuple[int, ...], Any] = {}
            subdict_key = bt_key if bt_key in compute_dict else None
            if subdict_key is None:
                continue

            for nbody, terms in compute_dict[subdict_key].items():
                for frag, bas in terms:
                    lbl = labeler(mc_level, frag, bas)
                    if lbl in component_results:
                        energy_by_frag[frag] = component_results[lbl]

            # VALIDATE COMPLETENESS (catches bugs in sub-cluster addition)
            fragment_tuples = set(energy_by_frag.keys())
            self._validate_hmbe_completeness(fragment_tuples, mc_level, bt_key)

            # Explicit Möbius inversion over available fragments
            contrib_by_frag: Dict[Tuple[int, ...], Any] = {}
            for n in range(1, max_nbody + 1):
                # Sort for determinism
                frags_of_size = sorted([f for f in energy_by_frag if len(f) == n])
                for frag in frags_of_size:
                    subtotal = copy_value(energy_by_frag[frag])

                    # Subtract contributions from all proper subsets we've already computed
                    for sub_frag, sub_contrib in contrib_by_frag.items():
                        if len(sub_frag) >= n:
                            continue
                        if set(sub_frag).issubset(frag):
                            subtotal -= sub_contrib

                    contrib_by_frag[frag] = subtotal
                    contribution_by_level[n] += subtotal

            # Build cumulative body dict (through-n values)
            body_dict: Dict[int, Any] = {}
            running = zero()
            for n in range(1, max_nbody + 1):
                running += contribution_by_level.get(n, zero())
                body_dict[n] = copy_value(running)

            results[f"{bt_key}_{property_label}_body_dict"] = body_dict

        # Overall return body dict & value for this property
        results[f"{property_label}_body_dict"] = results[
            f"{self.return_bsse_type.value}_{property_label}_body_dict"
        ]
        results[f"ret_{property_label}"] = copy_value(results[f"{property_label}_body_dict"][max_nbody])

        if not self.return_total_data:
            results[f"ret_{property_label}"] -= results[f"{property_label}_body_dict"][1]

        return results

    def _analyze(
        self,
        property_label: str,
        property_results: Dict[str, Union[float, np.ndarray]],  # Label to results
    ):
        # Initialize with zeros
        if not all_same_shape(property_results.values()):
            raise ValueError("All values in data dictionary must have the same shape.")

        # Use first data value to determine shape
        first_key = next(iter(property_results.keys()))
        property_shape = find_shape(property_results[first_key])

        property_result = shaped_zero(property_shape)
        property_body_dict = {bt.value: {} for bt in self.bsse_type}
        property_body_contribution = {bt.value: {} for bt in self.bsse_type}

        # results per model chemistry
        mc_results = {}
        species_results = {}

        # sort by nbody level, ignore supersystem
        sorted_nbodies = [(k, v) for k, v in self.nbodies_per_mc_level.items() if v != ["supersystem"]]
        sorted_nbodies = sorted(sorted_nbodies, reverse=True, key=lambda x: x[1])
        for mc_label, nbody_list in sorted_nbodies:
            # filter to only one model chemistry
            filtered_results = {k: v for k, v in property_results.items() if delabeler(k)[0] == mc_label}

            if not filtered_results:
                if nbody_list == [1]:
                    # Note A.2: Note A.1 holds, but for the special case of CP-only
                    #   and rtd=False and multilevel with a separate level for
                    #   1-body, the skipped tasks run afoul of sanity checks, so
                    #   we'll add a dummy result.
                    filtered_results = {labeler(mc_label, [1000], [1000]): shaped_zero(property_shape)}
                else:
                    raise RuntimeError(f"No data found for model chemistry {mc_label}")

            nb_component_results = self._assemble_nbody_components(property_label, filtered_results)
            mc_results[mc_label] = nb_component_results

            for n in nbody_list[::-1]:
                property_bsse_dict = {bt.value: shaped_zero(property_shape) for bt in self.bsse_type}

                for m in range(n - 1, n + 1):
                    if m == 0:
                        continue

                    # Subtract the (n-1)-body contribution from the n-body contribution to get the n-body effect
                    sign = (-1) ** (1 - m // n)
                    for bt in self.bsse_type:
                        property_bsse_dict[bt.value] += (
                            sign * mc_results[mc_label][f"{bt.value}_{property_label}_body_dict"][m]
                        )

                property_result += property_bsse_dict[self.return_bsse_type]
                for bt in self.bsse_type:
                    property_body_contribution[bt.value][n] = property_bsse_dict[bt.value]

        if self.has_supersystem:
            # Get the MC label for supersystem tasks
            supersystem_mc_level = self.levels["supersystem"]

            # Super system recovers higher order effects at a lower level
            frag_range = self._fragment_range_tuple

            ss_cresults = {k: v for k, v in property_results.items() if delabeler(k)[0] == supersystem_mc_level}
            ss_component_results = self._assemble_nbody_components(property_label, ss_cresults)
            mc_results[supersystem_mc_level] = ss_component_results

            # Compute components at supersystem level of theory
            ss_label = labeler(supersystem_mc_level, frag_range, frag_range)
            supersystem_result = property_results[ss_label]
            property_result += supersystem_result - ss_component_results[f"{property_label}_body_dict"][self.max_nbody]

            for bt in self.bsse_type:
                property_body_contribution[bt][self.nfragments] = (
                    supersystem_result - ss_component_results[f"{property_label}_body_dict"][self.max_nbody]
                )

        for bt in self.bsse_type:
            bstr = bt.value
            for n in property_body_contribution[bstr]:
                property_body_dict[bstr][n] = sum(
                    [
                        property_body_contribution[bstr][i]
                        for i in range(1, n + 1)
                        if i in property_body_contribution[bstr]
                    ]
                )

        if not self.return_total_data:
            # Remove monomer contribution for interaction data
            property_result -= property_body_dict[self.return_bsse_type][1]

        nbody_results = {
            f"ret_{property_label}": property_result,
            f"{property_label}_body_dict": property_body_dict,
            "mc_results": mc_results,
        }
        return nbody_results

    def analyze(
        self,
        component_results: Dict[str, Dict[str, Union[float, np.ndarray]]],
    ):
        """

        Parameters
        ----------
        component_results
            Nested dictionary with results from all individual molecular system
            calculations, including all subsystem/basis combinations, all model
            chemistries, and all properties (e.g., e/g/h).

            For example, the below is the format for a nocp gradient run on a
            helium dimer with 1-body at CCSD and 2-body at MP2. The outer string
            key can be generated with the ``qcmanybody.utils.labeler`` function.
            The inner string key is any property; QCManyBody presently knows how
            to process energy/gradient/Hessian.
            ```
            {'["ccsd", [1], [1]]': {'energy': -2.87, 'gradient': array([[0., 0., 0.]])},
             '["ccsd", [2], [2]]': {'energy': -2.87, 'gradient': array([[0., 0., 0.]])},
             '["mp2", [1], [1]]': {'energy': -2.86, 'gradient': array([[0., 0., 0.]])},
             '["mp2", [2], [2]]': {'energy': -2.86, 'gradient': array([[0., 0., 0.]])},
             '["mp2", [1, 2], [1, 2]]': {'energy': -5.73, 'gradient': array([[ 0., 0., 0.0053], [ 0., 0., -0.0053]])},
            }
            ```
        """

        # Convert AtomicResult objects to property dictionaries if needed
        # This handles both old format (Dict[str, Dict]) and new format (Dict[str, AtomicResult])
        converted_results = {}
        for label, result_data in component_results.items():
            if hasattr(result_data, 'properties'):
                # This is an AtomicResult object - extract properties
                properties_dict = {}

                # Always include the return_result (typically energy)
                if hasattr(result_data, 'return_result') and result_data.return_result is not None:
                    # Convert to float to ensure compatibility with analysis functions
                    properties_dict['energy'] = float(result_data.return_result)

                # Extract all non-None properties from the properties object
                if result_data.properties is not None:
                    props_dict = result_data.properties.dict()
                    for prop_name, prop_value in props_dict.items():
                        if prop_value is not None:
                            # Convert numeric values to appropriate types
                            if isinstance(prop_value, (int, float)):
                                prop_value = float(prop_value)

                            # Map common property names to expected format
                            if prop_name == 'return_energy':
                                properties_dict['energy'] = prop_value
                            elif prop_name == 'return_gradient':
                                properties_dict['gradient'] = prop_value
                            elif prop_name == 'return_hessian':
                                properties_dict['hessian'] = prop_value
                            else:
                                # Include other properties as-is
                                properties_dict[prop_name] = prop_value

                converted_results[label] = properties_dict
            else:
                # Already in the expected dictionary format
                converted_results[label] = result_data

        # Use converted results for the rest of the analysis
        component_results = converted_results

        # All properties that were passed to us
        # * seed with "energy" so free/no-op jobs can process
        available_properties: Set[str] = {"energy"}
        for property_data in component_results.values():
            available_properties.update(property_data.keys())

        # reorganize to component_results_inv[property][label] = 1.23
        component_results_inv = {k: {} for k in available_properties}

        for cluster_label, property_data in component_results.items():
            for property_label, property_value in property_data.items():
                component_results_inv[property_label][cluster_label] = property_value

        # Remove any missing data
        component_results_inv = {k: v for k, v in component_results_inv.items() if v}
        if not component_results_inv:
            # Note B: Rarely, "no results" is expected, like for CP-only,
            #   rtd=False, and max_nbody=1. We'll add a dummy entry so
            #   processing can continue.
            component_results_inv["energy"] = {'["dummy", [1000], [1000]]': 0.0}

        # Actually analyze
        is_embedded = bool(self.embedding_charges)
        component_properties = defaultdict(dict)
        all_results = {}
        nbody_dict = {}
        stdout = ""
        #        all_results["energy_body_dict"] = {"cp": {1: 0.0}}

        for property_label, property_results in component_results_inv.items():
            # Expand gradient and hessian
            if property_label == "gradient":
                property_results = {k: self.resize_gradient(v, delabeler(k)[2]) for k, v in property_results.items()}
            if property_label == "hessian":
                property_results = {k: self.resize_hessian(v, delabeler(k)[2]) for k, v in property_results.items()}

            r = self._analyze(property_label, property_results)
            for k, v in property_results.items():
                component_properties[k]["calcinfo_natom"] = len(self.molecule.symbols)
                component_properties[k][f"return_{property_label}"] = v
            all_results.update(r)

        for bt in self.bsse_type:
            stdout += print_nbody_energy(
                all_results["energy_body_dict"][bt],
                f"{bt.formal()} ({bt.abbr()})",
                self.nfragments,
                modelchem_labels(self.nbodies_per_mc_level, presorted=True),
                is_embedded,
                self.supersystem_ie_only,
                self.max_nbody if self.has_supersystem else None,
            )

        for property_label in available_properties:
            for bt in self.bsse_type:
                nbody_dict.update(
                    collect_vars(
                        bt,
                        property_label,
                        all_results[f"{property_label}_body_dict"][bt],
                        self.max_nbody,
                        is_embedded,
                        self.supersystem_ie_only,
                        self.has_supersystem,
                    )
                )

        all_results["results"] = nbody_dict
        all_results["component_properties"] = component_properties

        # Make dictionary with "1cp", "2cp", etc
        ebd = all_results["energy_body_dict"]
        all_results["energy_body_dict"] = {str(k) + bt: v for bt in ebd for k, v in ebd[bt].items()}
        all_results["stdout"] = stdout

        return all_results


class ManyBodyCalculator(ManyBodyCore):
    # retire after a grace period
    def __init__(
        self,
        molecule: Molecule,
        bsse_type: Sequence[BsseEnum],
        levels: Mapping[Union[int, Literal["supersystem"]], str],
        return_total_data: bool,
        supersystem_ie_only: bool,
        embedding_charges: Mapping[int, Sequence[float]],
    ):
        super().__init__(
            molecule,
            bsse_type,
            levels,
            return_total_data=return_total_data,
            supersystem_ie_only=supersystem_ie_only,
            embedding_charges=embedding_charges,
        )
