"""
HMBE preprocessor for converting hierarchical input to standard ManyBodyInput.

This module orchestrates the conversion from HMBE hierarchical fragment structure
to flat elementary fragments suitable for standard ManyBodyComputer execution.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

from qcmanybody.models.v1 import ManyBodyInput

from .fragment_builder import HMBEFragmentBuilder
from .hierarchy import FragmentHierarchy
from .tuple_generator import HMBETupleGenerator

logger = logging.getLogger(__name__)

__all__ = ["HMBEPreprocessor"]


class HMBEPreprocessor:
    """Converts HMBE input to standard ManyBodyInput with restricted compute list.

    This is the main integration point between HMBE and QCManyBody. It:
    1. Detects HMBE mode from input
    2. Parses hierarchical structure
    3. Builds flat elementary fragment molecule
    4. Generates restricted N-mer tuples based on HMBE constraints
    5. Stores HMBE metadata in ManyBodyInput.extras for builder.py to use

    The preprocessor operates transparently - users provide HMBE input and receive
    standard ManyBodyInput that works with existing ManyBodyComputer infrastructure.
    """

    @staticmethod
    def is_hmbe_input(input_data: ManyBodyInput) -> bool:
        """Check if input specifies HMBE calculation.

        HMBE mode is indicated by presence of 'hmbe_hierarchy' in keywords.

        Parameters
        ----------
        input_data : ManyBodyInput
            Input to check

        Returns
        -------
        bool
            True if HMBE mode is enabled
        """
        if input_data.specification is None:
            return False
        if input_data.specification.keywords is None:
            return False

        # Use getattr for safer attribute access (handles both old and new model versions)
        hmbe_hierarchy = getattr(input_data.specification.keywords, 'hmbe_hierarchy', None)
        return hmbe_hierarchy is not None

    @staticmethod
    def preprocess_hmbe_input(
        input_data: ManyBodyInput,
    ) -> ManyBodyInput:
        """Convert HMBE input to standard ManyBodyInput with flat fragments.

        Strategy:
        1. Parse hierarchical structure from input.specification.keywords.hmbe_hierarchy
        2. Build flat elementary fragment molecule via HMBEFragmentBuilder
        3. Generate restricted N-mer tuple list via HMBETupleGenerator
        4. Store HMBE metadata in input.extras for downstream use
        5. Return modified ManyBodyInput with flat molecule and HMBE metadata

        Parameters
        ----------
        input_data : ManyBodyInput
            Input with HMBE hierarchy specified in keywords

        Returns
        -------
        ManyBodyInput
            Modified input with:
            - molecule: flat elementary fragments
            - extras['hmbe']: metadata including valid tuples and hierarchy info

        Raises
        ------
        ValueError
            If HMBE hierarchy is malformed or incompatible with other settings

        Notes
        -----
        The returned ManyBodyInput can be used with standard ManyBodyComputer.
        The builder.py module will read the restricted tuples from extras['hmbe']
        to generate only valid HMBE N-mer calculations.
        """
        if not HMBEPreprocessor.is_hmbe_input(input_data):
            raise ValueError(
                "Input does not specify HMBE calculation. "
                "Set specification.keywords.hmbe_hierarchy to enable HMBE mode."
            )

        logger.info("Starting HMBE preprocessing")

        # Step 1: Parse hierarchical structure
        hierarchy_data = input_data.specification.keywords.hmbe_hierarchy
        hierarchy = FragmentHierarchy.from_dict(hierarchy_data)

        logger.info(
            f"Parsed HMBE hierarchy: {hierarchy.tiers} tiers, "
            f"{hierarchy.num_primary} primary fragments, "
            f"{hierarchy.num_elementary} elementary fragments"
        )

        # Step 2: Build flat molecule from elementary fragments
        builder = HMBEFragmentBuilder(hierarchy)
        flat_molecule = builder.build_flat_molecule(
            molecule_name=input_data.molecule.name or "hmbe_system"
        )

        logger.info(
            f"Built flat molecule: {len(flat_molecule.symbols)} atoms, "
            f"{len(flat_molecule.fragments)} fragments"
        )

        # Step 3: Generate restricted N-mer tuples
        max_nbody = input_data.specification.keywords.max_nbody
        if max_nbody is None:
            # Default to number of elementary fragments
            max_nbody = hierarchy.num_elementary
            logger.warning(
                f"max_nbody not specified, using num_elementary={max_nbody}"
            )

        generator = HMBETupleGenerator(hierarchy)
        valid_tuples = generator.generate_nmer_tuples(max_nbody)

        # Log tuple counts
        for nbody, tuples in valid_tuples.items():
            logger.info(f"Generated {len(tuples)} valid {nbody}-body tuples")

        # Step 4: Create HMBE metadata for builder.py to use
        hmbe_metadata = {
            "hmbe_mode": True,
            "tiers": hierarchy.tiers,
            "max_primary_per_nmer": hierarchy.max_primary_per_nmer,
            "num_primary": hierarchy.num_primary,
            "num_elementary": hierarchy.num_elementary,
            "fragment_mapping": builder.get_fragment_mapping(),
            "primary_mapping": builder.get_primary_mapping(),
            "valid_tuples": {
                str(nbody): [list(t) for t in tuples]
                for nbody, tuples in valid_tuples.items()
            },
            "hierarchy_data": hierarchy_data,  # Preserve original for reference
        }

        # Step 5: Create modified ManyBodyInput
        # Update molecule
        modified_input = input_data.copy(update={"molecule": flat_molecule})

        # Update extras with HMBE metadata
        existing_extras = modified_input.extras or {}
        updated_extras = {**existing_extras, "hmbe": hmbe_metadata}
        modified_input = modified_input.copy(update={"extras": updated_extras})

        # Remove hmbe_hierarchy from keywords (no longer needed, in extras now)
        modified_keywords = modified_input.specification.keywords.copy(
            update={"hmbe_hierarchy": None}
        )
        modified_spec = modified_input.specification.copy(
            update={"keywords": modified_keywords}
        )
        modified_input = modified_input.copy(update={"specification": modified_spec})

        logger.info("HMBE preprocessing complete")

        return modified_input

    @staticmethod
    def get_hmbe_metadata(input_data: ManyBodyInput) -> Optional[Dict[str, Any]]:
        """Extract HMBE metadata from preprocessed input.

        Parameters
        ----------
        input_data : ManyBodyInput
            Preprocessed input (output of preprocess_hmbe_input)

        Returns
        -------
        Optional[Dict[str, Any]]
            HMBE metadata dictionary if present, None otherwise
        """
        if input_data.extras is None:
            return None
        return input_data.extras.get("hmbe")

    @staticmethod
    def is_preprocessed(input_data: ManyBodyInput) -> bool:
        """Check if input has already been preprocessed for HMBE.

        Parameters
        ----------
        input_data : ManyBodyInput
            Input to check

        Returns
        -------
        bool
            True if input has HMBE metadata in extras
        """
        metadata = HMBEPreprocessor.get_hmbe_metadata(input_data)
        return metadata is not None and metadata.get("hmbe_mode", False)
