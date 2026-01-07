"""
Hierarchical Many-Body Expansion (HMBE) module for QCManyBody.

This module provides support for HMBE calculations, which restrict the many-body
expansion based on hierarchical fragment structure. HMBE dramatically reduces
computational cost for systems with natural hierarchical organization.

Key Components:
--------------
- HierarchicalFragment: Represents fragments with recursive nesting
- FragmentHierarchy: Manages complete N-tier hierarchical structure
- HMBEFragmentBuilder: Converts hierarchical to flat fragment representation
- HMBETupleGenerator: Generates valid N-mer tuples respecting hierarchy constraints
- HMBEPreprocessor: Orchestrates conversion from HMBE to standard ManyBodyInput
"""

from .fragment_builder import HMBEFragmentBuilder
from .hierarchy import FragmentHierarchy, HierarchicalFragment
from .preprocessor import HMBEPreprocessor
from .tuple_generator import HMBETupleGenerator

__all__ = [
    "HierarchicalFragment",
    "FragmentHierarchy",
    "HMBEFragmentBuilder",
    "HMBETupleGenerator",
    "HMBEPreprocessor",
]
