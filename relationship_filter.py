"""
Relationship filtering module for relatedness pruning workflow.

This module handles filtering of closely related individuals based on KING
relationship types and building relationship networks for pruning.
"""

import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict, deque
import networkx as nx


class RelationshipFilter:
    """Handles filtering and analysis of closely related individuals."""

    # Define closely related relationship types (up to 2nd degree)
    CLOSELY_RELATED_TYPES = {'Dup/MZ', 'PO', 'FS', '2nd'}

    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize RelationshipFilter.

        Args:
            logger: Optional logger instance. If None, creates a basic logger.
        """
        self.logger = logger or self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Set up basic logger if none provided."""
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger

    def filter_closely_related_pairs(self, kinship_df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter kinship data to include only closely related pairs.

        Args:
            kinship_df: DataFrame with kinship data from KING

        Returns:
            DataFrame with only closely related pairs
        """
        self.logger.info("Filtering closely related pairs...")

        # Filter for closely related types
        closely_related_mask = kinship_df['InfType'].isin(self.CLOSELY_RELATED_TYPES)
        filtered_df = kinship_df[closely_related_mask].copy()

        self.logger.info(f"Found {len(filtered_df)} closely related pairs out of {len(kinship_df)} total pairs")

        # Log distribution of relationship types
        if len(filtered_df) > 0:
            relationship_counts = filtered_df['InfType'].value_counts()
            self.logger.info("Distribution of closely related pairs:")
            for rel_type, count in relationship_counts.items():
                self.logger.info(f"  {rel_type}: {count} pairs")

        return filtered_df

    def build_relationship_graph(self, closely_related_df: pd.DataFrame) -> nx.Graph:
        """
        Build a network graph of closely related individuals.

        Args:
            closely_related_df: DataFrame with closely related pairs

        Returns:
            NetworkX graph where nodes are individuals and edges represent relationships
        """
        self.logger.info("Building relationship graph...")

        G = nx.Graph()

        # Add edges for each closely related pair
        for _, row in closely_related_df.iterrows():
            individual1 = (row['FID1'], row['ID1'])
            individual2 = (row['FID2'], row['ID2'])

            # Add nodes with attributes
            G.add_node(individual1, fid=row['FID1'], iid=row['ID1'])
            G.add_node(individual2, fid=row['FID2'], iid=row['ID2'])

            # Add edge with relationship type
            G.add_edge(individual1, individual2,
                      relationship_type=row['InfType'],
                      kinship=row.get('Kinship', np.nan))

        self.logger.info(f"Built relationship graph with {G.number_of_nodes()} individuals and {G.number_of_edges()} relationships")

        return G

    def find_connected_components(self, relationship_graph: nx.Graph) -> List[Set[Tuple[str, str]]]:
        """
        Find connected components in the relationship graph.

        Each connected component represents a group of related individuals
        where at least one individual must be removed.

        Args:
            relationship_graph: NetworkX graph of relationships

        Returns:
            List of sets, where each set contains tuples of (FID, IID) for related individuals
        """
        self.logger.info("Finding connected components...")

        connected_components = list(nx.connected_components(relationship_graph))

        # Filter out single-node components (unrelated individuals)
        multi_individual_components = [comp for comp in connected_components if len(comp) > 1]

        self.logger.info(f"Found {len(connected_components)} total components")
        self.logger.info(f"Found {len(multi_individual_components)} components with multiple individuals")

        # Log component size distribution
        if multi_individual_components:
            component_sizes = [len(comp) for comp in multi_individual_components]
            size_counts = pd.Series(component_sizes).value_counts().sort_index()
            self.logger.info("Component size distribution:")
            for size, count in size_counts.items():
                self.logger.info(f"  Size {size}: {count} components")

        return multi_individual_components

    def get_relationship_summary(self, kinship_df: pd.DataFrame,
                               closely_related_df: pd.DataFrame,
                               connected_components: List[Set[Tuple[str, str]]]) -> Dict:
        """
        Generate summary statistics about relationships.

        Args:
            kinship_df: Original kinship DataFrame
            closely_related_df: Filtered closely related pairs
            connected_components: List of connected components

        Returns:
            Dictionary with summary statistics
        """
        summary = {
            'total_pairs': len(kinship_df),
            'closely_related_pairs': len(closely_related_df),
            'total_individuals_in_kinship': len(set(
                [(row['FID1'], row['ID1']) for _, row in kinship_df.iterrows()] +
                [(row['FID2'], row['ID2']) for _, row in kinship_df.iterrows()]
            )),
            'total_components': len(connected_components),
            'individuals_in_components': sum(len(comp) for comp in connected_components),
            'relationship_type_distribution': closely_related_df['InfType'].value_counts().to_dict() if len(closely_related_df) > 0 else {}
        }

        # Add component size distribution
        if connected_components:
            component_sizes = [len(comp) for comp in connected_components]
            summary['component_size_distribution'] = pd.Series(component_sizes).value_counts().sort_index().to_dict()
        else:
            summary['component_size_distribution'] = {}

        return summary

    def validate_individuals_in_components(self, connected_components: List[Set[Tuple[str, str]]],
                                         trait_individuals: Set[Tuple[str, str]]) -> Dict:
        """
        Validate which individuals in components are present in trait data.

        Args:
            connected_components: List of connected components
            trait_individuals: Set of individuals present in trait data

        Returns:
            Dictionary with validation results
        """
        validation_results = {
            'components_with_trait_data': 0,
            'components_without_trait_data': 0,
            'individuals_in_components_with_trait_data': 0,
            'individuals_in_components_without_trait_data': 0,
            'missing_individuals': []
        }

        for component in connected_components:
            component_individuals_in_trait = component.intersection(trait_individuals)

            if len(component_individuals_in_trait) > 0:
                validation_results['components_with_trait_data'] += 1
                validation_results['individuals_in_components_with_trait_data'] += len(component_individuals_in_trait)
            else:
                validation_results['components_without_trait_data'] += 1

            # Track missing individuals
            missing_in_trait = component - trait_individuals
            validation_results['individuals_in_components_without_trait_data'] += len(missing_in_trait)
            validation_results['missing_individuals'].extend(list(missing_in_trait))

        self.logger.info(f"Components with trait data: {validation_results['components_with_trait_data']}")
        self.logger.info(f"Components without trait data: {validation_results['components_without_trait_data']}")
        self.logger.info(f"Individuals in components with trait data: {validation_results['individuals_in_components_with_trait_data']}")

        return validation_results

    def create_individual_lookup(self, trait_df: pd.DataFrame) -> Dict[Tuple[str, str], Dict]:
        """
        Create a lookup dictionary for individuals in trait data.

        Args:
            trait_df: DataFrame with trait data

        Returns:
            Dictionary mapping (FID, IID) to individual data
        """
        lookup = {}
        trait_column = trait_df.attrs.get('trait_column', trait_df.columns[2])

        for _, row in trait_df.iterrows():
            individual_id = (row['FID'], row['IID'])
            lookup[individual_id] = {
                'fid': row['FID'],
                'iid': row['IID'],
                'trait_value': row[trait_column],
                'has_trait_value': not pd.isna(row[trait_column])
            }

        return lookup


