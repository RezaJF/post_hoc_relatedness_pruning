"""
Pruning engine module for relatedness pruning workflow.

This module implements the core pruning algorithm that selects which individuals
to retain from closely related pairs based on trait values and prioritization rules.
"""

import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
from lab_measurement_analyzer import LabMeasurementAnalyzer


class PruningEngine:
    """Handles the core pruning algorithm for related individuals."""

    def __init__(self, logger: Optional[logging.Logger] = None,
                 lab_analyzer: Optional[LabMeasurementAnalyzer] = None):
        """
        Initialize PruningEngine.

        Args:
            logger: Optional logger instance. If None, creates a basic logger.
            lab_analyzer: Optional lab measurement analyzer instance.
        """
        self.logger = logger or self._setup_logger()
        self.lab_analyzer = lab_analyzer

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

    def prune_related_individuals(self, connected_components: List[Set[Tuple[str, str]]],
                                individual_lookup: Dict[Tuple[str, str], Dict],
                                measurement_lookup: Optional[Dict[str, Dict]] = None) -> Tuple[Set[Tuple[str, str]], Set[Tuple[str, str]]]:
        """
        Prune related individuals from connected components.

        For each connected component, retain one individual based on prioritization rules:
        If measurement_lookup is provided (measurement-based):
        1. Higher number of measurements > lower number of measurements
        2. Lower variance (tie-breaker) > higher variance
        3. Random selection (final tie-breaker)

        Otherwise (trait-based fallback):
        1. Non-missing trait value > missing trait value
        2. Larger absolute trait value > smaller absolute trait value
        3. Lexicographically smaller IID (tie-breaker)

        Args:
            connected_components: List of sets of related individuals
            individual_lookup: Dictionary mapping individuals to their data
            measurement_lookup: Optional dictionary with measurement statistics

        Returns:
            Tuple of (retained_individuals, removed_individuals)
        """
        self.logger.info(f"Starting pruning process for {len(connected_components)} connected components...")

        retained_individuals = set()
        removed_individuals = set()

        pruning_stats = {
            'components_processed': 0,
            'individuals_retained': 0,
            'individuals_removed': 0,
            'decisions_by_rule': {
                'non_missing_vs_missing': 0,
                'larger_absolute_value': 0,
                'lexicographic_tie_breaker': 0
            }
        }

        for i, component in enumerate(connected_components):
            self.logger.debug(f"Processing component {i+1}/{len(connected_components)} with {len(component)} individuals")

            # Filter component to only include individuals present in trait data
            component_in_trait_data = {ind for ind in component if ind in individual_lookup}

            if len(component_in_trait_data) == 0:
                self.logger.warning(f"Component {i+1} has no individuals in trait data, skipping")
                continue

            if len(component_in_trait_data) == 1:
                # Single individual, automatically retain
                retained_individual = list(component_in_trait_data)[0]
                retained_individuals.add(retained_individual)
                pruning_stats['individuals_retained'] += 1
                self.logger.debug(f"Single individual in component, retaining: {retained_individual}")
            else:
                # Multiple individuals, apply pruning rules
                retained_individual = self._select_individual_to_retain(
                    component_in_trait_data, individual_lookup, measurement_lookup
                )

                retained_individuals.add(retained_individual)
                removed_from_component = component_in_trait_data - {retained_individual}
                removed_individuals.update(removed_from_component)

                pruning_stats['individuals_retained'] += 1
                pruning_stats['individuals_removed'] += len(removed_from_component)

                self.logger.debug(f"Retained: {retained_individual}, Removed: {list(removed_from_component)}")

            pruning_stats['components_processed'] += 1

        self.logger.info(f"Pruning completed:")
        self.logger.info(f"  Components processed: {pruning_stats['components_processed']}")
        self.logger.info(f"  Individuals retained: {pruning_stats['individuals_retained']}")
        self.logger.info(f"  Individuals removed: {pruning_stats['individuals_removed']}")

        return retained_individuals, removed_individuals

    def _select_individual_to_retain(self, component: Set[Tuple[str, str]],
                                   individual_lookup: Dict[Tuple[str, str], Dict],
                                   measurement_lookup: Optional[Dict[str, Dict]] = None) -> Tuple[str, str]:
        """
        Select which individual to retain from a component based on prioritization rules.

        Args:
            component: Set of related individuals
            individual_lookup: Dictionary mapping individuals to their data
            measurement_lookup: Optional dictionary with measurement statistics

        Returns:
            Tuple of (FID, IID) for the individual to retain
        """
        # Convert to list for easier processing
        individuals = list(component)

        # Apply prioritization rules
        best_individual = individuals[0]
        best_decision_rule = "initial"

        for individual in individuals[1:]:
            decision_rule = self._compare_individuals(
                best_individual, individual, individual_lookup, measurement_lookup
            )

            if decision_rule == "individual2_better":
                best_individual = individual
                best_decision_rule = decision_rule
            elif decision_rule == "individual1_better":
                best_decision_rule = decision_rule
            # If decision_rule is "tie", keep the current best_individual

        self.logger.debug(f"Selected {best_individual} using rule: {best_decision_rule}")
        return best_individual

    def _compare_individuals(self, individual1: Tuple[str, str], individual2: Tuple[str, str],
                           individual_lookup: Dict[Tuple[str, str], Dict],
                           measurement_lookup: Optional[Dict[str, Dict]] = None) -> str:
        """
        Compare two individuals and determine which is better to retain.

        Uses measurement-based prioritization exclusively.

        Args:
            individual1: First individual (FID, IID)
            individual2: Second individual (FID, IID)
            individual_lookup: Dictionary mapping individuals to their data
            measurement_lookup: Dictionary with measurement statistics (required)

        Returns:
            String indicating which individual is better or if it's a tie
        """
        # Measurement-based comparison is mandatory
        if not self.lab_analyzer or measurement_lookup is None:
            raise ValueError("Measurement-based comparison is required but lab_analyzer or measurement_lookup is not available")

        return self.lab_analyzer.compare_individuals_by_measurements(
            individual1, individual2, measurement_lookup
        )

    def create_pruned_trait_data(self, trait_df: pd.DataFrame,
                               retained_individuals: Set[Tuple[str, str]],
                               rank_normalize: bool = False,
                               inverse_rank_normalize: bool = True) -> pd.DataFrame:
        """
        Create pruned trait data with only retained individuals.

        Args:
            trait_df: Original trait DataFrame
            retained_individuals: Set of individuals to retain
            rank_normalize: Whether to apply rank normalization to trait values
            inverse_rank_normalize: Whether to apply inverse rank normalization (default True)

        Returns:
            DataFrame with only retained individuals
        """
        self.logger.info("Creating pruned trait data...")

        # Create mask for retained individuals
        retained_mask = trait_df.apply(
            lambda row: (row['FID'], row['IID']) in retained_individuals, axis=1
        )

        pruned_df = trait_df[retained_mask].copy()

        # Apply rank normalization if requested
        if rank_normalize:
            self.logger.info("Applying rank normalization to trait values...")
            trait_column = pruned_df.attrs.get('trait_column', pruned_df.columns[2])

            # Get non-missing values for rank normalization
            non_missing_mask = pruned_df[trait_column].notna()
            if non_missing_mask.sum() > 0:
                # Apply rank normalization to non-missing values
                pruned_df.loc[non_missing_mask, trait_column] = pruned_df.loc[non_missing_mask, trait_column].rank(method='average')
                self.logger.info(f"Rank normalization applied to {non_missing_mask.sum()} non-missing values")
            else:
                self.logger.warning("No non-missing values found for rank normalization")

        # Apply inverse rank normalization if requested (default behavior for enhanced pipeline)
        elif inverse_rank_normalize:
            self.logger.info("Applying inverse rank normalization to trait values...")
            trait_column = pruned_df.attrs.get('trait_column', pruned_df.columns[2])

            # Get non-missing values for inverse rank normalization
            non_missing_mask = pruned_df[trait_column].notna()
            if non_missing_mask.sum() > 0:
                from scipy import stats
                # Convert to ranks
                ranks = pruned_df.loc[non_missing_mask, trait_column].rank(method='average')
                # Scale ranks to (0, 1) interval
                n = len(ranks)
                scaled_ranks = (ranks - 0.5) / n
                # Apply inverse normal transformation
                pruned_df.loc[non_missing_mask, trait_column] = stats.norm.ppf(scaled_ranks)
                self.logger.info(f"Inverse rank normalization applied to {non_missing_mask.sum()} non-missing values")
            else:
                self.logger.warning("No non-missing values found for inverse rank normalization")

        self.logger.info(f"Pruned trait data: {len(pruned_df)} individuals retained out of {len(trait_df)} total")

        return pruned_df

    def generate_pruning_summary(self, original_trait_df: pd.DataFrame,
                               pruned_trait_df: pd.DataFrame,
                               removed_individuals: Set[Tuple[str, str]],
                               individual_lookup: Dict[Tuple[str, str], Dict]) -> Dict:
        """
        Generate summary statistics for the pruning process.

        Args:
            original_trait_df: Original trait DataFrame
            pruned_trait_df: Pruned trait DataFrame
            removed_individuals: Set of removed individuals
            individual_lookup: Dictionary mapping individuals to their data

        Returns:
            Dictionary with pruning summary statistics
        """
        trait_column = original_trait_df.attrs.get('trait_column', original_trait_df.columns[2])

        # Calculate statistics for original data
        original_total = len(original_trait_df)
        original_non_missing = original_trait_df[trait_column].notna().sum()
        original_missing = original_total - original_non_missing

        # Calculate statistics for pruned data
        pruned_total = len(pruned_trait_df)
        pruned_non_missing = pruned_trait_df[trait_column].notna().sum()
        pruned_missing = pruned_total - pruned_non_missing

        # Calculate statistics for removed individuals
        removed_total = len(removed_individuals)
        removed_non_missing = sum(1 for ind in removed_individuals
                                if individual_lookup[ind]['has_trait_value'])
        removed_missing = removed_total - removed_non_missing

        summary = {
            'original_data': {
                'total_individuals': original_total,
                'non_missing_trait_values': original_non_missing,
                'missing_trait_values': original_missing,
                'proportion_non_missing': original_non_missing / original_total if original_total > 0 else 0
            },
            'pruned_data': {
                'total_individuals': pruned_total,
                'non_missing_trait_values': pruned_non_missing,
                'missing_trait_values': pruned_missing,
                'proportion_non_missing': pruned_non_missing / pruned_total if pruned_total > 0 else 0
            },
            'removed_individuals': {
                'total_individuals': removed_total,
                'non_missing_trait_values': removed_non_missing,
                'missing_trait_values': removed_missing,
                'proportion_non_missing': removed_non_missing / removed_total if removed_total > 0 else 0
            },
            'pruning_efficiency': {
                'individuals_removed': removed_total,
                'individuals_retained': pruned_total,
                'removal_rate': removed_total / original_total if original_total > 0 else 0,
                'retention_rate': pruned_total / original_total if original_total > 0 else 0
            }
        }

        # Add trait value statistics if available
        if original_non_missing > 0:
            original_trait_values = original_trait_df[trait_column].dropna()
            summary['original_data']['trait_statistics'] = {
                'mean': float(original_trait_values.mean()),
                'std': float(original_trait_values.std()),
                'min': float(original_trait_values.min()),
                'max': float(original_trait_values.max()),
                'median': float(original_trait_values.median())
            }

        if pruned_non_missing > 0:
            pruned_trait_values = pruned_trait_df[trait_column].dropna()
            summary['pruned_data']['trait_statistics'] = {
                'mean': float(pruned_trait_values.mean()),
                'std': float(pruned_trait_values.std()),
                'min': float(pruned_trait_values.min()),
                'max': float(pruned_trait_values.max()),
                'median': float(pruned_trait_values.median())
            }

        return summary
