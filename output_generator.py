"""
Output generator module for relatedness pruning workflow.

This module handles generating output files, logs, and reports for the
related individual pruning pipeline.
"""

import pandas as pd
import numpy as np
import logging
import json
from typing import Dict, List, Set, Tuple, Optional, Union
from pathlib import Path
from datetime import datetime


class OutputGenerator:
    """Handles generation of output files and reports."""

    def __init__(self, output_dir: Union[str, Path], logger: Optional[logging.Logger] = None):
        """
        Initialize OutputGenerator.

        Args:
            output_dir: Directory to save output files
            logger: Optional logger instance. If None, creates a basic logger.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
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

    def save_pruned_trait_data(self, pruned_df: pd.DataFrame,
                              original_filename: str) -> Path:
        """
        Save pruned trait data to file.

        Args:
            pruned_df: Pruned trait DataFrame
            original_filename: Original trait filename for naming

        Returns:
            Path to the saved pruned file
        """
        # Generate output filename
        original_path = Path(original_filename)
        output_filename = f"{original_path.stem}_related_pruned{original_path.suffix}"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving pruned trait data to: {output_path}")

        # Create a copy to avoid modifying the original
        output_df = pruned_df.copy()

        # Remove any empty rows (rows where all values are NaN, empty strings, or whitespace)
        # First, convert empty strings to NaN
        output_df = output_df.replace('', np.nan)
        # Then remove rows where all values are NaN
        output_df = output_df.dropna(how='all')

        # Ensure proper NA coding - replace any remaining NaN values with 'NA'
        output_df = output_df.fillna('NA')

        # Rename the third column based on the output filename
        if len(output_df.columns) >= 3:
            # Get the column name without the .tsv suffix
            new_column_name = output_filename.replace('.tsv', '')
            # Rename the third column (index 2)
            output_df.columns = list(output_df.columns[:2]) + [new_column_name] + list(output_df.columns[3:])
            self.logger.info(f"Renamed third column to: {new_column_name}")

        # Save the pruned data
        output_df.to_csv(output_path, sep='\t', index=False)

        self.logger.info(f"Pruned trait data saved: {len(output_df)} individuals")

        return output_path

    def save_description_file(self, trait_filename: str, omop_concept_id: str,
                            atc_codes: List[str], months_before: int,
                            unrelated_count: int, total_measurements: int) -> Path:
        """
        Save description file for the pruned trait data.

        Args:
            trait_filename: Original trait filename for naming
            omop_concept_id: OMOP concept ID for the lab test
            atc_codes: List of ATC codes for drug filtering
            months_before: Months before drug purchase
            unrelated_count: Number of unrelated individuals after pruning
            total_measurements: Total number of measurements used

        Returns:
            Path to the saved description file
        """
        # Generate output filename
        original_path = Path(trait_filename)
        output_filename = f"{original_path.stem}_related_pruned_descriptionfile.tsv"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving description file to: {output_path}")

        # Create description content
        atc_codes_str = ", ".join(atc_codes)
        description_content = f"{omop_concept_id} pre-medication (BLUPs, calculated from 1% Winsorised measurements and rank-normalised), OMOP_concept_id:{omop_concept_id}- Longitudinal measurements across {unrelated_count:,} unrelated individuals with >=3 measurements up to {months_before} days before {atc_codes_str} intake or no drug purchase history.\n"

        # Save description file
        with open(output_path, 'w') as f:
            f.write(description_content)

        self.logger.info(f"Description file saved: {output_path}")

        return output_path

    def save_removed_individuals(self, removed_individuals: Set[Tuple[str, str]],
                               trait_filename: str) -> Path:
        """
        Save list of removed individuals to file.

        Args:
            removed_individuals: Set of removed individuals
            trait_filename: Original trait filename for naming

        Returns:
            Path to the saved removed individuals file
        """
        # Generate output filename
        original_path = Path(trait_filename)
        output_filename = f"{original_path.stem}_removed_individuals.txt"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving removed individuals list to: {output_path}")

        # Convert to DataFrame for easier handling
        removed_data = []
        for fid, iid in sorted(removed_individuals):
            removed_data.append({'FID': fid, 'IID': iid})

        removed_df = pd.DataFrame(removed_data)
        removed_df.to_csv(output_path, sep='\t', index=False)

        self.logger.info(f"Removed individuals list saved: {len(removed_individuals)} individuals")

        return output_path

    def save_enhanced_pruning_log(self, trait_filename: str,
                                 relationship_summary: Dict,
                                 pruning_summary: Dict,
                                 measurement_summary: Dict,
                                 missing_in_kinship: List[Tuple[str, str]],
                                 validation_results: Dict) -> Path:
        """
        Save enhanced pruning log with measurement statistics.

        Args:
            trait_filename: Original trait filename for naming
            relationship_summary: Summary of relationship analysis
            pruning_summary: Summary of pruning process
            measurement_summary: Summary of measurement statistics
            missing_in_kinship: List of individuals missing in kinship data
            validation_results: Results of individual validation

        Returns:
            Path to the saved log file
        """
        # Generate output filename
        original_path = Path(trait_filename)
        output_filename = f"{original_path.stem}_enhanced_pruning_log.txt"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving enhanced pruning log to: {output_path}")

        # Create enhanced log content
        log_content = self._generate_enhanced_log_content(
            trait_filename, relationship_summary, pruning_summary, measurement_summary,
            missing_in_kinship, validation_results
        )

        # Save log file
        with open(output_path, 'w') as f:
            f.write(log_content)

        self.logger.info(f"Enhanced pruning log saved: {output_path}")

        return output_path

    def save_enhanced_summary_json(self, trait_filename: str,
                                  relationship_summary: Dict,
                                  pruning_summary: Dict,
                                  measurement_summary: Dict,
                                  missing_in_kinship: List[Tuple[str, str]],
                                  validation_results: Dict) -> Path:
        """
        Save enhanced summary statistics as JSON file.

        Args:
            trait_filename: Original trait filename for naming
            relationship_summary: Summary of relationship analysis
            pruning_summary: Summary of pruning process
            measurement_summary: Summary of measurement statistics
            missing_in_kinship: List of individuals missing in kinship data
            validation_results: Results of individual validation

        Returns:
            Path to the saved JSON file
        """
        # Generate output filename
        original_path = Path(trait_filename)
        output_filename = f"{original_path.stem}_enhanced_pruning_summary.json"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving enhanced summary JSON to: {output_path}")

        # Create enhanced summary dictionary
        summary_data = {
            'timestamp': datetime.now().isoformat(),
            'input_file': trait_filename,
            'relationship_summary': relationship_summary,
            'pruning_summary': pruning_summary,
            'measurement_summary': measurement_summary,
            'validation_results': validation_results,
            'missing_in_kinship': [f"{fid}_{iid}" for fid, iid in missing_in_kinship]
        }

        # Save JSON file
        with open(output_path, 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)

        self.logger.info(f"Enhanced summary JSON saved: {output_path}")

        return output_path

    def create_enhanced_batch_summary(self, results: List[Dict]) -> Path:
        """
        Create an enhanced summary of batch processing results.

        Args:
            results: List of result dictionaries from processing multiple files

        Returns:
            Path to the enhanced batch summary file
        """
        output_filename = "enhanced_batch_pruning_summary.txt"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Creating enhanced batch summary: {output_path}")

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        summary_content = f"""
ENHANCED BATCH RELATEDNESS PRUNING SUMMARY
==========================================

Generated: {timestamp}
Number of files processed: {len(results)}

DETAILED RESULTS
================

"""

        total_original = 0
        total_pruned = 0
        total_removed = 0
        total_measurements = 0

        for i, result in enumerate(results, 1):
            trait_filename = result['trait_filename']
            pruning_summary = result['pruning_summary']
            measurement_summary = result.get('measurement_summary', {})

            total_original += pruning_summary['original_data']['total_individuals']
            total_pruned += pruning_summary['pruned_data']['total_individuals']
            total_removed += pruning_summary['removed_individuals']['total_individuals']
            total_measurements += measurement_summary.get('total_individuals', 0)

            summary_content += f"{i}. {trait_filename}\n"
            summary_content += f"   Original: {pruning_summary['original_data']['total_individuals']:,} individuals\n"
            summary_content += f"   Pruned: {pruning_summary['pruned_data']['total_individuals']:,} individuals\n"
            summary_content += f"   Removed: {pruning_summary['removed_individuals']['total_individuals']:,} individuals\n"
            summary_content += f"   Retention rate: {pruning_summary['pruning_efficiency']['retention_rate']:.3f}\n"
            summary_content += f"   Measurements analyzed: {measurement_summary.get('total_individuals', 0):,}\n"
            summary_content += f"   Mean measurements per individual: {measurement_summary.get('mean_measurements_per_individual', 0):.2f}\n\n"

        summary_content += f"""
OVERALL SUMMARY
===============

Total individuals across all files:
- Original: {total_original:,}
- Pruned: {total_pruned:,}
- Removed: {total_removed:,}
- Overall retention rate: {total_pruned/total_original:.3f}
- Overall removal rate: {total_removed/total_original:.3f}
- Total measurements analyzed: {total_measurements:,}

END OF ENHANCED BATCH SUMMARY
============================
"""

        # Save enhanced batch summary
        with open(output_path, 'w') as f:
            f.write(summary_content)

        self.logger.info(f"Enhanced batch summary saved: {output_path}")

        return output_path

    def _generate_enhanced_log_content(self, trait_filename: str,
                                      relationship_summary: Dict,
                                      pruning_summary: Dict,
                                      measurement_summary: Dict,
                                      missing_in_kinship: List[Tuple[str, str]],
                                      validation_results: Dict) -> str:
        """
        Generate the content for the enhanced pruning log file.

        Args:
            trait_filename: Original trait filename
            relationship_summary: Summary of relationship analysis
            pruning_summary: Summary of pruning process
            measurement_summary: Summary of measurement statistics
            missing_in_kinship: List of individuals missing in kinship data
            validation_results: Results of individual validation

        Returns:
            String content for the enhanced log file
        """
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        log_content = f"""
ENHANCED RELATEDNESS PRUNING LOG
================================

Generated: {timestamp}
Input file: {trait_filename}

SUMMARY STATISTICS
==================

Original Data:
- Total individuals: {pruning_summary['original_data']['total_individuals']:,}
- Non-missing trait values: {pruning_summary['original_data']['non_missing_trait_values']:,}
- Missing trait values: {pruning_summary['original_data']['missing_trait_values']:,}
- Proportion with non-missing values: {pruning_summary['original_data']['proportion_non_missing']:.3f}

Pruned Data:
- Total individuals retained: {pruning_summary['pruned_data']['total_individuals']:,}
- Non-missing trait values: {pruning_summary['pruned_data']['non_missing_trait_values']:,}
- Missing trait values: {pruning_summary['pruned_data']['missing_trait_values']:,}
- Proportion with non-missing values: {pruning_summary['pruned_data']['proportion_non_missing']:.3f}
- Number of individuals with non-missing values in pruned data: {pruning_summary['pruned_data']['non_missing_trait_values']:,}

Removed Individuals:
- Total individuals removed: {pruning_summary['removed_individuals']['total_individuals']:,}
- Non-missing trait values: {pruning_summary['removed_individuals']['non_missing_trait_values']:,}
- Missing trait values: {pruning_summary['removed_individuals']['missing_trait_values']:,}
- Proportion with non-missing values: {pruning_summary['removed_individuals']['proportion_non_missing']:.3f}

Pruning Efficiency:
- Removal rate: {pruning_summary['pruning_efficiency']['removal_rate']:.3f}
- Retention rate: {pruning_summary['pruning_efficiency']['retention_rate']:.3f}

MEASUREMENT STATISTICS
======================

Total individuals with measurements: {measurement_summary.get('total_individuals', 0):,}
Mean measurements per individual: {measurement_summary.get('mean_measurements_per_individual', 0):.2f}
Standard deviation of measurements: {measurement_summary.get('std_measurements_per_individual', 0):.2f}
Mean variance across measurements: {measurement_summary.get('mean_variance', 0):.6f}
Standard deviation of variance: {measurement_summary.get('std_variance', 0):.6f}
Minimum measurements per individual: {measurement_summary.get('min_measurements', 0)}
Maximum measurements per individual: {measurement_summary.get('max_measurements', 0)}

RELATIONSHIP ANALYSIS
=====================

Total pairs in kinship data: {relationship_summary['total_pairs']:,}
Closely related pairs (up to 2nd degree): {relationship_summary['closely_related_pairs']:,}
Total individuals in kinship data: {relationship_summary['total_individuals_in_kinship']:,}
Connected kinship records with multiple individuals: {relationship_summary['total_components']:,}
Individuals in connected kinship records: {relationship_summary['individuals_in_components']:,}

Relationship Type Distribution:
"""

        for rel_type, count in relationship_summary['relationship_type_distribution'].items():
            log_content += f"- {rel_type}: {count:,} pairs\n"

        if relationship_summary['component_size_distribution']:
            log_content += "\nKinship Record Size Distribution:\n"
            for size, count in relationship_summary['component_size_distribution'].items():
                log_content += f"- Size {size}: {count:,} kinship records\n"

        # Add trait value statistics if available
        if 'trait_statistics' in pruning_summary['original_data']:
            log_content += "\nTRAIT VALUE STATISTICS\n"
            log_content += "=====================\n\n"

            log_content += "Original Data:\n"
            stats = pruning_summary['original_data']['trait_statistics']
            log_content += f"- Mean: {stats['mean']:.6f}\n"
            log_content += f"- Std: {stats['std']:.6f}\n"
            log_content += f"- Min: {stats['min']:.6f}\n"
            log_content += f"- Max: {stats['max']:.6f}\n"
            log_content += f"- Median: {stats['median']:.6f}\n"

            if 'trait_statistics' in pruning_summary['pruned_data']:
                log_content += "\nPruned Data:\n"
                stats = pruning_summary['pruned_data']['trait_statistics']
                log_content += f"- Mean: {stats['mean']:.6f}\n"
                log_content += f"- Std: {stats['std']:.6f}\n"
                log_content += f"- Min: {stats['min']:.6f}\n"
                log_content += f"- Max: {stats['max']:.6f}\n"
                log_content += f"- Median: {stats['median']:.6f}\n"

        # Add validation results
        log_content += "\nVALIDATION RESULTS\n"
        log_content += "==================\n\n"

        log_content += f"Kinship records with trait data: {validation_results['components_with_trait_data']:,}\n"
        log_content += f"Kinship records without trait data: {validation_results['components_without_trait_data']:,}\n"
        log_content += f"Individuals in kinship records with trait data: {validation_results['individuals_in_components_with_trait_data']:,}\n"
        log_content += f"Individuals in kinship records without trait data: {validation_results['individuals_in_components_without_trait_data']:,}\n"

        # Add missing individuals information
        if missing_in_kinship:
            log_content += f"\nIndividuals in trait data but not in kinship data: {len(missing_in_kinship):,}\n"
            if len(missing_in_kinship) <= 10:
                log_content += "Missing individual IDs:\n"
                for fid, iid in missing_in_kinship:
                    log_content += f"- {fid}_{iid}\n"
            else:
                log_content += "Sample missing individual IDs (first 10):\n"
                for fid, iid in missing_in_kinship[:10]:
                    log_content += f"- {fid}_{iid}\n"
                log_content += f"... and {len(missing_in_kinship) - 10} more\n"

        log_content += "\nEND OF ENHANCED LOG\n"
        log_content += "==================\n"

        return log_content

    def save_pruning_log(self, trait_filename: str,
                        relationship_summary: Dict,
                        pruning_summary: Dict,
                        missing_in_kinship: List[Tuple[str, str]],
                        validation_results: Dict) -> Path:
        """
        Save detailed pruning log to file.

        Args:
            trait_filename: Original trait filename for naming
            relationship_summary: Summary of relationship analysis
            pruning_summary: Summary of pruning process
            missing_in_kinship: List of individuals missing in kinship data
            validation_results: Results of individual validation

        Returns:
            Path to the saved log file
        """
        # Generate output filename
        original_path = Path(trait_filename)
        output_filename = f"{original_path.stem}_pruning_log.txt"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving pruning log to: {output_path}")

        # Create log content
        log_content = self._generate_log_content(
            trait_filename, relationship_summary, pruning_summary,
            missing_in_kinship, validation_results
        )

        # Save log file
        with open(output_path, 'w') as f:
            f.write(log_content)

        self.logger.info(f"Pruning log saved: {output_path}")

        return output_path

    def _generate_log_content(self, trait_filename: str,
                            relationship_summary: Dict,
                            pruning_summary: Dict,
                            missing_in_kinship: List[Tuple[str, str]],
                            validation_results: Dict) -> str:
        """
        Generate the content for the pruning log file.

        Args:
            trait_filename: Original trait filename
            relationship_summary: Summary of relationship analysis
            pruning_summary: Summary of pruning process
            missing_in_kinship: List of individuals missing in kinship data
            validation_results: Results of individual validation

        Returns:
            String content for the log file
        """
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        log_content = f"""
RELATEDNESS PRUNING LOG
=======================

Generated: {timestamp}
Input file: {trait_filename}

SUMMARY STATISTICS
==================

Original Data:
- Total individuals: {pruning_summary['original_data']['total_individuals']:,}
- Non-missing trait values: {pruning_summary['original_data']['non_missing_trait_values']:,}
- Missing trait values: {pruning_summary['original_data']['missing_trait_values']:,}
- Proportion with non-missing values: {pruning_summary['original_data']['proportion_non_missing']:.3f}

Pruned Data:
- Total individuals retained: {pruning_summary['pruned_data']['total_individuals']:,}
- Non-missing trait values: {pruning_summary['pruned_data']['non_missing_trait_values']:,}
- Missing trait values: {pruning_summary['pruned_data']['missing_trait_values']:,}
- Proportion with non-missing values: {pruning_summary['pruned_data']['proportion_non_missing']:.3f}
- Number of individuals with non-missing values in pruned data: {pruning_summary['pruned_data']['non_missing_trait_values']:,}

Removed Individuals:
- Total individuals removed: {pruning_summary['removed_individuals']['total_individuals']:,}
- Non-missing trait values: {pruning_summary['removed_individuals']['non_missing_trait_values']:,}
- Missing trait values: {pruning_summary['removed_individuals']['missing_trait_values']:,}
- Proportion with non-missing values: {pruning_summary['removed_individuals']['proportion_non_missing']:.3f}

Pruning Efficiency:
- Removal rate: {pruning_summary['pruning_efficiency']['removal_rate']:.3f}
- Retention rate: {pruning_summary['pruning_efficiency']['retention_rate']:.3f}

RELATIONSHIP ANALYSIS
=====================

Total pairs in kinship data: {relationship_summary['total_pairs']:,}
Closely related pairs (up to 2nd degree): {relationship_summary['closely_related_pairs']:,}
Total individuals in kinship data: {relationship_summary['total_individuals_in_kinship']:,}
Connected kinship records with multiple individuals: {relationship_summary['total_components']:,}
Individuals in connected kinship records: {relationship_summary['individuals_in_components']:,}

Relationship Type Distribution:
"""

        for rel_type, count in relationship_summary['relationship_type_distribution'].items():
            log_content += f"- {rel_type}: {count:,} pairs\n"

        if relationship_summary['component_size_distribution']:
            log_content += "\nKinship Record Size Distribution:\n"
            for size, count in relationship_summary['component_size_distribution'].items():
                log_content += f"- Size {size}: {count:,} kinship records\n"

        # Add trait value statistics if available
        if 'trait_statistics' in pruning_summary['original_data']:
            log_content += "\nTRAIT VALUE STATISTICS\n"
            log_content += "=====================\n\n"

            log_content += "Original Data:\n"
            stats = pruning_summary['original_data']['trait_statistics']
            log_content += f"- Mean: {stats['mean']:.6f}\n"
            log_content += f"- Std: {stats['std']:.6f}\n"
            log_content += f"- Min: {stats['min']:.6f}\n"
            log_content += f"- Max: {stats['max']:.6f}\n"
            log_content += f"- Median: {stats['median']:.6f}\n"

            if 'trait_statistics' in pruning_summary['pruned_data']:
                log_content += "\nPruned Data:\n"
                stats = pruning_summary['pruned_data']['trait_statistics']
                log_content += f"- Mean: {stats['mean']:.6f}\n"
                log_content += f"- Std: {stats['std']:.6f}\n"
                log_content += f"- Min: {stats['min']:.6f}\n"
                log_content += f"- Max: {stats['max']:.6f}\n"
                log_content += f"- Median: {stats['median']:.6f}\n"

        # Add validation results
        log_content += "\nVALIDATION RESULTS\n"
        log_content += "==================\n\n"

        log_content += f"Kinship records with trait data: {validation_results['components_with_trait_data']:,}\n"
        log_content += f"Kinship records without trait data: {validation_results['components_without_trait_data']:,}\n"
        log_content += f"Individuals in kinship records with trait data: {validation_results['individuals_in_components_with_trait_data']:,}\n"
        log_content += f"Individuals in kinship records without trait data: {validation_results['individuals_in_components_without_trait_data']:,}\n"

        # Add missing individuals information
        if missing_in_kinship:
            log_content += f"\nIndividuals in trait data but not in kinship data: {len(missing_in_kinship):,}\n"
            if len(missing_in_kinship) <= 10:
                log_content += "Missing individual IDs:\n"
                for fid, iid in missing_in_kinship:
                    log_content += f"- {fid}_{iid}\n"
            else:
                log_content += "Sample missing individual IDs (first 10):\n"
                for fid, iid in missing_in_kinship[:10]:
                    log_content += f"- {fid}_{iid}\n"
                log_content += f"... and {len(missing_in_kinship) - 10} more\n"

        log_content += "\nEND OF LOG\n"
        log_content += "==========\n"

        return log_content

    def save_summary_json(self, trait_filename: str,
                         relationship_summary: Dict,
                         pruning_summary: Dict,
                         missing_in_kinship: List[Tuple[str, str]],
                         validation_results: Dict) -> Path:
        """
        Save summary statistics as JSON file.

        Args:
            trait_filename: Original trait filename for naming
            relationship_summary: Summary of relationship analysis
            pruning_summary: Summary of pruning process
            missing_in_kinship: List of individuals missing in kinship data
            validation_results: Results of individual validation

        Returns:
            Path to the saved JSON file
        """
        # Generate output filename
        original_path = Path(trait_filename)
        output_filename = f"{original_path.stem}_pruning_summary.json"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Saving summary JSON to: {output_path}")

        # Create summary dictionary
        summary_data = {
            'timestamp': datetime.now().isoformat(),
            'input_file': trait_filename,
            'relationship_summary': relationship_summary,
            'pruning_summary': pruning_summary,
            'validation_results': validation_results,
            'missing_in_kinship': [f"{fid}_{iid}" for fid, iid in missing_in_kinship]
        }

        # Save JSON file
        with open(output_path, 'w') as f:
            json.dump(summary_data, f, indent=2)

        self.logger.info(f"Summary JSON saved: {output_path}")

        return output_path

    def create_batch_summary(self, results: List[Dict]) -> Path:
        """
        Create a summary of batch processing results.

        Args:
            results: List of result dictionaries from processing multiple files

        Returns:
            Path to the batch summary file
        """
        output_filename = "batch_pruning_summary.txt"
        output_path = self.output_dir / output_filename

        self.logger.info(f"Creating batch summary: {output_path}")

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        summary_content = f"""
BATCH RELATEDNESS PRUNING SUMMARY
=================================

Generated: {timestamp}
Number of files processed: {len(results)}

DETAILED RESULTS
================

"""

        total_original = 0
        total_pruned = 0
        total_removed = 0

        for i, result in enumerate(results, 1):
            trait_filename = result['trait_filename']
            pruning_summary = result['pruning_summary']

            total_original += pruning_summary['original_data']['total_individuals']
            total_pruned += pruning_summary['pruned_data']['total_individuals']
            total_removed += pruning_summary['removed_individuals']['total_individuals']

            summary_content += f"{i}. {trait_filename}\n"
            summary_content += f"   Original: {pruning_summary['original_data']['total_individuals']:,} individuals\n"
            summary_content += f"   Pruned: {pruning_summary['pruned_data']['total_individuals']:,} individuals\n"
            summary_content += f"   Removed: {pruning_summary['removed_individuals']['total_individuals']:,} individuals\n"
            summary_content += f"   Retention rate: {pruning_summary['pruning_efficiency']['retention_rate']:.3f}\n\n"

        summary_content += f"""
OVERALL SUMMARY
===============

Total individuals across all files:
- Original: {total_original:,}
- Pruned: {total_pruned:,}
- Removed: {total_removed:,}
- Overall retention rate: {total_pruned/total_original:.3f}
- Overall removal rate: {total_removed/total_original:.3f}

END OF BATCH SUMMARY
====================
"""

        # Save batch summary
        with open(output_path, 'w') as f:
            f.write(summary_content)

        self.logger.info(f"Batch summary saved: {output_path}")

        return output_path
