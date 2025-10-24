#!/usr/bin/env python3
"""
Enhanced main script for relatedness pruning workflow with lab measurement analysis.

This script performs related individual pruning from large cohorts for GWAS analysis
based on lab measurement statistics and prioritization criteria.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

# Import our custom modules
from data_loader import DataLoader
from relationship_filter import RelationshipFilter
from pruning_engine import PruningEngine
from output_generator import OutputGenerator
from lab_measurement_analyzer import LabMeasurementAnalyzer


def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None) -> logging.Logger:
    """
    Set up logging configuration.

    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional log file path

    Returns:
        Configured logger instance
    """
    # Create logger
    logger = logging.getLogger("enhanced_relatedness_pruning")
    logger.setLevel(getattr(logging, log_level.upper()))

    # Clear existing handlers
    logger.handlers.clear()

    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def process_single_trait_file_enhanced(config_entry: Dict, output_dir: str,
                                     logger: logging.Logger, inverse_rank_normalize: bool = True) -> Dict:
    """
    Process a single trait file for enhanced relatedness pruning.

    Args:
        config_entry: Configuration entry for the trait file
        output_dir: Output directory
        logger: Logger instance
        inverse_rank_normalize: Whether to apply inverse rank normalization after pruning

    Returns:
        Dictionary with processing results
    """
    trait_file = config_entry['trait_file_path']
    omop_concept_id = config_entry['omop_concept_id']
    atc_codes = config_entry['atc_codes']
    months_before = config_entry.get('months_before', 1)

    logger.info(f"Processing trait file: {trait_file}")
    logger.info(f"OMOP concept ID: {omop_concept_id}")
    logger.info(f"ATC codes: {atc_codes}")

    # Initialize components
    data_loader = DataLoader(logger)
    relationship_filter = RelationshipFilter(logger)
    output_generator = OutputGenerator(output_dir, logger)

    # Initialize lab measurement analyzer
    db_config_path = config_entry.get('db_config_path', '/mnt/longGWAS_disk_100GB/long_gwas/Github_clones/fganalysis-r/config/db_config_local.json')
    lab_analyzer = LabMeasurementAnalyzer(db_config_path, logger)
    pruning_engine = PruningEngine(logger, lab_analyzer)

    try:
        # Load kinship data
        logger.info("Loading kinship data...")
        kinship_file = config_entry.get('kinship_file', '/mnt/longGWAS_disk_100GB/long_gwas/analysis_covariates/finngen_R13.kin0')
        kinship_df = data_loader.load_kinship_data(kinship_file)

        # Load trait data
        logger.info("Loading trait data...")
        trait_df = data_loader.load_trait_data(trait_file)

        # Get measurement statistics
        logger.info("Getting lab measurement statistics...")
        measurement_stats = lab_analyzer.get_measurement_statistics(
            omop_concept_id, atc_codes, months_before
        )

        # Create measurement lookup
        measurement_lookup = lab_analyzer.create_measurement_lookup(measurement_stats)

        # Validate individual overlap
        logger.info("Validating individual overlap...")
        missing_in_kinship, missing_in_trait = data_loader.validate_individual_overlap(
            kinship_df, trait_df
        )

        # Filter closely related pairs
        logger.info("Filtering closely related pairs...")
        closely_related_df = relationship_filter.filter_closely_related_pairs(kinship_df)

        # Build relationship graph
        logger.info("Building relationship graph...")
        relationship_graph = relationship_filter.build_relationship_graph(closely_related_df)

        # Find connected components
        logger.info("Finding connected components...")
        connected_components = relationship_filter.find_connected_components(relationship_graph)

        # Create individual lookup
        individual_lookup = relationship_filter.create_individual_lookup(trait_df)

        # Validate individuals in components
        trait_individuals = set(individual_lookup.keys())
        validation_results = relationship_filter.validate_individuals_in_components(
            connected_components, trait_individuals
        )

        # Prune related individuals using measurement-based logic
        logger.info("Pruning related individuals using measurement-based prioritization...")
        retained_individuals, removed_individuals = pruning_engine.prune_related_individuals(
            connected_components, individual_lookup, measurement_lookup
        )

        # Create pruned trait data
        logger.info("Creating pruned trait data...")
        pruned_trait_df = pruning_engine.create_pruned_trait_data(trait_df, retained_individuals, False, inverse_rank_normalize)

        # Generate summaries
        logger.info("Generating summaries...")
        relationship_summary = relationship_filter.get_relationship_summary(
            kinship_df, closely_related_df, connected_components
        )

        pruning_summary = pruning_engine.generate_pruning_summary(
            trait_df, pruned_trait_df, removed_individuals, individual_lookup
        )

        # Get measurement summary
        measurement_summary = lab_analyzer.get_measurement_summary(measurement_lookup)

        # Save outputs
        logger.info("Saving outputs...")
        pruned_file_path = output_generator.save_pruned_trait_data(pruned_trait_df, trait_file)
        removed_file_path = output_generator.save_removed_individuals(removed_individuals, trait_file)

        # Save description file
        description_file_path = output_generator.save_description_file(
            trait_file, omop_concept_id, atc_codes, months_before,
            len(retained_individuals), len(measurement_stats)
        )

        # Generate enhanced log
        log_file_path = output_generator.save_enhanced_pruning_log(
            trait_file, relationship_summary, pruning_summary, measurement_summary,
            missing_in_kinship, validation_results
        )

        # Generate enhanced JSON summary
        summary_json_path = output_generator.save_enhanced_summary_json(
            trait_file, relationship_summary, pruning_summary, measurement_summary,
            missing_in_kinship, validation_results
        )

        logger.info("Processing completed successfully!")
        logger.info(f"Pruned data saved to: {pruned_file_path}")
        logger.info(f"Removed individuals saved to: {removed_file_path}")
        logger.info(f"Description file saved to: {description_file_path}")
        logger.info(f"Log saved to: {log_file_path}")
        logger.info(f"Summary JSON saved to: {summary_json_path}")

        return {
            'trait_filename': trait_file,
            'pruned_file_path': str(pruned_file_path),
            'removed_file_path': str(removed_file_path),
            'description_file_path': str(description_file_path),
            'log_file_path': str(log_file_path),
            'summary_json_path': str(summary_json_path),
            'relationship_summary': relationship_summary,
            'pruning_summary': pruning_summary,
            'measurement_summary': measurement_summary,
            'validation_results': validation_results,
            'missing_in_kinship': missing_in_kinship,
            'success': True
        }

    except Exception as e:
        logger.error(f"Error processing {trait_file}: {str(e)}")
        return {
            'trait_filename': trait_file,
            'error': str(e),
            'success': False
        }


def process_enhanced_config(config: Dict, output_dir: str, logger: logging.Logger, inverse_rank_normalize: bool = True) -> List[Dict]:
    """
    Process multiple trait files based on enhanced configuration.

    Args:
        config: Enhanced configuration dictionary
        output_dir: Output directory (can be overridden by config)
        logger: Logger instance
        inverse_rank_normalize: Whether to apply inverse rank normalization after pruning

    Returns:
        List of processing results
    """
    # Use output directory from config if specified, otherwise use provided output_dir
    final_output_dir = config.get('output_dir', output_dir)
    trait_files = config['trait_files']

    logger.info(f"Processing {len(trait_files)} trait files with enhanced measurement-based pruning...")

    results = []
    for i, trait_config in enumerate(trait_files, 1):
        logger.info(f"Processing file {i}/{len(trait_files)}: {trait_config['trait_file_path']}")
        result = process_single_trait_file_enhanced(trait_config, final_output_dir, logger, inverse_rank_normalize)
        results.append(result)

        if result['success']:
            logger.info(f"Successfully processed: {trait_config['trait_file_path']}")
        else:
            logger.error(f"Failed to process: {trait_config['trait_file_path']}")

    # Create batch summary
    successful_results = [r for r in results if r['success']]
    if successful_results:
        output_generator = OutputGenerator(final_output_dir, logger)
        batch_summary_path = output_generator.create_enhanced_batch_summary(successful_results)
        logger.info(f"Enhanced batch summary saved to: {batch_summary_path}")

    return results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Enhanced relatedness pruning workflow with lab measurement analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process with enhanced config
  python enhanced_relatedness_pruning.py --config enhanced_config.json --output-dir ./output

  # With custom log level
  python enhanced_relatedness_pruning.py --config enhanced_config.json --output-dir ./output --log-level DEBUG
        """
    )

    # Input options
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to enhanced JSON configuration file'
    )

    # Output options
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )

    # Logging options
    parser.add_argument(
        '--log-level',
        type=str,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level (default: INFO)'
    )

    parser.add_argument(
        '--log-file',
        type=str,
        help='Log file path (optional)'
    )

    parser.add_argument(
        '--no-inverse-rank-normalize',
        action='store_true',
        help='Disable inverse rank normalization after pruning (default: enabled)'
    )

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.log_level, args.log_file)

    try:
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load enhanced configuration
        logger.info(f"Loading enhanced configuration from: {args.config}")
        data_loader = DataLoader(logger)
        config = data_loader.load_enhanced_config(args.config)

        # Determine inverse rank normalization setting
        inverse_rank_normalize = not args.no_inverse_rank_normalize
        
        results = process_enhanced_config(config, str(output_dir), logger, inverse_rank_normalize)

        # Report results
        successful = sum(1 for r in results if r['success'])
        failed = len(results) - successful

        logger.info(f"Enhanced batch processing completed: {successful} successful, {failed} failed")

        if failed > 0:
            logger.warning("Some files failed to process. Check logs for details.")
            sys.exit(1)

    except KeyboardInterrupt:
        logger.info("Processing interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
