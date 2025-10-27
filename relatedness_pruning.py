#!/usr/bin/env python3
"""
Main script for relatedness pruning workflow.

This script performs related individual pruning from large cohorts for GWAS analysis
based on pre-calculated kinship coefficients and trait-based prioritization.

Usage:
    python relatedness_pruning.py --kinship-file <kinship_file> --trait-file <trait_file> --output-dir <output_dir>
    python relatedness_pruning.py --config <config_file> --output-dir <output_dir>
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
    logger = logging.getLogger("relatedness_pruning")
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


def process_single_trait_file(kinship_file: str, trait_file: str, output_dir: str,
                            logger: logging.Logger, rank_normalize: bool = False) -> Dict:
    """
    Process a single trait file for relatedness pruning.

    Args:
        kinship_file: Path to kinship file
        trait_file: Path to trait file
        output_dir: Output directory
        logger: Logger instance
        rank_normalize: Whether to apply rank normalization to trait values

    Returns:
        Dictionary with processing results
    """
    logger.info(f"Processing trait file: {trait_file}")

    # Initialize components
    data_loader = DataLoader(logger)
    relationship_filter = RelationshipFilter(logger)
    pruning_engine = PruningEngine(logger)
    output_generator = OutputGenerator(output_dir, logger)

    try:
        # Load data
        logger.info("Loading kinship data...")
        kinship_df = data_loader.load_kinship_data(kinship_file)

        logger.info("Loading trait data...")
        trait_df = data_loader.load_trait_data(trait_file)

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

        # Prune related individuals
        logger.info("Pruning related individuals...")
        retained_individuals, removed_individuals = pruning_engine.prune_related_individuals(
            connected_components, individual_lookup
        )

        # Create pruned trait data
        logger.info("Creating pruned trait data...")
        pruned_trait_df = pruning_engine.create_pruned_trait_data(trait_df, retained_individuals, rank_normalize)

        # Generate summaries
        logger.info("Generating summaries...")
        relationship_summary = relationship_filter.get_relationship_summary(
            kinship_df, closely_related_df, connected_components
        )

        pruning_summary = pruning_engine.generate_pruning_summary(
            trait_df, pruned_trait_df, removed_individuals, individual_lookup
        )

        # Save outputs
        logger.info("Saving outputs...")
        pruned_file_path = output_generator.save_pruned_trait_data(pruned_trait_df, trait_file)
        removed_file_path = output_generator.save_removed_individuals(removed_individuals, trait_file)
        log_file_path = output_generator.save_pruning_log(
            trait_file, relationship_summary, pruning_summary,
            missing_in_kinship, validation_results
        )
        summary_json_path = output_generator.save_summary_json(
            trait_file, relationship_summary, pruning_summary,
            missing_in_kinship, validation_results
        )

        logger.info("Processing completed successfully!")
        logger.info(f"Pruned data saved to: {pruned_file_path}")
        logger.info(f"Removed individuals saved to: {removed_file_path}")
        logger.info(f"Log saved to: {log_file_path}")
        logger.info(f"Summary JSON saved to: {summary_json_path}")

        return {
            'trait_filename': trait_file,
            'pruned_file_path': str(pruned_file_path),
            'removed_file_path': str(removed_file_path),
            'log_file_path': str(log_file_path),
            'summary_json_path': str(summary_json_path),
            'relationship_summary': relationship_summary,
            'pruning_summary': pruning_summary,
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


def process_multiple_trait_files(config: Dict, output_dir: str, logger: logging.Logger, rank_normalize: bool = False) -> List[Dict]:
    """
    Process multiple trait files based on configuration.

    Args:
        config: Configuration dictionary
        output_dir: Output directory
        logger: Logger instance
        rank_normalize: Whether to apply rank normalization to trait values

    Returns:
        List of processing results
    """
    kinship_file = config['kinship_file']
    trait_files = config['trait_files']

    logger.info(f"Processing {len(trait_files)} trait files...")

    results = []
    for i, trait_file in enumerate(trait_files, 1):
        logger.info(f"Processing file {i}/{len(trait_files)}: {trait_file}")
        result = process_single_trait_file(kinship_file, trait_file, output_dir, logger, rank_normalize)
        results.append(result)

        if result['success']:
            logger.info(f"Successfully processed: {trait_file}")
        else:
            logger.error(f"Failed to process: {trait_file}")

    # Create batch summary
    successful_results = [r for r in results if r['success']]
    if successful_results:
        output_generator = OutputGenerator(output_dir, logger)
        batch_summary_path = output_generator.create_batch_summary(successful_results)
        logger.info(f"Batch summary saved to: {batch_summary_path}")

    return results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Relatedness pruning workflow for GWAS analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single trait file
  python relatedness_pruning.py --kinship-file kinship.kin0 --trait-file trait.tsv --output-dir ./output

  # Process multiple trait files using config
  python relatedness_pruning.py --config config.json --output-dir ./output

  # With custom log level
  python relatedness_pruning.py --kinship-file kinship.kin0 --trait-file trait.tsv --output-dir ./output --log-level DEBUG

  # With rank normalization
  python relatedness_pruning.py --kinship-file kinship.kin0 --trait-file trait.tsv --output-dir ./output --rank-normalize
        """
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--kinship-file',
        type=str,
        help='Path to kinship file (KING output)'
    )
    input_group.add_argument(
        '--config',
        type=str,
        help='Path to JSON configuration file for batch processing'
    )

    # Trait file option (required if kinship-file is specified)
    parser.add_argument(
        '--trait-file',
        type=str,
        help='Path to trait file (required if --kinship-file is specified)'
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
        '--rank-normalize',
        action='store_true',
        help='Apply rank normalization to trait values after pruning'
    )

    args = parser.parse_args()

    # Validate arguments
    if args.kinship_file and not args.trait_file:
        parser.error("--trait-file is required when --kinship-file is specified")

    # Set up logging
    logger = setup_logging(args.log_level, args.log_file)

    try:
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if args.config:
            # Process multiple files using config
            logger.info(f"Loading configuration from: {args.config}")
            data_loader = DataLoader(logger)
            config = data_loader.load_config(args.config)

            results = process_multiple_trait_files(config, str(output_dir), logger, args.rank_normalize)

            # Report results
            successful = sum(1 for r in results if r['success'])
            failed = len(results) - successful

            logger.info(f"Batch processing completed: {successful} successful, {failed} failed")

            if failed > 0:
                logger.warning("Some files failed to process. Check logs for details.")
                sys.exit(1)

        else:
            # Process single file
            result = process_single_trait_file(
                args.kinship_file, args.trait_file, str(output_dir), logger, args.rank_normalize
            )

            if not result['success']:
                logger.error("Processing failed. Check logs for details.")
                sys.exit(1)

    except KeyboardInterrupt:
        logger.info("Processing interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
