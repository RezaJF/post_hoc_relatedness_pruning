"""
Data loading module for relatedness pruning workflow.

This module handles loading and validation of kinship data from KING software
and trait data files for the related individual pruning pipeline.
"""

import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
import json


class DataLoader:
    """Handles loading and validation of kinship and trait data."""

    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize DataLoader.

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

    def load_kinship_data(self, kinship_file: Union[str, Path]) -> pd.DataFrame:
        """
        Load kinship data from KING software output.

        Args:
            kinship_file: Path to the kinship file (e.g., finngen_R13.kin0)

        Returns:
            DataFrame with kinship data

        Raises:
            FileNotFoundError: If kinship file doesn't exist
            ValueError: If required columns are missing
        """
        kinship_path = Path(kinship_file)
        if not kinship_path.exists():
            raise FileNotFoundError(f"Kinship file not found: {kinship_path}")

        self.logger.info(f"Loading kinship data from: {kinship_path}")

        try:
            # Load the kinship data
            kinship_df = pd.read_csv(kinship_path, sep='\t')
            self.logger.info(f"Loaded kinship data with {len(kinship_df)} pairs")

            # Validate required columns
            required_columns = ['FID1', 'ID1', 'FID2', 'ID2', 'InfType']
            missing_columns = [col for col in required_columns if col not in kinship_df.columns]

            if missing_columns:
                raise ValueError(f"Missing required columns in kinship file: {missing_columns}")

            # Check for Kinship column (optional but recommended)
            if 'Kinship' not in kinship_df.columns:
                self.logger.warning("Kinship column not found. Only InfType will be used for filtering.")

            # Log basic statistics
            self.logger.info(f"Kinship data columns: {list(kinship_df.columns)}")
            if 'InfType' in kinship_df.columns:
                inf_type_counts = kinship_df['InfType'].value_counts()
                self.logger.info(f"Relationship type distribution:\n{inf_type_counts}")

            return kinship_df

        except Exception as e:
            self.logger.error(f"Error loading kinship data: {str(e)}")
            raise

    def load_trait_data(self, trait_file: Union[str, Path]) -> pd.DataFrame:
        """
        Load trait data file.

        Args:
            trait_file: Path to the trait data file

        Returns:
            DataFrame with trait data

        Raises:
            FileNotFoundError: If trait file doesn't exist
            ValueError: If required columns are missing
        """
        trait_path = Path(trait_file)
        if not trait_path.exists():
            raise FileNotFoundError(f"Trait file not found: {trait_path}")

        self.logger.info(f"Loading trait data from: {trait_path}")

        try:
            # Load the trait data
            trait_df = pd.read_csv(trait_path, sep='\t')
            self.logger.info(f"Loaded trait data with {len(trait_df)} individuals")

            # Validate required columns
            required_columns = ['FID', 'IID']
            missing_columns = [col for col in required_columns if col not in trait_df.columns]

            if missing_columns:
                raise ValueError(f"Missing required columns in trait file: {missing_columns}")

            # Identify the trait column (third column or first non-FID/IID column)
            trait_columns = [col for col in trait_df.columns if col not in ['FID', 'IID']]
            if not trait_columns:
                raise ValueError("No trait column found in trait file")

            # Use the first non-FID/IID column as trait column
            trait_column = trait_columns[0]
            self.logger.info(f"Using '{trait_column}' as trait column")

            # Add trait column name as metadata
            trait_df.attrs['trait_column'] = trait_column

            # Log basic statistics
            self.logger.info(f"Trait data columns: {list(trait_df.columns)}")

            # Check for missing values
            missing_count = trait_df[trait_column].isna().sum()
            non_missing_count = len(trait_df) - missing_count
            self.logger.info(f"Trait values: {non_missing_count} non-missing, {missing_count} missing")

            if non_missing_count > 0:
                trait_stats = trait_df[trait_column].describe()
                self.logger.info(f"Trait statistics:\n{trait_stats}")

            return trait_df

        except Exception as e:
            self.logger.error(f"Error loading trait data: {str(e)}")
            raise

    def load_config(self, config_file: Union[str, Path]) -> Dict:
        """
        Load configuration from JSON file.

        Args:
            config_file: Path to the JSON configuration file

        Returns:
            Dictionary with configuration parameters

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If config file is invalid JSON
        """
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        self.logger.info(f"Loading configuration from: {config_path}")

        try:
            with open(config_path, 'r') as f:
                config = json.load(f)

            # Validate required configuration keys
            required_keys = ['kinship_file', 'trait_files']
            missing_keys = [key for key in required_keys if key not in config]

            if missing_keys:
                raise ValueError(f"Missing required configuration keys: {missing_keys}")

            self.logger.info(f"Configuration loaded successfully")
            self.logger.info(f"Kinship file: {config['kinship_file']}")
            self.logger.info(f"Number of trait files: {len(config['trait_files'])}")

            return config

        except json.JSONDecodeError as e:
            self.logger.error(f"Invalid JSON in config file: {str(e)}")
            raise ValueError(f"Invalid JSON in config file: {str(e)}")
        except Exception as e:
            self.logger.error(f"Error loading configuration: {str(e)}")
            raise

    def load_enhanced_config(self, config_file: Union[str, Path]) -> Dict:
        """
        Load enhanced configuration from JSON file.

        Args:
            config_file: Path to the enhanced JSON configuration file

        Returns:
            Dictionary with enhanced configuration parameters

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If config file is invalid JSON
        """
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Enhanced config file not found: {config_path}")

        self.logger.info(f"Loading enhanced configuration from: {config_path}")

        try:
            with open(config_path, 'r') as f:
                config = json.load(f)

            # Validate required configuration keys
            required_keys = ['kinship_file', 'trait_files']
            missing_keys = [key for key in required_keys if key not in config]

            if missing_keys:
                raise ValueError(f"Missing required configuration keys: {missing_keys}")

            # Validate trait files structure
            for i, trait_file in enumerate(config['trait_files']):
                required_trait_keys = ['trait_file_path', 'omop_concept_id', 'atc_codes']
                missing_trait_keys = [key for key in required_trait_keys if key not in trait_file]

                if missing_trait_keys:
                    raise ValueError(f"Missing required keys in trait file {i+1}: {missing_trait_keys}")

            self.logger.info(f"Enhanced configuration loaded successfully")
            self.logger.info(f"Kinship file: {config['kinship_file']}")
            self.logger.info(f"Number of trait files: {len(config['trait_files'])}")

            return config

        except json.JSONDecodeError as e:
            self.logger.error(f"Invalid JSON in enhanced config file: {str(e)}")
            raise ValueError(f"Invalid JSON in enhanced config file: {str(e)}")
        except Exception as e:
            self.logger.error(f"Error loading enhanced configuration: {str(e)}")
            raise

    def validate_individual_overlap(self, kinship_df: pd.DataFrame, trait_df: pd.DataFrame) -> Tuple[List[str], List[str]]:
        """
        Validate overlap between kinship and trait data.

        Args:
            kinship_df: DataFrame with kinship data
            trait_df: DataFrame with trait data

        Returns:
            Tuple of (missing_in_kinship, missing_in_trait) individual lists
        """
        # Get unique individuals from both datasets
        kinship_individuals = set()
        for _, row in kinship_df.iterrows():
            kinship_individuals.add((row['FID1'], row['ID1']))
            kinship_individuals.add((row['FID2'], row['ID2']))

        trait_individuals = set()
        for _, row in trait_df.iterrows():
            trait_individuals.add((row['FID'], row['IID']))

        # Find individuals in trait data but not in kinship data
        missing_in_kinship = trait_individuals - kinship_individuals
        missing_in_trait = kinship_individuals - trait_individuals

        if missing_in_kinship:
            missing_ids = [f"{fid}_{iid}" for fid, iid in missing_in_kinship]
            self.logger.warning(f"Individuals in trait data but not in kinship data: {len(missing_ids)}")
            self.logger.warning(f"Sample missing IDs: {missing_ids[:5]}...")

        if missing_in_trait:
            missing_ids = [f"{fid}_{iid}" for fid, iid in missing_in_trait]
            self.logger.info(f"Individuals in kinship data but not in trait data: {len(missing_ids)}")

        return list(missing_in_kinship), list(missing_in_trait)

    def create_individual_id(self, fid: str, iid: str) -> str:
        """
        Create a unique individual identifier.

        Args:
            fid: Family ID
            iid: Individual ID

        Returns:
            Unique identifier string
        """
        return f"{fid}_{iid}"

    def parse_individual_id(self, individual_id: str) -> Tuple[str, str]:
        """
        Parse individual identifier back to FID and IID.

        Args:
            individual_id: Unique identifier string

        Returns:
            Tuple of (FID, IID)
        """
        parts = individual_id.split('_', 1)
        if len(parts) != 2:
            raise ValueError(f"Invalid individual ID format: {individual_id}")
        return parts[0], parts[1]
