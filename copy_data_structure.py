#!/usr/bin/env python3
"""
Script to copy SAXS and WAXS data files from atT/ directory structure 
to larger_test/2D/ with flattened organization.

Excludes: standard/, poni/, OneD_integrated* directories
Copies: .raw, .pdi files from SAXS/ and WAXS/ subdirectories, and .csv files from experiment root
"""

import os
import shutil
from pathlib import Path
import argparse
import logging

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('copy_data_structure.log'),
            logging.StreamHandler()
        ]
    )

def should_exclude_directory(dir_name):
    """
    Check if directory should be excluded based on naming rules
    
    Args:
        dir_name (str): Directory name to check
        
    Returns:
        bool: True if directory should be excluded
    """
    exclude_patterns = [
        'standard',
        'poni', 
        'OneD_integrated',
        'Run13',
        'Run12'
    ]
    
    for pattern in exclude_patterns:
        if dir_name.lower() == pattern.lower() or dir_name.startswith(pattern):
            return True
    return False

def copy_saxs_waxs_files(source_dir, target_dir, dry_run=False):
    """
    Copy SAXS and WAXS files from source directory structure to target
    
    Args:
        source_dir (Path): Source directory (atT/)
        target_dir (Path): Target directory (larger_test/2D/)
        dry_run (bool): If True, only show what would be copied
    """
    logger = logging.getLogger(__name__)
    
    # Create target directories
    saxs_target = target_dir / 'SAXS'
    waxs_target = target_dir / 'WAXS'
    
    if not dry_run:
        saxs_target.mkdir(parents=True, exist_ok=True)
        waxs_target.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created target directories: {saxs_target}, {waxs_target}")
    else:
        logger.info(f"DRY RUN: Would create directories: {saxs_target}, {waxs_target}")
    
    # Stats tracking
    stats = {
        'directories_processed': 0,
        'directories_excluded': 0,
        'saxs_files_copied': 0,
        'waxs_files_copied': 0,
        'errors': 0
    }
    
    # Iterate through all directories in source
    for item in source_dir.iterdir():
        if not item.is_dir():
            continue
            
        # Check if directory should be excluded
        if should_exclude_directory(item.name):
            logger.info(f"Excluding directory: {item.name}")
            stats['directories_excluded'] += 1
            continue
            
        logger.info(f"Processing directory: {item.name}")
        stats['directories_processed'] += 1
        
        # Copy CSV files from experiment root directory
        csv_files_copied = copy_csv_files(item, target_dir, dry_run)
        stats['saxs_files_copied'] += csv_files_copied  # Add to saxs count for summary
        
        # Process SAXS subdirectory
        saxs_dir = item / 'SAXS'
        if saxs_dir.exists():
            files_copied = copy_detector_files(saxs_dir, saxs_target, target_dir, 'SAXS', dry_run)
            stats['saxs_files_copied'] += files_copied
        else:
            logger.warning(f"No SAXS directory found in {item.name}")
            
        # Process WAXS subdirectory  
        waxs_dir = item / 'WAXS'
        if waxs_dir.exists():
            files_copied = copy_detector_files(waxs_dir, waxs_target, target_dir, 'WAXS', dry_run)
            stats['waxs_files_copied'] += files_copied
        else:
            logger.warning(f"No WAXS directory found in {item.name}")
    
    # Print summary statistics
    logger.info("="*50)
    logger.info("COPY OPERATION SUMMARY")
    logger.info("="*50)
    logger.info(f"Directories processed: {stats['directories_processed']}")
    logger.info(f"Directories excluded: {stats['directories_excluded']}")
    logger.info(f"SAXS files copied: {stats['saxs_files_copied']}")
    logger.info(f"WAXS files copied: {stats['waxs_files_copied']}")
    logger.info(f"Total files copied: {stats['saxs_files_copied'] + stats['waxs_files_copied']}")
    if stats['errors'] > 0:
        logger.error(f"Errors encountered: {stats['errors']}")

def copy_csv_files(experiment_dir, target_dir, dry_run):
    """
    Copy CSV files from experiment root directory to target 2D directory
    
    Args:
        experiment_dir (Path): Source experiment directory
        target_dir (Path): Target 2D directory
        dry_run (bool): If True, only show what would be copied
        
    Returns:
        int: Number of CSV files copied
    """
    logger = logging.getLogger(__name__)
    files_copied = 0
    
    # Find all .csv files in experiment root
    csv_files = list(experiment_dir.glob('*.csv'))
    
    for file_path in csv_files:
        target_file = target_dir / file_path.name
        
        if target_file.exists():
            logger.warning(f"CSV file already exists, skipping: {file_path.name}")
            continue
            
        try:
            if not dry_run:
                shutil.copy2(file_path, target_file)
                logger.debug(f"Copied CSV to 2D/: {file_path.name}")
            else:
                logger.info(f"DRY RUN: Would copy CSV to 2D/: {file_path.name}")
            files_copied += 1
            
        except Exception as e:
            logger.error(f"Failed to copy CSV {file_path}: {e}")
            
    if files_copied > 0:
        logger.info(f"Copied {files_copied} CSV files from {experiment_dir.name}")
        
    return files_copied

def copy_detector_files(source_detector_dir, target_detector_dir, csv_target_dir, detector_type, dry_run):
    """
    Copy .raw and .pdi files from detector directory to target
    
    Args:
        source_detector_dir (Path): Source SAXS or WAXS directory
        target_detector_dir (Path): Target SAXS or WAXS directory
        csv_target_dir (Path): Not used (kept for compatibility)
        detector_type (str): 'SAXS' or 'WAXS' for logging
        dry_run (bool): If True, only show what would be copied
        
    Returns:
        int: Number of files copied
    """
    logger = logging.getLogger(__name__)
    files_copied = 0
    
    # Find all .raw and .pdi files (CSV files handled separately)
    file_extensions = ['*.raw', '*.pdi']
    files_to_copy = []
    
    for pattern in file_extensions:
        files_to_copy.extend(source_detector_dir.glob(pattern))
    
    for file_path in files_to_copy:
        target_file = target_detector_dir / file_path.name
        
        if target_file.exists():
            logger.warning(f"File already exists, skipping: {file_path.name}")
            continue
            
        try:
            if not dry_run:
                shutil.copy2(file_path, target_file)
                logger.debug(f"Copied {detector_type}: {file_path.name}")
            else:
                logger.info(f"DRY RUN: Would copy {detector_type}: {file_path.name}")
            files_copied += 1
            
        except Exception as e:
            logger.error(f"Failed to copy {file_path}: {e}")
            
    if files_copied > 0:
        logger.info(f"Copied {files_copied} {detector_type} files from {source_detector_dir.parent.name}")
        
    return files_copied

def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(description='Copy SAXS/WAXS data with flattened structure')
    parser.add_argument('--source', required= True, 
                       help='Source directory')
    parser.add_argument('--target', required= True,
                       help='Target directory')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be copied without actually copying')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        
    # Convert to Path objects
    source_dir = Path(args.source)
    target_dir = Path(args.target)
    
    # Validate source directory
    if not source_dir.exists():
        logger.error(f"Source directory does not exist: {source_dir}")
        return 1
        
    logger.info(f"Source directory: {source_dir.absolute()}")
    logger.info(f"Target directory: {target_dir.absolute()}")
    
    if args.dry_run:
        logger.info("DRY RUN MODE - No files will be copied")
    
    # Perform the copy operation
    try:
        copy_saxs_waxs_files(source_dir, target_dir, args.dry_run)
        logger.info("Copy operation completed successfully")
        return 0
    except Exception as e:
        logger.error(f"Copy operation failed: {e}")
        return 1

if __name__ == '__main__':
    exit(main())