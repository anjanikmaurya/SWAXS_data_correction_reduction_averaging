#!/usr/bin/env python3
"""
SWAXS Data Correction and 1D Integration Pipeline
=================================================
Refactored implementation of Step1 notebook with improved modularity and configuration management.
Processes all .raw files in directory structure with essential data corrections and 1D integration.

Status: ACTIVE DEVELOPMENT - Core functionality implemented with accurate calculations, requires extension

Key Features:
- Experiment class for centralized configuration management from config.yml
- Multi-file processing (all .raw files in SAXS/WAXS subdirectories)
- Essential data corrections (transmission, normalization, material properties)
- Accurate SLD calculation using periodictable for precise material properties
- 1D radial integration with PyFAI
- Configurable experimental parameters and detector geometry
- Output in .dat format compatible with existing analysis workflows
- Uses PDI files instead of CSV

File Structure Requirements:
- Input: {data_directory}/[experiment_dirs]/SAXS/*.raw and WAXS/*.raw
- Metadata: Corresponding .csv files with same name as .raw files
- Output: 1D/[experiment_dirs]/*.dat (3-column: q, intensity, error)

Dependencies: pyFAI, fabio, xraydb, numpy, pandas, yaml
"""
import os
import glob
import numpy as np
import pandas as pd
import yaml
import fabio
import pyFAI
import periodictable as pt
import xraydb
from typing import Dict, List, Tuple
import re
import utils


class Experiment:
    """
    SWAXS Experiment class for managing configuration and data processing.
    
    Loads experimental parameters from config.yml and provides methods for
    processing SAXS/WAXS data with proper corrections and normalization.
    """
    
    def __init__(self, config_path: str = "config.yml"):
        """
        Initialize Experiment with configuration from YAML file.
        
        Parameters
        ----------
        config_path : str
            Path to the configuration file
        """
        self.config_path = config_path
        self.config = self._load_config()
        self._setup_directories()
        self._setup_parameters()
        self._load_pyfai_integrators()
        
    def _load_config(self) -> Dict:
        """Load configuration from YAML file."""
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _setup_directories(self):
        """Setup directory paths from config."""
        self.data_directory_2d = os.path.join(self.config['data_directory'], "2D")
        self.poni_directory = os.path.join(self.config['data_directory'], self.config['poni_directory'])
        self.output_directory_1d = os.path.join(self.config['data_directory'], "1D")
        os.makedirs(self.output_directory_1d, exist_ok = True)
        
        self.saxs_subdir = "SAXS" 
        self.waxs_subdir = "WAXS"
        self.detector_shapes = self.config['detector_shapes']
        self.poni_files = self.config['poni_files']
        self.mask_files = self.config['mask_files']
        
    def _setup_parameters(self):
        """Setup experimental parameters from config."""
        self.compound = self.config['compound']
        self.energy_keV = self.config['energy_keV']
        self.density_g_cm3 = self.config['density_g_cm3']
        self.i0_offset = self.config['i0_offset']
        self.bstop_offset = self.config['bstop_offset']
        self.i0_air = self.config['i0_air']
        self.bstop_air = self.config['bstop_air']
        self.thickness = self.config.get('thickness', None)
        self.npt_radial = self.config['npt_radial']
        self.error_model = self.config['error_model']
        
        self.metadata_function = self.config['read_metadata_function']
    
    def _load_pyfai_integrators(self):
        """Load PyFAI integrator objects for SAXS and WAXS."""
        self.saxs_poni_path = os.path.join(self.poni_directory, self.poni_files['saxs'])
        self.waxs_poni_path = os.path.join(self.poni_directory, self.poni_files['waxs'])
        
        if not os.path.exists(self.saxs_poni_path):
            raise FileNotFoundError(f"SAXS PONI file not found: {self.saxs_poni_path}")
        if not os.path.exists(self.waxs_poni_path):
            raise FileNotFoundError(f"WAXS PONI file not found: {self.waxs_poni_path}")
            
        self.ai_saxs = pyFAI.load(self.saxs_poni_path)
        self.ai_waxs = pyFAI.load(self.waxs_poni_path)
        
        # Load detector masks
        self._load_detector_masks()
    
    def _load_detector_masks(self):
        """Load detector mask files."""
        saxs_mask_path = os.path.join(self.poni_directory, self.mask_files['saxs'])
        waxs_mask_path = os.path.join(self.poni_directory, self.mask_files['waxs'])
        
        self.saxs_mask = fabio.open(saxs_mask_path).data
        self.waxs_mask = fabio.open(waxs_mask_path).data  

    def read_raw_detector_file(self, raw_file_path: str, detector_type: str) -> np.ndarray:
        """
        Read raw detector file and return detector image array.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the .raw detector file
        detector_type : str
            Type of detector ('SAXS' or 'WAXS')
            
        Returns
        -------
        np.ndarray
            Detector image array with proper shape
        """
        if detector_type.upper() == 'SAXS':
            shape = tuple(self.detector_shapes['saxs'])
        elif detector_type.upper() == 'WAXS':
            shape = tuple(self.detector_shapes['waxs'])
        else:
            raise ValueError("detector_type must be 'SAXS' or 'WAXS'")
        
        if not os.path.exists(raw_file_path):
            raise FileNotFoundError(f"Raw detector file not found: {raw_file_path}")
        
        data = np.fromfile(raw_file_path, dtype=np.int32).reshape(shape)
        return data

    def calculate_sld_mu_thickness(self, transmission: float) -> Dict[str, float]:
        """
        Calculate material properties from transmission using experiment parameters.
        
        Parameters
        ----------
        transmission : float
            Measured transmission value
            
        Returns
        -------
        Dict[str, float]
            Dictionary containing mu, sld, and thickness_m
        """
        if transmission <= 0 or transmission > 1:
            raise ValueError("Transmission must be in the range (0, 1]")
        
        # Use xraydb for absorption coefficient calculation
        mu_xraydb = xraydb.material_mu(self.compound, 
                                     energy=self.energy_keV * 1000, 
                                     density=self.density_g_cm3) * 100  # Convert to 1/m
        
        # Calculate scattering length density (SLD) using periodictable
        material = pt.formula(self.compound)
        sld, mu_pt = pt.xray_sld(material, energy=self.energy_keV, density=self.density_g_cm3)

        
        # Calculate thickness: T = exp(-mu * t) → t = -ln(T) / mu
        thickness = -np.log(transmission) / mu_xraydb
        
        return {
            "mu": mu_xraydb, 
            "sld": sld, 
            "thickness_m": thickness
        }

    def create_output_directory(self, raw_file_path: str, detector_type: str) -> str:
        """
        Create output directory structure with SAXS/Reduction or WAXS/Reduction subdirectories.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the raw file being processed
        detector_type : str
            Type of detector ('SAXS' or 'WAXS')
            
        Returns
        -------
        str
            Output directory path
        """
        # Create flat structure: run5_test/1D/SAXS/Reduction or run5_test/1D/WAXS/Reduction
        output_dir = os.path.join(self.output_directory_1d, detector_type.upper(), "Reduction")
        os.makedirs(output_dir, exist_ok=True)
        return output_dir
    
    def get_corrections_full(self, raw_file_path: str):
        """Processed metadata file and computes correction factors"""
        metadata_function = getattr(utils, self.metadata_function)
        i0, bstop, csv_path = metadata_function(raw_file_path)
        
        # Calculate correction factors
        transmission_factor_raw = bstop
        transmission_ratio = bstop / i0
        
        # Calculate material properties if thickness not provided
        if self.thickness is None:
            material_props = self.calculate_sld_mu_thickness(transmission_ratio)
            thickness = material_props['thickness_m']
        else:
            thickness = self.thickness
        
        # Calculate normalization factor: trans_factor * i0_corrected
        # Based on Step1 notebook: normfactor = trans_factor*i0_avg_all[0]
        normalization_factor = transmission_factor_raw * i0
        
        return {
            'i0_corrected': i0,
            'bstop_corrected': bstop,
            'transmission_factor': transmission_factor_raw,
            'transmission_ratio': transmission_ratio,
            'thickness': thickness,
            'normalization_factor': normalization_factor,
            'metadata_path': csv_path
        }
        
    def process_saxs_file(self, raw_file_path: str):
        """
        Process a single SAXS .raw file with corrections and 1D integration.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the SAXS .raw file
        """
        print(f"Processing SAXS file: {raw_file_path}")
        
        # Read raw detector data
        detector_data = self.read_raw_detector_file(raw_file_path, 'SAXS')
        
        # csv_file_path = raw_file_path.replace(".raw", ".csv")
        # assert os.path.exists(csv_file_path), "Metadata file not found: {csv_file_path}"
        # parameters = self.read_csv_parameters(csv_file_path)
        
        corrections = self.get_corrections_full(raw_file_path)
        
        print(f"  I0: {corrections['i0_corrected']:.3f}, "
              f"Bstop: {corrections['bstop_corrected']:.3f}, "
              f"Trans_factor: {corrections['transmission_factor']:.3f}, "
              f"Trans_ratio: {corrections['transmission_ratio']:.4f}, "
              f"Normalization: {corrections['normalization_factor']:.3f}")
        
        # Create output directory and filename
        output_dir = self.create_output_directory(raw_file_path, 'SAXS')
        output_filename = os.path.basename(raw_file_path).replace('.raw', '.dat')
        output_path = os.path.join(output_dir, output_filename)
        
        # Perform 1D radial integration with automatic file writing
        # PyFAI will automatically generate headers with detector config and column names
        q, intensity, error = self.ai_saxs.integrate1d(
            detector_data, 
            self.npt_radial,
            error_model=self.error_model,
            mask=self.saxs_mask,
            normalization_factor=corrections['normalization_factor'],
            filename=output_path
        )
        
        print(f"  Saved to: {output_path}")

    def process_waxs_file(self, raw_file_path: str):
        """
        Process a single WAXS .raw file with corrections and 1D integration.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the WAXS .raw file
        """
        print(f"Processing WAXS file: {raw_file_path}")
        
        detector_data = self.read_raw_detector_file(raw_file_path, 'WAXS')
        corrections = self.get_corrections_full(raw_file_path)

        # Create output directory and filename
        output_dir = self.create_output_directory(raw_file_path, 'WAXS')
        output_filename = os.path.basename(raw_file_path).replace('.raw', '.dat')
        output_path = os.path.join(output_dir, output_filename)
        
        # Perform 1D radial integration with automatic file writing
        # PyFAI will automatically generate headers with detector config and column names
        q, intensity, error = self.ai_waxs.integrate1d(
            detector_data,
            self.npt_radial,
            error_model=self.error_model,
            mask=self.waxs_mask,
            normalization_factor=corrections['normalization_factor'],
            filename=output_path
        )
        
        print(f"  Saved to: {output_path}")


def find_all_raw_files(data_directory_path: str) -> Tuple[List[str], List[str]]:
    """
    Find all .raw files in SAXS and WAXS directories.
    
    Parameters
    ----------
    data_directory_path : str
        Full path to the data directory (e.g., "run12_test/2D/")
        
    Returns
    -------
    Tuple[List[str], List[str]]
        Lists of SAXS and WAXS .raw file paths
    """
    # Find all SAXS .raw files
    saxs_pattern = os.path.join(data_directory_path, "**/SAXS/**/*.raw")
    saxs_files = [f for f in glob.glob(saxs_pattern, recursive=True) 
                  if not f.endswith('.raw.pdi')]
    
    # Find all WAXS .raw files  
    waxs_pattern = os.path.join(data_directory_path, "**/WAXS/**/*.raw")
    waxs_files = [f for f in glob.glob(waxs_pattern, recursive=True)
                  if not f.endswith('.raw.pdi')]
    
    return sorted(saxs_files), sorted(waxs_files)


def main():
    """
    Main function to process all .raw files in the dataset.
    """
    print("SWAXS Data Correction and 1D Integration Demo")
    print("=" * 50)
    
    # Initialize experiment from config
    experiment = Experiment("config.yml")
    print(f"Loaded configuration for compound: {experiment.compound}")
    print(f"Energy: {experiment.energy_keV} keV, Density: {experiment.density_g_cm3} g/cm³")
    print()
    
    # Find all .raw files to process
    saxs_files, waxs_files = find_all_raw_files(experiment.data_directory_2d)
    
    print(f"Found {len(saxs_files)} SAXS files and {len(waxs_files)} WAXS files to process")
    print()
    
    # Process all SAXS files
    if saxs_files:
        print("Processing SAXS files...")
        for saxs_file in saxs_files:
            experiment.process_saxs_file(saxs_file)
        print()
    
    # Process all WAXS files
    if waxs_files:
        print("Processing WAXS files...")
        for waxs_file in waxs_files:
            experiment.process_waxs_file(waxs_file)

        print()
    
    print("Processing complete!")
    
    return 0


if __name__ == "__main__":
    exit(main())
    
# Previous functions for experiment class
    # def get_counters_from_pdi(self, pdi_file: str):
    #     """Get counter names and values from PDI file
    #     Args:
    #         pdi_file (str): PDI file name with path
    #     Returns:
    #         dict: Dictionary containing counter names and values
    #     """
    #     print(f'pdi_file: {pdi_file}')
    #     with open(pdi_file, 'r') as f:
    #         data = f.read()
    #     if 'All Counters' not in data: 
    #         # Likely an empty file
    #         raise RuntimeError("Empty PDI not supported yet")
    #         return None
    #     data = data.replace('\n', ';')
        
    #     counters_match = re.search('All Counters;(.*);;# All Motors', data)
    #     assert counters_match is not None, f"Failed to parse PDI file {pdi_file}"
    #     counters = counters_match.group(1)
    #     cts = re.split(';|=', counters)
    #     counters = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
    #     return {
    #         'i0': counters['i0'],
    #         'bstop': counters['bstop'],
    #         'ctemp': counters['CTEMP'],
    #         'timer': counters['Timer']
    #     }
    
    # def read_csv_parameters(self, csv_file_path: str) -> Dict[str, float]:
    #     """
    #     Read and extract parameters from CSV metadata files.
        
    #     Parameters
    #     ----------
    #     csv_file_path : str
    #         Path to the CSV parameter file
            
    #     Returns
    #     -------
    #     Dict[str, float]
    #         Dictionary containing averaged parameters
    #     """
    #     if not os.path.exists(csv_file_path):
    #         raise FileNotFoundError(f"CSV file not found: {csv_file_path}")
            
    #     df = pd.read_csv(csv_file_path)
        
    #     # Extract columns as in Step1 notebook [2,3,6,10,29,30] 
    #     # Column indices: i0=3, bstop=6, ctemp=29, timer=30 (0-indexed)
    #     i0 = pd.to_numeric(df.iloc[:, 3], errors='coerce')
    #     bstop = pd.to_numeric(df.iloc[:, 6], errors='coerce')
    #     ctemp = pd.to_numeric(df.iloc[:, 29], errors='coerce')
    #     timer = pd.to_numeric(df.iloc[:, 30], errors='coerce')
        
    #     return {
    #         'i0': i0.mean(skipna=True),
    #         'bstop': bstop.mean(skipna=True),
    #         'ctemp': ctemp.mean(skipna=True),
    #         'timer': timer.mean(skipna=True)
    #     }
