#!/usr/bin/env python3
"""
SWAXS Data Correction and 1D Integration Demo
==============================================
Simplified version of Step1 notebook focusing on core data correction and 1D integration.
Processes individual .raw files without averaging or complex file searching logic.
All parameters read from config.yml file.

Key Features:
- Experiment class for configuration management
- Individual file processing (SAXS and WAXS)
- Essential data corrections (dark current, transmission, normalization)
- 1D radial integration with PyFAI
- Output in same format as existing .dat files
"""
import utils.py
import os
import glob
import numpy as np
import pandas as pd
import yaml
import fabio
import pyFAI
# import periodictable as pt  # Not needed for simplified demo
import xraydb
from typing import Dict, List, Tuple


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
        self._setup_detector_parameters()
        self._setup_experimental_parameters()
        self._setup_integration_parameters()
        self._load_pyfai_integrators()
        
    def _load_config(self) -> Dict:
        """Load configuration from YAML file."""
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _setup_directories(self):
        """Setup directory paths from config."""
        self.poni_directory = self.config['poni_directory']
        self.data_directory_2d = self.config['data_directory']  # Read from config
        self.output_directory_1d = "1D"
        
        self.saxs_subdir = "SAXS" 
        self.waxs_subdir = "WAXS"
        
    def _setup_detector_parameters(self):
        """Setup detector-specific parameters from config."""
        self.detector_shapes = self.config['detector_shapes']
        self.poni_files = self.config['poni_files']
        self.mask_files = self.config['mask_files']
        
    def _setup_experimental_parameters(self):
        """Setup experimental parameters from config."""
        self.compound = self.config['compound']
        self.energy_keV = self.config['energy_keV']
        self.density_g_cm3 = self.config['density_g_cm3']
        self.i0_offset = self.config['i0_offset']
        self.bstop_offset = self.config['bstop_offset']
        self.i0_air = self.config['i0_air']
        self.bstop_air = self.config['bstop_air']
        self.thickness = self.config.get('thickness', None)
        
    def _setup_integration_parameters(self):
        """Setup 1D integration parameters from config."""
        self.npt_radial = self.config['npt_radial']
        self.error_model = self.config['error_model']
    
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
        
        if os.path.exists(saxs_mask_path):
            self.saxs_mask = fabio.open(saxs_mask_path).data
        else:
            self.saxs_mask = None
            print(f"Warning: SAXS mask file not found: {saxs_mask_path}")
            
        if os.path.exists(waxs_mask_path):
            self.waxs_mask = fabio.open(waxs_mask_path).data  
        else:
            self.waxs_mask = None
            print(f"Warning: WAXS mask file not found: {waxs_mask_path}")

    def read_csv_parameters(self, csv_file_path: str) -> Dict[str, float]:
        """
        Read and extract parameters from CSV metadata files.
        
        Parameters
        ----------
        csv_file_path : str
            Path to the CSV parameter file
            
        Returns
        -------
        Dict[str, float]
            Dictionary containing averaged parameters
        """
        if not os.path.exists(csv_file_path):
            raise FileNotFoundError(f"CSV file not found: {csv_file_path}")
            
        df = pd.read_csv(csv_file_path)
        
        # Extract columns as in Step1 notebook [2,3,6,10,29,30] 
        # Column indices: i0=3, bstop=6, ctemp=29, timer=30 (0-indexed)
        i0 = pd.to_numeric(df.iloc[:, 3], errors='coerce')
        bstop = pd.to_numeric(df.iloc[:, 6], errors='coerce')
        ctemp = pd.to_numeric(df.iloc[:, 29], errors='coerce')
        timer = pd.to_numeric(df.iloc[:, 30], errors='coerce')
        
        return {
            'i0_avg': i0.mean(skipna=True),
            'bstop_avg': bstop.mean(skipna=True),
            'ctemp_avg': ctemp.mean(skipna=True),
            'timer_avg': timer.mean(skipna=True)
        }

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
        
        # Read raw binary file as int32 and reshape
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
        
        # Note: SLD calculation using periodictable - for simplified demo we'll use a default value
        sld = 2.0e-6  # Default SLD value for polymer materials
        
        # Calculate thickness: T = exp(-mu * t) → t = -ln(T) / mu
        thickness = -np.log(transmission) / mu_xraydb
        
        return {
            "mu": mu_xraydb, 
            "sld": sld, 
            "thickness_m": thickness
        }

    def calculate_correction_factors(self, parameters: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate transmission and normalization factors from raw parameters.
        IMPORTANT: Based on Step1 notebook analysis, offsets are set to 0!
        
        Parameters
        ----------
        parameters : Dict[str, float]
            Raw parameters from CSV file
            
        Returns
        -------
        Dict[str, float]
            Dictionary containing correction factors
        """
        # CRITICAL FIX: Step1 notebook uses offset=0, so NO dark current correction
        i0_corrected = parameters['i0_avg']  # - 0  (no offset applied)
        bstop_corrected = parameters['bstop_avg']  # - 0  (no offset applied)
        
        # Calculate transmission factor - this is the raw beamstop value for SAXS
        # Based on Step1 notebook: trans_factor = (bstop_avg_all[0]) for SAXS
        transmission_factor_raw = bstop_corrected  # Raw bstop value
        
        # For thickness calculation, need transmission ratio
        transmission_ratio = bstop_corrected / i0_corrected
        
        # Calculate material properties if thickness not provided
        if self.thickness is None:
            material_props = self.calculate_sld_mu_thickness(transmission_ratio)
            thickness = material_props['thickness_m']
        else:
            thickness = self.thickness
        
        # Calculate normalization factor: trans_factor * i0_corrected
        # Based on Step1 notebook: normfactor = trans_factor*i0_avg_all[0]
        normalization_factor = transmission_factor_raw * i0_corrected
        
        return {
            'i0_corrected': i0_corrected,
            'bstop_corrected': bstop_corrected,
            'transmission_factor': transmission_factor_raw,
            'transmission_ratio': transmission_ratio,
            'thickness': thickness,
            'normalization_factor': normalization_factor
        }

    def create_output_directory(self, raw_file_path: str) -> str:
        """
        Create output directory structure matching input file location.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the raw file being processed
            
        Returns
        -------
        str
            Output directory path
        """
        # Extract directory structure relative to 2D/
        rel_path = os.path.relpath(raw_file_path, self.data_directory_2d)
        dir_parts = os.path.dirname(rel_path).split(os.sep)
        
        # Remove SAXS/WAXS subdirectory from path
        if len(dir_parts) > 0 and dir_parts[-1] in ['SAXS', 'WAXS']:
            dir_parts = dir_parts[:-1]
            
        output_dir = os.path.join(self.output_directory_1d, *dir_parts)
        os.makedirs(output_dir, exist_ok=True)
        return output_dir

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
        
        # Find corresponding CSV file (same directory, same name, .csv extension)
        csv_file_path = raw_file_path.replace('.raw', '.csv')
        
        if not os.path.exists(csv_file_path):
            print(f"Warning: CSV file not found for {raw_file_path}, skipping")
            return
        
        # Read parameters from CSV
        parameters = self.read_csv_parameters(csv_file_path)
        
        # Calculate correction factors
        corrections = self.calculate_correction_factors(parameters)
        
        print(f"  I0: {corrections['i0_corrected']:.3f}, "
              f"Bstop: {corrections['bstop_corrected']:.3f}, "
              f"Trans_factor: {corrections['transmission_factor']:.3f}, "
              f"Trans_ratio: {corrections['transmission_ratio']:.4f}, "
              f"Normalization: {corrections['normalization_factor']:.3f}")
        
        # Perform 1D radial integration
        q, intensity, error = self.ai_saxs.integrate1d(
            detector_data, 
            self.npt_radial,
            error_model=self.error_model,
            mask=self.saxs_mask,
            normalization_factor=corrections['normalization_factor']
        )
        
        # Create output directory and filename
        output_dir = self.create_output_directory(raw_file_path)
        output_filename = os.path.basename(raw_file_path).replace('.raw', '.dat')
        output_path = os.path.join(output_dir, output_filename)
        
        # Save data in 3-column format (q, intensity, error)
        data_array = np.column_stack([q, intensity, error])
        np.savetxt(output_path, data_array, fmt='%.6e', delimiter='    ')
        
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
        
        # Read raw detector data
        detector_data = self.read_raw_detector_file(raw_file_path, 'WAXS')
        
        # Find corresponding CSV file (same directory, same name, .csv extension)
        csv_file_path = raw_file_path.replace('.raw', '.csv')
        
        if not os.path.exists(csv_file_path):
            print(f"Warning: CSV file not found for {raw_file_path}, skipping")
            return
        
        # Read parameters from CSV
        parameters = self.read_csv_parameters(csv_file_path)
        
        # Calculate correction factors
        corrections = self.calculate_correction_factors(parameters)
        
        print(f"  I0: {corrections['i0_corrected']:.3f}, "
              f"Bstop: {corrections['bstop_corrected']:.3f}, "
              f"Trans_factor: {corrections['transmission_factor']:.3f}, "
              f"Trans_ratio: {corrections['transmission_ratio']:.4f}, "
              f"Normalization: {corrections['normalization_factor']:.3f}")
        
        # Perform 1D radial integration
        q, intensity, error = self.ai_waxs.integrate1d(
            detector_data,
            self.npt_radial,
            error_model=self.error_model,
            mask=self.waxs_mask,
            normalization_factor=corrections['normalization_factor']
        )
        
        # Create output directory and filename
        output_dir = self.create_output_directory(raw_file_path)
        output_filename = os.path.basename(raw_file_path).replace('.raw', '.dat')
        output_path = os.path.join(output_dir, output_filename)
        
        # Save data in 3-column format (q, intensity, error)
        data_array = np.column_stack([q, intensity, error])
        np.savetxt(output_path, data_array, fmt='%.6e', delimiter='    ')
        
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
    
Experiment.