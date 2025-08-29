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

Logging rules:
1. Logs all config information just once — when Experiment._load_and_assign_config() is called
2. Also performs logging for all corrections called on the very first run of SAXS and the first run of WAXS (if both are run).

"""
# TODO: Split up into different files 
# TODO: Logging -- add PONI files, function by function
# TODO: Create averaged folder for SAXS and WAXS after doing correction and integration
# TODO: Plots: Use Anjani's function to create y offsets
# TODO: Remove prefix for SAXS and WAXS PDI files (sone_) 
# TODO: 1-d data should end with _SAXS and _WAXS
import os
import glob
import numpy as np
import yaml
import fabio
import pyFAI
import periodictable as pt
import xraydb
from typing import Dict, List, Tuple
import re
import utils
from pathlib import Path
import logging

# Start logger, overwrites existing logging files, includes timestamp in message
logging.basicConfig(filename = "corrections.log",
                    filemode = "w", 
                    level = logging.INFO,
                    format="%(asctime)s -  %(name)s - %(levelname)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
class Experiment:
    # Type hints for dynamically assigned config attributes
    data_directory: str
    poni_directory: str
    compound: str
    energy_keV: float
    density_g_cm3: float
    detector_shapes: Dict[str, List[int]]
    poni_files: Dict[str, str]
    mask_files: Dict[str, str]
    read_metadata_function: str
    thickness: float
    npt_radial: int
    error_model: str
    i0_offset: float
    bstop_offset: float
    i0_air: float
    bstop_air: float
    mode: str
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
        # Logging state tracking for first-run detailed logging
        self._saxs_log_pending = True
        self._waxs_log_pending = True
        
        self._load_and_assign_config()
        self._setup_directories()
        self._load_pyfai_integrators()
        
    def _load_and_assign_config(self):
        """Load configuration from YAML file and assign all values to instance variables."""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Log all configuration information once
        logging.info("=" * 60)
        logging.info("EXPERIMENT CONFIGURATION")
        logging.info("=" * 60)
        logging.info(f"Configuration file: {self.config_path}")
        logging.info("")
        
        # Loop through all config keys and assign them as instance variables
        for key, value in config.items():
            setattr(self, key, value)
            # Log each configuration parameter
            if isinstance(value, dict):
                logging.info(f"{key}:")
                for sub_key, sub_value in value.items():
                    logging.info(f"  {sub_key}: {sub_value}")
            elif isinstance(value, list):
                logging.info(f"{key}: {value}")
            else:
                logging.info(f"{key}: {value}")
        
        # Make mode upper-case (SAXS, SWAXS, WAXS) for case-insensitive string comparison
        self.mode = self.mode.upper()
        logging.info("")
        logging.info("=" * 60)
    
    def _setup_directories(self):
        """Setup computed directory paths."""
        self.data_directory_2d = os.path.join(self.data_directory, "2D")
        self.poni_directory = os.path.join(self.data_directory, self.poni_directory)
        self.output_directory_1d = os.path.join(self.data_directory, "1D")
        os.makedirs(self.output_directory_1d, exist_ok = True)
        
        self.saxs_subdir = "SAXS" 
        self.waxs_subdir = "WAXS"
    
    def _load_saxs_integrator(self):
        """Load PyFAI integrator and mask for SAXS detector."""
        self.saxs_poni_path = os.path.join(self.poni_directory, self.poni_files['saxs'])
        
        if not os.path.exists(self.saxs_poni_path):
            raise FileNotFoundError(f"SAXS PONI file not found: {self.saxs_poni_path}")
            
        self.ai_saxs = pyFAI.load(self.saxs_poni_path)
        logging.info(f"Loaded SAXS PyFAI integrator from: {self.saxs_poni_path}")
        
        # Load SAXS detector mask
        if self.mask_files['saxs'] is not None:
            saxs_mask_path = os.path.join(self.poni_directory, self.mask_files['saxs'])
            self.saxs_mask = fabio.open(saxs_mask_path).data
            logging.info(f"Loaded SAXS mask from: {saxs_mask_path}")
        else:
            self.saxs_mask = None
            logging.info("No SAXS mask file specified")
    
    def _load_waxs_integrator(self):
        """Load PyFAI integrator and mask for WAXS detector."""
        self.waxs_poni_path = os.path.join(self.poni_directory, self.poni_files['waxs'])
        
        if not os.path.exists(self.waxs_poni_path):
            raise FileNotFoundError(f"WAXS PONI file not found: {self.waxs_poni_path}")
            
        self.ai_waxs = pyFAI.load(self.waxs_poni_path)
        logging.info(f"Loaded WAXS PyFAI integrator from: {self.waxs_poni_path}")
        
        # Load WAXS detector mask
        if self.mask_files['waxs'] is not None:
            waxs_mask_path = os.path.join(self.poni_directory, self.mask_files['waxs'])
            self.waxs_mask = fabio.open(waxs_mask_path).data 
            logging.info(f"Loaded WAXS mask from: {waxs_mask_path}")
        else:
            self.waxs_mask = None
            logging.info("No WAXS mask file specified")

    def _load_pyfai_integrators(self):
        """Load PyFAI integrator objects based on processing mode."""
        
        if self.mode not in ['SAXS', 'WAXS', 'SWAXS']:
            raise ValueError(f"Invalid mode '{self.mode}'. Must be 'SAXS', 'WAXS', or 'SWAXS'")

        if self.mode == 'WAXS':
            self._load_waxs_integrator()
            self._saxs_log_pending = False
            
        elif self.mode == 'SAXS':
            self._load_saxs_integrator()
            self._waxs_log_pending = False

        else:
            self._load_saxs_integrator()
            self._load_waxs_integrator()


            

    

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
        # Log function call for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logging.info(f"Called read_raw_detector_file() for {detector_type} detector")
        
        if detector_type.upper() == 'SAXS':
            shape = tuple(self.detector_shapes['saxs'])
        elif detector_type.upper() == 'WAXS':
            shape = tuple(self.detector_shapes['waxs'])
        else:
            raise ValueError("detector_type must be 'SAXS' or 'WAXS'")
        
        if not os.path.exists(raw_file_path):
            raise FileNotFoundError(f"Raw detector file not found: {raw_file_path}")
        
        data = np.fromfile(raw_file_path, dtype=np.int32).reshape(shape)
        
        # Log data information for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logging.info(f"  Raw data loaded: shape {data.shape}, dtype {data.dtype}")
            logging.info(f"  Data range: min={data.min()}, max={data.max()}, mean={data.mean():.2f}")
        
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
        if (self._saxs_log_pending and self._waxs_log_pending):
            logging.info(f"Called calculate_sld_mu_thickness() with transmission = {transmission:.6f}")
        
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
        output_dir = os.path.join(self.output_directory_1d, detector_type.upper(), "Reduction")
        os.makedirs(output_dir, exist_ok=True)
        return output_dir
    
    def get_saxs_pdi_from_waxs(self, pdi_path) -> str:
        """Gets SAXS PDI from WAXS PDI. Also would work in vice-versa, but not used"""
        corresponding_pdi_regex = r"(Run\d{1,3}.*scan\d_\d\d\d\d)"
        match = re.search(corresponding_pdi_regex, pdi_path)
        
        if match is None: 
            raise RuntimeError("")
        shared_expression = match.group(1) 
        # Convert to pathlib to change directories
        search_directory = Path(pdi_path).parent.parent / "SAXS"
        for pdi_file in search_directory.glob("*.pdi"):
            if shared_expression in pdi_file.name:
                return str(pdi_file)
        # No file found
        raise FileNotFoundError(f"""Empty WAXS File: {shared_expression} and no 
                        corresponding SAXS found in directory {search_directory}""")
        

    def process_pdi_full(self, raw_file_path: str, detector_type: str):
        pdi_file_path = raw_file_path + ".pdi"
        if not Path(pdi_file_path).exists():
            raise FileNotFoundError(f"PDI file not found: {pdi_file_path}") 
        with open(pdi_file_path, 'r') as f:
            data = f.read()
    
        if (detector_type.upper() != "SAXS") and (detector_type.upper() != "WAXS"):
            raise RuntimeError("Detector Type must be either SAXS or WAXS")
        if ('All Counters' not in data) and (detector_type.upper() == "WAXS"):
            # Empty PDI file found in WAXS -- check SAXS to see if corresponding file exists
            pdi_file_path = self.get_saxs_pdi_from_waxs(pdi_file_path)
        counters, motors, extras = self.get_meta_from_pdi(pdi_file_path)
        if counters is None:
            raise RuntimeError(f"Error Parsing PDI File: {pdi_file_path} -- counters not found")
        return counters['i0'], counters['bstop'], pdi_file_path
          
    def get_meta_from_pdi(self, pdi_file: str):
        """Get motor and counter names and values from PDI file
        Function returns empty dictionary for counters if it cannot process successfully
        """
        with open(pdi_file, 'r') as f:
            data = f.read()
        data = data.replace('\n', ';')
        try:
            counters = re.search('All Counters;(.*);;# All Motors', data).group(1)
            cts = re.split(';|=', counters)
            Counters = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
            motors = re.search('All Motors;(.*);#', data).group(1)
            cts = re.split(';|=', motors)
            Motors = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
        except AttributeError:
            ss1 = '# Diffractometer Motor Positions for image;# '
            ss2 = ';# Calculated Detector Calibration Parameters for image:'
            motors = re.search(f'{ss1}(.*){ss2}', data).group(1)
            cts = re.split(';|=', motors)
            Motors = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
            Motors['TwoTheta'] = Motors['2Theta']
            Counters = {}
        Extras = {}
        if len(data[data.rindex(';') + 1:]) > 0:
            Extras['epoch'] = data[data.rindex(';') + 1:]
        return Counters, Motors, Extras


    def get_corrections_full(self, raw_file_path: str, detector_type: str):
        """Processed metadata file and computes correction factors"""
        # Log function call for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logging.info(f"Called get_corrections_full() for {detector_type} detector")
        
        # Use specified metadata function from utils if available, otherwise use PDI processing
        if self.read_metadata_function is not None and hasattr(utils, self.read_metadata_function):
            metadata_function = getattr(utils, self.read_metadata_function)
            i0, bstop, csv_path = metadata_function(raw_file_path, detector_type)
        else:
            # Using current PDI metadata processing
            i0, bstop, csv_path = self.process_pdi_full(raw_file_path, detector_type)
        
        # Log metadata processing for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logging.info(f"  Processing metadata from: {csv_path}")
            logging.info(f"  Raw I0: {i0:.6f}, Raw Bstop: {bstop:.6f}")
        
        # Calculate correction factors
        transmission_factor_raw = bstop
        transmission_ratio = bstop / i0
        
        # Log transmission calculations for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logging.info(f"  Calculating transmission: ratio = {transmission_ratio:.6f}")
        
        # Calculate material properties if thickness not provided
        if self.thickness is None:
            
            material_props = self.calculate_sld_mu_thickness(transmission_ratio)
            thickness = material_props['thickness_m']
            
            # Log material calculations for first runs only
            if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
                (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
                logging.info(f"  Material properties calculated:")
                logging.info(f"    Compound: {self.compound}")
                logging.info(f"    Energy: {self.energy_keV} keV")
                logging.info(f"    Density: {self.density_g_cm3} g/cm³")
                logging.info(f"    Absorption coefficient: {material_props['mu']:.6f} 1/m")
                logging.info(f"    SLD: {material_props['sld']:.6e}")
                logging.info(f"    Calculated thickness: {thickness:.6f} m")
        else:
            thickness = self.thickness
            if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
                (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
                logging.info(f"  Using provided thickness: {thickness:.6f} m")
        
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
            
        Returns
        -------
        Dict
            Dictionary containing plotting data: q, intensity, error, filename, corrections
        """
                
        # Log detailed information for first SAXS file only
        if self._saxs_log_pending:
            logging.info("")
            logging.info("=" * 60)
            logging.info("FIRST SAXS FILE PROCESSING")
            logging.info("=" * 60)
            logging.info(f"Raw file: {raw_file_path}")
            logging.info(f"SAXS PONI file: {self.saxs_poni_path}")
            logging.info(f"SAXS mask file: {os.path.join(self.poni_directory, self.mask_files['saxs']) if self.mask_files['saxs'] else 'None'}")
            logging.info("")
        
        # Read raw detector data
        detector_data = self.read_raw_detector_file(raw_file_path, 'SAXS')
        
        # csv_file_path = raw_file_path.replace(".raw", ".csv")
        # assert os.path.exists(csv_file_path), "Metadata file not found: {csv_file_path}"
        # parameters = self.read_csv_parameters(csv_file_path)
        
        corrections = self.get_corrections_full(raw_file_path, "SAXS")
        
        # Log corrections for first SAXS file only
        if self._saxs_log_pending:
            logging.info(f"  Metadata file: {corrections['metadata_path']}")

            self._saxs_log_pending = False
        
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
                
        # Return plotting data
        return {
            'q': q,
            'intensity': intensity,
            'error': error,
            'filename': output_filename,
            'corrections': corrections,
            'raw_file_path': raw_file_path
        }
    def process_waxs_file(self, raw_file_path: str):
        """
        Process a single WAXS .raw file with corrections and 1D integration.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the WAXS .raw file
            
        Returns
        -------
        Dict
            Dictionary containing plotting data: q, intensity, error, filename, corrections
        """        
        print(f"Processing WAXS file: {raw_file_path}")
        
        # Log detailed information for first WAXS file only
        if self._waxs_log_pending:
            logging.info("")
            logging.info("=" * 60)
            logging.info("FIRST WAXS FILE PROCESSING")
            logging.info("=" * 60)
            logging.info(f"Raw file: {raw_file_path}")
            logging.info(f"WAXS PONI file: {self.waxs_poni_path}")
            logging.info(f"WAXS mask file: {os.path.join(self.poni_directory, self.mask_files['waxs']) if self.mask_files['waxs'] else 'None'}")
            logging.info(f"Detector shape: {self.detector_shapes['waxs']}")
            logging.info("")
        
        detector_data = self.read_raw_detector_file(raw_file_path, 'WAXS')
        corrections = self.get_corrections_full(raw_file_path, "WAXS")

        # Log corrections for first WAXS file only
        if self._waxs_log_pending:
            logging.info(f"  Metadata file: {corrections['metadata_path']}")
            logging.info("=" * 60)
            self._waxs_log_pending = False

        # Create output directory and filename
        output_dir = self.create_output_directory(raw_file_path, 'WAXS')
        output_filename = os.path.basename(raw_file_path).replace('.raw', '.dat')
        output_path = os.path.join(output_dir, output_filename)
        
        # Perform D radial integration with automatic file writing
        # PyFAI will automatically generate headers with detector config and column names
        q, intensity, error = self.ai_waxs.integrate1d(
            detector_data,
            self.npt_radial,
            error_model=self.error_model,
            mask=self.waxs_mask,
            normalization_factor=corrections['normalization_factor'],
            filename=output_path
        )
                
        # Return plotting data
        return {
            'q': q,
            'intensity': intensity,
            'error': error,
            'filename': output_filename,
            'corrections': corrections,
            'raw_file_path': raw_file_path
        }


def find_all_raw_files(experiment: Experiment, data_directory_path: str) -> Tuple[List[str], List[str]]:
    
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
    if (experiment.mode != "WAXS"):
        saxs_pattern = os.path.join(data_directory_path, "**/SAXS/**/*.raw")
        saxs_files = [f for f in glob.glob(saxs_pattern, recursive=True) 
                    if not f.endswith('.raw.pdi')]
    else:
        # Don't include SAXS files in just WAXS
        saxs_files = []
    
    if (experiment.mode != "SAXS"):
        waxs_pattern = os.path.join(data_directory_path, "**/WAXS/**/*.raw")
        waxs_files = [f for f in glob.glob(waxs_pattern, recursive=True)
                    if not f.endswith('.raw.pdi')]
    else:
        # Don't include WAXS files in just SAXS
        waxs_files = []
        
    return sorted(saxs_files), sorted(waxs_files)

def full_correction_integration(plotting=False):
    """
    Main function to process all .raw files in the dataset.
    
    Parameters
    ----------
    plotting : bool, default False
        If True, return plotting data instead of just processing files
        
    Returns
    -------
    Dict or None
        If plotting=True, returns dictionary with plotting data.
        If plotting=False, returns None (just processes files).
    """
    print("SWAXS Data Correction and 1D Integration Demo")
    print("=" * 50)
    
    # Initialize experiment from config
    experiment = Experiment("config.yml")
    print(f"Loaded configuration for compound: {experiment.compound}")
    print(f"Energy: {experiment.energy_keV} keV, Density: {experiment.density_g_cm3} g/cm³")
    print()
    
    # Find all .raw files to process
    saxs_files, waxs_files = find_all_raw_files(experiment, experiment.data_directory_2d)
    
    print(f"Found {len(saxs_files)} SAXS files and {len(waxs_files)} WAXS files to process")
    print()
    
    # Initialize plotting data collection if requested
    if plotting:
        plotting_data = {
            'saxs_data': [],
            'waxs_data': [],
            'metadata': {
                'timers': [],
                'temperatures': [], 
                'i0_values': [],
                'bstop_values': [],
                'filenames': []
            }
        }
    
    # Process all SAXS files
    if saxs_files:
        print("Processing SAXS files...")
        for saxs_file in saxs_files:
            saxs_result = experiment.process_saxs_file(saxs_file)
            if plotting:
                plotting_data['saxs_data'].append(saxs_result)
                # Extract metadata for monitoring plots
                corrections = saxs_result['corrections']
                plotting_data['metadata']['i0_values'].append(corrections['i0_corrected'])
                plotting_data['metadata']['bstop_values'].append(corrections['bstop_corrected'])
                plotting_data['metadata']['filenames'].append(saxs_result['filename'])
        print()
    
    # Process all WAXS files
    if waxs_files:
        print("Processing WAXS files...")
        for waxs_file in waxs_files:
            waxs_result = experiment.process_waxs_file(waxs_file)
            if plotting:
                plotting_data['waxs_data'].append(waxs_result)
        print()
    
    print("Processing complete!")
    
    # Return plotting data if requested
    if plotting:
        return plotting_data
    else:
        return None
    
def main():
    return full_correction_integration()

if __name__ == "__main__":
    main()
