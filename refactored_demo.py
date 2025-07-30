#!/usr/bin/env python3
"""
Refactored SWAXS Data Processing Pipeline
==========================================
Objective:
Refactor the Step1_SWAXS_Reduction notebook and keep all corrections and azimuthal
1-D integration. Simplify complex file_path handling.
Processes a single run of data as a test-case. Read information
provided from config.yml such as data directory, poni file name, compound.
"""

import numpy as np
import os
import glob
import pandas as pd
import yaml
import fabio
import pyFAI
from typing import Dict, Optional, Tuple, List
import periodictable as pt
import xraydb




def read_csv_parameters(csv_file_path: str) -> Dict[str, float]:
    """
    Read and extract parameters from CSV files.
    
    Parameters
    ----------
    csv_file_path : str
        Path to the CSV parameter file
        
    Returns
    -------
    Dict[str, float]
        Dictionary containing averaged parameters:
        - 'i0_avg': Average I0 reading
        - 'bstop_avg': Average beamstop reading  
        - 'ctemp_avg': Average temperature reading
        - 'timer_avg': Average timer reading
    """
    df = pd.read_csv(csv_file_path)
    
    # Extract columns [2,3,6,10,29,30] as in notebook
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


def read_raw_detector_file(raw_file_path: str, detector_type: str) -> np.ndarray:
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
    # Load config to get detector shapes
    with open("config.yml", 'r') as f:
        config = yaml.safe_load(f)
    
    if detector_type.upper() == 'SAXS':
        shape = tuple(config['detector_shapes']['saxs'])
    elif detector_type.upper() == 'WAXS':
        shape = tuple(config['detector_shapes']['waxs'])
    else:
        raise ValueError("detector_type must be 'SAXS' or 'WAXS'")
    
    # Read raw binary file as int32 and reshape
    data = np.fromfile(raw_file_path, dtype=np.int32).reshape(shape)
    
    return data


def find_detector_files(data_directory: str, detector_type: str) -> List[str]:
    """
    Find detector files in the specified directory.
    
    Parameters
    ---------- 
    data_directory : str
        Path to the data directory
    detector_type : str
        Type of detector ('SAXS' or 'WAXS')
        
    Returns
    -------
    List[str]
        List of detector file paths
    """
    detector_dir = os.path.join(data_directory, detector_type.upper())
    raw_pattern = os.path.join(detector_dir, '*.raw')
    
    # Find all .raw files, excluding .raw.pdi files
    raw_files = [f for f in glob.glob(raw_pattern) if not f.endswith('.raw.pdi')]
    
    return sorted(raw_files)


def calculate_duration(timer_values: List[float]) -> List[float]:
    """
    Calculate relative duration in hours from absolute timer values.
    
    Parameters
    ----------
    timer_values : List[float]
        List of absolute timer values
        
    Returns
    -------
    List[float]
        List of relative durations in hours, starting from 0
    """
    if not timer_values:
        return []
    
    min_timer = min(timer_values)
    duration = [(timer - min_timer) / 3600 for timer in timer_values]
    
    return duration




def calculate_sld_mu_thickness(energy_keV: float, density: float, transmission: float, compound: str):
    """
    Calculate scattering length density (SLD), absorption coefficient (mu), and thickness
    from transmission data.
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'mu': absorption coefficient in 1/m
        - 'sld': scattering length density
        - 'thickness_m': sample thickness in meters
    """
    # Define PE composition correctly
    material = pt.formula(compound)

    # Calculate scattering length density (SLD) using periodictable
    # sld, _ = material.xray_sld(energy=energy_keV, density=density)
    sld, _ = material.xray_sld(energy=energy_keV)

    # Use xraydb for a more reliable absorption coefficient (mu) calculation
    mu_xraydb = xraydb.material_mu(compound, energy=energy_keV * 1000, density=density) * 100  # Convert to 1/m

    # Calculate thickness using transmission: T = exp(-mu * t) → t = -ln(T) / mu
    if transmission <= 0 or transmission > 1:
        raise ValueError("Transmission must be in the range (0, 1].")
    
    thickness = -np.log(transmission) / mu_xraydb  # Thickness in meters

    return {"mu": mu_xraydb, "sld": sld, "thickness_m": thickness}


def correct_detector_offsets(i0_raw: float, bstop_raw: float, 
                           i0_offset: float, bstop_offset: float) -> Tuple[float, float]:
    """
    Apply dark current corrections to detector readings.
    """
    i0_corrected = i0_raw - i0_offset
    bstop_corrected = bstop_raw - bstop_offset
    
    return i0_corrected, bstop_corrected


def calculate_transmission_factor(bstop_corrected: float, i0_corrected: float, 
                                bstop_air: Optional[float] = None, i0_air: Optional[float] = None) -> float:
    """
    Calculate transmission factor from corrected detector readings.
    
    Parameters
    ----------
    bstop_corrected : float
        Dark current corrected beamstop reading
    i0_corrected : float
        Dark current corrected I0 reading
    bstop_air : float, optional
        Beamstop reading for air/empty measurement
    i0_air : float, optional
        I0 reading for air/empty measurement
        
    Returns
    -------
    float
        Transmission factor
        
    Notes
    -----
    If air readings are provided, calculates normalized transmission.
    Otherwise returns simple beamstop reading.
    """
    if bstop_air is not None and i0_air is not None:
        # Normalized transmission calculation
        trans_factor = (bstop_corrected / i0_corrected) / (bstop_air / i0_air)
    else:
        # Simple transmission factor
        trans_factor = bstop_corrected
    
    return trans_factor


def calculate_normalization_factor(trans_factor: float, i0_corrected: float) -> float:
    """
    Notes
    -----
    The normalization factor accounts for transmission and incident beam intensity.
    """
    normfactor = trans_factor * i0_corrected
    
    # Optional thickness correction could be added here

    return float(normfactor)


def process_transmission_correction(i0_raw: float, bstop_raw: float,
                                  compound: str, config: Dict) -> Dict[str, float]:
    """
    Complete transmission and thickness correction workflow.
    
    Parameters
    ----------
    i0_raw : float
        Raw I0 detector reading
    bstop_raw : float
        Raw beamstop detector reading
    compound : str
        Chemical formula of the material (e.g., "C2H4")
    i0_offset : float, default I0_OFFSET
        I0 dark current offset
    bstop_offset : float, default BSTOP_OFFSET
        Beamstop dark current offset
    i0_air : Optional[float], default I0_AIR
        I0 reading for air/empty measurement
    bstop_air : Optional[float], default BSTOP_AIR
        Beamstop reading for air/empty measurement
    energy_keV : float, default ENERGY_KEV
        X-ray energy in keV
    density : float, default DENSITY_PE
        Material density in g/cm³
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing all correction parameters:
        - 'i0_corrected': Dark current corrected I0
        - 'bstop_corrected': Dark current corrected beamstop
        - 'trans_factor': Transmission factor
        - 'thickness_m': Sample thickness in meters
        - 'normfactor': Normalization factor
        - 'mu': Absorption coefficient
        - 'sld': Scattering length density
    """
    # Get parameters from config
    i0_offset = config['experimental']['i0_offset']
    bstop_offset = config['experimental']['bstop_offset']
    i0_air = config['experimental']['i0_air']
    bstop_air = config['experimental']['bstop_air']
    energy_keV = config['experimental']['energy_keV']
    density = config['experimental']['density_g_cm3']
    
    # Step 1: Apply dark current corrections
    i0_corrected, bstop_corrected = correct_detector_offsets(i0_raw, bstop_raw, 
                                                           i0_offset, bstop_offset)
    
    # Step 2: Apply air corrections if provided
    if i0_air is not None and bstop_air is not None:
        i0_air_corrected = i0_air - i0_offset
        bstop_air_corrected = bstop_air - bstop_offset
    else:
        i0_air_corrected = None
        bstop_air_corrected = None
    
    # Step 3: Calculate transmission factor
    trans_factor = calculate_transmission_factor(bstop_corrected, i0_corrected,
                                               bstop_air_corrected, i0_air_corrected)
    
    # Step 4: Calculate transmission for thickness calculation
    trans_factor_thickness = bstop_corrected / i0_corrected if i0_corrected != 0 else bstop_corrected
    
    # Step 5: Calculate material properties and thickness
    material_props = calculate_sld_mu_thickness(energy_keV, density, trans_factor_thickness, compound)
    
    # Step 6: Calculate normalization factor
    normfactor = float(trans_factor * i0_corrected)
    
    return {
        'i0_corrected': i0_corrected,
        'bstop_corrected': bstop_corrected,
        'trans_factor': trans_factor,
        'thickness_m': material_props['thickness_m'],
        'normfactor': normfactor,
        'mu': material_props['mu'],
        'sld': material_props['sld']
    }


def process_swaxs_data(base_directory: str) -> Dict[str, Dict]:
    """
    Complete SWAXS data processing workflow using configuration file.
    
    Parameters
    ----------
    base_directory : str
        Base directory path for the project
        
    Returns
    -------
    Dict[str, Dict]
        Dictionary containing processed results for SAXS and WAXS:
        - 'saxs': SAXS processing results
        - 'waxs': WAXS processing results
        - 'parameters': Experimental parameters
        - 'config': Configuration used
    """
    # Load configuration
    with open("config.yml", 'r') as f:
        config = yaml.safe_load(f)
    
    # Construct absolute paths using base directory
    data_dir = os.path.join(base_directory, config['data_directory'].lstrip('/'))
    poni_dir = os.path.join(base_directory, config['poni_directory'].lstrip('/'))
    
    # Find CSV parameter files
    csv_pattern = os.path.join(data_dir, "*_scan1.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        raise FileNotFoundError(f"No CSV parameter files found in {data_dir}")
    
    # Process first CSV file for parameters
    csv_file = csv_files[0]
    parameters = read_csv_parameters(csv_file)
    
    # Process transmission correction
    corrections = process_transmission_correction(
        i0_raw=parameters['i0_avg'],
        bstop_raw=parameters['bstop_avg'],
        compound=config['compound'],
        config=config
    )
    
    # Load PyFAI integrators from config
    saxs_poni = os.path.join(poni_dir, config['poni_files']['saxs'])
    waxs_poni = os.path.join(poni_dir, config['poni_files']['waxs'])
    
    if not os.path.exists(saxs_poni):
        raise FileNotFoundError(f"SAXS PONI file not found: {saxs_poni}")
    if not os.path.exists(waxs_poni):
        raise FileNotFoundError(f"WAXS PONI file not found: {waxs_poni}")
        
    ai_saxs = pyFAI.load(saxs_poni)
    ai_waxs = pyFAI.load(waxs_poni)
    
    # Load masks from config
    saxs_mask_path = os.path.join(poni_dir, config['mask_files']['saxs'])
    waxs_mask_path = os.path.join(poni_dir, config['mask_files']['waxs'])
    
    saxs_mask = fabio.open(saxs_mask_path).data if os.path.exists(saxs_mask_path) else None
    waxs_mask = fabio.open(waxs_mask_path).data if os.path.exists(waxs_mask_path) else None
    
    # Find and process detector files
    saxs_files = find_detector_files(data_dir, "SAXS")
    waxs_files = find_detector_files(data_dir, "WAXS")
    
    results = {
        'config': config,
        'parameters': parameters,
        'corrections': corrections,
        'saxs': {},
        'waxs': {}
    }
    
    # Process SAXS data
    if not saxs_files:
        raise RuntimeError("No SAXS files found in directory")

    saxs_data = read_raw_detector_file(saxs_files[0], "SAXS")
    
    # Perform 1D integration using config parameters
    q_saxs, I_saxs, error_saxs = ai_saxs.integrate1d(
        saxs_data, 
        config['integration']['npt_radial'], 
        error_model=config['integration']['error_model'],
        mask=saxs_mask,
        normalization_factor=corrections['normfactor']
    )
    
    results['saxs'] = {
        'q': q_saxs,
        'intensity': I_saxs,
        'error': error_saxs,
        'raw_file': saxs_files[0]
    }
    
    # Process WAXS data
    if not waxs_files:
        raise RuntimeError("No SAXS files found in directory")

    waxs_data = read_raw_detector_file(waxs_files[0], "WAXS")
    
    # Perform 1D integration using config parameters
    q_waxs, I_waxs, error_waxs = ai_waxs.integrate1d(
        waxs_data,
        config['integration']['npt_radial'],
        error_model=config['integration']['error_model'], 
        mask=waxs_mask,
        normalization_factor=corrections['normfactor']
    )
    
    results['waxs'] = {
        'q': q_waxs,
        'intensity': I_waxs,
        'error': error_waxs,
        'raw_file': waxs_files[0]
    }

    return results


def main():
    """
    Main execution function for single run SWAXS data processing.
    
    This function demonstrates how to use the refactored pipeline to process
    a single experimental run using the current working directory as base.
    """
    import sys
    
    # Get base directory (current working directory)
    base_directory = os.getcwd()
    
    # Allow config path to be specified as command line argument
    config_path = "config.yml"
    
    print(f"Processing SWAXS data from: {base_directory}")
    print(f"Using config file: {config_path}")
    
    try:
        # Process the data
        results = process_swaxs_data(base_directory)
        
        print("\n=== Processing Results ===")
        print(f"Configuration: {results['config']}")
        print(f"Parameters: {results['parameters']}")
        print(f"Corrections: {results['corrections']}")
        
        if results['saxs']:
            print(f"SAXS data processed: {len(results['saxs']['q'])} q-points")
            print(f"SAXS file: {results['saxs']['raw_file']}")
            
        if results['waxs']:
            print(f"WAXS data processed: {len(results['waxs']['q'])} q-points") 
            print(f"WAXS file: {results['waxs']['raw_file']}")
            
        print("\nProcessing completed successfully!")
        return results
        
    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()