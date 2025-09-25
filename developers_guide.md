# Developers Guide

## Notes
1. Path modifications are taken from the pathlib module. This [quickstart](https://realpython.com/python-pathlib/) is a very useful start guide.  

Paths are their own type with the pathlib module, but some functions, for instance pyfai's integrade 1d, require that a string is passed in so this path sometimes needs to be converted back to a string.
2. If the periodic table library doesn't support a compound, then you may need to change the code to allow the user to specify needed values in config.yml
3. If you add another attribute to correction_integration_demo.py, make sure to also add it to the type hints at the top.
4. If metadata has a string type that's not convertible to an integer, files_averaging may have issues averaging metadata

## How to add other features

## 1. Adding other beamlines (like 17-2)
Add the custom raw file and metadata file reading logic to the respective files in src/. Make sure to add new functions for these beamlines. Then, chain the main_correction script to call these functions when the configurations correspond. It is strongly recommended to make all raw file reading functions and all metadata reading functions have the same inputs and outputs.

## 2. Live data processing 
- Use asynchronous programming with watchdog to monitor for new files being created and call a wrapper function if so. If it is created, sleep the program until the file has the right contents uploaded and the metadata also exists. You can make it sleep for 5 seconds and do a check in to see if the data is done uploading for the file. After these steps are completed, run the main corrections and reductions to get the 1d data.

## 3. Extra Corrections
To add additional corrections, you can replace get_corrections_full with a different corrections function

## Possible Issues
If you get this kind of an error, this was caused by poni_directory not being a Path type.
File "/Users/vsoma_ogh6dsk/Documents/SLAC/SWAXS_data_correction_reduction_averaging/src/main_correction_reduction_v1.py", line 154, in _load_saxs_integrator
    self.saxs_poni_path = self.poni_directory / self.poni_files['saxs']
TypeError: unsupported operand type(s) for /: 'str' and 'str'

## API Documentation - Updated Method Signatures

### Experiment Class Methods

#### `calculate_sld_mu_thickness(self, transmission: float) -> Dict[str, float]`
Calculate material properties from transmission using experiment parameters.

**Parameters:**
- `transmission` (float): Measured transmission value

**Returns:**
- `Dict[str, float]`: Dictionary containing mu, sld, and thickness_m

#### `get_corrections_full(self, raw_file_path: Path, metadata_dict) -> Dict[str, float]`
Process metadata and compute correction factors for data normalization.

**Parameters:**
- `raw_file_path` (Path): Path to the raw file being processed
- `metadata_dict` (dict): Dictionary containing metadata with i0 and bstop values

**Returns:**
- `Dict[str, float]`: Dictionary containing correction factors including i0_corrected, bstop_corrected, transmission_factor, transmission_ratio, thickness, and normalization_factor

#### `process_saxs_file(self, raw_file_path: Path) -> Dict`
Process a single SAXS .raw file with corrections and 1D integration.

**Parameters:**
- `raw_file_path` (Path): Path to the SAXS .raw file

**Returns:**
- `Dict`: Dictionary containing q, intensity, error, filename, corrections, and raw_file_path

#### `process_waxs_file(self, raw_file_path: Path) -> Dict`
Process a single WAXS .raw file with corrections and 1D integration.

**Parameters:**
- `raw_file_path` (Path): Path to the WAXS .raw file

**Returns:**
- `Dict`: Dictionary containing q, intensity, error, filename, corrections, and raw_file_path

#### `_create_output_directory(self, detector_type: str) -> Path`
Create output directory structure with SAXS/Reduction or WAXS/Reduction subdirectories.

**Parameters:**
- `detector_type` (str): Type of detector ('SAXS' or 'WAXS')

**Returns:**
- `Path`: Output directory path