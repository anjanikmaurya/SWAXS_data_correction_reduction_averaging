# Developers Guide

## Notes
1. Path modifications use the pathlib module with the Path data type to provide more flexibility for path operations. This may be a source of issues early on. I strongly recommend reading this [quickstart](https://realpython.com/python-pathlib/) to get an overview of the module.  

Paths are their own type with the pathlib module, but some functions, for instance PyFAI's integrate1d, require that a string is passed in so this path sometimes needs to be converted back to a string.
2. If the periodic table library doesn't support a compound, then you may need to change the code to allow the user to specify the needed values in config.yml.
3. If you add another attribute to the main correction script, make sure to also add it to the type hints at the top.
4. If metadata has a string type that's not convertible to an integer, files_averaging may have issues averaging metadata.

## How to add other features

### 1. Adding Other Beamlines (like 17-2)
Add the custom raw file and metadata file reading logic to the respective files in src/. Make sure to add new functions for these beamlines. Then, chain the main correction script to call these functions when the configurations correspond. It is strongly recommended to make all raw file reading functions and all metadata reading functions have the same inputs and outputs as one another to enable future development. For instance, if the CSV metadata reading function takes in a raw file path and outputs a dictionary of metadata values, it is best if the PDI metadata reading function does the same (and it currently does).  

Note that you may need to perform additional steps, like extra corrections and PyFAI operations. 

### 2. Live Data Processing 
Use asynchronous programming with watchdog to monitor for new files being created and call a wrapper function when detected. If a file is created, pause the program until the file has the correct contents uploaded and the metadata also exists. You can make it sleep for 5 seconds and perform a check to see if the data is done uploading for the file. After these steps are completed, run the main corrections and reductions to get the 1D data.

### 3. Extra Corrections
To add additional corrections, you can replace `get_corrections_full` with a different correction function with the required steps and let the user specify which correction function they need in the YAML file. 

### 4. Plotting
`src/utils/read_dat_metadata` can help with plotting and data checking, since it returns both the data and the metadata. You can use this on a reduced file for individual dat and metadata or an averaged file for averaged dat and metadata. 

### 5. Adding more packages
All packages should be specified in the requirements file. If you have installed more packages and are using pip, then you can move all your packages to the requirements file with this command:

```bash
pip freeze > requirements.txt
```

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

### Module-Level Functions (main_correction_reduction_v1.py)

#### `find_all_raw_files(experiment: Experiment, data_directory_path: Path) -> Tuple[List[Path], List[Path]]`
Find all .raw files in SAXS and WAXS directories.

**Parameters:**
- `experiment` (Experiment): Experiment instance containing configuration
- `data_directory_path` (Path): Full path to the data directory (e.g., "run12_test/2D/")

**Returns:**
- `Tuple[List[Path], List[Path]]`: Lists of SAXS and WAXS .raw file paths

#### `full_correction_integration(config_file_path='config.yml', plotting=False)`
Main function to process all .raw files in the dataset.

**Parameters:**
- `config_file_path` (str): Path to configuration file (default: 'config.yml')
- `plotting` (bool): If True, return plotting data instead of just processing files (default: False)

**Returns:**
- `Dict or None`: If plotting=True, returns dictionary with plotting data. If plotting=False, returns None (just processes files).

#### `add_metadata_to_dat(dat_file_path: Path, metadata_dict: dict)`
Add metadata to the end of a .dat file as commented information.

**Parameters:**
- `dat_file_path` (Path): Path to the .dat file to modify
- `metadata_dict` (dict): Dictionary containing metadata with keys as column names and values as data

**Notes:**
- Appends metadata as comments (lines starting with #) to end of .dat file
- Writes dictionary keys and values in a structured format
- Preserves original .dat file structure and data columns

### Utility Functions (src/utils/read_dat_metadata.py)

#### `read_dat_data_metadata(file_path: Path) -> tuple`
Reads a .dat file and extracts header comments, data, and metadata.

**Parameters:**
- `file_path` (Path): Path to the .dat file to read

**Returns:**
- `tuple`: (header_lines, q_values, intensity_values, sigma_values, metadata_dict)
  - `header_lines` (list): List of header comment lines
  - `q_values` (np.ndarray): Q-vector values
  - `intensity_values` (np.ndarray): Intensity values
  - `sigma_values` (np.ndarray): Error/sigma values
  - `metadata_dict` (dict): Dictionary of metadata key-value pairs converted to floats