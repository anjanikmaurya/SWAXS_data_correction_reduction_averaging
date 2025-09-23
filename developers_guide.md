# Developers Guide

## Notes
1. If the periodic table library doesn't support a compound, then you may need to change the code
2. If you add another attribute to correction_integration_demo.py, make sure to also add it to the type hints at the top.
3. If metadata has a string type that's not convertible to an integer, files_averaging may have issues averaging metadata
4. Path conventions are taken from the pathlib module. This [quickstart](https://realpython.com/python-pathlib/) is a very useful start guide.

## How to add other features

## Adding other beamlines (like 17-2)
- Add raw file reading function and metadata file reading functions in src/
- Swap out the correction functions to a different one and return a list of kwargs
- Make sure to have the same inputs

## Live data processing steps




## Model of How Code works

## Possible Issues
If you get this kind of an error, this was caused by poni_directory not being a Path type.
File "/Users/vsoma_ogh6dsk/Documents/SLAC/SWAXS_data_correction_reduction_averaging/src/main_correction_reduction_v1.py", line 154, in _load_saxs_integrator
    self.saxs_poni_path = self.poni_directory / self.poni_files['saxs']
TypeError: unsupported operand type(s) for /: 'str' and 'str'