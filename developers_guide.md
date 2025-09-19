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

REVISIT: How important is it to specify a particular python version?
