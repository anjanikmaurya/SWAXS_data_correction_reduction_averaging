# SWAXS_data_reduction_correction_Analysis

This is a user-guide to use this repository, with detailed steps about how to use it. For a developer-guide, with details about adding code, use developer_guide.md. If you are using VSCode, using Cmd+Shift+V (or Ctrl+Shift+V on Windows) to preview the file in proper markdown.

## Purpose
This directory can run the corrections and reductions together. After that, it also provides file averaging functionality. This is likely to speed up these steps considerably and users will quickly be able to move on to the analysis for supported data.

## Instructions
## 1. Set up code

### 1a. Clone Repository
First, clone the repository.

Option 1: Cloning with Visual Studio Code (recommended if using vscode):
On the home screen, go to clone git repository
```bash
git clone https://github.com/anjanikmaurya/SWAXS_data_correction_reduction_averaging
```

### 1b. Create a virtual environment (optional but highly recommended)
Then, it is highly recommended to create a virtual environment for this project, but you may skip this step and the next . Right now, this is using pip and venv, but you can use a different package manager if you would like. 
```bash
python -m venv .venv  
```
Activate the virtual environment
```bash
source .venv/bin/activate
```

If it went correctly, you should now see the (.venv) on your terminal. 

### 1c. Install all packages
Next, within your virtual environment itself (if you created one), perform this command to install all dependencies in this repository. This command takes a while to run
```bash
pip install -r requirements.txt
```

## 2. Run Demo Example

Simple demo example configured in demo/ directory, with a built in config.yml file in this directory working. You can use these steps as guidance for what to do with your real experiment

### 2a. Look through demo structure
- Look at the proper folder structure in demo/2D
- Look through demo/config.yml file to see a configuration file. 
### 2b. Run Demo For Correction and Reduction
The follow code executes the script and performs corrections
```bash
python src/main_correction_reduction_v1.py demo/config.yml
```
src/main_correction_reduction_v1.py is the primary analysis script, which corrects and reduces data.  
demo/config.yml is a YML configuration file for the demo structure

### 2c. Look through corrections new files
You should see two new items created:
- Now, you should see a 1D directory get created in the demo directory.
- You should also see a logging file of the logs of the changes.


### 2d. Run Averaging Demo
Finally, run the demo with averaging. Use the src/files_averaging.ipynb notebook to perform the averaging of this example. 

## 3. Running a New Data Correction
Follow similar process from the demo, but make a copy of the config.yml file and set it up on your own to run a new data correction. You will not need to edit the code unless you would like to enable a new feature. After performing this, run
```bash
python src/main_correction_reduction_v1.py demo/config.yml
```

## Recommended File Structure Setup
File names may be different, but this is the recommended directory structure setup (and is in demo).
.
├── 2D
│   ├── Hor_scan_Run4_RampT20_ctr0_scan1.csv
│   ├── SAXS
│   │   ├── sone_Hor_scan_Run4_RampT20_ctr0_scan1_0000.raw
│   │   ├── sone_Hor_scan_Run4_RampT20_ctr0_scan1_0000.raw.pdi
│   └── WAXS
│       ├── b_tassone_Hor_scan_Run4_RampT20_ctr0_scan1_0000.raw
│       ├── b_tassone_Hor_scan_Run4_RampT20_ctr0_scan1_0000.raw.pdi
└── poni
    ├── atT_SAXS.poni
    ├── atT_WAXS.poni
    ├── RT_SAXS_mask_03.edf
    └── RT_WAXS_mask.edf

## Reformatting Directory Structure to correct one
Old directory structure. This structure is **not supported**, but can be transformed into the correct data structure using copy_data_structure.py. This is a simple example for demonstration.
.
├── poni
├── Run5
│   ├── Hor_scan_Run5_RampT20_ctr0_scan1.csv
│   └── SAXS
│       └── sone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.raw
├── Run6
│   ├── Hor_scan_Run6_RampT20_ctr0_scan1.csv
│   └── SAXS
│       └── sone_Hor_scan_Run6_RampT20_ctr0_scan1_0000.raw
└── Run7
    ├── Hor_scan_Run7_RampT20_ctr0_scan1.csv
    └── SAXS
        └── sone_Hor_scan_Run7_RampT20_ctr0_scan1_0000.raw

Note that this is segregated by run numbers, unlike the correct data structure. Use the following code to get it into the right structure.
```bash
python bin/copy_data_structure.py --source source_folder --target target_folder
```
Here, source_folder is the folder storing data in the old format and target folder is the new folder (now empty) which will store the data in the new format.

If you have other data organization structures, this script will not be able to convert it. You can create a program to change the structure, perhaps with LLM assistance.