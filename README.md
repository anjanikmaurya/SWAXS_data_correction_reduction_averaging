# SWAXS_data_reduction_correction_Analysis

This is a user guide for using this repository, with detailed steps on how to use it. For a developer guide with details about adding code, see `developers_guide.md`. If you are using VSCode, use Cmd+Shift+V (or Ctrl+Shift+V on Windows) to preview the file in proper markdown.

## Purpose
This directory can run the corrections and reductions together. After that, it also provides file averaging functionality. This will speed up these steps considerably and allow users to quickly move on to analysis for supported data.

## Instructions for.quickstart
Make sure you have python 3.9+ and git installed. If you would like to use an older version of python, minor code modifications are needed. Then, follow these steps
## 1. Set up code

### 1a. Clone Repository

First, you will need to clone the repository to get the code on your machine 

#### Option 1: Cloning with Visual Studio Code (recommended if using VSCode):
On the home screen of VSCode, click "Clone Git Repository" and type the repository url: https://github.com/anjanikmaurya/SWAXS_data_correction_reduction_averaging

#### Option 2: Clone from terminal
On the directory where you are active, you can type this command to clone the repository
```bash
git clone https://github.com/anjanikmaurya/SWAXS_data_correction_reduction_averaging
```

### 1b. Create a virtual environment (optional but highly recommended)
Then, it is highly recommended to create a virtual environment for this project, but you may skip this step and the next. Right now, this is using pip and venv, but **you can use a different package manager if you would like**. 
```bash
python -m venv .venv  
```

Also, activate the virtual environment
```bash
source .venv/bin/activate
```

If it went correctly, you should now see "(.venv)" in your terminal prompt. You will need to reactivate the virtual environment with the same command on each run 

### 1c. Install all packages
Next, within your virtual environment (if you created one), run this command to install all dependencies in this repository. This command takes a while to run:
```bash
pip install -r requirements.txt
```

## 2. Run Demo Example

Simple demo example configured in the demo/ directory, with a built-in config.yml file in this directory working. You can use these steps as guidance for what to do with your real experiment

### 2a. Look through demo structure
- Look at the proper folder structure in demo/2D. You will need to format your data in this way.
- Look through the demo/config.yml file to see an example configuration file. 
### 2b. Run Demo For Correction and Reduction
The following code executes the script and performs corrections
```bash
python src/main_correction_reduction_v1.py demo/config.yml
```
src/main_correction_reduction_v1.py is the primary analysis script, which corrects and reduces data.
demo/config.yml is a YAML configuration file for the demo structure

### 2c. Look through corrected files
You should see two new items created if the corrections were performed successfully:
- A 1D directory created in the demo directory includes a SAXS/Reduction and WAXS/Reduction directory with 5 processed files.
- A logging file with logs of the changes.
This is similar to how the processed data will look on a new run.

### 2d. Run Averaging Demo
Finally, run the demo with averaging. Use the src/files_averaging.ipynb notebook to perform averaging of this example. If this worked correctly, you should see a SAXS/Averaged and WAXS/Averaged directory with one file.

## 3. Running a New Data Correction
Follow a similar process from the demo, but make a copy of the config.yml file and set it up for your own data correction. You will not need to edit the code unless you would like to enable a new feature. After configuring this, run
```bash
python src/main_correction_reduction_v1.py your_config.yml
```

## Recommended File Structure Setup
File names may be different, but this is the recommended directory structure setup (as shown in the demo):

```
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
```

## Reformatting Directory Structure to the Correct One
If you have data with the legacy format (likely collected before this script was made) then you should reformat it to the correct format.

### Old directory structure
This structure is **not supported**, but can be transformed into the correct data structure using copy_data_structure.py. This is a simple example for demonstration:

```
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
```

Note that this is segregated by run numbers, unlike the correct data structure. 

### Correction Steps
Use the following code to convert it to the right structure
```bash
python bin/copy_data_structure.py --source source_folder --target target_folder
```
Here, source_folder is the folder storing data in the old format and target_folder is the new folder (initially empty) which will store the data in the new format.

If you have other data organization structures, this script will not be able to convert it automatically. One recommended solution is to create a program to change the structure from your specific format to this, perhaps with LLM assistance.