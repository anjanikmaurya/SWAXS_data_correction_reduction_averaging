# SWAXS_data_reduction_correction_Analysis

This is a user-guide to use this repository, with detailed steps about how to use it. For a developer-guide, with details about adding code, use developer_guide.md. If you are using VSCode, using Cmd+Shift+V (or Ctrl+Shift+V on Windows) to preview the file in proper markdown.

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

Basic demo example configured in demo/, with a built in config.yml file in this directory working.

### 2a. Look through demo structure
Just see the proper folder structure in demo/2D. Also, look through demo/config.yml file to see a configuration file. I recommend

### 2b. Run Demo For Correction and Reduction
To use the following script, run the following demo code.
```bash
python src/main_correction_reduction_v1.py demo/config.yml
```
src/main_correction_reduction_v1.py is the primary analysis script, which corrects and reduces data.
demo/config.yml is a YML configuration file for the demo structure

### 2c. Look through corrections
Now, you should see a 1D/Reductions get created in the demo directory along with a logging file of the logs of the changes. Look through these files to understand the corrections format

### 2d. Run Averaging Demo
Finally, run the demo with averaging. Use the src/files_averaging.ipynb notebook to perform the averaging of thi sexample

## Running a New Data Correction

After performing the data corrections
## Reformatting Directory Structure
If you would like to migrate an old directory structure into a new directory structure, then type:
```bash
python bin/copy_data_structure.py
```

## Notes About File Information
If you have an existing directory, then running this file will add

If you have an existing analyzed data using the same script in the same 1D directory, running this code will overwrite your existing files. 