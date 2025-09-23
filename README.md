# SWAXS_data_reduction_correction_Analysis

This is a user-guide to use this repository, with detailed steps about how to use it. For a developer-guide, with details about adding code, use developer_guide.md. If you are using VSCode, using Cmd+Shift+V (or Ctrl+Shift+V on Windows) to preview the file in proper markdown.

## 1. Set up code
First, clone the repository.
Option 1: Cloning with Visual Studio Code (recommended if using vscode):
On the home screen 
```bash
git clone https://github.com/anjanikmaurya/SWAXS_data_correction_reduction_averaging
```
Then, it is highly recommended to create a virtual environment for this project, but you may skip this step and the next . Right now, this is using pip and venv, but you can use a different package manager if you would like. 

Create this virtual 
```bash
python -m venv .venv  
```
Activate the virtual environment
```bash
source .venv/bin/activate
```
If it went correctly, you should now see the (.venv) on your terminal. Next, within your virtual environment itself (if you created one), perform this command to install all dependencies in this repository. This command takes a while to run
```bash
pip install -r requirements.txt
```

## 2. Run demo example

Basic demo example configured in demo/, with a built in config.yml file in this directory working. First, look through the demo directory to understand the proper file structure. Also, look through the config.yml file to see a simple configuration file.

To use the following script, run the following code:
```bash
python src/main_correction_reduction_v1.py demo/config.yml
```
src/main_correction_reduction_v1.py is the primary analysis script, which corrects and reduces data.
demo/config.yml is a YML configuration file for the demo structure

After you run this, you should see a logging file get created and 

## 3. Set up configuration YML file specific to this experiment
Create a configuration file (similar template to config.yml) specific to a single experiment data round. 

## Reformatting Directory Structure
If you would like to migrate an old directory structure into a new directory structure, then type:
```bash
python src/copy_data_structure.py -h
```

## Notes About File Information
If you have an existing directory, then running this file will add

If you have an existing analyzed data using the same script in the same 1D directory, running this code will overwrite your existing files. 