# MS_beta_fatigue
Code for Tatti Benelli et al., 2025 Brain communications

This repository contains MATLAB scripts used for the analysis and figure generation in:

Tatti E., et al. (2025). *Linking Movement-Related Beta Oscillations to Cortical Excitability, Structural Damage and Fatigue in Multiple Sclerosis* Brain Communications.

## Contents

- `code/main_analysis.m`: main pipeline to reproduce the computation of beta-modulation depth (ERD/ERS) and the correlation analyses with fatigue scales across ROIs.
- `code/compute_modulation_depth.m`: function to compute movement-related beta desynchronization (ERD) and synchronization (ERS) and modulation depth.
- `code/plot_*`: scripts for generating the main time–frequency and correlation figures.
- `example_data/` (optional): placeholder or toy data illustrating the expected data structure.

## Requirements

- MATLAB R2023 
- Signal Processing Toolbox and/or other relevant toolboxes 
- FieldTrip toolbox, see https://www.fieldtriptoolbox.org.
- EEGLAB Toolbox

## Usage

These scripts are provided to document the analysis pipeline. Due to privacy and consent constraints, the original EEG and clinical datasets cannot be shared publicly. To run the code end-to-end, users should adapt the input paths and replace our data folders with their own data.

For questions about the code, please contact etatti@med.cuny.edu.
