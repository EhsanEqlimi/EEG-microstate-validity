# EEG Microstate Analysis Pipeline

This MATLAB pipeline performs EEG microstate analysis for single subjects using **FieldTrip preprocessing** and **subject-specific k-means clustering**. Selected functions from the **Poulsen et al., 2018 Microstate EEGLAB Toolbox** are used programmatically (no GUI) for microstate backfitting, label smoothing, and statistics extraction.

## Features

- **Preprocessing** using FieldTrip:
  - Re-referencing, demeaning, detrending.
  - Band-pass filtering.
- **Global Field Power (GFP)** computation and peak detection.
- **Subject-specific microstate extraction** using modified k-means clustering.
- **Backfitting** of EEG data to individualized microstates using Poulsen et al., 2018 toolbox functions.
- **Smoothing of microstate labels** to remove noisy segments.
- **Microstate statistics** extracted from labels (duration, occurrence, transitions, GMD).
- Supports **testing a range of microstate numbers** (e.g., `k = 4:10`) for optimal clustering.

## Toolbox Requirements


This repository uses the **EEGLAB Microstate Toolbox** (Poulsen et al., 2018) for microstate clustering, backfitting, smoothing.

**Dependencies:**  
- [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) (for EEG data structures)  
- [FieldTrip](https://www.fieldtriptoolbox.org/) (for preprocessing and GFP extraction)  

**Notes:**  
- No GUI is required; all functions are called programmatically.  
- Microstate analysis is automated for a range of k (e.g., `4:10`).  

**Reference:**  
Poulsen, A. T., et al. (2018). *EEG microstate analysis with the EEGLAB Microstate Toolbox*. [Link](https://archive.compute.dtu.dk/files/public/users/atpo/Microstate)


## Functions

### `FnParamMaker4Microstate`

This function **creates a comprehensive parameter structure** for EEG microstate analysis. It defines all settings required for the pipeline, including preprocessing, GFP computation, k-means clustering, backfitting, smoothing, and visualization.  

It must be called **before** the main analysis function (`FnMicrostateOneSubject`) to ensure that all necessary parameters are available.

**Usage:**

```matlab
Param = FnParamMaker4Microstate();
MSResults=FnMicrostateOneSubject(EEGLabData,Param);
