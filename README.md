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

- **FieldTrip** – preprocessing and GFP computation.  
  [FieldTrip](https://www.fieldtriptoolbox.org/)  

- **EEGLAB** – loading `.mat` EEG datasets.  
  [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php)  

- **Microstate EEGLAB Toolbox (Poulsen et al., 2018)** – used programmatically for backfitting, smoothing, and statistics (no GUI).  
  [Microstate Toolbox](https://github.com/microstates/EEG-microstates)  

- **Signal Processing Toolbox** – filtering, peak detection, and signal operations.

## Functions

### `FnParamMaker4Microstate`

Creates a parameter struct for preprocessing, GFP computation, clustering, backfitting, and smoothing. Must be called **before** the main analysis function.

```matlab
Param = FnParamMaker4Microstate();
