# EEG Microstate Analysis Pipeline

This repository provides a MATLAB framework for EEG microstate analysis, supporting preprocessing, microstate extraction, backfitting, smoothing, and statistical evaluation for single-subject datasets. It integrates code compatible with **EEGLAB**, **FieldTrip**, and also includes routines inspired by **DTU EEG microstate toolbox**.

## Features

- Preprocessing: re-referencing, detrending, demeaning, and band-pass filtering.
- Global Field Power (GFP) computation and peak detection.
- Subject-specific microstate extraction using modified k-means clustering.
- Backfitting of EEG data to individualized microstates.
- Smoothing of microstate labels for noise reduction.
- Computation of microstate statistics, including duration, occurrence, transitions, and Global Map Dissimilarity (GMD).
- Storage of all relevant results: GFP, prototypes, labels, smoothed labels, clustering details, and statistics.
- Compatible with DTU EEG microstate toolbox routines for clustering and backfitting.

## Toolbox Requirements

This pipeline depends on the following MATLAB toolboxes and libraries:

- **EEGLAB**: For loading and handling EEG datasets in `.mat` format.  
  [EEGLAB Website](https://sccn.ucsd.edu/eeglab/index.php)

- **FieldTrip**: For EEG preprocessing, Global Field Power calculation, and trial handling.  
  [FieldTrip Website](https://www.fieldtriptoolbox.org/)

- **Signal Processing Toolbox**: For filtering, peak detection, and other signal operations.

- **EEG Microstate Toolbox**:  
  [EEG Microstate Toolbox](https://archive.compute.dtu.dk/files/public/users/atpo/Microstate)

Ensure that these toolboxes are installed and added to your MATLAB path before running the analysis.

## Functions

### `FnParamMaker4Microstate`

Creates a parameter structure for the EEG microstate pipeline. **This function should be called before `FnMicrostateOneSubject`.**

```matlab
Param = FnParamMaker4Microstate();

