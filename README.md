# EEG Microstate Analysis & Validity Toolbox

This MATLAB pipeline performs **EEG microstate analysis** for single subjects and extends the classical approach to **test the validity of common assumptions** in the field.  

Traditional microstate analysis assumes:
1. A **one-hot / winner-take-all representation** — only one microstate is active at a time.  
2. A strong **sparsity assumption** — EEG scalp patterns can be compressed into a small set of prototypical maps.  

This repository was created to provide not only the **classical microstate pipeline** (based on GFP peaks, modified k-means clustering, backfitting, smoothing, and statistics) but also to investigate these assumptions using:  

- **Dimensionality-based measures** (e.g., Rényi information dimension, compressibility, ARMA modeling).  
- **Spectral power analysis** (linking microstates to oscillatory dynamics).  
- **Comparison of GFP-peak vs. non-peak samples** to test whether restricting to GFP peaks is valid.  

In short: this toolbox aims to **bridge classical microstate analysis with modern information-theoretic and spectral approaches** to assess the validity of its core assumptions.

---

## Features

- **Preprocessing** using FieldTrip:
  - Re-referencing, demeaning, detrending  
  - Band-pass filtering  

- **Global Field Power (GFP)** computation and peak detection  

- **Subject-specific microstate extraction** using modified k-means clustering  

- **Backfitting** of EEG data to individualized microstates (Poulsen et al., 2018 toolbox functions)  

- **Smoothing of microstate labels** to remove noisy segments  

- **Microstate statistics** extracted from labels (duration, occurrence, transitions, GMD)  

- Supports **testing a range of microstate numbers** (e.g., `k = 4:10`)  

- **Validity extensions**:
  - **Dimensionality metrics** (Rényi information dimension, ARMA-based compressibility)  
  - **Spectral analysis** (oscillatory power linked to microstates)  
  - **Evaluation of one-hot assumption** by comparing GFP-peak vs. non-peak samples  

---

## Toolbox Requirements

This repository uses the **EEGLAB Microstate Toolbox** (Poulsen et al., 2018) for microstate clustering, backfitting, smoothing.

**Dependencies:**  
- [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) (for EEG data structures)  
- [FieldTrip](https://www.fieldtriptoolbox.org/) (for preprocessing and GFP extraction)  

**Notes:**  
- No GUI is required; all functions are called programmatically.  
- Microstate analysis is automated for a range of `k` (e.g., 4–10).  

**Reference:**  
Poulsen, A. T., et al. (2018). *EEG microstate analysis with the EEGLAB Microstate Toolbox*. [Link](https://archive.compute.dtu.dk/files/public/users/atpo/Microstate)

---

## Functions

### `FnParamMaker4Microstate`

This function **creates a comprehensive parameter structure** for EEG microstate analysis. It defines all settings required for the pipeline, including preprocessing, GFP computation, k-means clustering, backfitting, smoothing, spectral analysis, and dimensionality measures.  

It must be called **before** the main analysis function (`FnMicrostateOneSubject`) to ensure that all necessary parameters are available.

**Usage:**

```matlab
Param = FnParamMaker4Microstate();
MSResults = FnMicrostateOneSubject(EEGLabData, Param);
