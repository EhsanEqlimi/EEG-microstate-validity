
# EEG-microstate-validity

A MATLAB framework for evaluating the **reliability** and **validity** of EEG microstate segmentation and clustering. The toolbox provides tools for reproducibility testing, cluster stability evaluation, spectral and information-theoretic metrics, and guidelines for parameter selection.

---

## Motivation
EEG microstates are short-lived patterns of brain activity that reflect fundamental neural processes. While widely studied, the **validity and reliability** of their extraction methods often remain underexplored. This framework provides a standardized set of tools to evaluate:

- Stability of microstate segmentation
- Validity across subjects and sessions
- Effects of preprocessing and clustering parameters
- Robustness of temporal and spectral metrics

---

## Features
- **Segmentation & Clustering**: Multiple algorithms with flexible parameter settings.
- **Metrics**: Variance explained, reproducibility, cluster stability, information-theoretic, and spectral measures.
- **Validation**: Tools for testing reliability across sessions and resampling.
- **Visualization**: Plotting and reporting functions.
- **Extensible**: Modular design allows integration of new metrics or algorithms.

---

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/EhsanEqlimi/EEG-microstate-validity.git
````

2. Add the repository to your MATLAB path:

   ```matlab
   addpath(genpath('EEG-microstate-validity'))
   ```
3. Requirements:

   * MATLAB R2018b or later
   * Signal Processing Toolbox
   * Statistics and Machine Learning Toolbox

---

## Usage

### Basic example

```matlab
% Load sample data
data = load('examples/sample_data/subj01.mat');

% Run segmentation and clustering
results = run_full_pipeline(data);

% Visualize results
plot_microstates(results);
```

### Advanced example

See `examples/demo_advanced.m` for parameter tuning, reproducibility testing, and multi-subject validation.

---

## Input & Output

* **Input**: Preprocessed EEG data (MATLAB `.mat` format).

  * Ensure artifact removal and re-referencing are completed.
* **Output**:

  * Cluster maps and time series
  * Reliability and validity metrics
  * Figures and summary tables

---

## Parameters

Key parameters include:

* `nStates`: Number of microstates to extract (default = 4)
* `clusteringMethod`: K-means, modified K-means, or hierarchical
* `stabilityTest`: Bootstrapping, split-half, or session-based
* `sparsityConstraint`: Optional sparsity regularization

Detailed parameter descriptions: [docs/parameter\_tuning.md](docs/parameter_tuning.md)

---

## Validation & Metrics

Implemented measures:

* **Variance explained** (global and per-state)
* **Cluster stability** (via resampling)
* **Information-theoretic metrics** (entropy, mutual information)
* **Spectral features** of microstate sequences
* **Transition probabilities** and temporal dynamics

---

## Examples

* `examples/demo_basic.m`: Minimal pipeline demonstration
* `examples/demo_advanced.m`: Full reproducibility and robustness analysis

---

## Best Practices

* Always preprocess EEG (artifact removal, ICA, filtering) before using the toolbox.
* Validate the chosen number of states using multiple metrics.
* Compare across sessions to ensure reproducibility.
* Report both reliability and validity metrics in publications.

---

## Extending the Framework

To add a new metric or clustering algorithm:

1. Place your function in the relevant `src/` subdirectory.
2. Update the pipeline script to call your method.
3. Add documentation in `docs/`.

---

## Citation

If you use this toolbox in your research, please cite:

```
Eqlimi, E. (Year). EEG-microstate-validity: A framework for evaluating EEG microstate reliability and validity. GitHub. https://github.com/EhsanEqlimi/EEG-microstate-validity
```

---

## License

[MIT License](LICENSE)

---

## Contact

* Author: Ehsan Eqlimi
* GitHub: [@EhsanEqlimi](https://github.com/EhsanEqlimi)
* Issues: Please report via the [GitHub Issues page](https://github.com/EhsanEqlimi/EEG-microstate-validity/issues)

```
```

