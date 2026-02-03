# Contributing to scRNA-seq Analysis Pipeline

Thank you for your interest in contributing to this single-cell RNA-seq pipeline! This project was developed during a research internship at the Adebali Lab, Sabancı University, and aims to be a robust tool for the bioinformatics community.

## How Can You Contribute?

### 1. Reporting Bugs
If you find a bug (e.g., issues with Scanpy versions or ForceAtlas2 rendering), please open an "Issue" in the GitHub repository. Include the error log and your system specifications.

### 2. Suggesting New Features
We welcome suggestions for new analytical steps. For example:
* Adding new Marker Gene lists for different species.
* Implementing alternative batch-correction algorithms (e.g., Harmony, scVI).
* Enhancing the Trajectory visualization plots.

### 3. Submitting Pull Requests (PR)
If you've written code to improve the pipeline:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature/NewFeature`).
3. Commit your changes (`git commit -m 'Added scVI batch correction'`).
4. Push to the branch (`git push origin feature/NewFeature`).
5. Open a Pull Request!

**Note:** Please ensure any new Python code follows the PEP8 style guide and successfully runs on Python 3.10.
