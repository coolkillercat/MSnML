# Prerequisite
## Data Preparation
Mass spectra are stored in several CSV files.
* __mass_features.csv__ stores the learned mass features(m/z values) which are used by all spectra.
* __input.csv__ stores the peak intensities, each row represents one spectrum.
* __output.csv__ stores the ground truth for the relative a2'-3' to a2'-6' ratio.
* __metadata.csv__ stores the configuration and sample information for each mass spectrum.
We provide MATLAB code to process Thermo .raw files into the prepared .csv files.
## Software Requirements
Python 3, Pip, jupyter notebook
## Package Installation
Run the following code to install the packages
```
pip install -r requirements.txt
```
# Use guide
