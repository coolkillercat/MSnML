# Prerequisite
## Data Preparation
Mass spectra are stored in a folder which contains several CSV files.
* __folder name__ is in the format of ```test_{Product Name}_{MS level}``` or ```test_NMF_template_{Product Name}_{MS level}```. The first one indicates that the spectra are aligned to theoretical peaks. The latter one means spectra aligned with machine learned template.
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
* Put the preprocessed data into ```data``` folder(examples provided).
* Run notebook ```MS5_regression.ipynb``` for predicting 2-3 to 2-6 ratio on MS5.
* Run notebook ```MS234_classification.ipynb``` for predicting 2-3/2-6 linkage type using MS2/3/4.
* Modify ```datasets``` variables to use different Ms data.