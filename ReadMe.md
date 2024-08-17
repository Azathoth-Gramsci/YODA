The following is a guide to the code and data used for data analysis.

Timeseries bioreactor data is present in the form of CSVs for ODs, dilutions, and single cell fluorescence. The data can be loaded as a reactordata object.
All the functions used for analysis utilise this class.
For leakage and efficiency experiment raw cytometry data is provided.

lib folder includes the class definitions (cytoreactors) and methods (tools) required for preliminary treatment of the data like gating, doublet removal and deconvolution.
lib folder also contains functions for parsing the raw cytometry files (fcsparsing).

The jupyter notebooks consist of scripts and functions that were used to load and analyse the data to generate the figures reported in the paper.
These are separated by figures (check parentheses in the notebook name)
The code was developed in Windows 10 and verified by running in Windows 7.

Depenendencies for the code to run without hiccups,

Python 3.5+

Pandas 1.2.2+
Numpy 1.20.0+
Matplotlib 3.3.4+
Seaborn 0.11.1+
SciPy 1.6.0+
SciKit 0.24.2+
tkinter 8.6+
pickle 4.0
os
sys
stat

The raw_data directory available from (Zenodo- http://doi.org/10.5281/zenodo.4923833) must be placed in the same directory as this repository. Figures and processed data are exported to empty directories included in this repository (Plots and Processed_data respectively).

In addition to the analysis code, the two pickle files contained within this directory are requisite for deconvolution. These are the spectral signatures for different fluorescence proteins and the Autofluorescence values for WT strain.

Courage!
