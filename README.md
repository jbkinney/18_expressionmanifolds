# Data and analysis scripts for:
## *Precision dissection of cis-regulatory energetics in living cells.* 
## Forcier et al., bioRxiv preprint (in process)

Directory structure:
* ``code/``: Python analysis scripts and Jupyter notebooks
  * ``library_resampler.py``: Script for generating all bootstrap resamplings of the spacing libraries
  * ``plate_processor.py``: Script for generating all intermediate files from raw data
  * ``figures.ipynb``: iPython notebook for generating the files in ``figures/``
* ``data/``: Metadata, raw data, and processed data for this study
  * ``constructs/``: DNA sequence constructs used in this study
  * ``literature/``: Data extracted from the prior literature
  * ``plate_reader/``: Raw data from plate reader used for Miller assays
  * ``glycerol_stocks.xlsx``: List of glycerol stocks used in this study
* ``figures/``: Computationally rendered components of the figures in the manuscript.
* ``intermediate/``: Processed data
  * ``resamplings/``: Files of fits to bootstrap resamplings of the spacing libraries
* ``protocols/``: Primary protocols used in this study
