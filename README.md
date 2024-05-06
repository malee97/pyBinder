Welcome to pyBinder! Created by Michael Lee
If you run into any problems, please email me at mlee97@mit.edu :)

Features:
   - Centroids .mzML files 
   - Generates plots for data checks including TIC and profile maps
   - Outputs feature maps (.featureXML) files for every file input
   - Consolidates feature maps into a consensus map for all total peptides found across all replicates and proteins
   - Integrates feature peaks and performs enrichment and statistical analysis
   - Outputs top protein-enriched features in Orbi-usable inclusion list 
   - Plots user-defined amount of EICs for quality control

Installation and use:
1. Download repository to local disk - Note: OpenMS (package used to handle .mzML files) does not work well with files stored externally. Local storage is recommended.
2. While there is built-in functionality for converting Thermo .RAW files to .mzML files using MSConvert from the the ProteoWizard package, we recommend converting files separately and working only with .mzML files
3. Open Jupyter notebook (.ipynb) file, making sure the "Methods.py" file is in the same directory as the notebook
4. Give inputs as requested in the top cells
   Notes on inputs:
   - Separate protein names by commas, do not have any spaces
   - Ensure that your entered protein names are present in the names of the actual RAW files (i.e. if input 12ca5, have 12ca5.RAW)
   - If doing variable protein concentrations, have a consistent naming convention throughout all samples (i.e. include "High", "Medium", and "Low" in all samples where appropriate)
   - Additionally for variable protein concentrations, make sure to have a "beads only" control named as such
   - Selection name decides the name of the results folder
   - Length of library sets a baseline observed mass filter (peptides below mass of a sequence of X glycines are removed)
   - Enrichment score cutoff determines features listed on final inclusion list, perfectly non-specific is 1/X, where X = number of proteins 
   - P score cutoff is for statistical analysis, depending on number of replicates I would recommend using 0.1 or 0.15 
   - Minimum peak height for secondary peak finder is used for second round of feature finding, generally just use 1E5 unless you know your data
   - Keep number of replicates consistent across proteins (don't have 2 replicates for one but 3 for another)
   - Number of EICs can be as high as you want (usually do 100-150)
5. pyBinder can take a fairly long time to run (1-2 hours) depending on data quality, keep an eye on memory usage
