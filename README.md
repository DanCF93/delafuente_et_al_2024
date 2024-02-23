# delafuente_et_al_2024
Code associated to the publication from De La Fuente et al. for the cell cycle analysis.
It contains the following files or pipelines:
  1. CP_pipelines
  2. data_analysis

The pipelines have been created to process and analyse images taken on an Perkin Elmer Opera Phenix. If you have used a different instrument, the regex associated to each metadata element will need to be updated accordingly.

Use as follows:
1. Calculate illumination correction function from all individual channels. This is necessary to correct for variations in fluorescence signal that arise from the hardware setup.
2. Analyse the images to quantify EdU incorporation and any other associated markers (e.g. KI67, pH3).
3. Use the R script in data_analysis to process the resulting .csv files from your image analysis and calculate the S-phase length and TC.
