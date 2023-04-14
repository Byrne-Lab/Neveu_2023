# Neveu_etal_2023

The main folder contains all of the code used to collect and analyze the current traces of the various ion channels.  The channel is indicated in the file names.  Each .m file has a corresponding .mat file.  The .mat files are the precollected data and the .m file is the code to analyze the data file.  The files needed for the SNNAP model are in the Neuron_models folder.  

To build the SNNAP model, download the files in the makeSNNAP directory.  Run makesnnap('win','excelfilepath') or makesnnap('mac','excelfilepath') to populate the SNNAP model.  Run the varius simulations.  Then compare the results to the emprical data using the checksim function. 