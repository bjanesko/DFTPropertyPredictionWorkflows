# DFTPropertyPredictionWorkflows
Computational workflows for predicting pKa, distribution coefficients (logD), and related properties. Workflows use physics-based methods, typically automated conformational search plus DFT structure/energy refinement

The workflow for pKa prediction is as follows. We illustratre for fumaricacid.xyz, the fully protonated structure of fumaric acid. Run  as

python3 writeCREST.py fumaricacid.xyz --qfull 0 --nsolv 0 --prots 7 11

[run all CREST jobs]
[move into directory fumaricacid_0w.xyz]

perl convertCRESTToGaussian.prl

[run all Gaussian jobs]

python3 readGaussian.py

This should predict fumaric acid pKas 7.06 and 4.27. 

# DFT prediction workflow (clogP and logD)
The workflow consists of two main main scripts:
readPartitionCoefficient.pl – Calculates the octanol/water partition coefficient (clogP)
titra_curve_cal.py – Generates titration curves and computes logD values 
All from Gaussian output files

Usage: 
1. Calculate Partition Coefficient (clogP)
Run the Perl script to compute the partition coefficient between water and octanol:
perl readPartitionCoefficient.pl
This script processes the outputs the calculated clogP values.

2. Generate Titration Curves and Calculate logD
Run the Python script to compute pH-dependent logD values and generate titration curves:
python3 titra_curve_cal.py
By default, the script operates on the current working directory containing the required subdirectories (neutral, plusOne, twoplus).
To view all available command-line options:
python3 titra_curve_cal.py --help

3. Analyze Hbonds and calculate total surface area
Run these Python scripts
python3 hbondanalyzer.py and python3 total_surface_area.py

Sample of the Gaussian input and output files are in Gaussian_input_and_output.zip
