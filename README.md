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
