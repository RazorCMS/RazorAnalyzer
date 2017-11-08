Procedure for deriving QCD transfer factors:
1) Run `calculateMCScaleFactors.C` to add the MC scale factor weights to the background MC trees
2) Run `getInclusiveTransferFactors.py` to create the ROOT file containing the raw transfer factors
3) Run `QCDFits.py` to perform linear fits to the transfer factors for extrapolation to higher Rsq
4) Run `getBTagFractions.py` and copy the resulting output file to the scale factors folder for the extrapolation in nbtags
