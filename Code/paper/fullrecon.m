MR=MRecon
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.CoilCombination='yes' 
MR.Parameter.Recon.ACNrVirtualChannels=5;
% Produce k-space Data (using existing MRecon functions)
MR.ReadData;
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.GridData;
MR.RingingFilter;
MR.ZeroFill;
K=MR.Data;
%%
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Excercise'))
P=struct;
P.reconslices=200:220;
P.TVWeight=0.001;
P.xfmWeight = 0;
P.resultsfolder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/20161206_ICVW_VSNA/results' 
P=reconVarNSA(K,P) 