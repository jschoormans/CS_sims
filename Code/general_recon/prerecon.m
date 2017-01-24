% prerecon: saves K.mat from .raw data from scanner

cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/20161206_ICVW_VSNA')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code'))
addpath(genpath('/opt/amc/bart-0.3.00'));vars;
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Excercise'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/codefrombram/various/export_fig'))

clear all; close all; clc;
MR=MRecon()
MR.Parameter.Parameter2Read.typ=1;

MR.Parameter.Recon.ArrayCompression='Yes'
MR.Parameter.Recon.ACNrVirtualChannels=6;
% MR.Parameter.Parameter2Read.chan=[17]

MR.ReadData;
MR.RandomPhaseCorrection;
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
K=MR.Data;
%% find filename
startindex=regexp(MR.Parameter.Filename.Data,'/')
endindex=regexp(MR.Parameter.Filename.Data,'.raw')

filename=MR.Parameter.Filename.Data(startindex(end)+1:endindex-1)
if strcmp(MR.Parameter.Recon.ArrayCompression,'Yes')==1
    filenameK=['K_',filename,'_VC',num2str(MR.Parameter.Recon.ACNrVirtualChannels),'.mat']
else
    filenameK=['K_',filename,'_noVC.mat']
end
%% SAVE K
save(filenameK,'K')
disp('first part finished!')