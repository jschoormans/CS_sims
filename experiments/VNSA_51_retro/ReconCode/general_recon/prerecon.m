% prerecon: saves K.mat from .raw data from scanner

addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code'))
addpath(genpath('/opt/amc/bart-0.3.00'));vars;
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Excercise'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/codefrombram/various/export_fig'))

clear all; close all; clc;

inp=1; i=1
while inp==1 
FF{i} = uigetfile('*.raw')
inp = input('another one? type 1')


if inp==1
   i=i+1
   
else
    break;
end
end



for i=1:length(FF)
MR=Recon_varNSA(FF{i})
MR.Parameter.Parameter2Read.typ=1;

MR.Parameter.Recon.ArrayCompression='No'
MR.Parameter.Recon.ACNrVirtualChannels=6;
% MR.Parameter.Parameter2Read.chan=[17]
MR.ReadData;
MR.RandomPhaseCorrection;
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

if max(MR.Parameter.Labels.Index.dyn)~=0; % if dynamics are used, convert to aver;
    MR.Parameter.Labels.Index.aver=MR.Parameter.Labels.Index.dyn;
    MR.Parameter.Labels.Index.dyn=zeros(size(MR.Parameter.Labels.Index.dyn));
end
if max(MR.Parameter.Labels.Index.aver)~=0; % if dynamics are used, convert to aver;
% do nothing; 
else
MR.addAveragestoLabels
end

MR.Parameter.Recon.ImmediateAveraging='No'
MR.SortData;
K=squeeze(MR.Data);
%% find filename
startindex=regexp(MR.Parameter.Filename.Data,'/');
startindex=1;
endindex=regexp(MR.Parameter.Filename.Data,'.raw')

filename=MR.Parameter.Filename.Data(startindex(end)+1:endindex-1)
if strcmp(MR.Parameter.Recon.ArrayCompression,'Yes')==1
    filenameK=['K_',filename,'_VC',num2str(MR.Parameter.Recon.ACNrVirtualChannels),'.mat']
else
    filenameK=['K_',filename,'_noVC.mat']
end
%% SAVE K
% cd('/scratch/jschoormans/saved_K')
save(filenameK,'K','-v7.3')
disp('first part finished!')

end