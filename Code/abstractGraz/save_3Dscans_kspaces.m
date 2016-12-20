% LOAD AND SAVE KSPACE OF 3D SCANS BECAUSE THE SERVE IS TOO FUCKING SLOW :(

%%
clear all; close all; clc; 
disp('adding paths and setting up stuff...')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/opt/amc/matlab/toolbox/MRecon'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
addpath(genpath('/opt/amc/bart')); vars
% addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab'))


%% defining parameters 
clear dir;
disp('setting parameters')
P=struct;
P.chan=1;
P.folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/CS-3D-prospective-grapefruit/'
cd(P.folder)
files=dir('*.raw');



for filenumbers=[6,7,9,10,11,12]
    clear K MR
P.file=files(filenumbers).name

%% PHANTOM RESOLUTION

MR=MRecon([P.folder,P.file]); 
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Recon.ImmediateAveraging='no'
MR.ReadData;
disp('random phase corr')
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
disp('sorting')
MR.SortData;

% MR.Average;
MR.RemoveOversampling;
MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
 disp('data to k-space')
 %%
K=MR.Data;

cd(P.folder)

save(['K',num2str(filenumbers),'.mat'],'K','-v7.3')

end


%% FOR THE GOLD STANDARD


for filenumbers=[8]
    clear K MR
P.file=files(filenumbers).name

%% PHANTOM RESOLUTION

MR=MRecon([P.folder,P.file]); 
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Recon.ImmediateAveraging='yes'
MR.ReadData;
disp('random phase corr')
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
disp('sorting')
MR.SortData;

% MR.Average;
MR.RemoveOversampling;
MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
 disp('data to k-space')
 %%
K=MR.Data;

cd(P.folder)

save([files(filenumbers).name,'K'],'K')

end