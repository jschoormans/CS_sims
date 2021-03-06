%% KLADBLOK EXPERIMENT VNSA 51
clear all; close all; clc;

%% load k-space
if ispc()
    cd('L:\basic\divi\Ima\parrec\Jasper\VNSA\VNSA_51\VNSA_51')
    addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations'))
    addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\Wavelab850'))
    addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'))
    experimentfolder=['L:\basic\divi\Projects\cosart\CS_simulations\experiments\VNSA_51_retro\',date]
    mkdir(experimentfolder)
else
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
addpath(genpath('/opt/amc/bart')); vars
cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA_51/VNSA_51')
    
experimentfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_51_retro/',date]
mkdir(experimentfolder)

end
files=dir('*.mat')
load(files(1).name)

%% make sense maps
Kref=mean(K,4); %still use coils though
Kref=permute(Kref,[4 1 2 3]);

sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?

imagine(sens)
%% Imref
reg=0.01;
ImRef=fftshift(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Kref,sens)),2); %to compare image
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
figure(1); imshow(squeeze(abs(ImRef)),[])
%% TODO: HOW ABOUT NOISE CORRELATION (CHECK PRUESSMAN PAPER) 

%% make k-spaces I DONT NEED FUKCING NOISE ONLY UNDERSAMPLKING AND AVERAGING 
accvector=[1,2,3,4,5,6];

P.usedyns=5; %for example
for jj=1:3
for ii=1:length(accvector);
P.acc=accvector(ii);
P.jjj=jj+4 % 5 6 7 
P.noiselevel=0
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspacefromdynamics(K,P);
end
end

%% RECON HERE
%K should be [kx ky kz ncoils nNSA]
% FOR CERTAIN ACC FACTORS - LOOP OVERNIGHT

PR=struct;
PR.outeriter=4;
PR.Itnlim=10;
PR.noNSAcorr=false;
PR.TVWeight=(0);
PR.TGVfactor=0;
PR.xfmWeight=5e-3;
PR.reconslices=1;
PR.squareksp=true;
PR.resultsfolder=''
% PR.sensemaps=permute(se ns,[2 3 1 4]); 
PR.sensemapsprovided=0
for jj=1:3
for ii=1:length(accvector);
Ko=KD{3,ii,jj}.Ku_N2;
Ko=permute(Ko,[5 1 2 3 4]);
R{2,ii,jj}=reconVarNSA(Ko,PR)
end
end
%% SAVE R (WITHOUT SENSEMAPS)
for jj=1:3
for ii=1:length(accvector); 
    
R{2,ii,jj}.sensemaps=[];
end
end
filename='exp_3_51_11_usedyns5_acc1-6'
cd(experimentfolder)
save(filename,'R','P','PR','-v7.3')

