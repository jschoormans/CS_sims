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
    cd('linuxpath')
    
end
files=dir('*.mat')
load(files(2).name)

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

%% RECON HERE
%K should be [kx ky kz ncoils nNSA]
% FOR CERTAIN ACC FACTORS - LOOP OVERNIGHT

PR=struct;
PR.outeriter=4;
PR.Itnlim=10;
PR.noNSAcorr=false;
PR.TVWeight=1e-3;
PR.TGVfactor=0;
PR.xfmWeight=2e-3;
PR.reconslices=1;
PR.squareksp=true;
PR.resultsfolder=''
% PR.sensemaps=permute(se ns,[2 3 1 4]); 
PR.sensemapsprovided=0

accvector=[1,2,3,4,5,6];
nNSA=5;

for kk=1:nNSA
for jj=1:3
for ii=1:length(accvector);
    P.usedyns=kk; %for example
    P.acc=accvector(ii);
    P.jjj=jj+4 % 5 6 7
    P.noiselevel=0
    [~,~, KD]=makeNoisyKspacefromdynamics(K,P);
    
    Ko=KD.Ku_N2;
    Ko=permute(Ko,[5 1 2 3 4]);
    R=reconVarNSA(Ko,PR); 
    %append/save/clear to save memory??
    save(['R_',num2str(ii),'_',num2str(jj),'_',num2str(kk)],'R')
end
end
end

