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
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations'))
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'))
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
        addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'))

    experimentfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_51_retro/',date]
    mkdir(experimentfolder)
    cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA_51/VNSA_51')
    vars
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

%% make k-spaces I DONT NEED FUKCING NOISE ONLY UNDERSAMPLKING AND AVERAGING 
P.usedyns=5; %for example
P.acc=4;
P.jjj=5 % 5 6 7 
P.noiselevel=0
[KD{1}, KD{2}, KD{3}]=makeNoisyKspacefromdynamics(K,P);

%% RECON HERE
%K should be [kx ky kz ncoils nNSA]
% FOR CERTAIN ACC FACTORS - LOOP OVERNIGHT

TVvec=[0,1,2,5,10].*1e-3
xfmvec=TVvec; 

PR=struct;
PR.outeriter=4;
PR.Itnlim=10;
PR.noNSAcorr=false;
PR.TGVfactor=0;
PR.reconslices=1;
PR.squareksp=true;
PR.resultsfolder=''
% PR.sensemaps=permute(se ns,[2 3 1 4]); 
PR.sensemapsprovided=0
for jj=1:length(TVvec)
    for ii=1:length(xfmvec);
        PR.TVWeight=TVvec(jj);
        PR.xfmWeight=xfmvec(ii);
        
        
        Ko=KD{3}.Ku_N2;
        Ko=permute(Ko,[5 1 2 3 4]);
        R{2,ii,jj}=reconVarNSA(Ko,PR)
        R{2,ii,jj}.sensemaps=[];

    end
end
filename='exp_5_51_Weights.mat'
cd(experimentfolder)
save(filename,'R','P','PR','-v7.3')

%% postproc
ctr=1;
figure(1); 
for jj=1:length(TVvec)
    for ii=1:length(xfmvec);
        subplot(5,5,ctr);ctr=ctr+1;
imshow(fftshift(abs(R{2,ii,jj}.recon),1),[]); axis off; 
        
    end
end


%%
kx=200:300; ky=200:300;
ctr=1;
figure(2); 
for jj=1:length(TVvec)
    for ii=1:length(xfmvec);
        subplot(5,5,ctr);ctr=ctr+1;
imshow(abs(R{2,ii,jj}.recon(kx,ky)),[]); axis off; 
        
    end
end












