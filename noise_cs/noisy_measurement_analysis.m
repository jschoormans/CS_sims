%Noisy measurement analysis 30-5
clear all;close all; clc
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/noise_cs')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations'))
MR=MRecon('cs_highres_noisy/20_30052016_1952515_3_2_wipmprageadnisenseV4.raw')
%%
MR.Parameter.Parameter2Read.dyn=[0]'
MR.Parameter.Parameter2Read.typ=[1]
MR.ReadData;
MR.RandomPhaseCorrection;
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
% MR.K2IM;

K=MR.Data;
%%
%make sampling pattern
acc=3
[pdf,val] = genPDF([size(K,2) size(K,3)],5,1/acc,2,0,0);
Mfull=genEllipse(size(K,2),size(K,3));
Mfull=repmat(Mfull,[1 1 size(K,1)]);
Mfull=permute(Mfull,[3 1 2]);
Mfull=repmat(Mfull,[1 1 1 size(K,4) size(K,5)]); %add coils
M=permute(repmat(genSampling(pdf,10,100),[1 1 size(K,1) size(K,4) size(K,5)]),[3 1 2 4 5]).*Mfull;

%% GENERAL CS RECON  VS FULL RECON
%%
[calib emaps]=bart('ecalib -r 10',K)
sens=bart('slice 4 0',calib)
%%
IMCS=bart(['pics -RW:7:0:0.05 -S  -e -i20 -d5'],K.*M,sens);
IMFull=bart(['pics -RW:7:0:0.05 -S  -e -i20 -d5'],K,sens);

IMCS=fftshift(IMCS,2);
IMFull=fftshift(IMFull,2);

figure(1)
imshow([abs(squeeze(IMCS(:,:,50))) abs(squeeze(IMFull(:,:,50)))],[])
export_fig 
%%
MR2=MR.Copy;
MR2.K2I; MR2.CombineCoils;
figure(2)
imshow(abs(squeeze(MR2.Data(:,:,25))),[])



%% RECON CS
sens=ones(size(K,1),size(K,2),size(K,3),size(K,4));
[kspace, traj]=calctrajBART(K); 
Recon=bart(['pics -RW:7:0:0.05 -S  -e -i20 -d5 -t'],traj./2,permute(kspace,[3 1 2 4]),sens);

%%
figure(1)
imshow(abs(squeeze(Recon(:,:,10))),[])


%% RECON NUFFT
ReconNufft=bart('nufft -i -t',traj(:,:,1:100)./2,permute(kspace(:,1:100,1,1),[3 1 2 4]));

%%
figure(2)
imshow(abs(squeeze(ReconNufft(:,:,20))),[])













