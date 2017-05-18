%Noisy measurement analysis 30-5
clear all;close all; clc
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/noise_cs')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
MR=MRecon('cs_highres_noisy/20_30052016_1952515_3_2_wipmprageadnisenseV4.raw')
coil_survey='20_30052016_1952161_1000_8_wipcoilsurveyscanV4.raw'
sense_ref='20_30052016_1952367_1000_11_wipsenserefscanclearV4.raw'

%%
MR.Parameter.Parameter2Read.dyn=[0:1:11]'
MR.Parameter.Parameter2Read.typ=[1]
MR.ReadData;
MR.RandomPhaseCorrection;
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
% MR.K2IM;

KF=MR.Data;
%% estimate sense maps;;
Kfull=mean(KF,5); 
Kfull=bart('fftmod 2 -i',Kfull);
ImFFT=bart('fft 7 -i',Kfull);
ImFFT=bart('fftmod 2 -i',ImFFT);
figure(1); imshow(abs(squeeze(ImFFT(:,:,50,12))),[])

%%
%make sampling pattern
K=mean(KF(:,:,:,:,1),5);

acc=5
[pdf,val] = genPDF([size(K,2) size(K,3)],5,1/acc,2,0,0);
Mfull=genEllipse(size(K,2),size(K,3));
Mfull=repmat(Mfull,[1 1 size(K,1)]);
Mfull=permute(Mfull,[3 1 2]);
Mfull=repmat(Mfull,[1 1 1 size(K,4) size(K,5)]); %add coils
M=permute(repmat(genSampling(pdf,10,100),[1 1 size(K,1) size(K,4) size(K,5)]),[3 1 2 4 5]).*Mfull;

%% GENERAL CS RECON  VS FULL RECON

%%
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/noise_cs/results')

K1=mean(KF(:,:,:,:,1:10),5); %CS acc 5; NSA=10
K2=mean(KF(:,:,:,:,1:2),5); %full; NSA=2
sens=bart('caldir 50',K1);%make sense maps from dat

IMCS=bart(['pics -RW:7:0:0.01 -S  -e -i50 -d5'],K1.*M(:,:,:,:),sens);
IMFull=bart(['pics -RW:7:0:0.01 -S  -e -i50 -d5'],K2,sens);

IMCS=fftshift(IMCS,2);
IMFull=fftshift(IMFull,2);

figure(1)
imshow([abs(squeeze(IMCS(:,:,1))) abs(squeeze(IMFull(:,:,1)))],[])
title(['CS acc=',num2str(acc),' NSA=10, versus full recon NSA=2; all coils'])
export_fig -native recon.eps
