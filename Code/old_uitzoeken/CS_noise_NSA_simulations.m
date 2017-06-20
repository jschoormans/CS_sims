clear all; close all; clc; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/opt/amc/bart')); vars;

ny=405
nz=512

Im=imread('konijn.bmp');   %image in image domain
Im=Im(:,:,1);
K=fftshift(fft2(Im)); K=permute(K,[3 1 2]);


bartoptions=['poisson -v -e -Y',num2str(ny),' -Z',num2str(nz)]
M=bart(bartoptions);
% acceleration=sum(M(:)==1)./
Ku=M.*K;
%
ImCS=bart('pics -l1 -d5 -i10',K,ones(1,405,512,1));
Imlin=ifft2(squeeze(Ku));


figure(1);
subplot(221); imshow(Im); subplot(222); imshow(squeeze(M))
subplot(223); imshow(abs(squeeze(Imlin)),[]); subplot(224); imshow(abs(squeeze(ImCS)),[])
