%% recon of one image (one compound only)


clear all; close all; clc;
load('L:\basic\divi\Projects\cosart\CS_simulations\chemical_shift_fluor19_sims\PFCEandPFOB\PFCEandPFOB\kspaces_readout1&2.mat')
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
%%
rr1 = @(I) reshape(I,[128,128]);


%% preprocess data to reconstruct (there is only one shared phase-encoding direction...)
KC1=KC1./max(abs(KC1(:)));  %normalize
KC2=KC2./max(abs(KC2(:))); %normalize

imagine(abs(fftshift(fftn(KC1))))
I1=abs(fftshift(fftn(KC1)));
figure(1); imshow(squeeze(I1(51,:,:)),[]); title('good image to do the test: x-slice 51') % GOOD IMAGE

%% 
KR=fftshift(ifft(KC1,[],1),1);
KR51=squeeze(KR(128-51,:,:));
I51=abs(ifft2(KR51));

figure(2); imshow(abs(ifft2(KR51)),[])

%% spectrum 

[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(282,44642/128) %126.4??

Spectrum=zeros(1,128);
Spectrum(round(PFOB)+abs(min(round(PFOB)))+1)=PFOB_alpha./sum(PFOB_alpha(:))
Spectrum=flipdim(Spectrum,2);
Spectrum=circshift(Spectrum,[0 -(92+6)]) % shift such that CF3 peak is at -7 pixel shift

%% blind deconvolution

[J P] = deconvblind(I51,Spectrum,50)
figure(3); subplot(221); imshow(I51,[]); subplot(222); stem(P); title('result of BLIND DCV');
subplot(223); imshow(J,[]); subplot(224); stem(Spectrum); title('theory')

%% use same PSF in CG algo 
A=opConvolve(128,128,Spectrum,[0 0],'cyclic')
F=opDFT2(128,128,1); %fourier operator

CGK=nl_conjgrad_fluor(F*A,KR51(:),zeros(size(KR51(:))),50,zeros(size(KR51(:))),0,128,128); 
figure(4); subplot(224); imshow(rr1(abs(CGK)),[]); axis off; title('CG-k space')



%% 
KR51ABS=(fft2(I51))

CGK=nl_conjgrad_fluor(F*A,KR51ABS(:),zeros(size(KR51(:))),50,zeros(size(KR51(:))),0,128,128); 
figure(4); subplot(224); imshow(rr1(abs(CGK)),[]); axis off; title('CG-k space')





