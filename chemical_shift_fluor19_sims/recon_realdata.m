clear all; close all; clc;
load('L:\basic\divi\Projects\cosart\CS_simulations\chemical_shift_fluor19_sims\PFCEandPFOB\PFCEandPFOB\kspaces_readout1&2.mat')
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
%%

rr = @(I) reshape(I,[128,256]);
vec= @(I) reshape(I,[128^2, 1]);
rr1 = @(I) reshape(I,[128,128]);

%% preprocess data to reconstruct (there is only one shared phase-encoding direction...)
KC1=KC1./max(abs(KC1(:)));  %normalize
KC2=KC2./max(abs(KC2(:))); %normalize

k1=ifftshift(ifft(KC1,[],2),2);
k2=ifftshift(ifft(KC2,[],2),2);
sl=33
; % slice to recon
k1=squeeze(k1(:,sl,:));
k2=squeeze(k2(:,sl,:));

k1=fft2(abs(ifft2(k1))); %remove phase info
k2=fft2(abs(ifft2(k2))); % remove phase info
data=[vec(k1);vec(k2)]; 

F=opDFT2(128,128,1); %fourier operator
F2=opBlockDiag(F,F)  % 2 fourier ops - one for each image 
linear_recon=opInverse(F2)*data; %linear recon of data 
figure(1); imshow(abs(rr(linear_recon)),[]) %not totally sure about fftshift yet.... 
%% empirical spectrum 

if true
    % make Spectrum
%     [PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(282,44642/125.5) %126.4??
    [PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(282.56877,44642/128) %126.4??
    
    Spectrum=zeros(1,128);
    Spectrum(round(PFOB)+abs(min(round(PFOB)))+1)=PFOB_alpha./sum(PFOB_alpha(:))
    Spectrum=flipdim(Spectrum,2);
    Spectrum=circshift(Spectrum,[0 -(92+6)]) % shift such that CF3 peak is at -7 pixel shift
    
else
    %find PSF by blind deconvolution
    Spectrum_e=mean(abs(linear_recon_reshape(:,129)),2).';
    Spectrum_e=Spectrum_e./sum(Spectrum_e(:));
    Spectrum_e=Spectrum_e.*(Spectrum_e>0.03); % empirical spectrum
    
    [J P] = deconvblind(abs(linear_recon_reshape(1:128,129:133)),Spectrum_e.',10)
    offset=-min(round(PFOB));
    [J2 P2] = deconvblind(abs(linear_recon_reshape(10:20,1:128)),Spectrum_e,10)
    Spectrum=circshift(P.',[0 -17]);
end


figure(2); subplot(312); stem(Spectrum); xlabel('voxels')
ylabel('relative amplitude'); title('PFOB PSF')

% make operators (include fftshift in operators for simplicity)
A=opConvolve(128,128,Spectrum,[64 64],'cyclic') %perhaps 65?

Spectrum_vert=Spectrum.'; figure(99); imshow(Spectrum_vert,[])
A2=opConvolve(128,128,Spectrum_vert,[64 64],'cyclic')

Spectrum2=[1]; 
B=opConvolve(128,128,Spectrum2,[64 64],'cyclic') 
B2=opConvolve(128,128,Spectrum2.',[64 64],'cyclic')

% deconvolve both images
deconv_image1=opInverse(F*A)*data(1:128^2);
deconv_image2=opInverse(F*A2)*data(128^2+1:end);

figure(3); 
subplot(211);imshow(abs(rr1(deconv_image1)),[])
subplot(212);imshow(abs(rr1(deconv_image2)),[])

%% recon of one image (one compound only)
CGK=nl_conjgrad_fluor(F*A,data(1:128^2),zeros(size(data(1:128^2))),50,zeros(size(data(1:128^2))),0,128,128); 
figure(4); subplot(224); imshow(rr1(abs(CGK)),[]); axis off; title('CG-k space')


%%
M1=[F*A,F*B];
M2=[F*A2,F*B2];
M=[M1;M2] %measurement operator 

CGK=nl_conjgrad_fluor(M,data,zeros(size(data)),150,zeros(size(data)),0,128,256); 
figure(2); imshow(rr(abs(CGK)),[0 5e-2]); axis off; 
%%

figure(1); imshow(abs(fftshift(rr1(linear_recon(128^2+1:end,:)))),[])