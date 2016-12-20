%develop and test conjugate gradient algo
% purpose: to be able to recon 3d VARIABLE DENSITY NSA measurements

%% LOAD KSPACE

% cd('/home/jschoormans/lood_storage/divi/Projects/cosart/scans/CS-3D-prospective-grapefruit')

%% PC
clear all; close all; clc;
cd('L:\basic\divi\Projects\cosart\scans\CS-3D-prospective-grapefruit')
addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\Code'))
addpath(genpath('C:\Users\jschoormans\Dropbox\phD\bart-0.3.01')); 
addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\sparseMRI_v0.2'))

%%
tic; disp('loading data...')
filenumber=6
filename=['K',num2str(filenumber),'.mat']
load(filename); toc; 

%% FFT IN MEASUREMENT DIRECTIONS
tic; disp('FFT in measurement direction')
K=ifft(K,[],1);
K=K./max(K(:)); %normalize kspace
K=squeeze(K);
toc;

%% make mask
tic; disp('make mask and setting up data...')
[nx ny nz nc nNSA]=size(K);
Ks=squeeze(K(100,:,:,1,:)); %k-space for one channel and one slice
fullmask=Ks~=0;             %find mask used for scan (nx*ny*NSA)
MNSA=sum(fullmask,3);       %find NSA for all k-points in mask 

data=sum(K,5);              %sum of data over NSA
data=data./permute(...
    repmat(MNSA,[1 1 nx nc]),[3 1 2 4]);            %mean of data over NSA 2    
data(isnan(data))=0;        %clear up NaN values (due to /0)

data=squareksp(data);       %make k-spsace square and size a 2^n (bit buggy, not needed for nx now)
MNSA=squareksp(MNSA);       %make MNSA square and size of 2^n

mask=double(MNSA>0);        %2d mask (no NSA dimension)
pdf=ones(size(data));       %pdf is used for first guess; should be fixed!
toc


%% CONJUGATE GRADIENT RECONS
tic; disp('recondata')
sl=50;
recondata=data(sl,:,:,:); % data of one slice to be used in recon 
toc;
%% ESTIMATE SENSE MAPS
disp('Estimating sense maps');tic 
sensemaps=bart('ecalib -m1',bart('fftmod -i 7',recondata));
figure;imshow(squeeze(sensemaps(1,:,:,4)))
toc
%%
tic; disp('setting l1-recon parameters');
% L1 Recon Parameters 
    TVWeight = 0.000; 	% Weight for TV penalty
    xfmWeight = 0.1;	% Weight for Transform L1 penalty
    Itnlim = 20;         % Number of iterations
    %generate transform operator
    XFM = Wavelet('Daubechies',4,4);	% Wavelet
    
    % initialize Parameters for reconstruction
    param = init;
    param.XFM = XFM;
    param.TV = TVOP;
    param.TVWeight =0;     % TV penalty 
    param.xfmWeight = xfmWeight;  % L1 wavelet penalty
    param.Itnlim = Itnlim;
    param.lineSearchItnlim=50
    param.data=recondata;
    param.lineSearchItnlim=50;    
    param.V=(MNSA.*mask);
    param.Debug=0;
    param.xfmWeight=xfmWeight*(mean(param.V(mask~=0)));
toc;
%%
param.xfmWeight=100
param.lineSearchAlpha=1e-3
param.sensemaps=sensemaps;
param.data=recondata;

image=bart('fftmod 7',bart('fft 7',bart('fftmod 7',recondata)));
imagec=bart('fmac -C -s8',image, sensemaps);
res=bart('cdf97 6',imagec);
x0=imagec;
res=MultiCoil_CG(zeros(size(res)),param);
figure;subplot(121); imshow(abs(squeeze(imagec)),[])
subplot(122); imshow(abs(squeeze(bart('cdf97 -i 6',res))),[])

% 
%%
params.lineSearchAlpha=1e-5
param.data=double(squeeze(recondata(1,:,:,:)));
N=[256 256];
FT = MCp2DFT(mask, N, squeeze(conj(sensemaps)), 1, 2);
param.V=repmat((MNSA.*mask),[1 1 nc]);
param.FT = FT;

res=XFM*(FT'*param.data);
for n=1:6
    res = fnlCg_test(res,param);
end
recon = XFM'*res;
recon=recon./max(recon(:));
figure(789); imshow(abs(recon),[])

%% test MCFT
N=[256 256]
FT = MCp2DFT(mask, N,squeeze(sensemaps), 1, 2);
FTr=FT'*squeeze(recondata);

nch=1
figure(2); imshow(abs(FTr(:,:,nch)),[])
figure(3); imshow(squeeze(recondata(1,:,:,nch)),[])

%%
FT = MCp2DFT(mask, N,squeeze(sensemaps), 1, 2);
FTr=FT'*squeeze(recondata);

FTrr=FT*(FTr);
nch=2
figure(4); imshow(abs(squeeze(FTrr(:,:,nch))),[0 1e-12])

FTrrr=FT'*(FTrr);
figure(5); imshow(abs(squeeze(FTrrr(:,:))),[])


