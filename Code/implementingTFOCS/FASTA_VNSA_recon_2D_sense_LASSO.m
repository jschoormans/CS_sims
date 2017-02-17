
clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'))
addpath /home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/TFOCS-1.4
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'))
addpath(genpath('/opt/amc/bart/')); vars; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/solvers'))

%%
cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA12')
disp('loading k-space...')
load('K_scan13_nocc')
K_orig=K; 

%% 1: simple 2D example; iFFT in z-direction, one NSA

disp('setting up measurement vector')
sl=1;                                           %slice to reconstruct

Ks=squeeze(K(round(size(K,1)/2),:,:,1,:));      %k-space for one channel and one slice
fullmask=Ks~=0;                                 %find mask used for scan (nx*ny*NSA)
MNSA=sum(fullmask,3);                           %find NSA for all k-points in mask 

data=sum(K_orig,5);                             %sum of data over NSA
[nx,n1,n2,ncoils]=size(data);                   %get parameters

K=data./permute(...
    repmat(MNSA,[1 1 nx ncoils]),[3 1 2 4]);    % TAKE MEAN OF MEASURED VALUES
K(isnan(K))=0;                                  % set nans to zero

K=ifft(K,[],1);                                 % FOR NOW: iFFT along read dim
K=squeeze(K(sl,:,:,:));                         % select subset of k-space

mask=K(:,:,1)~=0;                               % find mask that was used for measurements
Ku=vec(K(repmat(mask,[1 1 ncoils])));                                % vectorize measurement data (excluding non-meas values)

figure(1); imshow(abs([K(:,:,1)./max(K(:)),mask]),[0 1]); axis off; title('kspace and mask')

%% make operators 
disp('setting up undersampled Fourier operator')

mat = @(x) reshape(x,n1,n2,ncoils);             % function that reshapes vector to matrix 
mat2d = @(x) reshape(x,n1,n2*ncoils);           % functions that reshapes in 2d matrix for visualization purposes 
matcc =@(x) reshape(x,n1,n2);                   % function that reshapes vector to matrix 

F=opDFT2(n1,n2,1);                              % DFT2 MATRIX

l=[1:n1*n2];
mask2=(l(~mask(:)));                            % indices of rows in FDT2 matrix that should be removed (were not sampled!)
R=opExcise(F,mask2,'rows');                     % remove non-sampled rows from Fourier matrix 
E1=opBlockDiag(R,R,R,R,R,R,R,R,R,R,R,R,R);      % make block-diagonal copies, one for each channel 

%% perform linear recon, sensitivity maps 
disp('linear reconstruction')
linear_recon=E1'*Ku;                            % linear recon (should change to least-squares??)

disp('estimating sense maps and setting up S operator')
sensmaps=bart('ecalib -r20 -m1 -k5',permute(K,[4 1 2 3]));  %% make sense maps with espirit method 
sensmaps=fftshift(fftshift(sensmaps,3),2);

figure(2)
imshow(abs([...
    mat2d(linear_recon./max(linear_recon(:)));...
    mat2d(sensmaps./max(sensmaps(:)))]),[]);  axis off;  
title('linear recon for all channnels')

for i=1:ncoils; 
    S{i}=opDiag(conj(sensmaps(1,:,:,i)));       % make individual coil multiplication operators (conj because it is actually the inverse op!)
end; 
S2= horzcat(S{i},S{2},S{3},S{4},S{5},...
    S{6},S{7},S{8},S{9},S{10},...
    S{11},S{12},S{13});                         % operator to multiply sens maps with image

E2=E1*S2';                                      % Encoding matrix ; DFT and inverse sense maps 

linear_recon_s=E2'*Ku;                          % linear recon including coil combination      
figure(3); 
imshow(abs(matcc(linear_recon_s)),[]); axis off; title('coil-combined linear recon')

%% ADD PHASE OPERATOR 
disp('estimating phase, setting up phase correction operator')

phase_est=angle(linear_recon_s);                % estimate phase image from linear recon
P=opDiag(exp(-1i.*phase_est));                  % define phase correction operator 

E3=E2*P';                                       % add phase corr to encoding matrix
linear_recon_sp=E3'*Ku;                         % phase corrected linear reconstruction 

figure(4); subplot(121);imshow(real(matcc(linear_recon_s)),[]); axis off;
title('linear recon without phase correction')
subplot(122);imshow(real(matcc(linear_recon_sp)),[]); axis off; 
title('linear recon with phase correction')

%% Wavelet operator
disp('setting up wavelet operator')

W=opWavelet2(n1,n2,'daubechies',4,4,0);         % define 2D wavelet operator (unsure about filter and length!)
E4=E3*W';                                       % add inverse wavelet op to encoding matrix

%% FASTA RECONSTRUCTION ALGROTIHM 
disp('FASTA reconstruction')

opts = [];
opts.recordObjective = true;                    % Record the objective function so we can plot it
opts.verbose = true;
opts.stringHeader='    ';                       % Append a tab to all text output from FISTA.  This option makes formatting look a bit nicer. 
% opts.accelerate = true;
opts.tol=1e-3; 
% opts.maxIters=50

A=@(x) E4*x;                                    % convert operator to function form 
AT=@(x) E4'*x;

mu=0.1                                          % control parameter

[sol, outs_adapt] = fasta_sparseLeastSquares(A,AT,Ku,mu,ones(size(W*linear_recon_s)), opts);

figure(5); imshow(abs(matcc(W'*sol)),[])

%% BART COMPARISON (BART PICS DOES l1-ESPIRIT; MORE ADVANCED)

imRecon=bart('pics -RW:7:0:0.05 -d5',permute(K,[4 1 2 3]),fftshift(fftshift(sensmaps,2),3));

imReconBART=abs(squeeze(fftshift(imRecon)./mean(abs(imRecon(:)))));
imReconFASTA=abs(matcc(W'*sol))./mean(abs(W'*sol(:)));

figure(6); imshow([imReconFASTA,imReconBART],[]); title('FASTA, BART')
figure(7); imshow([imReconBART-imReconFASTA],[]); title('BART-FASTA')
