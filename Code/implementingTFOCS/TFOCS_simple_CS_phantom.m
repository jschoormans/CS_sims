%% SIMPLES PHANTOM
clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'))
addpath /home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/TFOCS-1.4
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'))
addpath(genpath('/opt/amc/bart/')); vars; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/SALSA_v2.0'))
%%
n1=128; n2=n1;

 mat = @(x) reshape(x,n1,n2);


K=phantom([n1]);
K=fftshift(fftshift(fft2(K),1),2);
mask=rand(n1,n2)>0.7;
K(~mask)=0;
%K=phantom([100,100])
Kvec=vec(K);
Ku=vec(K(mask));

mat = @(x) reshape(x,n1,n2); %function that reshapes vector to matrix 


AA=opDFT2(n1,n1,1) %DFT2 MATRIX
l=[1:n1*n2];
mask3=(l(mask(:))); %indices of rows in FDT2 matrix that should be removed (were not sampled!)
Re=opRestriction(n1*n2,mask3);
E4=Re*AA;

linear_recon_s=E4'*Ku; % not useful; output should be undersampled length 

figure(3); imshow(abs(mat(linear_recon_s)),[]); axis off; 


%% Call the TFOCS solver

mu              = 5;
er              = @(x) norm(x(:)-linear_recon_s(:))/norm(linear_recon_s(:)); %NOT SURE ABOUT THIS: WE DONT KNOW SOLUTION YET...
opts = [];
opts.errFcn     = @(f,dual,primal) er(primal);
opts.maxIts     = 200;
opts.printEvery = 20;
opts.tol        = 1e-8;
opts.stopcrit   = 4;
opts.alg='GRA'  

x0 = linear_recon_s;  %first guess
z0  = [];           % we don't have a good guess for the dual (is this true)


a=0;
% build operators:
% A           = linop_handles([N,N], @(x)fft2(x), @(x)ifft2(x),'c2c');
A=linop_spot(E4,'c2c')
normA2      = 1;
W_wavelet   = linop_spot(opWavelet2(n1,n2,'daubechies'),'c2c');
normWavelet      = linop_normest( W_wavelet );

contOpts            = [];
contOpts.maxIts     = 4;


EPS=2e10
[x_wavelets,out_wave] = solver_sBPDN_W( A, W_wavelet, double(Ku), EPS, mu, ...
    x0(:), z0, opts, contOpts);


 %possible reweighting here
 
 figure(3); imshow(abs([mat(x_wavelets./max(x_wavelets(:))),mat(linear_recon_s./max(linear_recon_s(:)))]),[])
 