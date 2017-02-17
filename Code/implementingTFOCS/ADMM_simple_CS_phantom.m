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
mask=rand(n1,n2)>0.3;
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

%% phase operator
P=diag(exp(-1i*angle(linear_recon_s)));
% wavelet ops
W   = opWavelet2(n1,n2,'daubechies'),'c2c';

Wlin=W*P*linear_recon_s;
figure(4); imshow(mat(Wlin),[])


%% ADMM

A=E4*P'*W';
b=Ku; 
lambda_max = norm( A'*b, 'inf' );
lambda = 0.1*lambda_max*1e3
mu=1;
rho=1.0; 
[x history]=lassoADMM(A, b, lambda, mu,rho);
x=P'*W'*x;
figure(3); imshow(abs([mat(x./max(x(:))),mat(linear_recon_s./max(linear_recon_s(:))),mat((x./linear_recon_s)./max(x(:)./linear_recon_s(:)))]),[])
