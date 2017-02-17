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

%% FASTA 


A=@(x) E4*x;
AT=@(x) E4'*x;
mu=1

[sol, outs_adapt] = fasta_sparseLeastSquares(A,AT,Ku,mu,ones(size(linear_recon_s)), opts);

figure(1); imshow(abs(mat(sol)),[])