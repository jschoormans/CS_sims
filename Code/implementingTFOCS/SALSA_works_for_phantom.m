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

% vec = @(x) x(:);


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

linear_recon_s3=E4'*Ku; % not useful; output should be undersampled length 

figure(3); imshow(abs(mat(linear_recon_s3)),[]); axis off; 

%% 
y=Ku;
A=E4;
AT=A' ;
lambda=9e-2
mu=lambda/50;

inneriters = 1;
outeriters = 1000;
tol = 5e-6;

invLS = @(x) ((opEye(n1*n2) - (1/(1+mu))*E4'*E4)/mu)*vec(x);
Phi_TV = @(x) TVnorm(real(x));



[x_salsa, numA, numAt, objective, distance,  times]= ...
         SALSA_v2(y,A,lambda,...
         'AT', AT,...
         'Phi', Phi_TV, ...
         'LS', invLS, ...         
         'Phi', Phi_TV, ...
         'TVINITIALIZATION', 1, ...
         'StopCriterion', 1,...
         'MAXITERA', outeriters, ...
         'Verbose', 1);
figure(7);imshow(abs(mat(x_salsa)),[])





%%

mask;
A = @(x)  masked_FFT(x,mask);
AT = @(x) (masked_FFT_t(x,mask));
ATA = @(x) (ifft2c(mask.*fft2c(x))) ;

global calls;
calls = 0;
A = @(x) callcounter(A,x);
AT = @(x) callcounter(AT,x);
ATA = @(x) callcounter(ATA,x);

y=Ku;
% A=E4;
% AT=A' ;
%%%% algorithm parameters
lambda = 9e-2;
mu = lambda*1;
inneriters = 1;
outeriters = 1000;
tol = 5e-6;


invLS = @(x) (x - (1/(1+mu))*ATA(x) )/mu;
Phi_TV = @(x) TVnorm(real(x));
Psi_TV = @(x) W*(real(x));



[x_salsa, numA, numAt, objective, distance,  times]= ...
         SALSA_v2(y,A,lambda,...
         'AT', AT, ...
         'Mu', mu, ...
         'Phi', Phi_TV, ...
         'TVINITIALIZATION', 1, ...
         'StopCriterion', 1,...
         'MAXITERA', outeriters, ...
         'LS', invLS, ...
         'Verbose', 1);
figure(8);imshow(abs(mat(x_salsa)),[])






