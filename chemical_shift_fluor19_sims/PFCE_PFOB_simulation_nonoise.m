% phantom experiemnt Fluor imaging/unfolding 
clear; close all
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))

%% parameters
Nx=64;  %number of voxels
Gx=0.0005 %T/m gradient strength 
Lx=0.1; %FOV size 

rr = @(I) reshape(I,[Nx,Nx])

%% create phantom 

P1 = phantom([100,0.4,0.2,0.5,0.25,60], [64]);
P2 = phantom([80,0.2,0.4,-0.5,0.5,10], [64]);
P3 = phantom([20,0.25,0.25,-0.5,-0.5,50], [64]);
I=P1+P2+P3;


figure(1); subplot(311); imshow(I,[]);axis off;


%%
[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra(Gx,Nx,Lx)

% HOW TO DEAL WITH THIS DISCRETIZATION???
Spectrum=zeros(1,15);

Spectrum(round(PFOB)+abs(min(round(PFOB)))+1)=PFOB_alpha./sum(PFOB_alpha(:))
figure(1); subplot(312); stem(Spectrum); 

offset=-(min(round(PFOB)));
A=opConvolve(Nx,Nx,Spectrum,[0 offset],'cyclic') %cyclic/ truncated?
Ic=A*I(:);
figure(1); subplot(313);  imshow(rr(Ic),[]); axis off;

%% inverse im space
Ir=opInverse(A)*Ic(:);
figure(2); subplot(221); imshow(rr(Ir),[]); axis off; title('pseudo-inverse')

%% inverse; kspace

F=opDFT2(64,64,1);
K=F*A*I(:);
A2=F*A; 
A2I=opInverse(A2)

Klin=opInverse(A2)*K;
figure(2); subplot(222); imshow(rr(Klin),[]); axis off; title('inverse (linear recon) -kspace')


%% recon
CGI=nl_conjgrad_fluor(A,Ic,zeros(size(Ir)),200,I(:),1e-1,64,64); 
figure(2); subplot(223); imshow(rr(abs(CGI)),[]); axis off; title('CG-image space')

%% recon in k-space (do spectrum shift in k-space, no need for FFT all the time??) 
% 
% CGK=nl_conjgrad_fluor(A2,K,zeros(size(K)),200,I(:),0,64,64); 
% figure(2); subplot(224); imshow(rr(abs(CGK)),[]); axis off; title('CG-k space')
% 
