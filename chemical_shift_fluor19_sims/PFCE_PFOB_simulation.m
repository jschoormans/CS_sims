% phantom experiemnt Fluor imaging/unfolding 
clear; close all
 if ispc()
     addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'))
 else
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/chemical_shift_fluor19_sims')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'));
 end
%% parameters
Nx=64;  %number of voxels
Gx=0.0004 %T/m gradient strength 
Lx=0.1; %FOV size 

rr = @(I) reshape(I,[Nx,Nx])

%% create phantom 

P1 = phantom([100,0.4,0.2,0.5,0.25,60], [64]);
P2 = phantom([80,0.2,0.4,-0.5,0.5,10], [64]);
% P3 = phantom([30,0.25,0.25,-0.5,-0.5,50], [64]);
I=P1+P2%+P3;

figure(1); subplot(311); imshow(I,[]);axis off;


%%
[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra(Gx,Nx,Lx)

% HOW TO DEAL WITH THIS DISCRETIZATION???
Spectrum=zeros(1,15);

Spectrum(round(PFOB)+abs(min(round(PFOB)))+1)=PFOB_alpha./sum(PFOB_alpha(:))
figure(1); subplot(312); stem(Spectrum); 

offset=-(min(round(PFOB)));
A=opConvolve(Nx,Nx,Spectrum,[1 offset],'cyclic') %cyclic/ truncated?
Ic=A*I(:);
figure(1); subplot(313);  imshow(rr(Ic),[]); axis off;

%% 
figure(2); subplot(221); imshow(rr(I),[]); axis off; title('phantom - no noise')

%% inverse; kspace

F=opDFT2(64,64,1);
K=F*A*I(:);

sigma2=4 %4/12
K=K+randn(size(K))*sigma2; %add noise to K

A2=F*A; 
A2I=opInverse(A2)

Klin=opInverse(A2)*K;

figure(2); subplot(222); imshow(rr(Klin),[]); axis off; title('deconvolution')


 figure(2); subplot(223); imshow(rr(abs(opInverse(F)*K)),[]); axis off; title('linear recon of noisy kspace')

%% recon in k-space (do spectrum shift in k-space, no need for FFT all the time??) 

CGK=nl_conjgrad_fluor(A2,K,zeros(size(K)),50,I(:),1e0,64,64); 
figure(2); subplot(224); imshow(rr(abs(CGK)),[]); axis off; title('CG-k space')
 if ispc()
 cd('L:\basic\divi\Projects\cosart\CS_simulations\chemical_shift_fluor19_sims')    
 else
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/chemical_shift_fluor19_sims/figs')
 end
%%
schaalmin=7
schaalmax=120
close all
% figure(3); imshow(rr(I),[0 schaalmax]); axis off;
% export_fig phantom.pdf -native
% figure(4); imshow(rr(Klin),[]); axis off
% export_fig deconvolution.pdf -native
figure(5); imshow(rr(abs(opInverse(F)*K)),[schaalmin schaalmax]); axis off;
export_fig linear_recon.pdf -native
figure(6); imshow(rr(abs(CGK)),[schaalmin schaalmax]); axis off;
export_fig CG.pdf -native