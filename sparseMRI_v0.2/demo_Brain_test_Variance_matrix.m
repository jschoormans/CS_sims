% This script demonstare a CS reconstruction from 
% Randomly undersmapled phase encodes of 2D FSE
% of a brain image.


addpath(strcat(pwd,'/utils'));

if exist('FWT2_PO') <2
	error('must have Wavelab installed and in the path');
end

load brain512

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(data); 	% image Size
DN = size(data); 	% data Size
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.005;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

acc=0.3;
[pdf,val] = genPDF(N,4,acc,1,0.1,0);
[mask,stat,actpctg] = genSampling(pdf,10,10);

%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);
data=FT*phantom(N(1)); %replace data with undersampled phantom

% scale data
im_dc = FT'*(data.*mask./pdf);
%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet


% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =0;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


%% %%%%%%%%%%% EXP 1
%{
% ADD NOISE TO THE DATA
noiselevel=2.5e-2
Noise=(randn(N).*noiselevel+randn(N)*noiselevel*1j);
NoiseM=Noise.*mask; %add noise to sample points
param.data=data+Noise;

param.data = param.data/max(abs(im_dc(:)));
im_dc2=FT'*(param.data);
im_dc2 = im_dc2/max(abs(im_dc2(:)));


figure(100), imshow(abs(im_dc2),[]);drawnow;

res = XFM*im_dc2;

% do iterations
param.V=ones(N) %variance matrix
param.V(N/4:3*N/4,N/4:3*N/4)=1e-1

for n=1:2
	res = fnlCg_test(res,param);
	im_res = XFM'*res;
	figure(100), imshow(abs(im_res),[]), drawnow
end



figure(101), imshow(abs(cat(2,im_dc2,im_res,(im_res-im_dc))),[]);
figure(102), imshow(abs(cat(2,im_dc2(155:304,110:209), im_res(155:304,110:209))),[0,1],'InitialMagnification',200);
title(' zf-w/dc              l_1 Wavelet');

figure(103); imshow(log(mask.*abs((FT*im_res)-param.data)),[])
%}
%%
%% %%%%%%%%%%% EXP 2
% ADD VARIABLE DENSITY NOISE TO THE DATA
close all
NSA=ceil(10.*pdf).^2;
%NSA=ones(size(pdf))
NoiseStd=1./sqrt(NSA);
noiselevel=50e-2.*NoiseStd;
Noise=(randn(N).*noiselevel+randn(N).*noiselevel.*1j);
NoiseM=Noise.*mask; %add noise to sample points
param.data=data+NoiseM;

param.data = param.data/max(abs(im_dc(:)));
im_dc2=FT'*(param.data);
im_dc2 = im_dc2/max(abs(im_dc2(:)));

figure(99); imshow(abs(NoiseStd))
figure(100), imshow(abs(im_dc2),[]);drawnow;

res = XFM*im_dc2;
res_orig = XFM*im_dc2;

% do iterations
param.V=1./(NSA);
% param.V=ones(N)
for n=1:2
	res = fnlCg_test(res,param);
	im_res = XFM'*res;
    
    res_orig = fnlCg(res_orig,param);
	im_res_orig = XFM'*res_orig;
	figure(100), imshow(cat(2,abs(im_dc2),abs(im_res),abs(im_res_orig)),[]), drawnow
end
im_res_pics=bart('pics -RW:7:0:0.02 -S',permute(param.data,[3 1 2]),ones([512 512 1]));



figure(101), imshow(abs(cat(2,im_dc2,im_res,(im_res-im_dc),im_res_orig,(im_res_orig-im_dc),im_res_pics,im_res_pics-im_dc)),[]);
figure(102), imshow(abs(cat(2,im_dc2(155:304,110:209), im_res(155:304,110:209),im_res_orig(155:304,110:209),im_res_pics(155:304,110:209))),[0,1],'InitialMagnification',200);
title(' zf-w/dc              l_1 Wavelet VAR MATRIX     l1-wavelet no VAR MATRIX    BART');

figure(103); subplot(2,2,1); imshow(log(mask.*abs((FT*im_res)-param.data)),[]); subplot(2,2,2);imshow(NSA,[])

subplot(2,2,3); 
imshow(fft2c(abs((FT*im_res)-param.data)),[]);

subplot(2,2,4); 
imshow(fft2c(abs((FT*im_res_orig)-param.data)),[]);