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


%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

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


%% %%%%%%%%%%%
% ADD VARIABLE NOISE TO THE DATA
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


