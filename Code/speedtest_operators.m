% TEST WAVELET SPEEDS


n1=400
mask=rand(n1,n1)>0.5;
XFM=Wavelet_SQKSP(size(mask),[1,2]);
A=rand(n1,n1);

niter=100
disp('@Wavelet operator')
tic
for iter=1:niter
B=XFM*A;
end
timeW1=toc

% SPOT OPERATOR

XFM2=opWavelet2(n1,n1,'Daubechies');
disp('SPOT operator')
tic
for iter=1:niter
B=(XFM2*A(:));
end;
timeW2=toc;

fprintf('time @wavelet: %d, time SPOT: %d',timeW1,timeW2)
%conclusion: SPOT OPERATOR ~20 x faster!
%% make spot operator for total variation
close all
A=zeros(n1,n1); A(5,5)=4; A(5,5)=4;  A(4,5)=4; 
TV1=TVOP
TVx=opConvolve(n1,n1,[-1 1].',[0 1],'cyclic');
TVy=opConvolve(n1,n1, [1 -1],[1 2],'cyclic');
TV2=vertcat(TVx,TVy);

B1=TV1*A;
B2=reshape(TV2*A(:),[n1,n1,2]);

fprintf('TVSPOT X: %d TVOP_X: %d',B3(5,5,1),B1(5,5,1))
figure(1); plot(B2(5,1:9,1),'or'); hold on; plot(B1(5,1:9,1),'k-');hold off
figure(2); plot(B2(5,1:10,2),'or'); hold on; plot(B1(5,1:10,2),'k-');hold off

%% compare efficiency
niter=100;

A=rand(n1,n1);

tic
for i=1:niter
TV1*A;
end
timeTV1=toc;


tic
for i=1:niter
TV2*A(:);
end
timeTV2=toc;
fprintf('TVOPERATOR: %d; SPOT OPERATPOR: %d', timeTV1,timeTV2)
 %conclusion: TVOP is faster (just a matrix multiplication instead of
 %convolution)

%% FFT operator 

niter=1000
% FT1=
FT2=opDFT2(n1,n1,1)
A=rand(n1,n1);

tic
for i=1:niter
ifft2c(A);
end
timeFT1=toc;


tic
for i=1:niter
FT2'*A(:);
end
timeFT2=toc;
fprintf(' FT OPERATOR: %d; SPOT OPERATOR: %d', timeFT1,timeFT2)

%CONCLUSION: SPOT OPERATOR A LITTLE BIT SLOWER (10%)
% IFFT: 2 x slower !
% FT OPERATOR: 2.810379e+00; SPOT OPERATOR: 3.219021e+00



