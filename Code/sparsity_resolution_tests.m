clear all; close all; clc; 
addpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.506'));
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab'));
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab2'));
addpath(genpath('/opt/amc/bart')); vars; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/')

%%
I=imread('brain.png');
I=double(I(:,:,1))./255;

imshow(abs(I),[])
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab Collection/CS Excercise'))
W=Wavelet;



%% zerofill to 1024x1024

I1024=padarray(I,[268,267]);
I1024=I1024(1:1024,1:1024);



%%

res=[1024,512,256,128,64,32,16];
for ii=1:length(res)
D{ii}.I=I1024(1:2^(ii-1):end,1:2^(ii-1):end)
D{ii}.K=ifftshift(ifftshift(ifft2(D{ii}.I),2),1);
D{ii}.W=W*D{ii}.I;
end
% I512=I1024(1:2:end,1:2:end);
% I256=I1024(1:4:end,1:4:end);
% I128=I1024(1:8:end,1:8:end);
% I64=I1024(1:16:end,1:16:end);
% %%
% W1024=W*I1024;
% W512=W*I512;
% W256=W*I256;
% W128=W*I128;
% W64=W*I64;
% %%
% K1024=ifft2(I1024);
% K512=ifft(I512);








%% plot imags and wavelet
figure(1);
subplot(2,2,1)
imshow(abs(I1024),[]);
axis off;
subplot(2,2,2)
imshow(abs(W1024),[])
axis off;


%% are high reeees images mpre sparsse?
figure(2)
for ii=1:length(res)
l0(ii)=sum(D{ii}.W(:)~=0)/res(ii)^2
end
plot(res,l0)

%% how do the centers of k-space compare to the low res vers?

ii=1
K1=D{ii}.K;
K1=K1(res(ii+2)+1:1:end-res(ii+2),res(ii+2)+1:1:end-res(ii+2));
K2=D{ii+1}.K;

ssim(real(K1),real(K2))
ssim(imag(K1),imag(K2))

% so if the higher resolution k-space is the entire previous k-space plus
% some things in the borders, the detail sparsity (marginal sparsity??)
%is quite dramatic...
%% lets try to calculate the marginal sparsity....
%i.e. the sparsity of the k-space you add by going to a higher
%res....


C=0.1;
for ii=1:length(res)-1
l0_marg(ii)=((l0(ii)*res(ii)^2)-(l0(ii+1)*res(ii+1)^2))/((3/4)*res(ii)^2)
end

for ii=1:length(res) %theoretical unders. bound 
   Slogn(ii)=C*log(res(ii).^2)*l0(ii) 
end

figure(3); hold on; plot(res,l0,'k+-');
plot(res,Slogn,'g')
plot(res(1:length(res)-1),l0_marg,'r+-'); hold off





