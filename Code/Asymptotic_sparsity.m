clear all; close all; clc; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
%%
% I=imread('brain.png');
I=imread('earth.jpg');
I=double(I(:,:,1))./255;

figure(1)
imshow(abs(I),[])
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab Collection/CS Excercise'))
W=Wavelet;
%%
% I1024=padarray(I,[218,217]);
% I1024=I1024(1:1024,:);
% WI=W*I1024;

I2048=I(1:2048,1:2048);
WI=W*I2048;

figure(2); imshow(abs(WI),[0 1])
%% SORT AND ORDER DATA
nmax=11
for  ii=2:nmax
    if ii==2
        F{ii}.W=reshape((WI(1:2^ii,1:2^ii)),1,[]);
        F{ii}.S=sort(abs(F{ii}.W),'descend'); % all l1 values of this level sorted
        F{ii}.C=abs(cumsum(F{ii}.S))./abs(sum(F{ii}.S)); %should be sum from largest to i
    else
        F{ii}.W=reshape([(WI(2^(ii-1)+1:2^ii,1:2^ii)),(WI(1:2^(ii-1),1+2^(ii-1):2^ii))],1,[]);
        F{ii}.S=sort(abs(F{ii}.W),'descend'); % all l1 values of this level sorted
        F{ii}.C=abs(cumsum(F{ii}.S))./abs(sum(F{ii}.S)); %should be sum from largest to i
    end
end
%% PLOT SPARSITY DATA

colormap('jet')
figure(3)
hold on
for ii=2:10
   x=linspace(0,1,length(F{ii}.C));
   plot(x,F{ii}.C) 
end
xlabel('fraction of ordered coeffs (low-high)')
ylabel('l1-norm of smallest x coeffs (normalized to l1 of all coeff)')

hold off
legend('level 1','2','3','4','5','6','7','8','9')
title('Asymptotic sparsity')

%% generate PDF from this

% 
% frac=0.1
% F{ii}.C(round(frac*length(F{ii}.C)))

% PDF=linspace(1,1,2^nmax);
% epsilon=0.1;
% 
% for ii=2:nmax
%     sparsity(ii)=1-(find(F{ii}.C>epsilon,1,'first')./length(F{ii}.C)) ;
%     s(ii)=sparsity(ii)*length(F{ii}.C);
%     % PDF(2^(ii-1):2^ii)=sparsity(ii);
% end


%{
%% make matrix pdf
pdf=ones(2^nmax,2^nmax);
for ii=1:(2^nmax)
    for jj=1:(2^nmax)
   pdf(ii,jj)=PDF(1+round(sqrt((ii-(2^nmax)/2)^2+(jj-(2^nmax)/2)^2)));
    end
end

figure(8); subplot(121);plot(PDF);
subplot(122);imshow(pdf,[]);

S=genSampling(pdf,10,1);

%%
pdf2=genPDF(1024,4,0.34,2,0.1)
S2=genSampling(pdf2,10,1);


%% TRY CS WITH NEW SAMPLING PATTERN; 


K=ifftshift(ifftshift(ifft2(I1024),2),1); 
Ku=S.*K; 

R=bart('pics -RW:3:0:0.01',Ku,ones(1024,1024));
R=ifftshift(ifftshift(R,2),1);

figure(9); imshow(abs(R),[])
%}