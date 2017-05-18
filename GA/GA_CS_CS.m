function [E,CS,MNSA]=GA_CS_CS(Si,K,ImRef,sens,Noise,XFM)

% calculate k-space for individual 
   MNSA=zeros(size(K)) ;
for ii=1:size(Si,3)
    MNSA(Si(1,1,ii),Si(1,2,ii))=MNSA(Si(1,1,ii),Si(1,2,ii))+1;
end

Ki=K+Noise.*(1./sqrt(MNSA));
Ki(MNSA==0)=0; %remove zeroes; 

% do CS recon 


%%
N=size(K);

if false
CS_GA_initCG;
param.data=Ki;
FT = p2DFT(MNSA~=0, N, 1, 2);
param.FT = FT;
im_dc2 = FT'*(param.data); %linear recon; scale data to prevent low-pass filtering
res = XFM*im_dc2;
param.V=(MNSA);
param.Debug=0;
param.xfmWeight=param.xfmWeight*(mean(param.V(MNSA~=0)));

for n=1:param.niter
    res = fnlCg_test(res,param);
    CS = XFM'*res;
end
CS=fftshift(abs(CS),2);
else
    CS=bart('pics -l1 -r0.02 -i50',Ki,ones(size(Ki)));
    CS=fftshift(abs(CS),2);
end

CS=CS./max(CS(:));

%% evaluate
% E=psnr(abs(CS./max(CS(:))),abs(ImRef./max(ImRef(:))));
% E=ssim(abs(CS./max(CS(:))),abs(ImRef./max(ImRef(:))));
E=1/(1e-8+immse(abs(CS./max(CS(:))),abs(ImRef./max(ImRef(:)))))
end