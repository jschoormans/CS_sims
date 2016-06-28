function [Recon,D,D2,D3,acc,Msave,MSV]=analyse_CSSNRf(Im,K,sens,NoiseLevels,regs,iters,restrictellipse)
%function to analyse the reconstructions of undersampled and fully sampled
%Im for different NoiseLevels, regularization factors, and undersampling
%ratios. 
%Now only for fixed scan times!!

%TO DO/TO THINK ABOUT:
%How to calculate profiles??
%Sense maps??
%which kind of images?
%which error measures?

unders=[1:6]

Msave=zeros(size(K,1),size(K,2),length(unders),iters);
for iter=1:iters
    M=zeros(size(K,1),size(K,2),length(unders));
    for kk=1:length(unders)
        % M(:,:,kk)=bart(['poisson -Y',num2str(size(K,1)),' -Z',num2str(size(K,2)),' -v -e -C50 -y',num2str(unders(kk)),' -z',num2str(unders(kk))]);
        [pdf,val] = genPDF(size(K),5,1/unders(kk),2,0.05,0);
        Mfull=genEllipse(size(K,1),size(K,2));
        
        
        if  ~exist('restrictellipse')
            M(:,:,kk)=genSampling(pdf,10,1).*Mfull;
        elseif restrictellipse(iter)==1
            M(:,:,kk)=genSampling(pdf,10,1).*Mfull;
        else
            M(:,:,kk)=genSampling(pdf,10,1);
        end
        Msave(:,:,:,iter)=M(:,:,:);
        acc(kk,iter)=sum(sum(M(:,:,kk)==1))/((size(M,1)/2)*(size(M,2)/2)*(pi));
    end

for ll=1:length(regs)
reg=regs(ll)
for kk=1:length(unders)
waitbar(((iter-1)/(iters))+(1/(iters))*((ll-1)/length(regs))+(1/(iters))*((1/length(regs)))*(kk/length(unders)));
for jj=1:length(NoiseLevels);
clear K_N Ku_N


for iii=1:unders(kk)
K_N=addNoise(K,NoiseLevels(jj));
Ku_N1=squeeze(M(:,:,kk)).*K_N;
K_N=Mfull.*K_N; %set outside ellipse to zeroes
K_N=permute(K_N,[3 2 1]);
Ku_N2(:,:,:,:,:,iii)=permute(Ku_N1,[3 2 1,4,5]);
end
Ku_N=mean(Ku_N2,6);

MSV(jj,kk)=mean(Ku_N(Ku_N~=0)) %mean signal value in us to get an idea how the mean SNR would change!!!!

Recon(:,:,1,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -e -i100'],Ku_N,sens);
Recon(:,:,2,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -e -i100'],K_N,sens);
Recon(:,:,3,jj,kk,ll,iter)=ifftshift(ifft2((squeeze(Ku_N))));


%HERE: ADD FUNCTIONS THAT DO THE RECON ANALYSIS

D(1,jj,kk,ll,iter)=mean(mean(abs(squeeze(Recon(:,:,1,jj,kk,ll,iter))-Im').^2,1),2);
D(2,jj,kk,ll,iter)=mean(mean(abs(squeeze(Recon(:,:,2,jj,kk,ll,iter))-Im').^2,1),2);

D2(1,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,1,jj,kk,ll,iter))),abs(Im)');
D2(2,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,2,jj,kk,ll,iter))),abs(Im)');
D2(3,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,3,jj,kk,ll,iter))),abs(Im)');

I1=abs(squeeze(Recon(50:70,100:200,1,jj,kk,ll,iter)));
I2=abs(squeeze(Recon(50:70,100:200,2,jj,kk,ll,iter)));

D3(1,jj,kk,ll,iter)=ssim(I1,Im(50:70,100:200)); %SSIM for a highly detailed region
D3(2,jj,kk,ll,iter)=ssim(I2,Im(50:70,100:200)); 

end
end
end
end