function [Recon,D,D2,D3,acc,Msave,MSV,Ns,MNSA]=analyse_CSSNR_varNSA(Im,K,sens,NoiseLevels,regs,unders,iters,restrictellipse)
%function to analyse the reconstructions of undersampled and fully sampled
%Im for different NoiseLevels, regularization factors, and undersampling
%ratios. 
%Now only for fixed scan times!!

%TO DO/TO THINK ABOUT:
%How to calculate profiles??
%Sense maps??
%which kind of images?
%which error measures?



Msave=zeros(size(K,1),size(K,2),length(unders),iters);
for iter=1:iters
    M=zeros(size(K,1),size(K,2),length(unders));
    for kk=1:length(unders)
        % M(:,:,kk)=bart(['poisson -Y',num2str(size(K,1)),' -Z',num2str(size(K,2)),' -v -e -C50 -y',num2str(unders(kk)),' -z',num2str(unders(kk))]);
        [pdf,val] = genPDF(size(K),5,unders(kk),2,0,0);
        Mfull=genEllipse(size(K,1),size(K,2));
        MNSA(:,:,kk)=1./pdf;
        MNSA(:,:,kk)=round(MNSA(:,:,kk));
        MNSA2(:,:,kk)=pdf.*(1./unders(kk));
        MNSA2(:,:,kk)=round(MNSA2(:,:,kk));

        
        if  ~exist('restrictellipse')
            M(:,:,kk)=genSampling(pdf,10,100).*Mfull;
        elseif restrictellipse(iter)==1
            M(:,:,kk)=genSampling(pdf,10,1).*Mfull;
        else
            M(:,:,kk)=genSampling(pdf,10,1);
        end
        Msave(:,:,:,iter)=M(:,:,:);
        acc(kk,iter)=sum(sum(M(:,:,kk)==1))/((size(M,1)/2)*(size(M,2)/2)*(pi));
    end
Ns(:,iter)=sum(sum(MNSA.*M,1),2)
    
    
    
for ll=1:length(regs)
reg=regs(ll)
ImRef=squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],permute(K,[3 2 1]),sens)); %to compare image


for kk=1:length(unders)
waitbar(((iter-1)/(iters))+(1/(iters))*((ll-1)/length(regs))+(1/(iters))*((1/length(regs)))*(kk/length(unders)));
for jj=1:length(NoiseLevels);
clear K_N Ku_N



for iii=1:max(MNSA(:,:,kk)) %Matrix of NSA values
K_N=addNoise(K,NoiseLevels(jj));
Ku_N1=squeeze(M(:,:,kk).*(MNSA(:,:,kk)>=iii)).*K_N;
Ku_N2(1,:,:,1,1,iii)=permute(Ku_N1,[3,2,1,4,5]);

Ku_N3=squeeze(M(:,:,kk).*(MNSA2(:,:,kk)>=iii)).*K_N;
Ku_N4(1,:,:,1,1,iii)=permute(Ku_N3,[3,2,1,4,5]);


Ku_Nnonvar=squeeze(M(:,:,kk)).*K_N; %tempporary: just two times averaging for all points
Ku_Nnonvar2(1,:,:,1,1,iii)=permute(Ku_Nnonvar,[3,2,1,4,5]);


K_N=Mfull.*K_N; %set outside ellipse to zeroes
K_N=permute(K_N,[3 2 1]);

end


Ku_Nvar1=sum(Ku_N2,6)./permute(MNSA(:,:,kk),[3 2 1]); %mean of measures....
Ku_NnonvarF=sum(Ku_Nnonvar2(:,:,:,:,:,1:(1./unders)),6)./(1./unders); % %tempporary: just 4 times averaging for all points
Ku_NnoNSA=sum(Ku_Nnonvar2(:,:,:,:,:,1),6); % %tempporary: just 4 times averaging for all points
Ku_Nvar2=sum(Ku_N4,6)./permute(MNSA2(:,:,kk),[3 2 1]); %mean of measures....



MSV(jj,kk)=mean(Ku_Nvar1(Ku_Nvar1~=0)); %mean signal value in us to get an idea how the mean SNR would change!!!!

Recon(:,:,1,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Ku_Nvar1,sens);
Recon(:,:,2,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],K_N,sens);
Recon(:,:,3,jj,kk,ll,iter)=ifftshift(ifft2((squeeze(Ku_Nvar1))));
Recon(:,:,4,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Ku_NnonvarF,sens);
Recon(:,:,5,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Ku_NnoNSA,sens);
Recon(:,:,6,jj,kk,ll,iter)=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Ku_Nvar2,sens);


%HERE: ADD FUNCTIONS THAT DO THE RECON ANALYSIS

D(1,jj,kk,ll,iter)=mean(mean(abs(squeeze(Recon(:,:,1,jj,kk,ll,iter))-ImRef).^2,1),2);
D(2,jj,kk,ll,iter)=mean(mean(abs(squeeze(Recon(:,:,2,jj,kk,ll,iter))-ImRef).^2,1),2);

D2(1,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,1,jj,kk,ll,iter))),abs(ImRef));
D2(2,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,2,jj,kk,ll,iter))),abs(ImRef));
D2(3,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,3,jj,kk,ll,iter))),abs(ImRef));
D2(4,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,4,jj,kk,ll,iter))),abs(ImRef));
D2(5,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,5,jj,kk,ll,iter))),abs(ImRef));
D2(6,jj,kk,ll,iter)=ssim(abs(squeeze(Recon(:,:,6,jj,kk,ll,iter))),abs(ImRef));


I1=abs(squeeze(Recon(50:70,100:200,1,jj,kk,ll,iter)));
I2=abs(squeeze(Recon(50:70,100:200,2,jj,kk,ll,iter)));

D3(1,jj,kk,ll,iter)=ssim(I1,abs(ImRef(50:70,100:200))); %SSIM for a highly detailed region
D3(2,jj,kk,ll,iter)=ssim(I2,abs(ImRef(50:70,100:200))); 

end
end
end
end