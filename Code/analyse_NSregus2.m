function [Recon,D,D2,D3,acc]=analyse_NSregus2(Im,K,sens,NoiseLevels,reg,unders)

for kk=1:length(unders)
[pdf,val] = genPDF(size(K),8,unders(kk),2,0.05,0);
% M(:,:,kk)=bart(['poisson -Y',num2str(size(K,1)),' -Z',num2str(size(K,2)),' -v -e -C50 -y',num2str(unders(kk)),' -z',num2str(unders(kk))]);
M(:,:,kk)=genSampling(pdf,10,1);
acc(kk)=sum(sum(M(:,:,kk)==1))/((size(M,1)/2)*(size(M,2)/2)*(pi))
figure;
imshow(M(:,:,kk))
end

for kk=1:length(unders)
for jj=1:length(NoiseLevels);
clear K_N Ku_N
for i=1:4  %MAKE  four noisy undersampled meas (and one fully sampled noisy measurement)
K_N=addNoise(K,NoiseLevels(jj));
Ku_N(:,:,:,:,i)=squeeze(M(:,:,kk)).*K_N;
end
K_N=permute(K_N,[3 2 1]);
Ku_N=permute(Ku_N,[3 2 1,4,5]);
% Ku_N4=mean(Ku_N,5);

Recon(jj,1,:,:,kk)=bart(['pics -RW:6:0:',num2str(reg),' -e -i100'],Ku_N(:,:,:,1),sens);
% Recon(jj,2,:,:,kk)=bart(['pics -RW:6:0:',num2str(reg),' -e -i100'],Ku_N4,sens);
Recon(jj,2,:,:,kk)=bart(['pics -RW:6:0:',num2str(reg),' -e -i100'],K_N,sens);

D(1,jj,kk)=mean(mean(abs(squeeze(Recon(jj,1,:,:,kk))-Im').^2,1),2);
D(2,jj,kk)=mean(mean(abs(squeeze(Recon(jj,2,:,:,kk))-Im').^2,1),2);

D2(1,jj,kk)=ssim(abs(squeeze(Recon(jj,1,:,:,kk))),abs(Im)');
D2(2,jj,kk)=ssim(abs(squeeze(Recon(jj,2,:,:,kk))),abs(Im)');

I1=abs(squeeze(Recon(jj,1,150:end,:,kk)));
I2=abs(squeeze(Recon(jj,2,150:end,:,kk)));

D3(1,jj,kk)=std(I1(:))
D3(2,jj,kk)=std(I2(:))

end
end