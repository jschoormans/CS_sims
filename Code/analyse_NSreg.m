function [Recon,D,D2,D3]=analyse_NSreg(Im,K,M,sens,NoiseLevels,regs)

for kk=1:length(regs)
for jj=1:length(NoiseLevels);
clear K_N Ku_N
for i=1:4  %MAKE  four noisy undersampled meas (and one fully sampled noisy measurement)
K_N=addNoise(K,NoiseLevels(jj));
Ku_N(:,:,:,:,i)=squeeze(M).*K_N;
end
K_N=permute(K_N,[3 2 1]);
Ku_N=permute(Ku_N,[3 2 1,4,5]);
Ku_N4=mean(Ku_N,5);


reg=regs(kk)
Recon(jj,1,:,:,kk)=bart(['pics -RW:6:0:',num2str(reg),' -e -i100'],Ku_N(:,:,:,1),sens);
Recon(jj,2,:,:,kk)=bart(['pics -RW:6:0:',num2str(reg),' -e -i100'],Ku_N4,sens);
Recon(jj,3,:,:,kk)=bart(['pics -RW:6:0:',num2str(reg),' -e -i100'],K_N,sens);

D(1,jj,kk)=mean(mean(abs(squeeze(Recon(jj,1,:,:,kk))-Im').^2,1),2);
D(2,jj,kk)=mean(mean(abs(squeeze(Recon(jj,2,:,:,kk))-Im').^2,1),2);
D(3,jj,kk)=mean(mean(abs(squeeze(Recon(jj,3,:,:,kk))-Im').^2,1),2);

D2(1,jj,kk)=ssim(abs(squeeze(Recon(jj,1,:,:,kk))),abs(Im)');
D2(2,jj,kk)=ssim(abs(squeeze(Recon(jj,2,:,:,kk))),abs(Im)');
D2(3,jj,kk)=ssim(abs(squeeze(Recon(jj,3,:,:,kk))),abs(Im)');

I1=abs(squeeze(Recon(jj,1,150:end,:,kk)));
I2=abs(squeeze(Recon(jj,2,150:end,:,kk)));
I3=abs(squeeze(Recon(jj,3,150:end,:,kk)));

D3(1,jj,kk)=std(I1(:))
D3(2,jj,kk)=std(I2(:))
D3(3,jj,kk)=std(I3(:))

end
end