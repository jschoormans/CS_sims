function [Recon,D]=analyse_NS(Im,K,M,sens,NoiseLevels)


for jj=1:length(NoiseLevels);
clear K_N Ku_N
for i=1:4  %MAKE  four noisy undersampled meas (and one fully sampled noisy measurement)
K_N=addNoise(K,NoiseLevels(jj));
Ku_N(:,:,:,:,i)=squeeze(M).*K_N;
end
K_N=permute(K_N,[3 2 1]);
Ku_N=permute(Ku_N,[3 2 1,4,5]);
Ku_N4=mean(Ku_N,5);

Recon(jj,1,:,:)=bart('pics -RW:6:0:0.05 -e -i100',Ku_N(:,:,:,1),sens);
Recon(jj,2,:,:)=bart('pics -RW:6:0:0.05 -e -i100',Ku_N4,sens);
Recon(jj,3,:,:)=bart('pics -RW:6:0:0.05 -e -i100',K_N,sens);

D(1,jj,:)=sum(sum(abs(squeeze(Recon(jj,1,:,:))-Im'),1),2);
D(2,jj,:)=sum(sum(abs(squeeze(Recon(jj,2,:,:))-Im'),1),2);
D(3,jj,:)=sum(sum(abs(squeeze(Recon(jj,3,:,:))-Im'),1),2);

end