clear all; close all; clc;

folder='/home/jschoormans/lood_storage/divi/Temp/jasper/20151104_3DimSDE_ISMRM/v2/'
file='20_22102015_1512550_9_2_wip3d05mm35mmacc4nsa4senseV4.raw'
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
MR=CS_Recon([folder,file])


for ii=0:3
MR.Parameter.Parameter2Read.aver=ii
MR.Perform1;
MR.Sort;
A{ii+1}=MR.Data;
end

M=cat(5,A{1},A{2},A{3},A{4}); %5D measuremtn matrix
for nc=1:8 %separately for coils
S(:,:,:,nc)=mean(M(:,:,:,nc,:),5);
V(:,:,:,nc,:)=std(real(M(:,:,:,nc,:)),[],5);
v=squeeze(V(:,:,:,nc,:));
N(nc)=(mean((v(v(:)~=0))))
end

figure(1);
title('SNR per coil for one slice in kspace')
sl=40;
for ii=1:8
    subplot(4,2,ii)
    imshow(S(:,:,sl,ii)/N(ii),[])
end

%%
for ii=1:4
E=M(:,:,:,:,ii)-S; end
%%
figure(2)

hold on 
x=[1:1:500]
[f,x]=hist(real(E(E~=0)),1000) %histogram of all estimated k-space errors
norm=normpdf(x,0,mean(N))
bar(x,f/trapz(x,f))
plot(x,norm);
hold off
%%

nc=1
E_k=(1)/(2*pi).*sum(sum(sum(abs(M(:,:,:,nc,1)).^2,1),2),3)
E_k_mean=mean(mean(mean(abs(M(:,:,:,nc,1)).^2)))

I=ifftn(M(:,:,:,nc,1));
E_I=sum(sum(sum(abs(I(:,:,:)).^2,1),2),3)

SNR_estimator=mean(N)/E_k_mean 
%%
%2D DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder='/home/jschoormans/lood_storage/divi/Temp/jasper/20151104_3DimSDE_ISMRM/v2/'
file='20_22102015_1512550_9_2_wip3d05mm35mmacc4nsa4senseV4.raw'
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
MR=CS_Recon([folder,file])


for ii=0:3
MR.Parameter.Parameter2Read.aver=ii
MR.Perform1;
MR.Sort;
MR.K2IM
A{ii+1}=MR.Data;
end

M=cat(5,A{1},A{2},A{3},A{4}); %5D measuremtn matrix
for nc=1:8 %separately for coils
S(:,:,:,nc)=mean(M(:,:,:,nc,:),5);
V(:,:,:,nc,:)=std(real(M(:,:,:,nc,:)),[],5);
v=squeeze(V(:,:,:,nc,:));
N(nc)=(mean((v(v(:)~=0))))
end

figure(1);
title('SNR per coil for one slice in kspace')
sl=43;
for ii=1:8
    subplot(4,2,ii)
    imshow(S(:,:,sl,ii)/N(ii),[])
end

%%
for ii=1:4
E=M(:,:,:,:,ii)-S; end
%%
figure(2)
hold on 
x=[1:1:500];
[f,x]=hist(real(E(E~=0)),1000); %histogram of all estimated k-space errors
norm=normpdf(x,0,mean(N));
bar(x,f/trapz(x,f));
plot(x,norm);
hold off
%%

nc=1
E_k=(1)/(2*pi).*sum(sum(sum(abs(M(:,:,:,nc,1)).^2,1),2),3)
E_k_mean=mean(mean(mean(abs(M(280,:,:,nc,1)).^2)))

I=ifft(ifft(M(:,:,:,nc,1),[],2),[],3);
E_I=sum(sum(sum(abs(I(:,:,:)).^2,1),2),3)

SNR_estimator=mean(N)/E_k_mean ;

%% for one slice:
E_k=(1)/(2*pi).*sum(sum(sum(abs(M(:,:,30,nc,1)).^2,1),2),3)
E_I=sum(sum(sum(abs(I(:,:,30)).^2,1),2),3)


SNR(:,:,:)=S(:,:,:,ii)/N(ii);
[Dc, Ac]=MeanValuefromCenter(SNR(:,:,43));


Mu=squeeze(M(280,:,:,nc,1)); %kspace of one 2D image right...
figure;
hist(abs(Mu(abs(Mu)~=0)))
figure; hist(abs(K(abs(K)~=0)),1000) % WHERE K IS ..... FROM SIMULATIONS

%%
figure(3)
plot(Dc,abs(Ac),'.')


mask=bart('poisson -Y576 -Z288 -y2 -z2 -v -C30');
SNRmask=SNR.*squeeze(mask);

meanSNR=mean(abs(SNR(SNR~=0)))

meanSNRmask=mean(abs(SNRmask(SNRmask~=0)))