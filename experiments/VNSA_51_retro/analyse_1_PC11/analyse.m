cd('L:\basic\divi\Projects\cosart\CS_simulations\experiments\VNSA_51_retro\analyse_1_PC11')
clear all; close all; clc;
load('ImRef_11.mat')
ImRef=bart('fft 7',ImRef);
ImRef=bart('resize 0 1024 1 1024',ImRef);
size(ImRef)
ImRef=bart('fft -i 7',ImRef);
ImRef=abs(ImRef)./max(abs(ImRef(:)));

accvector=[1,2,3,4,5,6];
nNSA=5;

MSE=NaN(6,3,3)
SSIM=NaN(6,3,3)

for kk=1:3 %NSA
    for jj=1:3  %vartypes
        for ii=1:length(accvector);
            ii
            load(['R_',num2str(ii),'_',num2str(jj),'_',num2str(kk)],'R')
            Recon{ii,jj,kk}=fftshift(abs(R.recon)./max(abs(R.recon(:))),1);
            MSE(ii,jj,kk)=mean(mean(abs(squeeze(Recon{ii,jj,kk})-ImRef).^2,1),2);
            SSIM(ii,jj,kk)=ssim(abs(squeeze(Recon{ii,jj,kk})),abs(ImRef))
            % ADD RECON HERE AS WELL FOR IMAGE PRINT
        end
end
end


%% Visualize MSE
figure(101)
subplot(311)
plot(1./accvector,squeeze(MSE(:,:,1)),'.-'); title('NSA1'); legend('uni','peri','center');
xlabel('sampling fraction'); ylabel('MSE')
subplot(312)
plot(1./accvector,squeeze(MSE(:,:,2)),'.-'); title('NSA2');
xlabel('sampling fraction'); ylabel('MSE')
subplot(313)
plot(1./accvector,squeeze(MSE(:,:,3)),'.-'); title('NSA3');xlabel('sampling fraction'); ylabel('MSE')

%% Visualize SSIM
figure(102)
subplot(311)
plot(1./accvector,squeeze(SSIM(:,:,1)),'.-'); title('NSA1'); legend('uni','peri','center');
xlabel('sampling fraction'); ylabel('SSIM')
subplot(312)
plot(1./accvector,squeeze(SSIM(:,:,2)),'.-'); title('NSA2');
xlabel('sampling fraction'); ylabel('SSIM')
subplot(313)
plot(1./accvector,squeeze(SSIM(:,:,3)),'.-'); title('NSA3');xlabel('sampling fraction'); ylabel('SSIM')

%%
I4=[]; I4temp=[];
for ii=1:length(accvector); 
for jj=1:3;
I4temp=[I4temp,abs(Recon{ii,jj,2})];
end;
I4=[I4;I4temp];
clear I4temp; I4temp=[];
end;
figure(104);
imshow(I4,[])
export_fig -native '4_RECONS.eps'
export_fig -native '4_RECONS.png'

xv=170:250; yv=xv;
I5=[]; I5temp=[];
for ii=1:length(accvector); 
for jj=1:3;
I5temp=[I5temp,abs(recon{ii,jj}(xv,yv))];
end;
I5=[I5;I5temp];
clear I5temp; I5temp=[];
end;
figure(105);
imshow(I5,[])
% export_fig -native '5_RECONS.eps'
% export_fig -native '5_RECONS.png'

figure(107)
imshow(abs(ImRef),[]);
% export_fig -native '7_imRef.eps'
% export_fig -native '7_imRef.png'

I5=[]; I5temp=[];
for ii=1:length(accvector); 
for jj=1:3;
I5temp=[I5temp,KD{3,ii,jj}.M.*KD{3,ii,jj}.MNSA];
end;
I5=[I5;I5temp];
clear I5temp; I5temp=[]
end;
figure(108)
cmap=hot(200); cmap=cmap([1,20:200],:);
imshow(I5,[]);colormap(cmap)
% export_fig -native '8_MASKS.eps'
% export_fig -native '8_MASKS.png'

figure(109);
imshow(abs(ImRef(xv,yv)),[]); axis off
% export_fig -native '9_ImRefZ.eps'
% export_fig -native '9_ImRefZ.png'