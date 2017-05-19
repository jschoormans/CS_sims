%% Estimate errors
clear MSE SSIM hfen
ImRef=bart('fft -i 7',ImRef);
ImRef=bart('resize -c 0 1024 1 1024',ImRef);
ImRef=bart('fft 7',ImRef);
ImRef=abs(ImRef)./max(abs(ImRef(:)));


for ii=1:length(accvector);
    for jj=1:3
        recon{ii,jj}=ifftshift(R{2,ii,jj}.recon,1);
        MSE(2,ii,jj)=mean(mean(abs(squeeze(recon{ii,jj})-ImRef).^2,1),2)
        SSIM(2,ii,jj)=ssim(abs(squeeze(recon{ii,jj})),abs(ImRef))
%         hfen(2,ii,jj)=HFEN(abs(squeeze(recon)),abs(ImRef))
    end
end


%% Visualize
mkdir([experimentfolder,'/',filename])
cd([experimentfolder,'/',filename])
figure(101)
plot(1./accvector,squeeze(MSE(2,:,:)),'.-')
title('MSE')
xlabel('sampling fraction')
ylabel('MSE')
% export_fig -native '1_MSE.eps'
% export_fig -native '1_MSE.png'

figure(102)
plot(1./accvector,squeeze(SSIM(2,:,:)),'.-')
title('SSIM')
xlabel('acceleration')
ylabel('SSIM')
% export_fig -native '2_SSIM.eps'
% export_fig -native '2_SSIM.png'

% figure(103)
% plot(1./accvector,squeeze(hfen(2,:,:)),'.-')
% title('HFEN')
% xlabel('acceleration')
% ylabel('HFEN')
% export_fig -native '3_HFEN.eps'
% export_fig -native '3_HFEN.png'

I4=[]; I4temp=[];
for ii=1:length(accvector); 
for jj=1:3;
I4temp=[I4temp,abs(recon{ii,jj})];
end;
I4=[I4;I4temp];
clear I4temp; I4temp=[];
end;
figure(104);
imshow(I4,[])
% export_fig -native '4_RECONS.eps'
% export_fig -native '4_RECONS.png'

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


figure(106)
plot([1:3],squeeze(SSIM(2,:,:)).','.-')
title('SSIM')
xlabel('averaging type')
ylabel('SSIM')
% export_fig -native '6_SSIM.eps'
% export_fig -native '6_SSIM.png'

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