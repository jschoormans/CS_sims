clear all; close all; clc; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/opt/amc/bart-0.3.00'));% vars;
setenv('TOOLBOX_PATH', '/opt/amc/bart-0.3.00/bin/');
setenv('OMP_NUM_THREADS','4');
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab Collection/utils'));

%VARS
experiment_number=3;

Im=imread('imageblok.png');   %image in image domain

% Im=imread('image.bmp');   %image in image domain
% Im=imread('imagebrain.png');   %image in image domain
ny=size(Im,2)
nz=size(Im,1)
Im=double(Im(:,:,1))./255;

if false %make example image more sparse in W-domainby soft-thresholding --> probably only good for multilevel wavelets!!
wavcoeff=bart('cdf97 3', Im); %wavelet coefficietns;
wavcoeff= bart('threshold 0.1',wavcoeff);
Im=bart('cdf97 3 -i', wavcoeff);
imshow(Im)
sparsity=sum(wavcoeff(:)>0)/(ny*nz)
end


%MAKE SENS MAPS AND K_SPACE (SIMULATION DATA)
sens=ones(1,ny,nz,1);
K=bart('fft 3',Im);

load('mask');
% M=mask';
M=bart(['poisson -Z',num2str(ny),' -Y',num2str(nz),' -z2 -y2 -C50 -v'])
acceleration=sum(M(:)==1)/((ny/2)*(nz/2)*(pi))
% MAKE ELLIPSE MASK 
%
Mfull=genEllipse(ny,nz);

K=squeeze(Mfull').*K;
Ku=squeeze(M).*K;      %undersampled k-space%
Ku=permute(Ku,[3 2 1]);
%CS
ImCS6=bart('pics -RW:7:0:0.05 -d5 -i100',Ku,sens);
ImCS7=bart('pics -RW:7:0:0.05 -i100 ',Ku,sens);

Imlin=ifftshift(ifft2((squeeze(Ku))));

figure(1);
subplot(221); imshow(Im);axis off; title('original');subplot(222); imshow(squeeze(M)'); axis off; title('mask')
subplot(223); imshow(abs(squeeze(Imlin)'),[]);axis off; title('linear recon'); subplot(224); imshow([abs(squeeze(ImCS6))' abs(squeeze(ImCS7))'],[]);axis off; title('CS Recon')
%%

%SIMULATION VARIABLES
% NoiseLevels=[[1:1:9],[10:2:30],[35:5:100]]
NoiseLevels=[1,5,10,20,40]

unders=[0.1:0.1:1]
reg=[1e-2]
iter=1

[Recon,D,D2,D3,acc,masks,MSV]=analyse_CSSNR(Im,K,sens,NoiseLevels,reg,unders,iter);

%%
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures')
dir=mkdir(['7-4-',num2str(experiment_number)])
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/','7-4-',num2str(experiment_number)]
cd(figfolder)

% save variables
save vars.mat
%%
figure(2) %PLOT SNR on x-axis
m=1;n=1;u=1;r=1;i=1;
hold on 
m=1
for u=1:10
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'*-')
end
m=2
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'--')
hold off
title('SSIM: us 1 to 10')
xlabel('SNR')
export_fig -native '2_SSIM.eps'
export_fig -native '2_SSIM.png'


%%
colors=prism(70);
close all
n6=5
figure(3) % PLOT US-RATE ON X AXIS
for r=1:2
subplot(2,2,r)
hold on 
for n=1:n6
plot(acc(:,1), squeeze(mean(D2(1,n,:,r,:),5)),'-*','Color',colors(n,:))
plot(acc(:,1), ones(length(acc),1).*mean(squeeze(mean(D2(2,n,:,r,:),5))),'--','Color',colors(n,:))
end
hold off
title(['SSIM: SNR ',num2str(NoiseLevels(1)),' to ',num2str(NoiseLevels(6))',' reg is: ',num2str(reg(r))])
xlabel('SNR')
ylabel('SSIM')
end

for r=1:2

subplot(2,2,r+2)
hold on 
for n=1:n6
plot(acc(:,1), squeeze(mean(D(1,n,:,r,:),5)),'-*','Color',colors(n,:))
plot(acc(:,1), ones(length(acc),1).*mean(squeeze((mean(D(2,n,:,r,:),5)))),'--','Color',colors(n,:))
end
hold off
title(['MSE: SNR ',num2str(NoiseLevels(1)),' to ',num2str(NoiseLevels(6))',' reg is: ',num2str(reg(r))])
xlabel('SNR')
ylabel('MSE')
end

export_fig -native '3_SNRdep.eps'
export_fig -native '3_SNRdep.png'

%%
n1=10; n6=15;
figure(4) % PLOT US-RATE ON X AXIS
for r=1:2
subplot(2,2,r)
hold on 
for n=n1:n6
plot(acc(:,1), squeeze(mean(D2(1,n,:,r,:),5)),'-*','Color',colors(n,:))
plot(acc(:,1), ones(length(acc),1).*mean(squeeze(mean(D2(2,n,:,r,:),5))),'--','Color',colors(n,:))
end
hold off
title(['SSIM: SNR ',num2str(NoiseLevels(n1)),' to ',num2str(NoiseLevels(n6)),' reg is: ',num2str(reg(r))])
xlabel('SNR')
ylabel('SSIM')
end

for r=1:2

subplot(2,2,r+2)
hold on 
for n=n1:n6
plot(acc(:,1), squeeze(mean(D(1,n,:,r,:),5)),'-*','Color',colors(n,:))
plot(acc(:,1), ones(length(acc),1).*mean(squeeze((mean(D(2,n,:,r,:),5)))),'--','Color',colors(n,:))
end
hold off
title(['MSE: SNR ',num2str(NoiseLevels(n1)),' to ',num2str(NoiseLevels(n6)),' reg is: ',num2str(reg(r))])
xlabel('SNR')
ylabel('MSE')
end

export_fig -native '4_SNRdep.eps'
export_fig -native '4_SNRdep.png'




%%
figure(5) %show the masks used; 
imshow(reshape(masks,[size(masks,1),size(masks,2)*size(masks,3)]))

export_fig -native '5_masks.eps'
export_fig -native '5_masks.png'

%%  for restrictellipse=[1 0], to visualize if this makes a difference!!!
% NEEDS MORE ITERATIONS!!!
% figure(6) % PLOT US-RATE ON X AXIS
% hold on 
% for n=1:6
% plot(acc(:,1), squeeze((D2(1,n,:,r,1))),'-*','Color',colors(n,:))
% plot(acc(:,2), squeeze((D2(1,n,:,r,2))),'-+','Color',colors(n,:))
% plot(acc, ones(length(acc),1).*mean(squeeze(mean(D2(2,n,:,r,:),5))),'--','Color',colors(n,:))
% end
% hold off
% title('SSIM-detail')

%%
 % show images next to each other
for n=5;
close all

FigHandle = figure(7)
set(FigHandle, 'Position', [0,0, size(Recon,1)*10, 2*size(Recon,2)]);

for r=1:length(reg)
    subplot(length(reg),1,r)
    iter=1;
    PLOT=[]
    for u=1:10
        PLOT=[PLOT,squeeze(Recon(:,:,1,n,u,r,iter))'];
        x=10; y=10;
    end
    hold on
    imshow(PLOT,[])
    for u=1:10
    text(x+(u-1)*size(Recon,1),y,[num2str(unders(u)),' SSIM:',num2str(D2(1,n,u,r,iter))],'Color','white')
    end
    hold off
    axis off
    title(['SNR is: ',num2str(NoiseLevels(n)),' reg is: ',num2str(reg(r))])
end

name=['7_reconsSNR',num2str(NoiseLevels(n)),'.eps']
name2=['7_reconsSNR',num2str(NoiseLevels(n)),'.png']

export_fig(name)
export_fig(name2)

end
%%
figure(99)
plot(unders,abs((MSV)),'.b'); %mean absolute value of the k-space signal
title('mean SNR values per measurement')
xlabel('undersampling value')

export_fig -native '99_MSV.eps'
export_fig -native '99_MSV.png'