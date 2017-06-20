clear all; close all; clc; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))

addpath(genpath('/opt/amc/bart-0.3.00'));% vars;
setenv('TOOLBOX_PATH', '/opt/amc/bart-0.3.00/bin/');
setenv('OMP_NUM_THREADS','4');
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab Collection/utils'));

CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/',CCC]
dir=mkdir(figfolder);cd(figfolder)


%VARS

% Im=imread('imageblok.png');   %image in image domain

Im=imread('image.bmp');   %image in image domain
% Im=imread('imagebrain.png');   %image in image domain
Im=imread('brain.png');   %image in image domain
% Im=imread('lelijk.png');   %image in image domain
% Im=imread('earth.jpg');

Im=double(Im(:,:,1))./255;
ny=size(Im,2)
nz=size(Im,1)


%% MAKE SENS MAPS AND K_SPACE (SIMULATION DATA)
K=bart('fft 3',Im);


[K]=genPhantomKspace(256,1)
ny=size(K,2)
nz=size(K,1)
sens=ones(ny,nz); %change to sense based on phantoms...

%%
% load('mask');
% M=mask';
% M=bart(['poisson -Z',num2str(ny),' -Y',num2str(nz),' -z2 -y2 -C20 -v'])
[pdf,val] = genPDF(size(K),5,1/4,2,0.05,0);
M=double(genSampling(pdf,10,10));
            
acceleration=sum(M(:)==1)/((ny/2)*(nz/2)*(pi))
% MAKE ELLIPSE MASK 
%
Mfull=genEllipse(ny,nz);

K=squeeze(Mfull').*K;
Ku=squeeze(M).*K;      %undersampled k-space%
Ku=permute(Ku,[3 2 1]);
%CS
ImCS6=bart('pics -RW:7:0:0.05 -d5 -i100',Ku,sens);
ImCS6=ImCS6./max(ImCS6(:))
mu=50e-1
ImCS7=mrics(squeeze(M)',squeeze(Ku),mu,mu,mu/100,30,5);

ImCS7=ifftshift(permute(ImCS7,[3 1 2])./max(ImCS7(:)));

Imlin=ifftshift(ifft2((squeeze(Ku))));

cd(figfolder)
figure(1);
subplot(221); imshow(Im);axis off; title('original');subplot(222); imshow(squeeze(M)'); axis off; title('mask')
subplot(223); imshow(abs(squeeze(Imlin)'),[]);axis off; title('linear recon'); subplot(224); imshow([abs(squeeze(ImCS6))' abs(squeeze(ImCS7))'],[]);axis off; title('CS Recon')

export_fig -native '1.eps'
export_fig -native '1.png'


%%

%SIMULATION VARIABLES
% NoiseLevels=[[1:1:9],[10:2:30],[35:5:100]]
% NoiseLevels=[10]
NoiseLevels=[1e-5 5e-5 1e-4 2.5e-4 5e-4 7.5e-4 1e-3]
reg=[1e-2]
iter=1

% unders=[1 2 3 5 10]
% [Recon,D,D2,D3,acc,masks,MSV]=analyse_CSSNRf(Im,K,sens,NoiseLevels,reg,iter);%DOES
% NSA=R oversampling
unders=[6]
% [Recon,D,D2,D3,acc,masks,MSV]=analyse_CSSNR(Im,K,sens,NoiseLevels,reg,1./unders,iter); 

[Recon,D,D2,D3,acc,masks,MSV,Ns,MNSA]=analyse_CSSNR_varNSA(Im,K,sens,NoiseLevels,reg,1./unders,iter); 

% save variables
save vars.mat
%%
figure(2) %PLOT SNR on x-axis
m=1;n=1;u=1;r=1;i=1;
hold on 
m=1
for u=1:length(unders)
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'k*-')
end
m=2
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'b*-')
m=4
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'r*-')
m=5
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'g.-')
m=6
plot(NoiseLevels, squeeze(mean(D2(m,:,u,r,:),5)),'y*-')


hold off
title('SSIM: us 1 to 10')
xlabel('SNR')
legend('undersVar','full','undersUni','no NSA, R=4','var reversed')
export_fig -native '2_SSIM.eps'
export_fig -native '2_SSIM.png'


%%

%{
colors=jet(5);
close all
n6=2
figure(3) % PLOT US-RATE ON X AXIS
for r=1:length(reg)
subplot(2,length(reg),r)
hold on 
for n=1:n6
plot(acc(:,1), squeeze(mean(D2(1,n,:,r,:),5)),'-*','Color',colors(n,:))
plot(acc(:,1), ones(length(acc),1).*mean(squeeze(mean(D2(2,n,:,r,:),5))),'--','Color',colors(n,:))
% plot(acc(:,1), ones(length(acc),1).*mean(squeeze(mean(D2(3,n,:,r,:),5))),'.-','Color',colors(n,:))
end
hold off
title(['SSIM: SNR ',num2str(NoiseLevels(1)),' to ',num2str(NoiseLevels(n6)),' reg is: ',num2str(reg(r))])
xlabel('acc')
ylabel('SSIM')
end

for r=1:length(reg)
subplot(2,length(reg),r+length(reg))
hold on 
for n=1:n6
plot(acc(:,1), squeeze(mean(D(1,n,:,r,:),5)),'-*','Color',colors(n,:))
plot(acc(:,1), ones(length(acc),1).*mean(squeeze((mean(D(2,n,:,r,:),5)))),'--','Color',colors(n,:))
end
hold off
title(['MSE: SNR ',num2str(NoiseLevels(1)),' to ',num2str(NoiseLevels(n6)),' reg is: ',num2str(reg(r))])
xlabel('acc')
ylabel('MSE')
end

export_fig -native '3_SNRdep.eps'
export_fig -native '3_SNRdep.png'
%}
%%
%{
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
%}



%%
figure(5) %show the masks used; 
imshow(MNSA.*masks,[]); colormap('hot');
colorbar
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
 %{
for n=1:length(NoiseLevels);
close all

FigHandle = figure(7)
% set(FigHandle, 'Position', [1,1, 1000,100]);

for r=1:length(reg)
    iter=1;
    PLOT=[]
    for u=1:length(unders)
        PLOT=[PLOT,squeeze(Recon(:,:,1,n,u,r,iter))'];
        x=10; y=10;
    end
    hold on
    imshow(abs(PLOT),[],'InitialMagnification',50)
    for u=1:length(unders)
    text(x+(u/2-1/2)*size(Recon,1),y,[num2str(unders(u)),' SSIM:',num2str(D2(1,n,u,r,iter))],'Color','white')
    end
    hold off
    axis off
    title(['N is: ',num2str(NoiseLevels(n)),' reg is: ',num2str(reg(r))])
end

name=['7_reconsSNR',num2str(NoiseLevels(n)),'.eps']
name2=['7_reconsSNR',num2str(NoiseLevels(n)),'.tiff']

export_fig(name,'-native')
% export_fig(name2,'-native')
imwrite(PLOT./max(PLOT(:)),name2);

end
%}
%%

for n=1:length(NoiseLevels);
close all

FigHandle = figure(8)
set(FigHandle, 'Position', [0,0, size(Recon,1)*2, 2*size(Recon,2)]);

for r=1:length(reg)
    iter=1;
    PLOT=[]
        PLOT=[squeeze(Recon(:,:,1,n,u,r,iter))',squeeze(Recon(:,:,2,n,u,r,iter))';squeeze(Recon(:,:,4,n,u,r,iter))',squeeze(Recon(:,:,6,n,u,r,iter))'];
        x=10; y=10;
    hold on
    imshow(abs(PLOT),[])
    for u=1
    text(x,y,['R=',num2str(unders(u)),'varNSA SSIM:',num2str(D2(1,n,u,r,iter))],'Color','white')
        text(x+size(Recon,1),y,['full SSIM:',num2str(D2(2,n,u,r,iter))],'Color','white')
                text(x,y+size(Recon,2),['R=',num2str(unders(u)),' uniNSA SSIM:',num2str(D2(4,n,u,r,iter))],'Color','white')

%         text(x+size(Recon,1),y+size(Recon,2),['R=',num2str(unders(u)),' NSA=1 SSIM:',num2str(D2(5,n,u,r,iter))],'Color','white')
        text(x+size(Recon,1),y+size(Recon,2),['R=',num2str(unders(u)),' NSA=var rev SSIM:',num2str(D2(6,n,u,r,iter))],'Color','white')

    end
    hold off
    axis off
    title(['varNSA versus full versus uniNSA; N is: ',num2str(NoiseLevels(n)),' reg is: ',num2str(reg(r))])
end

name=['8_reconsSNR',num2str(NoiseLevels(n)),'.eps']
name2=['8_reconsSNR',num2str(NoiseLevels(n)),'.png']

export_fig(name)
export_fig(name2)

end


%%
figure(10);
% MNSA

%%
figure(99)
plot(unders,abs((MSV)),'.b'); %mean absolute value of the k-space signal
title('mean SNR values per measurement')
xlabel('undersampling value')

export_fig -native '99_MSV.eps'
export_fig -native '99_MSV.png'

%% autocorrelation and noise frequencies
clear nps_2Dr_measured
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/nps_package'))
for m=[1,2,3,4,5,6]
    n=34;
RECON=abs(squeeze(Recon(101:500,101:500,m,n,u,r,iter)).');
RECON=RECON./max(RECON(:));
% RECON=randn(500,500); %end result should be WN !
figure(98); subplot(2,3,m);imshow(RECON,[]); pause(1)
[nps_2D_measured,fx]=calc_digital_nps(RECON,2,1,1,0);
[nps_2Dr_measured(m,:), f_2Dr] = cart2rad(nps_2D_measured, fx, 2, uniquify);
nps_2Dr_measured(m,1:5)=nan(1,5)
end

export_fig -native '98_NPS_ROI.eps'
export_fig -native '98_NPS_ROI.png'
%%
figure(10)
plot(f_2Dr,log(nps_2Dr_measured))
legend('var','full','linear','uniNSA','noNSA','varRev')
export_fig -native '10_NPS.eps'
export_fig -native '10_NPS.png'



%%
clear nps_2Dr_measured_conv
for m=1:6
nps_2Dr_measured_conv(m,:)=conv(nps_2Dr_measured(m,:),ones(1,30),'same');
end
figure(11)
plot(f_2Dr,log(nps_2Dr_measured_conv))
legend('var','full','linear','uniNSA','noNSA','varRev')

export_fig -native '11_NPS_sm.eps'
export_fig -native '11_NPS_sm.png'

%%
figure(12);
hold on
plot(f_2Dr,nps_2Dr_measured_conv(1,:)./nps_2Dr_measured_conv(2,:),'r'); 
plot(f_2Dr,nps_2Dr_measured_conv(4,:)./nps_2Dr_measured_conv(2,:),'k');
plot(f_2Dr,nps_2Dr_measured_conv(5,:)./nps_2Dr_measured_conv(2,:),'y');
plot(f_2Dr,nps_2Dr_measured_conv(6,:)./nps_2Dr_measured_conv(2,:),'b');
hold off
xlabel('spatial frequency')
ylabel('noise power relative to full sampling')
legend('varNSA / noNSA','uniNSA/nonNSA','no NSA','varRevNSA / noNSA')

export_fig -native '12_NPS_relative.eps'
export_fig -native '12_NPS_relative.png'

