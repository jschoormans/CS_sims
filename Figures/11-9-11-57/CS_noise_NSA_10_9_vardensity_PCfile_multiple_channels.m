%% 10-9 scans: PC MATLAB FILE

clear all; close all; clc; 
disp('adding paths and setting up stuff...')
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection/exportfig'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection/tightfig/'))
% addpath(genpath('/opt/amc/matlab/toolbox\MRecon-3.0.519'))
addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations/'))
% addpath(genpath('C:\Users\jschoormans\Dropbox\phD\bart-0.3.01\bart-0.3.01')); vars
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab'))
CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['L:\basic\divi\Projects\cosart\CS_simulations\Figures\',CCC]
dir=mkdir(figfolder);cd(figfolder)

%% defining parameters 
clear dir;
disp('setting parameters')
P=struct;
P.chan=1;
P.folder='L:\basic\divi\Projects\cosart\scans\10-9-grapefruit-CS-Graz\'
cd(P.folder)
files=dir('*.raw');
P.file=files(16).name;
load('K15_ch1_ringing.mat')
cd(figfolder) %CD figfolder to save figures in good folder


%% make ref image (should be of all channels!)
K1=K(:,:,:,1,:,:,:,:,:,:,:,:,:);
disp('making ref image')
reg=0.02;
Kref=mean(K1,12); %still use coils though
ImRef=fftshift(fft2(Kref(:,:,1,1)),1);
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
figure(1000); imshow(abs(ImRef(:,:,1,1)),[])




%% loop over channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0.005; 	% Weight for TV penalty
xfmWeight = 0.00;	% Weight for Transform L1 penalty (0.001) (0.003 too much)
Itnlim = 8;		% Number of iterations (8)
innerloopiter=3 %(3)

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

for ch=[1:8];
K1=K(:,:,:,ch,:,:,:,:,:,:,:,:,:);

%% make undersampled k-spaces (TO CHANGE: SHOULD BE THE SAME FOR ALL CHANNELS!!)
disp('making the undersampled k-spaces')
accvector=[1,2,3,4,5]
% accvector=[1,4,8]
noisevector=[5,6,7] % in this case the different W 
for jj=1:length(noisevector)
for ii=1:length(accvector);
P.jjj=noisevector(jj); %weighting scheme (3=normal)
P.usedyns=1; %for example
P.acc=accvector(ii);
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspacefromdynamics(K1,P);
end
end



%% CONJUGATE GRADIENT RECONS

% initialize Parameters for reconstruction
param = init;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =0;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
N=[size(K1,1) size(K1,2)]

for ii=1:length(accvector);
    for jj=1:length(noisevector)
%         if KD{3,ii,jj}.jjj==3 || KD{3,ii,jj}.jjj==4 ||KD{3,ii,jj}.jjj==5
%             R{2,ii,jj}=fftshift(bart(['pics -RW:7:0:',num2str(reg),' -d5 -S -e -i100'],KD{2,ii,jj},sens),2)
%             R{2,ii,jj}=abs(R{2,ii,jj});
%             R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:)); %scaling!
%         else
            %generate Fourier sampling operator
            param.data=KD{2,ii,jj};
            FT = p2DFT(KD{3,ii,jj}.M, N, 1, 2);
            param.FT = FT;
            param.xfmWeight=xfmWeight;
            im_dc2 = FT'*(param.data.*KD{3,ii,jj}.M./KD{3,ii,jj}.pdf); %linear recon; scale data to prevent low-pass filtering
            res = XFM*im_dc2;
            param.V=(KD{3,ii,jj}.MNSA.*KD{3,ii,jj}.M);
            param.Debug=0;
            param.xfmWeight=xfmWeight*(mean(param.V(KD{3,ii,jj}.M~=0)));
            
            for n=1:innerloopiter
                res = fnlCg_test(res,param);
                R{2,ii,jj} = XFM'*res;
            end
            R{2,ii,jj}=fftshift(abs(R{2,ii,jj}),2);
            R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:));
            RT(ii,jj,ch,:,:)=R{2,ii,jj}; %has all channels, accelerations and sampling patterns

    end
end

end
%% combine coils 
for ii=1:length(accvector);
    for jj=1:length(noisevector)
            R{2,ii,jj}=squeeze(sqrt(sum(RT(ii,jj,:,:,:).^2,3)))
    end
end


%% Estimate errors
clear MSE SSIM hfen S Center

if false %one line
for ii=1:length(accvector);
    for jj=1:length(noisevector)
        MSE(2,ii,jj)=mean(mean(abs(squeeze(R{2,ii,jj})-ImRef).^2,1),2);
%       SSIM(2,ii,jj)=ssim(abs(squeeze(R{2,ii,jj})),abs(ImRef));
        hfen(2,ii,jj)=HFEN(abs(squeeze(R{2,ii,jj})),abs(ImRef));
%       PSNR(2,ii,jj)=psnr(abs(squeeze(R{2,ii,jj})),abs(ImRef))

        coordsx=95:105; coordsy=250
        Edge=mean(abs(squeeze(R{2,ii,jj}(coordsx,coordsy))),2)
        S(2,ii,jj)=SigmoidFitting(Edge.',Center,1,0)
        
    end
end

else %for all coordsy 
    coordsyvec=[240:246];
    coordsx=85:125;


    for coordsy=coordsyvec
        
        Y=mean(abs(squeeze(R{2,5,3}(coordsy,coordsx))),1)
        Yn=(Y-min(Y))/(max(Y)-min(Y));
        [ind1]=find((Yn-0.5)>0,1,'first');
        ind0=ind1-1;
        val1=Yn(ind1)-0.5;
        val0=Yn(ind0)-0.5;
        Center(coordsy)=ind0+abs(val0)/(val1-val0)
    
        
        for ii=1:length(accvector);
            for jj=1:length(noisevector)
                MSE(2,ii,jj)=mean(mean(abs(squeeze(R{2,ii,jj})-ImRef).^2,1),2);
                %       SSIM(2,ii,jj)=ssim(abs(squeeze(R{2,ii,jj})),abs(ImRef));
                hfen(2,ii,jj)=HFEN(abs(squeeze(R{2,ii,jj})),abs(ImRef));
                %       PSNR(2,ii,jj)=psnr(abs(squeeze(R{2,ii,jj})),abs(ImRef))
                
                Edge=abs(squeeze(R{2,ii,jj}(coordsy,coordsx)))
                S(2,ii,jj,coordsy)=SigmoidFitting(Edge,Center(coordsy),1,2)
                
            end
        end
    end
mean(S(:,:,:,coordsyvec),4)
    
end


figure(106)
plot(accvector(1:end),squeeze(median(S(2,1:end,:,coordsyvec),4)),'.-')
title('sigmoid')
xlabel('NSA')
ylabel('width (pix)')
legend('uni','more outside','center','Location','NorthWest')
% export_fig -native '6_SSIM.eps'
% export_fig -native '6_SSIM.png'

figure(199)
imshow(R{2,ii,jj}(coordsyvec,coordsx))

%% Visualization code 
%{
cd(figfolder)
figure(101)
set(gcf,'color','w');
subplot(221)
plot(1./accvector,squeeze(MSE(2,:,:)),'.-')
title('MSE')
xlabel('sampling fraction')
ylabel('MSE')
legend('5','6','7')

subplot(222)
plot(1./accvector,squeeze(SSIM(2,:,:)),'.-')
title('SSIM')
xlabel('sampling fraction')
ylabel('SSIM')

subplot(223)
plot(1./accvector,squeeze(hfen(2,:,:)),'.-')
title('HFEN')
xlabel('sampling fraction')
ylabel('HFEN')

subplot(224)
plot(1./accvector,squeeze(PSNR(2,:,:)),'.-')
title('PSNR')
xlabel('sampling fraction')
ylabel('PSNR')

legend('sampling 1','sampling 2','sampling 3')
export_fig -native '1_errors.eps'
export_fig -native '1_errors.png'

figure(102)
subplot(221)
plot(noisevector,squeeze(MSE(2,:,:)).','.')
title('MSE')
xlabel('sampling pattern')
ylabel('MSE')

subplot(222)
plot(noisevector,squeeze(SSIM(2,:,:)).','.')
title('SSIM')
xlabel('sampling pattern')
ylabel('SSIM')

subplot(223)
plot(noisevector,squeeze(hfen(2,:,:)).','.')
title('HFEN')
xlabel('sampling pattern')
ylabel('HFEN')

subplot(224)
plot(noisevector,squeeze(PSNR(2,:,:)).','.')
title('PSNR')
xlabel('sampling pattern')
ylabel('PSNR')

export_fig -native '2_errors.eps'
export_fig -native '2_errors.png'

figure(103)
% QQ PLOT
for ii=1:length(accvector); 
for jj=1:length(noisevector);
    subplot(length(accvector),length(noisevector),(ii-1)*length(noisevector)+jj)
temp=real(squeeze(R{2,ii,jj})-ImRef);
qqplot(temp(:))
end
end
export_fig -native '3_QQ.eps'
export_fig -native '3_QQ.png'
%}
%%

I4=[]; I4temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I4temp=[I4temp,abs(R{2,ii,jj})./max(abs(R{2,ii,jj}(:)))];
end;
I4=[I4;I4temp];
clear I4temp; I4temp=[];
end;
figure(104);
imshow(I4,[])
% export_fig -native '4_RECONS.eps'
% export_fig -native '4_RECONS.png'

%%
xv=[200:300];
yv=[200:300];
I5=[]; I5temp=[];

for ii=1:length(accvector); 
for jj=1:length(noisevector);
I5temp=[I5temp,abs(R{2,ii,jj}(xv,yv))./max(abs(R{2,ii,jj}(:)))];
end;
I5=[I5;I5temp];
clear I5temp; I5temp=[];
end;

figure(105);
imshow(I5,[])
axis off

% export_fig -native '5_RECONS.eps'
% export_fig -native '5_RECONS.png'

%%
figure(107) %REFERENCE IMAGE
imshow(abs(ImRef),[]);
% export_fig -native '7_imRef.eps'
% export_fig -native '7_imRef.png'

% PLOT MASKS
I5=[]; I5temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I5temp=[I5temp,KD{3,ii,jj}.M.*KD{3,ii,jj}.MNSA];
end;
I5=[I5;I5temp];
clear I5temp; I5temp=[]
end;
figure(108)
cmap=hot(200); cmap=cmap([1,20:200],:);
imshow(I5,[]);colormap(cmap)
colorbar
% export_fig -native '8_MASKS.eps'
% export_fig -native '8_MASKS.png'

figure(109); % ZOOMED IN OF REFERENCE IMAGE
imshow(abs(ImRef(xv,yv)),[]); axis off
% export_fig -native '9_ImRefZ.eps'
% export_fig -native '9_ImRefZ.png'
 %%
 I10=[]; I10temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I10temp=[I10temp,abs(R{2,ii,jj})./max(abs(R{2,ii,jj}(:)))-flipdim(flipdim(abs(ImRef)./max(abs(ImRef(:))),2),1)];
end;
I10=[I10;I10temp];
clear I4temp; I4temp=[];
end;
figure(110);
imshow(I10,[])

%% SHOW 1/3 OF IMAGE NEXT TO EACH OTHER
xv=[1:512];
yv=[50:200];

ii=4
I5temp=[abs(R{2,1,1}(xv,yv))./(0.5.*max(abs(R{2,1,1}(:))))]
for jj=1:length(noisevector);
I5temp=[I5temp,abs(R{2,ii,jj}(xv,yv))./max(abs(R{2,ii,jj}(:)))];
end;
I5=[I5temp];
clear I5temp; 

figure(111);
imshow(I5,[])
axis off




