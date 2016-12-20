% 14-9-2016
%GOAL PLOT NSAxAXX (x -axis) vs sigmoid width 
%


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
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.0002;	% Weight for Transform L1 penalty (0.001) (0.003 too much)
Itnlim = 8;		% Number of iterations (8)
innerloopiter=3 %(3)

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

for ch=[1:8];
K1=K(:,:,:,ch,:,:,:,:,:,:,:,:,:);

%% make undersampled k-spaces (TO CHANGE: SHOULD BE THE SAME FOR ALL CHANNELS!!)
disp('making the undersampled k-spaces')
accvector=[1,2,3,4,5,6,7]
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
clear MSE SSIM hfen S

coordsx=[210:220]; coordsy=[80:150]

if false %one line
for ii=1:length(accvector);
    for jj=1:length(noisevector)   
        MSE(2,ii,jj)=mean(mean(abs(squeeze(R{2,ii,jj})-ImRef).^2,1),2);
%       SSIM(2,ii,jj)=ssim(abs(squeeze(R{2,ii,jj})),abs(ImRef));
        hfen(2,ii,jj)=HFEN(abs(squeeze(R{2,ii,jj})),abs(ImRef));
%       PSNR(2,ii,jj)=psnr(abs(squeeze(R{2,ii,jj})),abs(ImRef))

        Edge=mean(abs(squeeze(R{2,ii,jj}(coordsx,coordsy))),2)
        S(2,ii,jj)=SigmoidFitting(Edge.',[],1,0)
    end
end
else %for all coordsy 
    for coordsxi=coordsx
        
% FIND CENTER FOR THIS X
            Y=abs(squeeze(R{2,5,3}(coordsxi,coordsy))) %best one
            Yn=(Y-min(Y))/(max(Y)-min(Y));
            [val1,ind1]=find((Yn-0.5)>0,1,'first');
            ind0=ind1-1;
            val1=Yn(ind1)-0.5;
            val0=Yn(ind0)-0.5;
            Center(coordsxi)=ind0+abs(val0)/(val1-val0);
        
        for ii=1:length(accvector);
            for jj=1:length(noisevector)
                MSE(2,ii,jj)=mean(mean(abs(squeeze(R{2,ii,jj})-ImRef).^2,1),2);
                %       SSIM(2,ii,jj)=ssim(abs(squeeze(R{2,ii,jj})),abs(ImRef));
                hfen(2,ii,jj)=HFEN(abs(squeeze(R{2,ii,jj})),abs(ImRef));
                %       PSNR(2,ii,jj)=psnr(abs(squeeze(R{2,ii,jj})),abs(ImRef))
                Edge=abs(squeeze(R{2,ii,jj}(coordsxi,coordsy)));
                S(2,ii,jj,coordsxi)=SigmoidFitting(Edge,Center(coordsxi),1,0);
                
            end
        end
    end
    
end
mean(S(:,:,:,coordsx),4);
    
%% MAKE PLOTS
figure(2)
imshow(abs(squeeze(R{2,5,3}(coordsx,coordsy))))
figure(4)
plot(coordsx,Center(coordsx),'.')
figure(3)
 plot(squeeze(mean(S(2,:,:,coordsx),4)),'.-')
 legend('uni','more-outside','more-inside')
 xlabel('NSA-acc')
 ylabel('width of sigmoid (pix)')
 %% PLOT RECON
figure(107) %REFERENCE IMAGE
imshow(abs(R{2,5,3}),[]);  
 %% PLOT RECONS FOR ALL THINGS (5x3)
 nx=100:200
 ny=100:300
 
I4=[]; I4temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I4temp=[I4temp,abs(R{2,ii,jj}(nx,ny))./max(abs(R{2,ii,jj}(:)))];
end;
I4=[I4;I4temp];
clear I4temp; I4temp=[];
end;
figure(104);
imshow(I4,[])
%% add other error functions again
 
 
